#==============================================================================
# LOSS IN EXPECTATION OF LIFE (LEL) FOR 5-YEAR SURVIVORS OF CHILDHOOD CANCER
# Flexible parametric relative survival model (rstpm2::stpm2) on the attained-age
# timescale with delayed entry (age at 5-yr survival) and expected hazard from
# England & Wales all-cause rates ($rates/allcauses.dta).
#
# Inputs : $temp/x-mort3-prepforstset.dta  (one row per survivor, n=34,488)
#          $rates/allcauses.dta            (sex x ageband x yeargrp, newrate)
# Outputs: $output/table_lel.csv            (by cancer type, a0 = 20/30/40/50)
#          $output/table_lel_bydecade.csv   (by decade of diagnosis, a0 = 20/40)
#          $output/table_lel_bysex.csv      (sex-specific, overall and by diag)
#          $output/fig_lel_bydiag.png        (LEL at age 20 by cancer type)
#          $output/fig_le_bydecade.png       (survivor vs expected LE at age 20)
#          $output/lel_model_checks.txt      (expected-deaths check, SMR, fit stats)
#          scratch/x-mort3-lel-split.rds     (split dataset, reused by nothing else)
#
# Method
#  - Follow-up split by 5-y age band and single calendar year (as the Stata
#    pipeline mort3-02a), merged to newrate (annual E&W all-cause rate).
#  - stpm2: log cumulative excess hazard = spline(log age, df=5) + diag + sex +
#    decade (+ time-dependent diag effect if it converges / improves AIC).
#    With left truncation the likelihood only needs the expected hazard at each
#    record's exit, so fitting on the split records is equivalent to fitting on
#    one record per person with newrate at the exit age/year; the latter is used
#    for the final fit (identical likelihood, 30x faster) and both are compared.
#  - Expected residual LE from a0: 2022 E&W period life table (sex-specific),
#    5-y band rates interpolated log-linearly between band midpoints to single
#    ages, 85+ band placed at age 90, log-linear extrapolation beyond 90,
#    survival integrated to age 100.
#  - Survivor residual LE from a0: integral_{a0}^{100} S_exp(a|a0) R(a)/R(a0) da
#    where R(a) = exp(-exp(eta(a))) is the model relative survival; the model is
#    extrapolated beyond the oldest observed ages (~94; sparse above ~75).
#  - LEL = expected - survivor LE; % lost = LEL / expected LE.
#  - 95% CIs by parametric simulation: 500 draws from MVN(coef, vcov).
#  - Overall / by-diag / by-decade estimates are standardised to the cohort:
#    average of profile-specific LE over the (diag, sex, decade) distribution of
#    survivors under follow-up at age a0 (falling back to the group's whole
#    distribution if nobody in the group is under follow-up at a0).
#
# RCR (Claude-assisted), 23 Oct 2026
#==============================================================================
suppressPackageStartupMessages({
  library(haven); library(data.table); library(rstpm2); library(MASS); library(ggplot2)
})
set.seed(20261023)

od      <- "/Users/raoulreulen/Library/CloudStorage/OneDrive-UniversityofBirmingham/mort3"
rates   <- file.path(od, "rates")
output  <- file.path(od, "output")
temp    <- "/Users/raoulreulen/Documents/TEMPWORK"
scratch <- "/private/tmp/claude-501/-Users-raoulreulen-Library-CloudStorage-OneDrive-UniversityofBirmingham-mort3/62a57f68-e958-4b8b-adf8-b5fdd0bba6df/scratchpad"
checkfile <- file.path(output, "lel_model_checks.txt")
cat("", file = checkfile)
note <- function(...) { msg <- sprintf(...); cat(msg, "\n"); cat(msg, "\n", file = checkfile, append = TRUE) }

NSIM   <- 500
A0S    <- c(20, 30, 40, 50)
AMAX   <- 100
STEP   <- 0.1          # integration step (years)
DIAGLAB <- c("leukaemia","Hodgkin lymphoma","non-Hodgkin lymphoma","CNS","neuroblastoma",
             "heritable retinoblastoma","non-heritable retinoblastoma","Wilms","bone","soft tissue sarcoma","other")
DECLAB  <- c("1900"="<1970","1970"="1970-79","1980"="1980-89","1990"="1990-99","2000"="2000-06")

#------------------------------------------------------------------------------
# 1. Data and split by age band x calendar year
#------------------------------------------------------------------------------
d <- as.data.table(read_dta(file.path(temp, "x-mort3-prepforstset.dta")))
d <- d[, .(indexno, sex = as.integer(sex), diag = as.integer(diag), decdxcat = as.integer(decdxcat),
           dob, doe, dox, event = as.integer(!is.na(allcauses) & allcauses == 1))]
d[, age0 := as.numeric(doe - dob) / 365.25]
d[, age1 := as.numeric(dox - dob) / 365.25]
stopifnot(all(d$age1 > d$age0))
# calendar time in fractional years (Stata: after(time=d(1/1/1900)) at(45(1)122) => year at record start)
yrfrac <- function(dt) { y <- as.integer(format(dt, "%Y")); y + as.numeric(dt - as.Date(paste0(y, "-01-01"))) / (365 + (y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0))) }
d[, cal0 := yrfrac(doe)]
d[, age_at_cal0 := age0]  # age at doe

rt <- as.data.table(read_dta(file.path(rates, "allcauses.dta")))[, .(sex = as.integer(sex), ageband, yeargrp, newrate)]
setkey(rt, sex, ageband, yeargrp)

split_fu <- function(d) {
  # cut points in age: age-band boundaries 5,10,...,85 and the age at each 1 Jan
  d[, pid := .I]
  cuts <- d[, {
    bands <- seq(5, 85, 5)
    firstjan <- ceiling(cal0)                          # next 1 Jan after entry
    jans <- seq(firstjan, floor(cal0 + (age1 - age0)) + 1, by = 1)
    ages_jan <- age0 + (jans - cal0)
    cp <- sort(unique(c(age0, bands[bands > age0 & bands < age1], ages_jan[ages_jan > age0 & ages_jan < age1], age1)))
    .(age_start = cp[-length(cp)], age_end = cp[-1])
  }, by = .(pid, indexno, sex, diag, decdxcat, age0, age1, cal0, event)]
  cuts[, last := age_end == age1]
  cuts[, event := as.integer(event == 1 & last)]
  cuts[, ageband := pmin(85, 5 * floor(age_start / 5))]
  cuts[, yeargrp := floor(cal0 + (age_start - age0) + 1e-9)]
  cuts[, yeargrp := pmin(2024, pmax(1950, yeargrp))]
  cuts <- merge(cuts, rt, by = c("sex", "ageband", "yeargrp"), all.x = TRUE)
  stopifnot(!anyNA(cuts$newrate))
  setorder(cuts, pid, age_start)
  cuts
}
t0 <- Sys.time()
s <- split_fu(copy(d))
note("Split dataset: %d records, %d persons, %.1f s", nrow(s), uniqueN(s$pid), as.numeric(Sys.time() - t0, units = "secs"))
s[, py := age_end - age_start]
s[, e := py * newrate]
O <- sum(s$event); E <- sum(s$e); PY <- sum(s$py)
note("Observed deaths = %d, expected deaths = %.1f (target ~861.5), person-years = %.0f, SMR = %.2f",
     O, E, PY, O / E)
stopifnot(abs(E - 861.5) / 861.5 < 0.03)
saveRDS(s, file.path(scratch, "x-mort3-lel-split.rds"))

# One record per person with newrate at exit age/year (likelihood-equivalent)
p1 <- s[last == TRUE, .(pid, indexno, sex, diag, decdxcat, age_start = age0, age_end = age1, event, newrate)]
stopifnot(nrow(p1) == nrow(d))
p1[, `:=`(fdiag = factor(diag, levels = 1:11, labels = DIAGLAB),
          fsex  = factor(sex, levels = 1:2, labels = c("male", "female")),
          fdec  = factor(decdxcat, levels = c(1900, 1970, 1980, 1990, 2000), labels = DECLAB))]
s[, `:=`(fdiag = factor(diag, levels = 1:11, labels = DIAGLAB),
         fsex  = factor(sex, levels = 1:2, labels = c("male", "female")),
         fdec  = factor(decdxcat, levels = c(1900, 1970, 1980, 1990, 2000), labels = DECLAB))]

#------------------------------------------------------------------------------
# 2. Models
#------------------------------------------------------------------------------
fitfile <- file.path(scratch, "x-mort3-lel-fits.rds")
if (file.exists(fitfile) && Sys.getenv("LEL_REFIT") != "1") {
  fits <- readRDS(fitfile); m1_split <- fits$m1_split; m2 <- fits$m2; m3 <- fits$m3
  note("Model fits loaded from cache %s (set LEL_REFIT=1 to refit)", fitfile)
} else {
  t0 <- Sys.time()
  m1_split <- stpm2(Surv(age_start, age_end, event) ~ fdiag + fsex, data = s, bhazard = s$newrate, df = 5)
  note("Model 1 (diag+sex) on split records: %.0f s, loglik %.2f", as.numeric(Sys.time() - t0, units = "secs"), -m1_split@min)
  # likelihood-equivalence check on one record per person (rstpm2's default starting values are
  # unstable for the wide person-level intervals, so start from the split-data solution)
  t0 <- Sys.time()
  m1 <- tryCatch(stpm2(Surv(age_start, age_end, event) ~ fdiag + fsex, data = p1, bhazard = p1$newrate, df = 5,
                       init = coef(m1_split)), error = function(e) NULL)
  if (!is.null(m1)) {
    note("Model 1 (diag+sex) on person records (init from split fit): %.0f s, loglik %.2f (should equal split loglik)",
         as.numeric(Sys.time() - t0, units = "secs"), -m1@min)
    note("Max abs coef difference split vs person-level fit: %.2e", max(abs(coef(m1_split) - coef(m1))))
  } else note("Person-level fit failed; split-record fits used throughout")
  t0 <- Sys.time()
  m2 <- stpm2(Surv(age_start, age_end, event) ~ fdiag + fsex + fdec, data = s, bhazard = s$newrate, df = 5)
  note("Model 2 (diag+sex+decade): %.0f s, AIC %.1f (model 1 AIC %.1f)", as.numeric(Sys.time() - t0, units = "secs"), AIC(m2), AIC(m1_split))
  t0 <- Sys.time()
  m3 <- tryCatch(
    stpm2(Surv(age_start, age_end, event) ~ fdiag + fsex + fdec, data = s, bhazard = s$newrate, df = 5,
          tvc = list(fdiag = 2)),
    error = function(e) { note("Model 3 (tvc diag) failed: %s", conditionMessage(e)); NULL })
  if (!is.null(m3)) note("Model 3 (diag+sex+decade, tvc diag df=2): %.0f s, AIC %.1f", as.numeric(Sys.time() - t0, units = "secs"), AIC(m3))
  saveRDS(list(m1_split = m1_split, m2 = m2, m3 = m3), fitfile)
}
conv_ok <- FALSE
if (!is.null(m3)) {
  se3 <- sqrt(diag(vcov(m3)))
  conv_ok <- all(is.finite(se3)) && max(se3) < 10
  if (!conv_ok) note("Model 3 vcov not well behaved (max SE %.1f); using proportional model 2", max(se3))
}
use_m3 <- !is.null(m3) && conv_ok && AIC(m3) < AIC(m2) - 2
fit <- if (use_m3) m3 else m2
note("Final model: %s", if (use_m3) "model 3 (time-dependent diag effect)" else "model 2 (proportional)")
cf <- coef(fit); V <- vcov(fit)
note("Coefficients (excess hazard ratios):")
hr <- exp(cf[!grepl("nsx", names(cf))])
for (nm in names(hr)) note("  %-45s %.2f", nm, hr[nm])

#------------------------------------------------------------------------------
# 3. Expected (general population) mortality: 2022 E&W, single ages 0-100
#------------------------------------------------------------------------------
pop22 <- rt[yeargrp == 2022]
pop_single <- function(sx) {
  r <- pop22[sex == sx][order(ageband)]
  mid <- ifelse(r$ageband == 85, 90, r$ageband + 2.5)
  lr  <- log(r$newrate)
  ages <- seq(0, AMAX, STEP)
  # log-linear interpolation of midpoints; extrapolate with last slope beyond 90; flat below 2.5
  slope_end <- (lr[length(lr)] - lr[length(lr) - 1]) / (mid[length(mid)] - mid[length(mid) - 1])
  lh <- approx(mid, lr, xout = ages, rule = 2)$y
  lh[ages > max(mid)] <- lr[length(lr)] + slope_end * (ages[ages > max(mid)] - max(mid))
  data.table(sex = sx, age = ages, hexp = exp(lh))
}
pop <- rbind(pop_single(1), pop_single(2))
setkey(pop, sex, age)
# cumulative expected hazard on the grid (trapezoid)
pop[, Hexp := c(0, cumsum((hexp[-1] + hexp[-.N]) / 2 * diff(age))), by = sex]

grid <- pop$age[pop$sex == 1]
S_exp_cond <- function(sx, a0) {      # S_exp(a | a0) on grid (0 for a < a0)
  H <- pop[sex == sx]$Hexp
  H0 <- approx(grid, H, a0)$y
  out <- exp(-(H - H0)); out[grid < a0] <- NA; out
}
integ <- function(y, x = grid) { ok <- !is.na(y); x <- x[ok]; y <- y[ok]; sum((y[-1] + y[-length(y)]) / 2 * diff(x)) }
exp_le <- CJ(sex = 1:2, a0 = A0S)[, .(exp_le = integ(S_exp_cond(sex, a0))), by = .(sex, a0)]
note("Expected residual life expectancy (2022 E&W period table):")
for (i in seq_len(nrow(exp_le))) note("  sex %d age %d: %.1f y", exp_le$sex[i], exp_le$a0[i], exp_le$exp_le[i])

#------------------------------------------------------------------------------
# 4. Model relative survival on the age grid for each covariate profile; LE by
#    parametric simulation
#------------------------------------------------------------------------------
profiles <- CJ(diag = 1:11, sex = 1:2, decdxcat = c(1900, 1970, 1980, 1990, 2000))
profiles[, `:=`(fdiag = factor(diag, levels = 1:11, labels = DIAGLAB),
                fsex  = factor(sex, levels = 1:2, labels = c("male", "female")),
                fdec  = factor(decdxcat, levels = c(1900, 1970, 1980, 1990, 2000), labels = DECLAB))]
profiles[, prof := .I]
gpos <- grid[grid >= 5]               # model defined for age >= 5 (log age spline)
B <- t(mvrnorm(NSIM, cf, V))          # p x NSIM
B <- cbind(cf, B)                     # column 1 = point estimate

# design matrix per profile (ngrid x p), eta = X %*% B, R = exp(-exp(eta))
lel_profile <- function(pr) {
  nd <- data.frame(age_start = 5, age_end = gpos, event = 0, newrate = 0,
                   fdiag = pr$fdiag, fsex = pr$fsex, fdec = pr$fdec)
  X <- predict(fit, newdata = nd, type = "lpmatrix")
  eta <- X %*% B                       # ngrid x (NSIM+1)
  cumexc <- exp(eta)                   # cumulative excess hazard Lambda_e(a)
  res <- list()
  for (a0 in A0S) {
    i0 <- which.min(abs(gpos - a0))
    Rc <- exp(-(sweep(cumexc, 2, cumexc[i0, ], "-")))   # R(a)/R(a0), matrix
    Rc[gpos < a0, ] <- NA
    Se <- S_exp_cond(pr$sex, a0)[grid >= 5]
    # survivor LE for each column
    y <- Rc * Se
    ok <- gpos >= a0
    w <- diff(gpos[ok]); yy <- y[ok, , drop = FALSE]
    surv_le <- colSums((yy[-1, , drop = FALSE] + yy[-nrow(yy), , drop = FALSE]) / 2 * w)
    # all-cause survival to 60 given alive at a0 (sanity)
    i60 <- which.min(abs(gpos - 60))
    res[[length(res) + 1]] <- data.table(prof = pr$prof, a0 = a0, draw = 0:NSIM,
                                         surv_le = surv_le, R60 = Rc[i60, ] * Se[i60])
  }
  rbindlist(res)
}
t0 <- Sys.time()
sim <- rbindlist(lapply(seq_len(nrow(profiles)), function(i) lel_profile(profiles[i])))
note("Simulation over %d profiles x %d a0 x %d draws: %.0f s", nrow(profiles), length(A0S), NSIM, as.numeric(Sys.time() - t0, units = "secs"))
sim <- merge(sim, profiles[, .(prof, diag, sex, decdxcat)], by = "prof")
sim <- merge(sim, exp_le, by = c("sex", "a0"))

#------------------------------------------------------------------------------
# 5. Standardisation weights: distribution of (diag, sex, decade) among survivors
#    under follow-up at age a0, within each reporting group
#------------------------------------------------------------------------------
# a0-specific weight tables, grouped by: overall, diag, decade, sex, diag x sex
run_group <- function(a0v, expr_str) {
  at <- d[age0 <= a0v & age1 > a0v]
  gexpr <- parse(text = expr_str)[[1]]
  at[, grp := eval(gexpr)]
  w <- at[, .N, by = .(diag, sex, decdxcat, grp)]
  # groups with nobody under follow-up at a0: fall back to whole-group distribution
  allg <- copy(d)[, grp := eval(gexpr)]
  missing <- setdiff(unique(allg$grp), unique(w$grp))
  if (length(missing)) {
    w2 <- allg[grp %in% missing, .N, by = .(diag, sex, decdxcat, grp)]
    w <- rbind(w, w2)
    note("  a0=%d, grouping %s: groups %s have no survivors under follow-up at a0; whole-group mix used (extrapolation)",
         a0v, expr_str, paste(missing, collapse = ","))
  }
  w[, wt := N / sum(N), by = grp]
  natrisk <- at[, .(n_atrisk = .N), by = grp]
  x <- merge(sim[a0 == a0v], w, by = c("diag", "sex", "decdxcat"), allow.cartesian = TRUE)
  z <- x[, .(surv_le = sum(wt * surv_le), exp_le = sum(wt * exp_le), cum60 = 1 - sum(wt * R60)), by = .(grp, draw)]
  z[, `:=`(lel = exp_le - surv_le, plost = 100 * (exp_le - surv_le) / exp_le)]
  out <- z[, .(a0 = a0v,
               exp_le = exp_le[draw == 0], surv_le = surv_le[draw == 0],
               surv_le_lo = quantile(surv_le[draw > 0], .025), surv_le_hi = quantile(surv_le[draw > 0], .975),
               lel = lel[draw == 0], lel_lo = quantile(lel[draw > 0], .025), lel_hi = quantile(lel[draw > 0], .975),
               plost = plost[draw == 0], plost_lo = quantile(plost[draw > 0], .025), plost_hi = quantile(plost[draw > 0], .975),
               cum60 = cum60[draw == 0]), by = grp]
  merge(out, natrisk, by = "grp", all.x = TRUE)
}
res_overall <- rbindlist(lapply(A0S, run_group, expr_str = "0L"))[, `:=`(group = "overall", level = "overall")]
res_diag    <- rbindlist(lapply(A0S, run_group, expr_str = "diag"))[, `:=`(group = "diag", level = DIAGLAB[grp])]
res_dec     <- rbindlist(lapply(A0S, run_group, expr_str = "decdxcat"))[, `:=`(group = "decade", level = DECLAB[as.character(grp)])]
res_sex     <- rbindlist(lapply(A0S, run_group, expr_str = "sex"))[, `:=`(group = "sex", level = c("male", "female")[grp])]
res_diagsex <- rbindlist(lapply(A0S, run_group, expr_str = "diag*10+sex"))[, `:=`(group = "diag x sex",
                 level = paste(DIAGLAB[grp %/% 10], c("male", "female")[grp %% 10]))]
res <- rbind(res_overall, res_diag, res_dec, res_sex, res_diagsex)
res[, grp := NULL]
setcolorder(res, c("group", "level", "a0", "n_atrisk"))
num <- setdiff(names(res), c("group", "level", "a0", "n_atrisk"))
res[, (num) := lapply(.SD, round, 2), .SDcols = num]

fwrite(res[group %in% c("overall", "diag")], file.path(output, "table_lel.csv"))
fwrite(res[group == "decade" & a0 %in% c(20, 40)], file.path(output, "table_lel_bydecade.csv"))
fwrite(res[group %in% c("sex", "diag x sex")], file.path(output, "table_lel_bysex.csv"))
fwrite(res, file.path(output, "table_lel_all.csv"))

#------------------------------------------------------------------------------
# 6. Sanity checks: observed vs model cumulative mortality, implied SMR
#------------------------------------------------------------------------------
library(survival)
km <- survfit(Surv(age0, age1, event) ~ 1, data = d[age1 > 20][, age0 := pmax(age0, 20)])
km60 <- 1 - summary(km, times = 60)$surv
note("Observed (KM, age scale, conditional on under follow-up at 20) cumulative all-cause mortality by age 60: %.1f%%", 100 * km60)
note("Model-based (standardised to mix at 20) cumulative all-cause mortality by age 60 from age 20: %.1f%%", 100 * res_overall[a0 == 20]$cum60)
# like-for-like: survivors diagnosed before 1980 (the only ones contributing much follow-up beyond 50)
km80 <- survfit(Surv(age0, age1, event) ~ 1, data = d[age1 > 20 & decdxcat < 1980][, age0 := pmax(age0, 20)])
note("Observed KM cumulative mortality by 60 from age 20, diagnosed <1980: %.1f%%; model (decades <1970 / 1970-79): %.1f%% / %.1f%%",
     100 * (1 - summary(km80, times = 60)$surv), 100 * res_dec[a0 == 20 & grp == 1900]$cum60, 100 * res_dec[a0 == 20 & grp == 1970]$cum60)
km5 <- survfit(Surv(age0, age1, event) ~ 1, data = d)
note("Observed KM cumulative all-cause mortality by age 40/50/60 from entry: %.1f%% / %.1f%% / %.1f%%",
     100 * (1 - summary(km5, times = 40)$surv), 100 * (1 - summary(km5, times = 50)$surv), 100 * (1 - summary(km5, times = 60)$surv))
# model-implied SMR: sum over split records of (newrate + excess hazard at midpoint) * py vs expected
s[, mid := (age_start + age_end) / 2]
nd <- data.frame(age_start = 5, age_end = s$mid, event = 0, newrate = 0, fdiag = s$fdiag, fsex = s$fsex, fdec = s$fdec)
s[, hexc := as.numeric(predict(fit, newdata = nd, type = "hazard"))]   # excess hazard (bhazard=0 in newdata)
Opred <- sum((s$newrate + s$hexc) * s$py)
note("Model-implied deaths = %.0f vs observed %d; implied SMR = %.2f vs observed %.2f", Opred, O, Opred / E, O / E)
smr_age <- s[, .(agegrp = cut(mid, c(0, 20, 30, 40, 50, 60, 110), right = FALSE))][, `:=`(O = s$event, E = s$e, Opred = (s$newrate + s$hexc) * s$py)][
  , .(obs_smr = sum(O) / sum(E), model_smr = sum(Opred) / sum(E), O = sum(O)), by = agegrp][order(agegrp)]
note("SMR by attained age, observed vs model-implied:")
for (i in seq_len(nrow(smr_age))) note("  %-9s O=%5d  obs SMR %.2f  model SMR %.2f", as.character(smr_age$agegrp[i]), smr_age$O[i], smr_age$obs_smr[i], smr_age$model_smr[i])
note("Oldest exit age in cohort: %.1f; person-years above age 70: %.0f (%.1f%% of total); above 80: %.0f",
     max(d$age1), s[mid >= 70, sum(py)], 100 * s[mid >= 70, sum(py)] / PY, s[mid >= 80, sum(py)])

#------------------------------------------------------------------------------
# 7. Figures
#------------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()))
f1 <- res_diag[a0 == 20][order(lel)]
f1[, level := factor(level, levels = level)]
ov <- res_overall[a0 == 20]
g1 <- ggplot(f1, aes(x = lel, y = level)) +
  geom_vline(xintercept = ov$lel, colour = "grey55", linetype = 2) +
  geom_errorbarh(aes(xmin = lel_lo, xmax = lel_hi), height = 0, colour = "#1f5fa0") +
  geom_point(size = 2.5, colour = "#1f5fa0") +
  annotate("text", x = ov$lel, y = 0.6, label = sprintf("overall %.1f y", ov$lel), hjust = -0.05, size = 3, colour = "grey30") +
  labs(x = "Loss in expectation of life from age 20 (years, 95% CI)", y = NULL,
       title = "Loss in expectation of life by childhood cancer type",
       subtitle = "5-year survivors alive at age 20, standardised to each type's sex and decade-of-diagnosis mix") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  theme(panel.grid.major.y = element_line(colour = "grey92"), plot.title.position = "plot",
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9, colour = "grey30"))
ggsave(file.path(output, "fig_lel_bydiag.png"), g1, width = 7, height = 4.5, dpi = 300, bg = "white")

f2 <- res_dec[a0 == 20]
f2[, level := factor(level, levels = unname(DECLAB))]
f2l <- rbind(f2[, .(level, type = "survivors", le = surv_le, lo = surv_le_lo, hi = surv_le_hi)],
             f2[, .(level, type = "general population", le = exp_le, lo = NA_real_, hi = NA_real_)])
g2 <- ggplot(f2l, aes(x = level, y = le, colour = type, group = type)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, na.rm = TRUE) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_colour_manual(values = c("general population" = "#c8562b", "survivors" = "#1f5fa0"), name = NULL) +
  scale_y_continuous(limits = c(40, 70)) +
  labs(x = "Decade of diagnosis", y = "Residual life expectancy at age 20 (years)",
       title = "Residual life expectancy at age 20 by decade of diagnosis",
       subtitle = "Survivors (model-based, 95% CI) vs England & Wales 2022 population; standardised to each decade's cancer type/sex mix") +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey92"), plot.title.position = "plot",
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 8.5, colour = "grey30"))
ggsave(file.path(output, "fig_le_bydecade.png"), g2, width = 7, height = 4.5, dpi = 300, bg = "white")

note("Done.")
print(res[group %in% c("overall", "decade")])
print(res[group == "diag" & a0 == 20])
