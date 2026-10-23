#===============================================================================
# mort3-12a-rr-lyl-bycause.R
#
# PURPOSE
#   Life years lost (LYL) by 5-year survivors of childhood cancer (BCCSS),
#   decomposed by cause of death using Andersen's decomposition as implemented
#   in R package 'lillies' (Plana-Ripoll et al. 2020, PLoS ONE 15:e0228073),
#   and compared with the general population (England & Wales, 2022 period
#   life table built from $rates/allcauses.dta).
#
#   Quantities (all "restricted" to a maximum age tau, from a given age a):
#     (i)   LE_surv   = survivors' remaining life expectancy from age a to tau
#                       (area under all-cause KM, delayed entry handled by lillies)
#           LYL_surv  = tau - a - LE_surv = survivors' total LYL before tau,
#                       split by cause (lillies: area under each cause-specific
#                       cumulative incidence function); causes sum to LYL_surv.
#     (ii)  LE_pop    = general-population remaining LE from a to tau
#                       (2022 E&W period rates; sex-weighted by the sex
#                       distribution of survivors alive at age a)
#           LYL_pop   = tau - a - LE_pop  (population LYL before tau, all causes)
#     (iii) Excess LYL = LE_pop - LE_surv = LYL_surv - LYL_pop
#   95% CIs: non-parametric bootstrap (lillies::lyl_ci, NITER resamples of
#   survivors; percentile CI). The population life table is treated as fixed,
#   so CI(excess) = CI(LYL_surv) shifted by LYL_pop (as in lillies::lyl_diff_ref).
#
#   tau: 80 years (main). Survivors at risk beyond age 70 = 583, beyond 80 = 28,
#   beyond 85 = 2, so tau = 90 would extrapolate the survivors' curve from
#   almost no data (lillies carries the last CIF value flat to tau).
#   Sensitivity at tau = 90 for the overall cohort only.
#   Decade of diagnosis: ages 20 and 30 only, with tau = 60 and tau = 80
#   (later decades cannot have reached older ages; 2000-06 reach at most ~37).
#
# INPUTS
#   $temp/x-mort3-prepforstset.dta     one row per 5-year survivor
#   $rates/allcauses.dta               E&W all-cause rates by sex x 5-y ageband x year
#
# OUTPUTS ($output)
#   table_lyl_bycause.csv        one row per stratum x age x tau (wide; numbers + formatted strings)
#   table_lyl_bycause_long.csv   long format (one row per stratum x age x tau x cause)
#   fig_lyl_bycause.png          stacked bars of LYL by cause at age 20
#   fig_lyl_bycause.csv          numbers behind the figure
#   Log of run: mort3-12a-rr-lyl-bycause.log (in $output)
#
# Date: 23 Oct 2026   Author: RCR (Claude-assisted)
#===============================================================================

suppressMessages({
  library(haven); library(dplyr); library(tidyr); library(ggplot2)
  library(lillies); library(parallel)
})

t_start <- Sys.time()

temp   <- "/Users/raoulreulen/Documents/TEMPWORK"
od     <- "/Users/raoulreulen/Library/CloudStorage/OneDrive-UniversityofBirmingham/mort3"
rates  <- file.path(od, "rates")
output <- file.path(od, "output")

NITER    <- as.integer(Sys.getenv("LYL_NITER", "500"))   # bootstrap iterations (env override for testing)
TAU_MAIN <- 80
TAU_SENS <- 90             # overall cohort only
TAU_DEC  <- 60             # decade of diagnosis
AGES     <- c(20, 30, 40, 50)
NCORES   <- max(1, min(8, detectCores() - 2))

#-------------------------------------------------------------------------------
# 1. Survivor data
#-------------------------------------------------------------------------------
d <- read_dta(file.path(temp, "x-mort3-prepforstset.dta"))
diag_lab <- c("1"="leukaemia", "2"="Hodgkin lymphoma", "3"="non-Hodgkin lymphoma",
              "4"="CNS", "5"="neuroblastoma", "6"="heritable retinoblastoma",
              "7"="non-heritable retinoblastoma", "8"="Wilms tumour", "9"="bone",
              "10"="soft tissue sarcoma", "11"="other")
dec_lab  <- c("1900"="<1970", "1970"="1970-79", "1980"="1980-89", "1990"="1990-99", "2000"="2000-06")
cause_levels <- c("Alive", "recurrence/progression", "SPN", "circulatory",
                  "respiratory", "external", "infection", "other causes")

d <- d %>%
  mutate(ageent  = as.numeric(doe - dob) / 365.25,
         ageexit = as.numeric(dox - dob) / 365.25,
         dead    = !is.na(allcauses) & allcauses == 1,
         cause   = case_when(!dead ~ "Alive",
                             recur == 1       ~ "recurrence/progression",
                             spn == 1         ~ "SPN",
                             circulation == 1 ~ "circulatory",
                             respiratory == 1 ~ "respiratory",
                             external == 1    ~ "external",
                             infection == 1   ~ "infection",
                             TRUE             ~ "other causes"),   # incl. unknown cause
         cause   = factor(cause, levels = cause_levels),
         sexlab  = factor(sex, levels = 1:2, labels = c("male", "female")),
         diaglab = factor(diag_lab[as.character(diag)], levels = diag_lab),
         declab  = factor(dec_lab[as.character(decdxcat)], levels = dec_lab)) %>%
  select(ageent, ageexit, cause, sex, sexlab, diaglab, declab) %>%
  as.data.frame()

stopifnot(all(d$ageexit > d$ageent), !any(is.na(d$cause)))
cat("Survivors:", nrow(d), " deaths:", sum(d$cause != "Alive"), "\n")
print(table(d$cause))
cat("At risk beyond age: ")
for (a in c(20, 30, 40, 50, 60, 70, 80, 85, 90)) cat(a, "=", sum(d$ageexit > a), "; ")
cat("\n")

#-------------------------------------------------------------------------------
# 2. General-population reference: 2022 E&W period life table, single years
#    S(a+1) = S(a) * exp(-r_a); LE(a..tau) = sum_{x=a}^{tau-1} S(x)/S(a) * (1-exp(-r_x))/r_x
#    (exact integral of a piecewise-constant hazard; equals integral of S to tau)
#-------------------------------------------------------------------------------
r22 <- read_dta(file.path(rates, "allcauses.dta")) %>%
  filter(yeargrp == 2022) %>% select(sex, ageband, newrate)
stopifnot(nrow(r22) == 36)
lt <- expand.grid(sex = 1:2, age = 0:99) %>%
  mutate(ageband = pmin(5 * floor(age / 5), 85)) %>%
  left_join(r22, by = c("sex", "ageband")) %>%
  arrange(sex, age) %>%
  group_by(sex) %>%
  mutate(S = cumprod(c(1, exp(-newrate)))[-(n() + 1)],      # S at start of year 'age'
         Ly = S * (1 - exp(-newrate)) / newrate) %>%        # person-years lived in [age, age+1)
  ungroup()
stopifnot(!any(is.na(lt$newrate)))

pop_rle <- function(a, tau, sex1) {   # restricted LE from a to tau for one sex
  x <- lt[lt$sex == sex1 & lt$age >= a & lt$age < tau, ]
  sum(x$Ly) / x$S[1]
}
pop_rle_w <- function(a, tau, pmale) pmale * pop_rle(a, tau, 1) + (1 - pmale) * pop_rle(a, tau, 2)

cat("\nPopulation restricted LE (2022 E&W) from age 20, tau 80: male",
    round(pop_rle(20, 80, 1), 2), " female", round(pop_rle(20, 80, 2), 2), "\n")

#-------------------------------------------------------------------------------
# 3. Task list (stratum x age x tau)
#-------------------------------------------------------------------------------
strata <- bind_rows(
  tibble(type = "overall", stratum = "all survivors", var = NA_character_, value = NA_character_),
  tibble(type = "sex", stratum = levels(d$sexlab), var = "sexlab", value = levels(d$sexlab)),
  tibble(type = "cancer type", stratum = levels(d$diaglab), var = "diaglab", value = levels(d$diaglab)),
  tibble(type = "decade of diagnosis", stratum = levels(d$declab), var = "declab", value = levels(d$declab))
)
tasks <- bind_rows(
  strata %>% filter(type != "decade of diagnosis") %>% tidyr::crossing(age = AGES, tau = TAU_MAIN),
  strata %>% filter(type == "overall") %>% tidyr::crossing(age = AGES, tau = TAU_SENS),
  strata %>% filter(type == "decade of diagnosis") %>% tidyr::crossing(age = c(20, 30), tau = c(TAU_DEC, TAU_MAIN))
) %>% mutate(task = row_number())
cat("Number of lyl/lyl_ci runs:", nrow(tasks), "\n")

#-------------------------------------------------------------------------------
# 4. One run: lillies point estimate + bootstrap, plus reference and excess
#-------------------------------------------------------------------------------
run_task <- function(i) {
  tk <- tasks[i, ]
  dd <- if (is.na(tk$var)) d else d[d[[tk$var]] == tk$value, ]
  atrisk <- dd[dd$ageexit > tk$age, ]
  n_atrisk <- nrow(atrisk)
  pmale <- mean(atrisk$sex == 1)
  deaths_after <- table(atrisk$cause[atrisk$cause != "Alive"])
  if (n_atrisk < 20) return(NULL)

  est <- suppressMessages(lyl(dd, t0 = ageent, t = ageexit, status = cause,
                              age_specific = tk$age, tau = tk$tau))
  boot <- suppressMessages(lyl_ci(est, niter = NITER))
  causes <- est$death_labels
  pt <- est$LYL
  bs <- boot$LYL_ci
  lo <- function(v) unname(quantile(v, 0.025, na.rm = TRUE))
  hi <- function(v) unname(quantile(v, 0.975, na.rm = TRUE))

  le_pop   <- pop_rle_w(tk$age, tk$tau, pmale)
  lyl_pop  <- tk$tau - tk$age - le_pop
  lyl_tot  <- tk$tau - tk$age - pt$life_exp
  lyl_tot_bs <- tk$tau - tk$age - bs$life_exp

  out <- tibble(
    type = tk$type, stratum = tk$stratum, age = tk$age, tau = tk$tau,
    n_atrisk = n_atrisk, deaths_after_age = sum(deaths_after),
    max_age_exit = max(dd$ageexit), extrapolated = as.integer(max(dd$ageexit) < tk$tau), prop_male = pmale,
    le_surv = pt$life_exp, le_surv_lo = lo(bs$life_exp), le_surv_hi = hi(bs$life_exp),
    le_pop = le_pop,
    lyl_surv_total = lyl_tot, lyl_surv_total_lo = lo(lyl_tot_bs), lyl_surv_total_hi = hi(lyl_tot_bs),
    lyl_pop = lyl_pop,
    excess_lyl = lyl_tot - lyl_pop, excess_lyl_lo = lo(lyl_tot_bs) - lyl_pop, excess_lyl_hi = hi(lyl_tot_bs) - lyl_pop
  )
  for (cz in causes) {
    nm <- gsub("[^a-z]", "", tolower(cz)); nm <- sub("recurrenceprogression", "recur", nm)
    nm <- sub("othercauses", "other", nm); nm <- sub("circulatory", "circ", nm); nm <- sub("respiratory", "resp", nm)
    col <- gsub(" ", "", cz)                      # lillies strips spaces from column names
    if (is.null(pt[[col]])) { pt[[col]] <- 0; bs[[col]] <- 0 }   # cause absent in this stratum
    out[[paste0("lyl_", nm)]]       <- pt[[col]]
    out[[paste0("lyl_", nm, "_lo")]] <- lo(bs[[col]])
    out[[paste0("lyl_", nm, "_hi")]] <- hi(bs[[col]])
    out[[paste0("deaths_", nm)]]    <- as.integer(deaths_after[cz])
  }
  out
}

#-------------------------------------------------------------------------------
# 5. Run (parallel over tasks; reproducible seed)
#-------------------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG"); set.seed(20261023)
cat("Running", nrow(tasks), "tasks on", NCORES, "cores,", NITER, "bootstrap iterations each ...\n")
res_list <- mclapply(seq_len(nrow(tasks)), function(i) tryCatch(run_task(i),
                     error = function(e) { message("task ", i, " failed: ", conditionMessage(e)); NULL }),
                     mc.cores = NCORES, mc.set.seed = TRUE)
failed <- which(sapply(res_list, is.null))
if (length(failed)) { cat("Tasks skipped/failed (n<20 or error):\n"); print(tasks[failed, ]) }
res <- bind_rows(res_list)

#-------------------------------------------------------------------------------
# 6. Formatted strings and CSV outputs
#-------------------------------------------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
ci1 <- function(e, l, h) paste0(f1(e), " (", f1(l), " to ", f1(h), ")")
ci2 <- function(e, l, h) paste0(f2(e), " (", f2(l), " to ", f2(h), ")")
cause_vars <- c("recur", "spn", "circ", "resp", "external", "infection", "other")
cause_names <- cause_levels[-1]

res <- res %>%
  mutate(str_le_surv = ci1(le_surv, le_surv_lo, le_surv_hi),
         str_le_pop = f1(le_pop),
         str_lyl_surv_total = ci1(lyl_surv_total, lyl_surv_total_lo, lyl_surv_total_hi),
         str_lyl_pop = f1(lyl_pop),
         str_excess = ci1(excess_lyl, excess_lyl_lo, excess_lyl_hi))
for (v in cause_vars) res[[paste0("str_", v)]] <- ci2(res[[paste0("lyl_", v)]], res[[paste0("lyl_", v, "_lo")]], res[[paste0("lyl_", v, "_hi")]])

res <- res %>% arrange(tau != TAU_MAIN, match(type, c("overall", "sex", "cancer type", "decade of diagnosis")),
                       match(stratum, c("all survivors", levels(d$sexlab), levels(d$diaglab), levels(d$declab))), tau, age)

readme <- c(
  "# mort3 LYL by cause (lillies). Columns:",
  "#  n_atrisk = survivors alive and under follow-up beyond 'age'; deaths_after_age = deaths among them",
  "#  extrapolated = 1 if no survivor in the stratum was followed to tau: lillies carries the last cumulative incidence flat to tau, which UNDERSTATES LYL (relevant for decade of diagnosis >=1980)",
  "#  le_surv = survivors' remaining LE from 'age' restricted to tau (lillies, 95% bootstrap CI)",
  "#  le_pop = general-population (E&W 2022 period rates) remaining LE from 'age' restricted to tau, sex-weighted (prop_male)",
  "#  lyl_surv_total = tau - age - le_surv = survivors' LYL before tau, split into lyl_<cause> (Andersen decomposition; causes sum to total)",
  "#  lyl_pop = tau - age - le_pop = population LYL before tau (all causes)",
  "#  excess_lyl = le_pop - le_surv = lyl_surv_total - lyl_pop (CI = bootstrap CI of lyl_surv_total shifted; life table treated as fixed)",
  paste0("#  bootstrap iterations = ", NITER, "; 'other causes' includes unknown cause of death")
)
csv_wide <- file.path(output, "table_lyl_bycause.csv")
writeLines(readme, csv_wide)
suppressWarnings(write.table(res, csv_wide, sep = ",", row.names = FALSE, append = TRUE, qmethod = "double"))

res_long <- res %>%
  select(type, stratum, age, tau, n_atrisk, starts_with("lyl_"), -lyl_pop) %>%
  pivot_longer(cols = matches("^lyl_(recur|spn|circ|resp|external|infection|other)(_lo|_hi)?$"),
               names_to = c("cause", "stat"), names_pattern = "^lyl_([a-z]+)(_lo|_hi)?$") %>%
  mutate(stat = ifelse(stat == "", "est", sub("_", "", stat))) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(cause = factor(cause_names[match(cause, cause_vars)], levels = cause_names))
write.csv(res_long, file.path(output, "table_lyl_bycause_long.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# 7. Figure: stacked bars of LYL by cause at age 20
#    Panel A: overall and by cancer type (tau = 80); Panel B: by decade (tau = 60)
#-------------------------------------------------------------------------------
figdat <- bind_rows(
  res_long %>% filter(age == 20, tau == TAU_MAIN, type %in% c("overall", "cancer type")) %>%
    mutate(panel = paste0("A. Overall and by cancer type (tau = ", TAU_MAIN, ")")),
  res_long %>% filter(age == 20, tau == TAU_DEC, type == "decade of diagnosis") %>%
    mutate(panel = paste0("B. By decade of diagnosis (tau = ", TAU_DEC, ")"))
) %>%
  mutate(stratum = factor(stratum, levels = c("all survivors", levels(d$diaglab), levels(d$declab))))
totals <- figdat %>% group_by(panel, stratum) %>% summarise(total = sum(est), .groups = "drop")
excess_lab <- res %>% filter(age == 20, (tau == TAU_MAIN & type %in% c("overall", "cancer type")) | (tau == TAU_DEC & type == "decade of diagnosis")) %>%
  select(stratum, tau, excess_lyl, lyl_pop) %>% mutate(stratum = factor(stratum, levels = levels(figdat$stratum)))
totals <- totals %>% left_join(excess_lab %>% select(-tau), by = "stratum")
write.csv(figdat %>% select(panel, stratum, cause, est, lo, hi) %>% left_join(totals, by = c("panel", "stratum")),
          file.path(output, "fig_lyl_bycause.csv"), row.names = FALSE)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(okabe) <- cause_names
p <- ggplot(figdat, aes(x = stratum, y = est, fill = cause)) +
  geom_col(width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_point(data = totals, aes(x = stratum, y = lyl_pop), inherit.aes = FALSE,
             shape = 23, size = 2.2, fill = "white", colour = "black") +
  geom_text(data = totals, aes(x = stratum, y = total, label = f1(total)), inherit.aes = FALSE,
            vjust = -0.4, size = 2.8) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = okabe, name = "Cause of death") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL,
       y = "Life years lost from age 20 (years)",
       caption = paste0("Bars: survivors' life years lost before age tau, decomposed by cause of death (lillies). ",
                        "\nDiamond: life years lost before tau by the general population (England & Wales 2022 rates, sex-weighted); ",
                        "bar height minus diamond = excess life years lost.\nPanel B: later decades have not reached tau (cumulative incidence carried flat to tau, so LYL are understated for 1990-99 and 2000-06).")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
        panel.grid.major.x = element_blank(), plot.caption = element_text(hjust = 0, size = 7.5),
        strip.background = element_rect(fill = "grey92"))
ggsave(file.path(output, "fig_lyl_bycause.png"), p, width = 11, height = 6.5, dpi = 300, bg = "white")

#-------------------------------------------------------------------------------
# 8. Console summary
#-------------------------------------------------------------------------------
cat("\n==== Overall (tau =", TAU_MAIN, ") ====\n")
print(as.data.frame(res %>% filter(type == "overall") %>%
      select(age, tau, n_atrisk, str_le_surv, str_le_pop, str_excess, str_recur, str_spn, str_circ, str_resp, str_external, str_infection, str_other)))
cat("\n==== Age 20, excess LYL by stratum ====\n")
print(as.data.frame(res %>% filter(age == 20) %>% select(type, stratum, tau, n_atrisk, max_age_exit, str_le_surv, str_le_pop, str_excess)))
cat("\nRuntime:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
