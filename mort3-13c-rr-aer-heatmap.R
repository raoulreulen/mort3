# =============================================================================
# MORT3-13c: Lexis-style heat maps of AER and SMR by attained age x decade of
#            diagnosis (BCCSS late mortality, v3)
#
# INPUT : $output/aer_heatmap_data.csv (from mort3-13b-rr-aer-heatmap.do)
# OUTPUT: $output/fig_aer_heatmap.png  (AER per 100,000 PY; recurrence = crude rate)
#         $output/fig_smr_heatmap.png  (SMR; all causes except recurrence)
# Colour: viridis (colour-blind safe), AER on a pseudo-log scale
#         (scales::pseudo_log_trans, sigma = 20) because the AER spans ~0 to
#         several thousand and can be negative; SMR on a log10 scale.
# Cells with <1000 person-years are blank (grey, no text).
# Date: 23 Oct 2026   Author: RCR (Claude-assisted)
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(scales) })
outdir <- "/Users/raoulreulen/Library/CloudStorage/OneDrive-UniversityofBirmingham/mort3/output"
d <- read.csv(file.path(outdir, "aer_heatmap_data.csv"), stringsAsFactors = FALSE)

cause_lv <- c("all causes", "recurrence/progression", "SPN", "circulatory", "respiratory", "external")
cause_lb <- c("All causes", "Recurrence/progression (crude rate)", "SPN", "Circulatory", "Respiratory", "External")
d <- d %>%
  mutate(cause  = factor(cause, levels = cause_lv, labels = cause_lb),
         agegrp = factor(agegrp, levels = c("5-19", "20-29", "30-39", "40-49", "50-59", "60+")),
         decade = factor(decade, levels = c("<1970", "1970-79", "1980-89", "1990-99", "2000-06")),
         blank  = blank == 1,
         lab_aer = ifelse(blank, "", paste0(formatC(round(aer), format = "d", big.mark = ","), "\n(", obs, ")")),
         lab_smr = ifelse(blank | is.na(smr), "", paste0(formatC(smr, format = "f", digits = 1), "\n(", obs, ")")))

use_viridis <- requireNamespace("viridisLite", quietly = TRUE)
fill_scale <- function(name, trans, breaks, labels) {
  if (use_viridis) scale_fill_viridis_c(name = name, option = "viridis", trans = trans, breaks = breaks,
                                        labels = labels, na.value = "grey90", direction = -1)
  else scale_fill_gradient(name = name, low = "#deebf7", high = "#08306b", trans = trans, breaks = breaks,
                           labels = labels, na.value = "grey90")
}
base_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 8), legend.position = "right",
        plot.caption = element_text(size = 7, hjust = 0), plot.background = element_rect(fill = "white", colour = NA))

# --- AER heat map -------------------------------------------------------------
aer_breaks <- c(10, 50, 200, 1000)
p1 <- ggplot(d, aes(x = agegrp, y = decade, fill = aer)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab_aer, colour = ifelse(!is.na(aer) & aer > 400, "w", "b")), size = 2.3, lineheight = 0.85) +
  scale_colour_manual(values = c(w = "white", b = "black"), guide = "none") +
  fill_scale("AER per\n100,000 PY", trans = pseudo_log_trans(sigma = 20, base = 10), breaks = aer_breaks,
             labels = comma(aer_breaks)) +
  facet_wrap(~cause, ncol = 3) +
  labs(x = "Attained age (years)", y = "Decade of diagnosis",
       caption = paste0("Absolute excess risk (AER = observed - expected deaths per 100,000 person-years; recurrence/progression: crude rate, no population comparator).\n",
                        "Cell text: AER (observed deaths). Colour on a pseudo-log scale. Grey cells: <1,000 person-years. BCCSS five-year survivors, follow-up to 31 Dec 2022.")) +
  base_theme
ggsave(file.path(outdir, "fig_aer_heatmap.png"), p1, width = 10, height = 7, dpi = 300)

# --- SMR heat map (all except recurrence) --------------------------------------
ds <- d %>% filter(cause != "Recurrence/progression (crude rate)") %>% mutate(cause = droplevels(cause),
        smr_fill = ifelse(!is.na(smr) & smr > 0, smr, NA_real_))
smr_breaks <- c(0.5, 1, 2, 5, 10, 20, 50)
p2 <- ggplot(ds, aes(x = agegrp, y = decade, fill = smr_fill)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab_smr, colour = ifelse(!is.na(smr_fill) & smr_fill > 6, "w", "b")), size = 2.3, lineheight = 0.85) +
  scale_colour_manual(values = c(w = "white", b = "black"), guide = "none") +
  fill_scale("SMR", trans = "log10", breaks = smr_breaks, labels = smr_breaks) +
  facet_wrap(~cause, ncol = 3) +
  labs(x = "Attained age (years)", y = "Decade of diagnosis",
       caption = paste0("Standardised mortality ratio (observed/expected deaths). Cell text: SMR (observed deaths). Colour on a log10 scale.\n",
                        "Grey cells: <1,000 person-years or no observed deaths. BCCSS five-year survivors, follow-up to 31 Dec 2022.")) +
  base_theme
ggsave(file.path(outdir, "fig_smr_heatmap.png"), p2, width = 10, height = 5, dpi = 300)
cat("done\n")
