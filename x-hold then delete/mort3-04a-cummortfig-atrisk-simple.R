# ==============================================================================
# FIGURE 1 - CUMULATIVE MORTALITY WITH NUMBER AT RISK TABLE (SIMPLIFIED R VERSION)
# Uses base R and minimal dependencies to avoid package conflicts
# ==============================================================================

# Load only essential packages
library(survival)
library(cmprsk)
library(haven)

# Set parameters
xtime <- 70
cod <- c("spn", "recur", "circulation", "respiratory", "external")

# ==============================================================================
# READ AND PREPARE DATA
# ==============================================================================

# Read the prepared data
# Adjust these paths based on your temp directory location
temp_dir <- Sys.getenv("temp")
if (temp_dir == "") {
  temp_dir <- "C:/Users/raoul/AppData/Local/Temp"  # Default Windows temp
}

# Read stset data for all causes
data_allcauses <- read_dta(file.path(temp_dir, "x-mort3-stset-allcauses.dta"))

# ==============================================================================
# CALCULATE NUMBER AT RISK BY 10-YEAR AGE BANDS
# ==============================================================================

ages <- c(5, 15, 25, 35, 45, 55, 65)
midpoints <- ages + 5
n_at_risk <- numeric(length(ages))

for (i in seq_along(ages)) {
  n_at_risk[i] <- sum(data_allcauses$`_t0` <= ages[i] & data_allcauses$`_t` > ages[i], na.rm = TRUE)
}

# ==============================================================================
# CALCULATE CUMULATIVE INCIDENCE FOR EACH CAUSE OF DEATH
# ==============================================================================

# Read the original file
data <- read_dta(file.path(temp_dir, "x-mort3-prepforstset.dta"))

# Storage for CI results
ci_list <- list()

for (cause in cod) {
  # Create event indicator: 1 = event of interest, 2 = competing risk, 0 = censored
  event_status <- ifelse(data[[cause]] == 1, 1,
                         ifelse(data$allcauses == 1, 2, 0))

  # Calculate follow-up time in years
  time <- as.numeric(difftime(data$dox, data$dob, units = "days")) / 365.25

  # Use cuminc for competing risks analysis
  fit <- cuminc(ftime = time, fstatus = event_status, cencode = 0)

  # Extract cumulative incidence for event of interest (first element)
  ci_time <- fit[[1]]$time
  ci_est <- fit[[1]]$est * 100  # Convert to percentage

  # Filter to age range 5-70
  keep_idx <- ci_time >= 5 & ci_time <= 70
  ci_list[[cause]] <- list(time = ci_time[keep_idx], ci = ci_est[keep_idx])
}

# ==============================================================================
# CALCULATE EXPECTED MORTALITY
# ==============================================================================

expected_data <- read_dta(file.path(temp_dir, "x-mort3-expected-allcause.dta"))
expected_data$expected <- (1 - expected_data$conditional) * 100

# Filter to age range
keep_exp <- expected_data$t_exp >= 5 & expected_data$t_exp <= 70
exp_time <- expected_data$t_exp[keep_exp]
exp_values <- expected_data$expected[keep_exp]

# Apply lowess smoothing (span = 0.3 to match Stata)
exp_smooth <- lowess(exp_time, exp_values, f = 0.3)

# ==============================================================================
# CREATE PLOT WITH NUMBER AT RISK TABLE
# ==============================================================================

# Set up PDF output
pdf("mort3-04a-cummortfig-atrisk.pdf", width = 8, height = 7)

# Create layout: main plot (80%) and risk table (20%)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.85, 0.15))

# ---- MAIN PLOT ----
par(mar = c(0.5, 4, 2, 2))

# Initialize plot
plot(NULL, xlim = c(5, 70), ylim = c(0, 16),
     xlab = "", ylab = "Cumulative mortality, %",
     xaxt = "n", yaxt = "n", bty = "n")

# Add grid
axis(2, at = 0:16, las = 1, cex.axis = 0.9)
abline(h = 0:16, col = "gray90", lty = 1, lwd = 0.5)

# Add cumulative incidence lines
lines(ci_list$recur$time, ci_list$recur$ci, col = "black", lwd = 2)
lines(ci_list$spn$time, ci_list$spn$ci, col = "blue", lwd = 2)
lines(ci_list$circulation$time, ci_list$circulation$ci, col = "red", lwd = 2)
lines(ci_list$respiratory$time, ci_list$respiratory$ci, col = "green", lwd = 2)
lines(ci_list$external$time, ci_list$external$ci, col = "black", lwd = 2)

# Add expected mortality (dashed line)
lines(exp_smooth$x, exp_smooth$y, col = "cyan", lwd = 2, lty = 2)

# Add x-axis at top
axis(1, at = seq(5, 70, 5), labels = FALSE)

# Add legend
legend("topleft",
       legend = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality"),
       col = c("black", "blue", "red", "green", "black", "cyan"),
       lty = c(1, 1, 1, 1, 1, 2),
       lwd = 2,
       bty = "n",
       cex = 0.8)

# ---- RISK TABLE ----
par(mar = c(4, 4, 0, 2))

# Create empty plot for risk table
plot(NULL, xlim = c(5, 70), ylim = c(0, 1),
     xlab = "Attained age, years", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")

# Add x-axis
axis(1, at = seq(5, 70, 5), cex.axis = 0.9)

# Add "Number at risk" label
text(5, 0.5, "Number at risk", adj = 1, cex = 0.9)

# Add risk numbers at midpoints
for (i in seq_along(midpoints)) {
  text(midpoints[i], 0.5, n_at_risk[i], cex = 0.9)
}

dev.off()

# Also create PNG version
png("mort3-04a-cummortfig-atrisk.png", width = 8, height = 7, units = "in", res = 300)

layout(matrix(c(1, 2), nrow = 2), heights = c(0.85, 0.15))

# ---- MAIN PLOT ----
par(mar = c(0.5, 4, 2, 2))

plot(NULL, xlim = c(5, 70), ylim = c(0, 16),
     xlab = "", ylab = "Cumulative mortality, %",
     xaxt = "n", yaxt = "n", bty = "n")

axis(2, at = 0:16, las = 1, cex.axis = 0.9)
abline(h = 0:16, col = "gray90", lty = 1, lwd = 0.5)

lines(ci_list$recur$time, ci_list$recur$ci, col = "black", lwd = 2)
lines(ci_list$spn$time, ci_list$spn$ci, col = "blue", lwd = 2)
lines(ci_list$circulation$time, ci_list$circulation$ci, col = "red", lwd = 2)
lines(ci_list$respiratory$time, ci_list$respiratory$ci, col = "green", lwd = 2)
lines(ci_list$external$time, ci_list$external$ci, col = "black", lwd = 2)
lines(exp_smooth$x, exp_smooth$y, col = "cyan", lwd = 2, lty = 2)

axis(1, at = seq(5, 70, 5), labels = FALSE)

legend("topleft",
       legend = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality"),
       col = c("black", "blue", "red", "green", "black", "cyan"),
       lty = c(1, 1, 1, 1, 1, 2),
       lwd = 2,
       bty = "n",
       cex = 0.8)

# ---- RISK TABLE ----
par(mar = c(4, 4, 0, 2))

plot(NULL, xlim = c(5, 70), ylim = c(0, 1),
     xlab = "Attained age, years", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")

axis(1, at = seq(5, 70, 5), cex.axis = 0.9)

text(5, 0.5, "Number at risk", adj = 1, cex = 0.9)

for (i in seq_along(midpoints)) {
  text(midpoints[i], 0.5, n_at_risk[i], cex = 0.9)
}

dev.off()

cat("Plots saved as:\n")
cat("  - mort3-04a-cummortfig-atrisk.pdf\n")
cat("  - mort3-04a-cummortfig-atrisk.png\n")
