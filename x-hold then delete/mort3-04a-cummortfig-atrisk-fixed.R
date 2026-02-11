# ==============================================================================
# FIGURE 1 - CUMULATIVE MORTALITY WITH NUMBER AT RISK TABLE
# Fixed version with better error handling
# ==============================================================================

# Load only essential packages
library(survival)
library(cmprsk)
library(foreign)  # Alternative to haven for reading Stata files

# Set parameters
xtime <- 70
cod <- c("spn", "recur", "circulation", "respiratory", "external")

# ==============================================================================
# SET UP PATHS - ADJUST THESE TO YOUR ACTUAL TEMP DIRECTORY
# ==============================================================================

# Set temp directory to match your Stata $temp location
temp_dir <- "D:/TEMPWORK"

# Set working directory to save plots in the same location
setwd(temp_dir)

cat("Looking for data files in:", temp_dir, "\n")

# Check if files exist
required_files <- c("x-mort3-stset-allcauses.dta",
                   "x-mort3-prepforstset.dta",
                   "x-mort3-expected-allcause.dta")

for (f in required_files) {
  filepath <- file.path(temp_dir, f)
  if (file.exists(filepath)) {
    cat("  Found:", f, "\n")
  } else {
    cat("  MISSING:", f, "\n")
    cat("  Please check your temp directory path!\n")
  }
}

# ==============================================================================
# READ DATA WITH ERROR HANDLING
# ==============================================================================

# Try reading with foreign package (more robust for older Stata files)
cat("\nReading stset-allcauses data...\n")
data_allcauses <- tryCatch({
  read.dta(file.path(temp_dir, "x-mort3-stset-allcauses.dta"),
           convert.factors = FALSE)
}, error = function(e) {
  cat("Error reading with foreign::read.dta, trying haven::read_dta...\n")
  library(haven)
  read_dta(file.path(temp_dir, "x-mort3-stset-allcauses.dta"))
})

cat("Data loaded. Rows:", nrow(data_allcauses), "Cols:", ncol(data_allcauses), "\n")
cat("Column names:", paste(head(names(data_allcauses), 10), collapse = ", "), "...\n")

# ==============================================================================
# CALCULATE NUMBER AT RISK BY 10-YEAR AGE BANDS
# ==============================================================================

ages <- c(5, 15, 25, 35, 45, 55, 65)
midpoints <- ages + 5
n_at_risk <- numeric(length(ages))

# Check for the survival time variables
if ("_t0" %in% names(data_allcauses) && "_t" %in% names(data_allcauses)) {
  for (i in seq_along(ages)) {
    n_at_risk[i] <- sum(data_allcauses$`_t0` <= ages[i] &
                        data_allcauses$`_t` > ages[i], na.rm = TRUE)
  }
  cat("\nNumber at risk calculated:\n")
  for (i in seq_along(midpoints)) {
    cat("  Age", midpoints[i], ":", n_at_risk[i], "\n")
  }
} else {
  cat("WARNING: Could not find _t0 and _t variables\n")
  cat("Available columns:", paste(names(data_allcauses), collapse = ", "), "\n")
  # Use dummy values for testing
  n_at_risk <- rep(1000, length(ages))
}

# ==============================================================================
# CALCULATE CUMULATIVE INCIDENCE FOR EACH CAUSE OF DEATH
# ==============================================================================

cat("\nReading prepforstset data...\n")
data <- tryCatch({
  read.dta(file.path(temp_dir, "x-mort3-prepforstset.dta"),
           convert.factors = FALSE)
}, error = function(e) {
  library(haven)
  read_dta(file.path(temp_dir, "x-mort3-prepforstset.dta"))
})

cat("Data loaded. Rows:", nrow(data), "\n")

# Storage for CI results
ci_list <- list()

cat("\nCalculating cumulative incidence for each cause...\n")
for (cause in cod) {
  cat("  Processing:", cause, "...")

  # Create event indicator
  event_status <- ifelse(data[[cause]] == 1, 1,
                         ifelse(data$allcauses == 1, 2, 0))

  # Calculate follow-up time in years
  time <- as.numeric(difftime(data$dox, data$dob, units = "days")) / 365.25

  # Use cuminc for competing risks analysis
  fit <- cuminc(ftime = time, fstatus = event_status, cencode = 0)

  # Extract cumulative incidence for event of interest
  ci_time <- fit[[1]]$time
  ci_est <- fit[[1]]$est * 100

  # Filter to age range 5-70
  keep_idx <- ci_time >= 5 & ci_time <= 70
  ci_list[[cause]] <- list(time = ci_time[keep_idx], ci = ci_est[keep_idx])

  cat(" Done (", sum(keep_idx), "time points)\n")
}

# ==============================================================================
# CALCULATE EXPECTED MORTALITY
# ==============================================================================

cat("\nReading expected mortality data...\n")
expected_data <- tryCatch({
  read.dta(file.path(temp_dir, "x-mort3-expected-allcause.dta"),
           convert.factors = FALSE)
}, error = function(e) {
  library(haven)
  read_dta(file.path(temp_dir, "x-mort3-expected-allcause.dta"))
})

expected_data$expected <- (1 - expected_data$conditional) * 100

# Filter to age range
keep_exp <- expected_data$t_exp >= 5 & expected_data$t_exp <= 70
exp_time <- expected_data$t_exp[keep_exp]
exp_values <- expected_data$expected[keep_exp]

# Apply lowess smoothing
exp_smooth <- lowess(exp_time, exp_values, f = 0.3)

cat("Expected mortality calculated (", length(exp_smooth$x), "time points)\n")

# ==============================================================================
# CREATE PLOT WITH NUMBER AT RISK TABLE
# ==============================================================================

cat("\nCreating plots...\n")

# Set up PDF output
pdf("mort3-04a-cummortfig-atrisk.pdf", width = 8, height = 7)

# Create layout: main plot (80%) and risk table (20%)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.80, 0.20))

# ---- MAIN PLOT ----
par(mar = c(0, 4, 2, 2))

# Find the maximum y value across all curves
max_y <- max(
  max(ci_list$recur$ci, na.rm = TRUE),
  max(ci_list$spn$ci, na.rm = TRUE),
  max(ci_list$circulation$ci, na.rm = TRUE),
  max(ci_list$respiratory$ci, na.rm = TRUE),
  max(ci_list$external$ci, na.rm = TRUE),
  max(exp_smooth$y, na.rm = TRUE)
)
# Round up to nearest integer
ylim_max <- ceiling(max_y)

# Initialize plot with dynamic y-axis
plot(NULL, xlim = c(5, 70), ylim = c(0, ylim_max),
     xlab = "", ylab = "Cumulative mortality, %",
     xaxt = "n", yaxt = "n", bty = "n")

# Add grid
axis(2, at = 0:ylim_max, las = 1, cex.axis = 0.9)
abline(h = 0:ylim_max, col = "gray90", lty = 1, lwd = 0.5)

# Add cumulative incidence lines
lines(ci_list$recur$time, ci_list$recur$ci, col = "black", lwd = 2)
lines(ci_list$spn$time, ci_list$spn$ci, col = "blue", lwd = 2)
lines(ci_list$circulation$time, ci_list$circulation$ci, col = "red", lwd = 2)
lines(ci_list$respiratory$time, ci_list$respiratory$ci, col = "green", lwd = 2)
lines(ci_list$external$time, ci_list$external$ci, col = "black", lwd = 2)

# Add expected mortality (dashed line)
lines(exp_smooth$x, exp_smooth$y, col = "cyan", lwd = 2, lty = 2)

# Add legend
legend("topleft",
       legend = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality"),
       col = c("black", "blue", "red", "green", "black", "cyan"),
       lty = c(1, 1, 1, 1, 1, 2),
       lwd = 2,
       bty = "n",
       cex = 0.8)

# ---- RISK TABLE ----
par(mar = c(4.5, 4, 0, 2))

# Create empty plot for risk table with more vertical space
plot(NULL, xlim = c(5, 70), ylim = c(-0.5, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")

# Add x-axis labels at the top (y=1)
axis(1, at = seq(5, 70, 5), cex.axis = 0.9, pos = 1, lwd.ticks = 1, lwd = 0)

# Add x-axis title
mtext("Attained age, years", side = 1, line = 2, cex = 0.9)

# Add "Number at risk" label much lower (y=0)
text(5, 0, "Number at risk", adj = 1, cex = 0.9)

# Add risk numbers at midpoints much lower (y=0)
for (i in seq_along(midpoints)) {
  text(midpoints[i], 0, n_at_risk[i], cex = 0.9)
}

dev.off()

# Also create PNG version
png("mort3-04a-cummortfig-atrisk.png", width = 8, height = 7, units = "in", res = 300)

layout(matrix(c(1, 2), nrow = 2), heights = c(0.80, 0.20))

par(mar = c(0, 4, 2, 2))
plot(NULL, xlim = c(5, 70), ylim = c(0, ylim_max),
     xlab = "", ylab = "Cumulative mortality, %",
     xaxt = "n", yaxt = "n", bty = "n")
axis(2, at = 0:ylim_max, las = 1, cex.axis = 0.9)
abline(h = 0:ylim_max, col = "gray90", lty = 1, lwd = 0.5)
lines(ci_list$recur$time, ci_list$recur$ci, col = "black", lwd = 2)
lines(ci_list$spn$time, ci_list$spn$ci, col = "blue", lwd = 2)
lines(ci_list$circulation$time, ci_list$circulation$ci, col = "red", lwd = 2)
lines(ci_list$respiratory$time, ci_list$respiratory$ci, col = "green", lwd = 2)
lines(ci_list$external$time, ci_list$external$ci, col = "black", lwd = 2)
lines(exp_smooth$x, exp_smooth$y, col = "cyan", lwd = 2, lty = 2)
legend("topleft",
       legend = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality"),
       col = c("black", "blue", "red", "green", "black", "cyan"),
       lty = c(1, 1, 1, 1, 1, 2),
       lwd = 2,
       bty = "n",
       cex = 0.8)

par(mar = c(4.5, 4, 0, 2))
plot(NULL, xlim = c(5, 70), ylim = c(-0.5, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")
axis(1, at = seq(5, 70, 5), cex.axis = 0.9, pos = 1, lwd.ticks = 1, lwd = 0)
mtext("Attained age, years", side = 1, line = 2, cex = 0.9)
text(5, 0, "Number at risk", adj = 1, cex = 0.9)
for (i in seq_along(midpoints)) {
  text(midpoints[i], 0, n_at_risk[i], cex = 0.9)
}

dev.off()

cat("\n=== DONE ===\n")
cat("Plots saved in:\n")
cat("  ", getwd(), "\n")
cat("Files:\n")
cat("  - mort3-04a-cummortfig-atrisk.pdf\n")
cat("  - mort3-04a-cummortfig-atrisk.png\n")

# Display plot in RStudio viewer
cat("\nDisplaying plot in RStudio...\n")

# Recreate plot for display in RStudio
layout(matrix(c(1, 2), nrow = 2), heights = c(0.80, 0.20))

par(mar = c(0, 4, 2, 2))
plot(NULL, xlim = c(5, 70), ylim = c(0, ylim_max),
     xlab = "", ylab = "Cumulative mortality, %",
     xaxt = "n", yaxt = "n", bty = "n")
axis(2, at = 0:ylim_max, las = 1, cex.axis = 0.9)
abline(h = 0:ylim_max, col = "gray90", lty = 1, lwd = 0.5)
lines(ci_list$recur$time, ci_list$recur$ci, col = "black", lwd = 2)
lines(ci_list$spn$time, ci_list$spn$ci, col = "blue", lwd = 2)
lines(ci_list$circulation$time, ci_list$circulation$ci, col = "red", lwd = 2)
lines(ci_list$respiratory$time, ci_list$respiratory$ci, col = "green", lwd = 2)
lines(ci_list$external$time, ci_list$external$ci, col = "black", lwd = 2)
lines(exp_smooth$x, exp_smooth$y, col = "cyan", lwd = 2, lty = 2)
legend("topleft",
       legend = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality"),
       col = c("black", "blue", "red", "green", "black", "cyan"),
       lty = c(1, 1, 1, 1, 1, 2),
       lwd = 2,
       bty = "n",
       cex = 0.8)

par(mar = c(4.5, 4, 0, 2))
plot(NULL, xlim = c(5, 70), ylim = c(-0.5, 1),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")
axis(1, at = seq(5, 70, 5), cex.axis = 0.9, pos = 1, lwd.ticks = 1, lwd = 0)
mtext("Attained age, years", side = 1, line = 2, cex = 0.9)
text(5, 0, "Number at risk", adj = 1, cex = 0.9)
for (i in seq_along(midpoints)) {
  text(midpoints[i], 0, n_at_risk[i], cex = 0.9)
}
