# ==============================================================================
# FIGURE 1 - CUMULATIVE MORTALITY WITH NUMBER AT RISK TABLE (R VERSION)
# ==============================================================================

# Load required packages
library(tidyverse)
library(survival)
library(cmprsk)
library(haven)
library(gridExtra)
library(grid)
library(scales)

# Set parameters
xtime <- 70
cod <- c("spn", "recur", "circulation", "respiratory", "external")

# ==============================================================================
# READ AND PREPARE DATA
# ==============================================================================

# Read the prepared data (assuming it's been saved as .dta or converted to RDS)
# You'll need to adjust the path based on where your temp files are stored
data <- read_dta(file.path(Sys.getenv("temp"), "x-mort3-prepforstset.dta"))

# ==============================================================================
# CALCULATE NUMBER AT RISK BY 10-YEAR AGE BANDS
# ==============================================================================

# Read stset data for all causes
data_allcauses <- read_dta(file.path(Sys.getenv("temp"), "x-mort3-stset-allcauses.dta"))

# Calculate number at risk at start of each 10-year age band
# Age bands: 5-14, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74
# Display at midpoints: 10, 20, 30, 40, 50, 60, 70

n_at_risk <- tibble(
  age = c(5, 15, 25, 35, 45, 55, 65),
  midpoint = age + 5
) %>%
  rowwise() %>%
  mutate(
    n = sum(data_allcauses$`_t0` <= age & data_allcauses$`_t` > age, na.rm = TRUE)
  ) %>%
  ungroup()

# ==============================================================================
# CALCULATE CUMULATIVE INCIDENCE FOR EACH CAUSE OF DEATH
# ==============================================================================

ci_results <- list()

for (cause in cod) {
  # Read the original file for each cause
  temp_data <- data

  # Create event indicator: 1 = event of interest, 2 = competing risk, 0 = censored
  temp_data <- temp_data %>%
    mutate(
      event_status = case_when(
        !!sym(cause) == 1 ~ 1,
        allcauses == 1 & !!sym(cause) == 0 ~ 2,
        TRUE ~ 0
      ),
      # Calculate follow-up time in years (from entry to exit, origin at birth)
      time = as.numeric(difftime(dox, dob, units = "days")) / 365.25,
      entry_time = as.numeric(difftime(doe, dob, units = "days")) / 365.25
    )

  # Use cuminc from cmprsk package for competing risks analysis
  # Need to create appropriate time and status variables
  fit <- with(temp_data,
              cuminc(ftime = time,
                     fstatus = event_status,
                     cencode = 0))

  # Extract cumulative incidence for event of interest (code 1)
  ci_data <- tibble(
    time = fit[[1]]$time,
    ci = fit[[1]]$est * 100  # Convert to percentage
  ) %>%
    filter(time >= 5 & time <= 70)

  ci_results[[cause]] <- ci_data %>% mutate(cause = cause)
}

# Combine all CI results
ci_combined <- bind_rows(ci_results)

# ==============================================================================
# CALCULATE EXPECTED MORTALITY
# ==============================================================================

# Read expected mortality data
expected_data <- read_dta(file.path(Sys.getenv("temp"), "x-mort3-expected-allcause.dta"))

expected_data <- expected_data %>%
  mutate(expected = (1 - conditional) * 100) %>%
  filter(t_exp >= 5 & t_exp <= 70)

# Apply lowess smoothing to match Stata's lowess with bandwidth 0.3
expected_smooth <- loess(expected ~ t_exp, data = expected_data, span = 0.3)
expected_fitted <- tibble(
  t_exp = expected_data$t_exp,
  expected = predict(expected_smooth)
)

# ==============================================================================
# CREATE PLOT WITH NUMBER AT RISK TABLE
# ==============================================================================

# Define colors to match Stata output
colors <- c("recur" = "black", "spn" = "blue", "circulation" = "red",
            "respiratory" = "green", "external" = "black", "expected" = "cyan")

linetypes <- c("recur" = "solid", "spn" = "solid", "circulation" = "solid",
               "respiratory" = "solid", "external" = "solid", "expected" = "dashed")

# Create main plot
main_plot <- ggplot() +
  # Add cumulative incidence lines for each cause
  geom_line(data = ci_combined %>% filter(cause == "recur"),
            aes(x = time, y = ci, color = "recurrence"), linewidth = 0.7) +
  geom_line(data = ci_combined %>% filter(cause == "spn"),
            aes(x = time, y = ci, color = "spn"), linewidth = 0.7) +
  geom_line(data = ci_combined %>% filter(cause == "circulation"),
            aes(x = time, y = ci, color = "circulation"), linewidth = 0.7) +
  geom_line(data = ci_combined %>% filter(cause == "respiratory"),
            aes(x = time, y = ci, color = "respiratory"), linewidth = 0.7) +
  geom_line(data = ci_combined %>% filter(cause == "external"),
            aes(x = time, y = ci, color = "external"), linewidth = 0.7) +
  # Add expected mortality line
  geom_line(data = expected_fitted,
            aes(x = t_exp, y = expected, color = "expected mortality"),
            linetype = "dashed", linewidth = 0.5) +
  # Scales and labels
  scale_x_continuous(breaks = seq(5, 70, 5), limits = c(5, 70)) +
  scale_y_continuous(breaks = seq(0, 16, 1), limits = c(0, 16)) +
  scale_color_manual(
    values = c("recurrence" = "black", "spn" = "blue", "circulation" = "red",
               "respiratory" = "green", "external" = "black", "expected mortality" = "cyan"),
    breaks = c("recurrence", "spn", "circulation", "respiratory", "external", "expected mortality")
  ) +
  labs(
    x = "",
    y = "Cumulative mortality, %",
    color = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

# Create number at risk table
risk_table <- ggplot(n_at_risk, aes(x = midpoint, y = 1)) +
  geom_text(aes(label = n), size = 3) +
  annotate("text", x = 5, y = 1, label = "Number at risk",
           hjust = 1, size = 3) +
  scale_x_continuous(breaks = seq(5, 70, 5), limits = c(5, 70)) +
  scale_y_continuous(limits = c(0.8, 1.2)) +
  labs(x = "Attained age, years", y = NULL) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, margin = margin(t = 5)),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

# Combine plots
combined_plot <- grid.arrange(
  main_plot,
  risk_table,
  ncol = 1,
  heights = c(5, 0.8)
)

# Save the plot
ggsave("mort3-04a-cummortfig-atrisk.pdf", combined_plot,
       width = 8, height = 6, units = "in")
ggsave("mort3-04a-cummortfig-atrisk.png", combined_plot,
       width = 8, height = 6, units = "in", dpi = 300)

# Display the plot
print(combined_plot)
