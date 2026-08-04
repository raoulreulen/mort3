# ========================================================================================
# BCCSS MORTALITY DATA - MENTAL HEALTH PROJECT (HAMA)
# READ DATA (N=5,922 DEATHS)
# DATE: 12 June 2025
# Converted from Stata to R
# ========================================================================================

# ========================================================================================
# TUTORIAL NOTES: GETTING STARTED WITH R
# ========================================================================================
# This script demonstrates how to read, merge, and manipulate data in R.
#
# KEY DIFFERENCES FROM STATA:
# - Stata uses "$data" for global macros; R uses variables or file paths directly
# - Stata's "use" command = R's read_dta() or read.csv()
# - Stata's "merge" command = R's merge() or dplyr's left_join(), inner_join(), etc.
# - Stata's "gen" command = R's <- assignment or mutate()
# - Stata's "label" command = R's factor() with labels or attributes
# - Stata's "." for missing values = R's NA
# ========================================================================================

# -----------------------------------------------------------------------------------------
# STEP 1: LOAD REQUIRED PACKAGES
# -----------------------------------------------------------------------------------------
# In R, you need to load packages (libraries) that contain functions you'll use.
# The 'haven' package reads Stata .dta files
# The 'dplyr' package provides data manipulation functions
# The 'lubridate' package makes working with dates easier

# Install packages if you haven't already (run these once):
# install.packages("haven")
# install.packages("dplyr")
# install.packages("lubridate")

# Load the packages for this session:
library(haven)      # For reading Stata .dta files
library(dplyr)      # For data manipulation (mutate, select, etc.)
library(lubridate)  # For working with dates

# -----------------------------------------------------------------------------------------
# STEP 2: SET UP FILE PATHS
# -----------------------------------------------------------------------------------------
# STATA EQUIVALENT: global data "path/to/data"
# In R, we create variables to store our file paths.
# This makes it easy to change paths in one place if needed.

# IMPORTANT: Update these paths to match your actual folder structure!
data_path <- "D:/DATAWORK/mort3"
temp_path <- "D:/TEMPWORK"

# Note: R uses forward slashes (/) or double backslashes (\\) in file paths,
# even on Windows. Single backslashes (\) won't work!

# -----------------------------------------------------------------------------------------
# STEP 3: READ THE MAIN DATA FILE
# -----------------------------------------------------------------------------------------
# STATA COMMAND: use "$data/BCCSSMort_Apr25.dta", clear
# R EQUIVALENT: read_dta()
#
# read_dta() reads Stata .dta files and stores them as a data frame (R's table format)
# The '<-' symbol is R's assignment operator (like '=' in Stata's gen command)

mort_data <- read_dta(file.path(data_path, "BCCSSMort_Apr25.dta"))

# TIP: file.path() automatically creates the correct file path for your operating system
# Alternative: mort_data <- read_dta(paste0(data_path, "/BCCSSMort_Apr25.dta"))

cat("Main data loaded. Number of rows:", nrow(mort_data), "\n")
# cat() prints messages to the console (like 'display' in Stata)
# nrow() counts the number of rows (observations)

# -----------------------------------------------------------------------------------------
# STEP 4: FIRST MERGE - ADD NEOPLASM CLASSIFICATION DATA
# -----------------------------------------------------------------------------------------
# STATA COMMAND: merge 1:1 lngpk_bccss o_indexno using "$data/BCCSSMort_Apr25_NEOP.dta"
# R EQUIVALENT: merge() or left_join()
#
# MERGE TYPES:
# - Stata "merge 1:1" = R's merge(..., all=FALSE) for inner join
# - Stata "merge 1:m" or "m:1" = same in R, just match multiple records
# - Stata "_merge" variable = R doesn't create this automatically, but we can check

neop_data <- read_dta(file.path(data_path, "BCCSSMort_Apr25_NEOP.dta"))

# Perform the merge on matching variables (keys)
# all.x = TRUE means keep all records from mort_data (left join)
# all.y = FALSE means don't keep unmatched records from neop_data
mort_data <- merge(mort_data, neop_data,
                   by = c("lngpk_bccss", "o_indexno"),
                   all.x = TRUE,
                   all.y = FALSE)

# STATA EQUIVALENT: assert _merge==3 | _merge==1
# In Stata, _merge==3 means matched, _merge==1 means only in master (left) data
# In R, we check this manually if needed (but merge handles it automatically)

cat("After NEOP merge. Number of rows:", nrow(mort_data), "\n")

# -----------------------------------------------------------------------------------------
# STEP 5: SECOND MERGE - ADD CAUSE OF DEATH CATEGORIES
# -----------------------------------------------------------------------------------------
# STATA COMMAND: merge 1:1 ... using "..." , keepusing(cat suicide circcat neop_recurr)
# R EQUIVALENT: merge() and then select only certain columns
#
# The 'keepusing' option in Stata means "only bring these variables from the merge file"
# In R, we can select columns before or after merging

cod_data <- read_dta(file.path(data_path, "BCCSSMortCODCat.dta"))

# Select only the columns we want to keep from cod_data
# (lngpk_bccss and o_indexno are needed for merging, plus the 4 variables we want)
cod_data_selected <- cod_data %>%
  select(lngpk_bccss, o_indexno, cat, suicide, circcat, neop_recurr)

# TUTORIAL NOTE: THE PIPE OPERATOR %>%
# The %>% operator (called "pipe") takes the output from the left and feeds it to the right
# Think of it like a pipeline: data flows from left to right
# cod_data %>% select(...) means "take cod_data, THEN select these columns"
# This is similar to chaining commands in Stata with '///'

# Perform inner join (only keep matches, equivalent to assert _merge==3)
mort_data <- merge(mort_data, cod_data_selected,
                   by = c("lngpk_bccss", "o_indexno"),
                   all.x = FALSE,  # Don't keep unmatched from left
                   all.y = FALSE)  # Don't keep unmatched from right (inner join)

cat("After COD categories merge. Number of rows:", nrow(mort_data), "\n")

# -----------------------------------------------------------------------------------------
# STEP 6: THIRD MERGE - ADD OFFSET DATA
# -----------------------------------------------------------------------------------------
# STATA COMMAND: merge 1:1 lngpk_bccss using "$temp/OFFSET_DATA_BCCSS.dta"
# Then: keep if _merge == 3 (only keep matched records)

offset_data <- read_dta(file.path(temp_path, "OFFSET_DATA_BCCSS.dta"))

# Merge and keep only matched records (inner join)
mort_data <- merge(mort_data, offset_data,
                   by = "lngpk_bccss",
                   all.x = FALSE,
                   all.y = FALSE)

cat("After offset merge. Number of rows:", nrow(mort_data), "\n")

# -----------------------------------------------------------------------------------------
# STEP 7: APPLY OFFSETS TO CREATE NATIVE VARIABLES
# -----------------------------------------------------------------------------------------
# STATA COMMAND: gen indexno = o_indexno - indexno_offset
# R EQUIVALENT: mutate() or direct assignment
#
# UNDERSTANDING OFFSETS:
# Offsets are used to decode/de-identify data. The original values (o_*) have been
# shifted by an offset amount. We subtract the offset to get the true values.
#
# TWO WAYS TO CREATE VARIABLES IN R:
# Method 1: Direct assignment with $
#   mort_data$indexno <- mort_data$o_indexno - mort_data$indexno_offset
# Method 2: Using dplyr's mutate() [RECOMMENDED for multiple variables]

mort_data <- mort_data %>%
  mutate(
    # Subtract offsets to get true values
    indexno = o_indexno - indexno_offset,
    cohort  = o_cohort - cohort_offset,
    country = o_country - country_offset,
    icdver  = o_icdver - icdver_offset,

    # Create separate day, month, year variables for date of death
    dod_d = o_dod_day - dod_offset,
    dod_m = o_dod_mth - dod_offset,
    dod_y = o_dod_yr - dod_offset
  )

# TUTORIAL NOTE: mutate()
# mutate() creates new variables or modifies existing ones
# You can create multiple variables at once (separated by commas)
# The new variables are immediately available within the same mutate() call

# -----------------------------------------------------------------------------------------
# STEP 8: CREATE DATE VARIABLE FROM DAY, MONTH, YEAR
# -----------------------------------------------------------------------------------------
# STATA COMMAND: gen dod = mdy(dod_m, dod_d, dod_y)
#                format dod %td
# R EQUIVALENT: make_date() from lubridate package
#
# In Stata, mdy() creates a date from month, day, year
# In R, make_date() does the same but takes year, month, day (different order!)

mort_data <- mort_data %>%
  mutate(
    dod = make_date(year = dod_y, month = dod_m, day = dod_d)
  ) %>%
  select(-dod_d, -dod_m, -dod_y)  # Drop the temporary day/month/year variables

# TUTORIAL NOTE: select() with minus sign
# select(-variable) means "remove this variable" (like 'drop' in Stata)
# select(variable1, variable2) means "keep only these variables" (like 'keep' in Stata)

# -----------------------------------------------------------------------------------------
# STEP 9: TIDY UP - REMOVE OFFSET AND ORIGINAL VARIABLES
# -----------------------------------------------------------------------------------------
# STATA COMMAND: drop o_* *offset
# R EQUIVALENT: select() with -starts_with() or -ends_with()
#
# We want to remove:
# - All variables starting with "o_" (original offset values)
# - All variables ending with "offset"

mort_data <- mort_data %>%
  select(-starts_with("o_"), -ends_with("offset"))

# TUTORIAL NOTE: HELPER FUNCTIONS IN select()
# - starts_with("prefix"): selects columns starting with "prefix"
# - ends_with("suffix"): selects columns ending with "suffix"
# - contains("text"): selects columns containing "text"
# - matches("regex"): selects columns matching a regular expression
# - Using minus (-) before any of these removes those columns

# -----------------------------------------------------------------------------------------
# STEP 10: LABEL VARIABLES (ADD DESCRIPTIONS)
# -----------------------------------------------------------------------------------------
# STATA COMMAND: label var indexno "unique identifier"
# R EQUIVALENT: Setting variable labels as attributes
#
# In Stata, variable labels are built-in and show up in output
# In R, labels are less commonly used, but we can add them as attributes
# They'll be preserved when saving to Stata format with write_dta()

# Method 1: Using attr() for individual variables
attr(mort_data$indexno, "label") <- "unique identifier"
attr(mort_data$ucause, "label") <- "underlying cause of death"
attr(mort_data$icdver, "label") <- "ICD version applicable at time of death"
attr(mort_data$country, "label") <- "country of cancer diagnosis"
attr(mort_data$cohort, "label") <- "cohort up to 1991 (1); cohort 1992-2006 (2)"

# TUTORIAL NOTE: VARIABLE LABELS IN R
# Variable labels in R are mainly useful if you plan to export to Stata
# For R-only analysis, descriptive variable names are often preferred
# Example: Instead of 'dod' with label "date of death", just name it 'date_of_death'

# -----------------------------------------------------------------------------------------
# STEP 11: CREATE VALUE LABELS (CATEGORICAL VARIABLES)
# -----------------------------------------------------------------------------------------
# STATA COMMAND: label define lblcountry 1 "England (1)" 2 "Scotland (2)" 3 "Wales (3)"
#                label values country lblcountry
# R EQUIVALENT: Convert to factor with labels
#
# In Stata, you create label definitions and apply them to variables
# In R, we convert numeric variables to factors with labeled levels

# Convert country to a factor with labels
mort_data <- mort_data %>%
  mutate(
    country = factor(country,
                     levels = c(1, 2, 3),
                     labels = c("England", "Scotland", "Wales"))
  )

# TUTORIAL NOTE: FACTORS IN R
# Factors are R's way of handling categorical data
# - levels: the underlying numeric values
# - labels: the text descriptions shown to users
# - factor() converts a variable to a categorical type
# Example: country = 1 displays as "England"

# Alternative approach for yes/no variables (like in Stata's lblyesno):
# If you had a yes/no variable, you'd do:
# mort_data$some_yesno_var <- factor(mort_data$some_yesno_var,
#                                    levels = c(0, 1),
#                                    labels = c("no", "yes"))

# -----------------------------------------------------------------------------------------
# STEP 12: REORDER COLUMNS
# -----------------------------------------------------------------------------------------
# STATA COMMAND: order indexno lngpk_bccss cohort country dod icdver
# R EQUIVALENT: select() in the order you want
#
# In Stata, 'order' moves specified variables to the front
# In R, we use select() to specify the exact order we want

mort_data <- mort_data %>%
  select(indexno, lngpk_bccss, cohort, country, dod, icdver, everything())

# TUTORIAL NOTE: everything()
# everything() is a helper function that means "all other variables not yet mentioned"
# So this line says: "Put these 6 variables first, then all the rest"

# -----------------------------------------------------------------------------------------
# STEP 13: SAVE THE PROCESSED DATA
# -----------------------------------------------------------------------------------------
# STATA COMMAND: save "$temp/x-mort3-readandlabel", replace
# R EQUIVALENT: write_dta() for Stata format, or saveRDS() for R format
#
# You have two main options for saving in R:

# OPTION 1: Save as Stata .dta file (to maintain compatibility with Stata)
write_dta(mort_data, file.path(temp_path, "x-mort3-readandlabel.dta"))
cat("Data saved as Stata .dta file\n")

# OPTION 2: Save as R data file (.rds format) [RECOMMENDED for R-only workflows]
# .rds files are faster to read/write and preserve all R data types perfectly
saveRDS(mort_data, file.path(temp_path, "x-mort3-readandlabel.rds"))
cat("Data saved as R .rds file\n")

# To read the .rds file later, use:
# mort_data <- readRDS(file.path(temp_path, "x-mort3-readandlabel.rds"))

# OPTION 3: Save as CSV (for maximum compatibility, but loses labels and date formats)
# write.csv(mort_data, file.path(temp_path, "x-mort3-readandlabel.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------------------
# STEP 14: VIEW YOUR DATA
# -----------------------------------------------------------------------------------------
# STATA EQUIVALENT: browse, describe, summarize, list
# R EQUIVALENTS: Multiple options for viewing data

# View basic information about the dataset:
cat("\n=== DATA SUMMARY ===\n")
cat("Number of observations (rows):", nrow(mort_data), "\n")
cat("Number of variables (columns):", ncol(mort_data), "\n")

# Show structure of the data (similar to Stata's 'describe'):
cat("\n=== DATA STRUCTURE ===\n")
str(mort_data)

# Show first few rows (similar to Stata's 'list in 1/6'):
cat("\n=== FIRST 6 ROWS ===\n")
print(head(mort_data))

# Show summary statistics (similar to Stata's 'summarize'):
cat("\n=== SUMMARY STATISTICS ===\n")
print(summary(mort_data))

# Open data viewer in RStudio (similar to Stata's 'browse'):
# View(mort_data)  # Uncomment this line if running in RStudio

# ========================================================================================
# QUICK REFERENCE: STATA TO R COMMAND TRANSLATION
# ========================================================================================
#
# DATA IMPORT/EXPORT:
# Stata: use "file.dta", clear          →  R: data <- read_dta("file.dta")
# Stata: save "file.dta", replace       →  R: write_dta(data, "file.dta")
# Stata: import delimited "file.csv"    →  R: data <- read.csv("file.csv")
#
# DATA MANIPULATION:
# Stata: gen newvar = oldvar * 2        →  R: data$newvar <- data$oldvar * 2
#                                            or: data <- data %>% mutate(newvar = oldvar * 2)
# Stata: replace var = 0 if var < 0     →  R: data$var[data$var < 0] <- 0
#                                            or: data <- data %>% mutate(var = ifelse(var < 0, 0, var))
# Stata: drop varname                   →  R: data$varname <- NULL
#                                            or: data <- data %>% select(-varname)
# Stata: keep var1 var2                 →  R: data <- data %>% select(var1, var2)
# Stata: rename oldname newname         →  R: data <- data %>% rename(newname = oldname)
#
# MERGING DATA:
# Stata: merge 1:1 id using "file.dta"  →  R: data <- merge(data, other_data, by = "id")
#                                            or: data <- data %>% inner_join(other_data, by = "id")
#
# SUBSETTING:
# Stata: keep if var > 10               →  R: data <- data %>% filter(var > 10)
# Stata: drop if var == .               →  R: data <- data %>% filter(!is.na(var))
#
# SORTING:
# Stata: sort var1 var2                 →  R: data <- data %>% arrange(var1, var2)
# Stata: gsort -var1                    →  R: data <- data %>% arrange(desc(var1))
#
# VIEWING DATA:
# Stata: browse                         →  R: View(data)
# Stata: list in 1/10                   →  R: head(data, 10)
# Stata: describe                       →  R: str(data) or glimpse(data)
# Stata: summarize                      →  R: summary(data)
# Stata: tab var                        →  R: table(data$var)
#
# ========================================================================================

cat("\n=== SCRIPT COMPLETED SUCCESSFULLY ===\n")
cat("Processed mortality data is ready for analysis!\n")

# NEXT STEPS:
# 1. Check that your file paths in STEP 2 are correct
# 2. Ensure all required input files exist in the data and temp folders
# 3. Run this script line by line in RStudio to understand each step
# 4. Modify the code as needed for your specific analysis
# 5. Refer to the Quick Reference section for common Stata-to-R translations
