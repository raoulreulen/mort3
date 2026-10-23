# Additional analyses for late mortality paper v3 (23 Oct 2026)

Compiled document: `mort3/manuscript/late death v3 - additional analyses v1.0 23Oct2026.docx`
(built by `mort3-19a-rr-compile-additional.do` from `mort3/output/`).

## Pipeline on the Mac
- `mort3-00c-rr-setwd-mac.do` globals; `mort3-00d-rr-runpipeline-mac.do` runs 01a-rrhb, 01b-rr-covardata, 01b-rrhb-gencodcats, then
  `mort3-02a-rr-stsetdata-rates2024.do` (= Hama's stset do-file but population rates used to 2024 instead of capped at 2019;
  all rate files were updated to 1950-2024 on 26 Aug 2026). Follow-up ends 31 Dec 2022.
- Observed deaths unchanged (5,922); expected all-cause deaths now 861.5 (Jan-2026 covariate file + 2020-22 rates) vs 803.3 in tables v1.4.

## Analyses (script -> output)
| # | Script | Output | Method |
|---|---|---|---|
| 10a | mort3-10a-rr-era-poisson.do | table_era_poisson.docx/.csv | Poisson model of SMR by decade of diagnosis adjusted for attained age, sex, cancer type, age at dx (expected deaths offset); recurrence as crude-rate RR |
| 10b | mort3-10b-rr-sex-smraer.do | table_sex_smraer.docx | SMR/AER by sex, p-heterogeneity |
| 10c | mort3-10c-rr-period-analysis.do | table_period_2010_2022.docx | person-time 2010-2022 only ("contemporary") |
| 10d | mort3-10d-rr-sensitivity.do | table_sensitivity.docx | truncate 2019, exclude embarkations (vitstat==1), exclude unknown cause, ICD era |
| 11a-b | mort3-11a-rr-lel-stpm2.R, 11b .do | table_lel.docx, fig_lel_bydiag.png, fig_le_bydecade.png, lel_model_checks.txt | flexible parametric relative survival (rstpm2, attained-age timescale, delayed entry, bhazard = E&W rates), LE/LEL conditional on age 20/30/40/50, 500 MVN draws for CIs, 2022 period expected rates, extrapolated to age 100 |
| 11c-d | mort3-11c-rr-mortality-age-gap.do, 11d .do | table_mortality_age_gap.docx, fig_mortality_age_gap.png | smoothed survivor rate by age vs 2022 E&W rates: equivalent population age |
| 12a-b | mort3-12a-rr-lyl-bycause.R, 12b .do | table_lyl_bycause.docx, fig_lyl_bycause.png | lillies (Andersen decomposition): life years lost before age 80 by cause, bootstrap CIs |
| 13a | mort3-13a-rr-multiple-cause.do | table_multiple_cause.docx, fig_multiple_cause.png | any-mention vs underlying cause; contributory causes only exist from 1993 (n=4,097 deaths) |
| 13b-c | mort3-13b-rr-aer-heatmap.do, 13c .R | fig_aer_heatmap.png, fig_smr_heatmap.png, aer_heatmap_data.csv | AER/SMR by attained age x decade of diagnosis |

## Two life-expectancy numbers, two estimands
- stpm2 LEL (11a): model-based, integrates to age 100 -> overall LEL at age 20 = 7.6 y (12% of expected life).
- lillies (12a): non-parametric, restricted to age 80 (tau) -> excess LYL at age 20 = 8.4 y; its purpose is the split by cause.
Quote 11a for "years of life lost" and 12a for "of which x years to SPN, y to recurrence...".

## Caveats to carry into the text
- Era models: adjusting for attained age reverses the crude trend; excess mortality from causes other than recurrence is ~1.3x higher
  (SPN ~1.7x) in later decades than <1970 (more intensive chemotherapy era / closer to diagnosis). Decade and time since diagnosis
  are near-collinear at ages 40+.
- LEL by decade for 1990s/2000s is extrapolated beyond ~age 40-50 (no person-time).
- lillies carries the cumulative incidence flat beyond the last observed age: strata not followed to 80 are flagged (*).
- smraer.ado AER CI uses invgammap(|O-E|) which is not a proper CI when O<E (negative AER cells).
- 28 deaths coded U07.1 (COVID-19) are codcat 0 "unknown"; consider a COVID row/footnote.
- Multiple-cause: ICD-9 era (1993-2000) certificates carry fewer contributory causes than ICD-10.
