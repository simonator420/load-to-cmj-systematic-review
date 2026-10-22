# Post hoc exploratory sensitivity analysis: coarser strata (definitions)

Written: 2026-10-18 22:44:05 CEST, BEFORE any stratum cohort counts or pooled estimates were computed.
Disclosure: the 59 harmonised rows had been listed (fields only) to design the definitions, so the family below is not blind to the data. It is fixed here and not altered after results.
Status: POST HOC, NOT PRESPECIFIED. Eligibility is NOT relaxed: only the 59 rows with poolable_candidate == "yes" in final/08_classified_correlation_rows.json.

## Common rules
- Outcome: jump_height_change or jump_height_absolute as stated per stratum.
- Lag uses numeric lag_hours (rows with missing lag_hours are excluded from every lag-bounded stratum). 24-72 h means 24 <= lag_hours <= 72; 0-<30 h means 0 <= lag_hours < 30; "<=168 h" means 0 <= lag_hours <= 168.
- External load = exposure_domain in {external_match, external_training}. Internal load = exposure_domain == internal_load.
- Independent cohort = cohort token. PRIMARY coding = "cohort" field. SENSITIVITY coding = reports flagged possible_overlap_unresolved in final/07_cohort_linkage.csv (Ishida, Shearer/Russell, Males/Males) get one token per report; confirmed_same_cohort stays merged.
- Gate: >= 3 independent cohorts (protocol analysis_plan.md). Gate is evaluated on cohorts contributing rows (reports and rows also counted).

## Declared strata (6)
| ID | Exposure | Outcome | Lag | Scientific defensibility |
|---|---|---|---|---|
| S1 (user-requested) | any external load | change OR absolute jump height | 24-72 h | Moderately defensible: same broad exposure class and recovery window; STRETCH on mixing change and absolute jump height (different estimands) |
| S1c | any external load | jump_height_change only | 24-72 h | Most defensible coarse stratum (single estimand, single window) |
| S2 | any load (external or internal) | change OR absolute | 24-72 h | STRETCH: internal and external load are different constructs, plus change/absolute mixing |
| S3 | any external load | change OR absolute | 0-<30 h | Defensible as "acute/next-day"; STRETCH on mixing immediate post-exposure (0 h) with next-morning tests, match with small-sided-game training, and change with absolute |
| S4 | any external load | change OR absolute | all known lags <=168 h | Weakest: one week of recovery kinetics collapsed; change/absolute and match/training mixed |
| S5 | any load | change OR absolute | all known lags <=168 h | Weakest; broadest bound, included only as a ceiling on how coarse the gate could possibly be |

Conclusion sensitivity is graded: S1c = strong evidence if it passes; S1 = moderate; S2-S5 = weak (stretch).

## Exploratory synthesis (only for strata passing the gate)
- One effect per cohort. Tier (analysis_plan hierarchy): 1 = within_athlete_modelled; 2 = repeated-measures/athlete-clustered (none identified); 3 = athlete-level change-score association (outcome change, dependence_status independent_athletes); 4 = ordinary correlation (unclear / repeated_observations_pooled / absolute).
- Within the best tier of a cohort, deterministic tie-break NOT based on effect size: construct priority total_distance > player_load > high_speed_running > sprint_distance > sprint_number/count > accel_decel variants > other; then lag closest to the stratum midpoint; then lowest row_id.
- Fisher z; variance 1/(n-3), n = independent athlete-level n only. If only pooled observations are available, the effect is excluded from pooling (cohort still counts toward the gate and this is reported).
- REML tau^2, Hartung-Knapp CI (t, k-1 df), I^2, 95% PI (t, k-2 df) when k>=3, leave-one-cohort-out. Back-transform via tanh.
- Selection-sensitivity: (a) cohort mean of z across all rows of the stratum; (b) alternative extremes (most negative / most positive row per cohort within the best tier) to show the range attributable to the arbitrary pick.
