# Data dictionary

| File | Unit | Purpose |
|---|---|---|
| `data/01_study_selection_and_PRISMA_ledger.csv` | selection stage | Audited PRISMA counts and their provenance. |
| `data/02_report_decisions.csv` | report | Final inclusion or exclusion decision and exclusion reason. |
| `data/03_effect_level_data.csv` | extracted effect | Exposure, outcome, timing, estimand and reported numerical result, with synthesis-eligibility flags. |
| `data/04_study_characterisation.csv` | included report | Citation, population, design, exposure, CMJ outcome, timing and evidence stream. |
| `data/05_cohort_linkage.csv` | report to cohort link | Conservative cohort-independence coding and sensitivity action. |
| `data/06_risk_of_bias_QUIPS.csv` | result unit by QUIPS domain | Domain judgement, assessment status and supporting rationale. |
| `data/06b_quips_source_gaps.csv` | included report | Reports for which complete QUIPS assessment was unavailable. |
| `data/07_prediction_model_assessment.csv` | candidate report | Classification against the prospective-prediction definition. |
| `data/08_synthesis_strata.csv` | compatibility stratum | Primary cohort-independence coding used for the synthesis gate (strata S01 to S41). |
| `data/08b_synthesis_strata_sensitivity.csv` | compatibility stratum | Alternative coding for potentially overlapping cohorts. |
| `data/08c_primary_effects_71.csv` | extracted effect | The 71 primary jump-height effects, with orientation of the correlation and cohort token. |
| `data/08d_gate_candidate_stratum_review.csv` | extracted effect | Rows of stratum S01, the only stratum that reached three independent cohorts. |
| `data/08_classified_correlation_rows.json` | extracted effect | Classification of effect rows by exposure construct, lag window and outcome class. |
| `data/09_meta_analysis_gate.json` | review | Machine-readable result of the prespecified pooling gate. |
| `data/10_fulltext_provenance.csv` | included report | Archived-source and full-text-based-assessment availability categories. |
| `data/11_candidate_report_stage_decisions.csv` | candidate report | Retrieval outcome and final decision for all 1,784 candidate reports. |
| `results/cascade/` | various | Effect-row attrition cascade tables. |
| `results/sensitivity/` | analysis | Results of the exploratory sensitivity analysis in coarser strata. |

Missing, unclear and unavailable values are retained explicitly and must not be interpreted as zero. Report counts are not equivalent to independent study, cohort, athlete or effect counts.
