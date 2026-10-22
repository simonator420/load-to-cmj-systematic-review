# Load to CMJ systematic review

Data and reproducibility package for the systematic review:

> Training and Match Load and Subsequent Countermovement Jump Performance in Field-Based Invasion Team Sports: A Systematic Review

## Author

Šimon Salaj

## Repository contents

```text
data/                Derived review data supporting the reported results
results/             Machine-readable outputs of the cascade and sensitivity scripts
scripts/             Python scripts for the effect-row cascade and the exploratory sensitivity analysis
search_strategies/   Exact database queries, dates and search documentation
```

The release contains the material needed to audit the review:

- the PRISMA selection ledger and the stage-by-stage decisions for all 1,784 candidate reports;
- characteristics for all 379 included reports (192 primary, 39 secondary, 148 contextual);
- 2,093 effect rows from 345 reports, of which 71 form the primary jump-height effect set from 12 independent cohorts;
- cohort-linkage decisions;
- 2,484 QUIPS domain judgements across 414 result units;
- prediction-model assessment and 41 synthesis strata;
- full-text availability categories; and
- executable strategies for PubMed, Scopus, Web of Science, Embase and SPORTDiscus.

No pooled meta-analysis was conducted. One stratum (S01, acceleration and deceleration measures, 0 to <6 h, change in jump height) reached three independent cohorts, but the exposure definitions were incompatible, so the prespecified gate stayed closed. A post hoc exploratory sensitivity analysis in coarser strata is provided in `results/sensitivity/`.

## Reproduce the derived results

Python 3.12 was used for the archived results.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python scripts/build_effect_cascade.py
python scripts/run_exploratory_sensitivity.py
```

Outputs are written to `results/cascade/` and `results/sensitivity/`.

## Exclusions from the public release

Copyrighted article full texts, licensed database exports, figures, manuscripts and local working files are not redistributed. Stable report identifiers, derived data, availability classifications and complete search strategies are provided instead.

## Licence

The derived data and documentation are released under the Creative Commons Attribution 4.0 International licence. This licence does not apply to third-party publications, article full texts or licensed database records.
