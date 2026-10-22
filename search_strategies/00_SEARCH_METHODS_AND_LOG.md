# Final search methods and database log

## Status

The searches documented in this directory are the final database searches used for screening and analysis in this systematic review.

## Search scope

The search combined three concepts:

1. countermovement jump, vertical jump or neuromuscular performance;
2. quantified training or match load, athlete monitoring, fatigue, recovery or temporally relevant state information;
3. eligible field-based invasion team sports.

The final searches did not encode date or language restrictions in the supplied query strings.

## Databases and results

| Database | Platform | Search date | Raw records | Exact strategy | Result export |
| --- | --- | ---: | ---: | --- | --- |
| MEDLINE/PubMed and PubMed-not-MEDLINE | PubMed via NCBI E-utilities | 2026-09-22 | 1,241 | `01_PubMed_final_query.txt` | `06_PubMed_results.csv`; `07_PubMed_results.xml` |
| Scopus | Scopus | 2026-09-23 | 1,419 | `02_Scopus_final_query.txt` | `09_Scopus_results.csv` |
| Web of Science Core Collection | Web of Science | 2026-09-23 | 1,302 | `03_Web_of_Science_final_query.txt` | `10_WebOfScience_results_part1.xls` (1,000); `11_WebOfScience_results_part2.xls` (302) |
| Embase | Ovid | 2026-09-23 | 1,107 | `04_Embase_Ovid_final_query.txt` | `12_Embase_results.xls` |
| SPORTDiscus | EBSCOhost | 2026-09-23 | 995 | `05_SPORTDiscus_EBSCO_final_query.txt` | `13_SPORTDiscus_results.csv` |
| **Total before deduplication** |  |  | **6,064** |  |  |

The PubMed run started at 2026-09-22 21:29:09 UTC and completed at 21:29:51 UTC. Exact execution times were not retained in the consolidated records for the other four databases; their dated strategy files and exports record the search date.

## Search validation

The final PubMed strategy retrieved all 15 known candidate reports in the expanded seed set. Search quality-control status is recorded in `14_SEARCH_QUALITY_CONTROL.md`.

## File provenance

- `01` to `05` are the exact database-specific query strings.
- `06` to `13` are the original or normalized exports used in the project.
- `08_PubMed_final_run_manifest.json` records the final status, execution timestamp, result count and provenance of the PubMed search.
- `14_SEARCH_QUALITY_CONTROL.md` records known-item validation and the status of independent search review.
- `SHA256SUMS.txt` provides a checksum for every file in this directory except itself.

## Use in reporting

For the Methods section, report that the final searches covered PubMed, Embase via Ovid, SPORTDiscus via EBSCOhost, Scopus and Web of Science Core Collection from database inception through 22-23 September 2026. Report 6,064 records before deduplication. Deduplicated screening and inclusion counts must be taken from the final screening and PRISMA datasets rather than calculated by summing database exports.
