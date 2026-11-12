# Data dictionary — Rv1830 (ResR) analysis

## analysis_set_40276.csv
One row per isolate; 40,276 isolates. This is the cleaned, harmonised,
analysis-ready table from which every figure panel and table is derived.

| Column | Type | Description |
|---|---|---|
| isolate_id | string | Isolate identifier within its source collection (C1: run accession; C2: CRyPTIC UNIQUEID; C3: SRA id). |
| collection | C1/C2/C3 | C1 genotypic (TB-Profiler); C2 CRyPTIC MIC compendium; C3 WGS + laboratory DST (TB Portals). |
| severity_level | 0–4 | Harmonised ordinal severity (see mapping below). |
| severity_label | string | 0 Susceptible, 1 Other resistance, 2 MDR, 3 Pre-XDR, 4 XDR. |
| resistant | 0/1 | 1 if severity_level ≥ 1. |
| rv1830_carrier | 0/1 | 1 if the isolate carries ≥1 protein-changing Rv1830 variant (missense, nonsense or in-frame indel). Promoter and synonymous changes do not set this flag. |
| rv1830_protein_substitutions | string | Semicolon-separated protein-changing substitutions; blank for non-carriers. |
| D144N_carrier | 0/1 | 1 if D144N is among the carried substitutions. |

### Severity harmonisation (native scheme → level)
- **C1 (TB-Profiler Phenotype):** Sensitive→0; HR-TB, RR-TB→1; MDR-TB→2; Pre-XDR-TB→3; XDR-TB→4.
- **C2 (CRyPTIC resistance_class):** no_resistance_called→0; resistant_other→1; MDR→2; MDR+FQ(pre-XDR)→3; MDR+FQ+inj(XDR-like)→4.
- **C3 (type_of_resistance):** Sensitive→0; Mono DR, Poly DR→1; MDR non XDR→2; Pre-XDR→3; XDR→4.

### Exclusion rule (raw → analysis set)
| Collection | Raw records | Excluded | Analysed |
|---|---|---|---|
| C1 (Data_set_1.xlsx) | 17,719 | 1,079 "other/uninterpretable" | 16,640 |
| C2 (Data_set_2.csv) | 21,685 | 6,525 no Rv1830 call; 86 untested | 15,074 |
| C3 (Data_set_3.csv) | 8,620 | 58 flagged Rv1830 sequence | 8,562 |
| **Total** | | **7,748** | **40,276** |

The three collections share no isolate accession.

## source/ (raw input data, unmodified)
- **Data_set_1.xlsx** — 17,719 genotypically classified isolates (sheet `Combined_17719_isolates`; includes its own `Data_dictionary`).
- **Data_set_2.csv** — 21,685 CRyPTIC-format records; per-drug `<DRUG>_MIC`, `<DRUG>_phenotype`, `<DRUG>_quality` columns for 13 drugs; `Rv1830` column holds the locus variants.
- **Data_set_3.csv** — 8,620 WGS records with Rv1830 QC fields and `type_of_resistance`.

## Reproduction
`code/verify_rv1830.py` rebuilds this analysis set from the three source
files and recomputes every statistic:
```
python verify_rv1830.py Data_set_1.xlsx Data_set_2.csv Data_set_3.csv <manuscript.docx>
```
`code/fig.py` regenerates Figure 1 from the values in `Figure1_source_data.xlsx`.
