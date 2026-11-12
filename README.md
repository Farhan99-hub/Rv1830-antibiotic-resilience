# Rv1830 (ResR) and drug resistance in *Mycobacterium tuberculosis* — data and code

Supporting data and analysis code for the correspondence
**"Protein-changing variation in Rv1830 (ResR) is associated with drug resistance
in *Mycobacterium tuberculosis*: an integrated analysis of 40,276 clinical genomes."**

This repository is self-contained and reproducible: raw source data → a documented
cleaning/harmonisation step → the analysis set → every figure, table and statistic
in the manuscript, together with the code that links them.

## Contents

```
manuscript/
  Manuscript_07112026.docx          Final manuscript (as submitted)
  Manuscript_07112026.pdf           PDF rendering
figures/
  Figure1.png                       Figure 1, 300 dpi
  Figure1.pdf                       Figure 1, vector
  Figure1_source_data.xlsx          Exact numbers behind panels a–e
tables/
  Table1_main.xlsx                  Table 1 — carriage by resistance category
  Supplementary_Table_S1_MIC_by_drug.xlsx     Eight-drug MIC comparison
  Supplementary_Table_S2_substitutions.xlsx   All 202 distinct substitutions (Fig 1c)
  Supplementary_Table_S3_robustness.xlsx      Robustness analyses (4 sheets)
data/
  raw/                              Unmodified source data (see provenance below)
    Data_set_1.xlsx                 17,719 genotypically classified isolates (TB-Profiler)
    Data_set_2.csv                  21,685 CRyPTIC-format records with MICs (13 drugs)
    Data_set_3.csv                  8,620 WGS records with phenotypic DST and outcomes
  processed/
    analysis_set_40276.csv          One row per isolate (the analysis-ready table)
  DATA_DICTIONARY.md                Columns, severity mapping, exclusion rule
code/
  build_outputs.py                  Rebuilds analysis set, Table 1, Supp. Table S2
  verify_rv1830.py                  Independent re-verification of every headline statistic
  robustness_rv1830.py              The four robustness analyses
  make_figure1.py                   Regenerates Figure 1
  requirements.txt                  Python dependencies
  results/
    verification_results.json       Output of verify_rv1830.py
    robustness_results.json         Output of robustness_rv1830.py
LICENSE.txt                         Recommended licensing (code MIT; data CC BY 4.0)
CITATION.cff                        Citation metadata for the repository
```

## Analysis set (40,276 isolates)

Three independent collections, harmonised to a shared severity scale, with no shared
isolate accession:

| Collection | Raw records | Excluded | Analysed |
|---|---|---|---|
| C1 (TB-Profiler) | 17,719 | 1,079 "other/uninterpretable" | 16,640 |
| C2 (CRyPTIC MICs) | 21,685 | 6,525 no Rv1830 call; 86 untested | 15,074 |
| C3 (WGS + DST + outcomes) | 8,620 | 58 flagged Rv1830 sequence | 8,562 |
| **Total** | | **7,748** | **40,276** |

A carrier has ≥1 protein-changing Rv1830 variant (missense, nonsense, in-frame indel);
non-coding (promoter) and synonymous changes are set aside and those isolates kept as
non-carriers. Full column definitions and the per-collection severity mapping are in
`data/DATA_DICTIONARY.md`.

## Reproducing the results

```bash
cd code
pip install -r requirements.txt
python build_outputs.py   ../data/raw/Data_set_1.xlsx ../data/raw/Data_set_2.csv ../data/raw/Data_set_3.csv .
python verify_rv1830.py   ../data/raw/Data_set_1.xlsx ../data/raw/Data_set_2.csv ../data/raw/Data_set_3.csv ../manuscript/Manuscript_07112026.docx
python robustness_rv1830.py ../data/raw/Data_set_1.xlsx ../data/raw/Data_set_2.csv ../data/raw/Data_set_3.csv
python make_figure1.py
```

- `build_outputs.py` → `analysis_set_40276.csv`, `Table1_main.xlsx`, `Supplementary_Table_S2_substitutions.xlsx`.
- `verify_rv1830.py` → rebuilds the set and checks each quoted value against the manuscript; writes `verification_results.json`. The S1 MIC numbers are taken from this output.
- `robustness_rv1830.py` → the four robustness analyses; writes `robustness_results.json` (formatted in Supplementary Table S3).
- `make_figure1.py` → `Figure1.png` / `Figure1.pdf`.

## Headline results

839 of 40,276 isolates (2.08%) carry a protein-changing Rv1830 variant; carriers are
resistant more often than non-carriers (63.8% vs 50.9%); crude OR 1.69 (1.47–1.95),
Mantel–Haenszel 1.82 (1.58–2.11). Carriage rises across severity (1.55/2.17/2.50/3.44/
3.00%; trend z = 8.14). D144N (151 isolates, all sublineage 4.1.1.1) is inversely
associated crudely (OR 0.62) but null after lineage adjustment. 44.7% of variant
observations fall in residues 65–137. No rifampicin MIC shift (p = 0.21); no drug
significant after Bonferroni. Robustness analyses (lineage adjustment, random-effects
heterogeneity, target-gene-stratified MIC, treatment outcomes) are in Supplementary
Table S3 and `robustness_results.json`.

## Data provenance and third-party terms

Collection 2 derives from the CRyPTIC Consortium dataset (Zenodo v3.4.1; manuscript
refs 8–9) and Collection 3 from the NIH TB Portals programme. Redistribution of the
raw files is subject to those sources' terms; cite them accordingly. See `LICENSE.txt`.

## Notes

- `verify_rv1830.py` is the authors' own re-verification script (29/29 checks).
- The c-87t promoter variant is deliberately excluded from the analysis as a
  non-coding, lineage-associated change (isolates retained as non-carriers).
- Definitive control of population structure (a genome-wide tree/kinship test) is
  identified in the manuscript as the confirmatory next step.
