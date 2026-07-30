# A. oryzae Secretory Pathway Atlas

A current, traceable map of the genes that *Aspergillus oryzae* uses to fold,
modify, and secrete proteins — rebuilt from the 2014 reference dataset and
resolved against today's databases.

## Why this exists

*A. oryzae* (the koji mold behind soy sauce and miso) is one of the best
protein secretors in industry. But its genome is thinly characterized, and the
standard reference list of its secretion machinery was published in 2014 using
gene identifiers that predate a decade of annotation changes.

This project takes that list, connects it to current databases, keeps its
experimental evidence intact, and outputs something you can draw as a diagram.

## What this project produces

After running all three notebooks:

| File | What it is |
|---|---|
| `aoryzae_secretory_components.csv` | The secretion machinery — 369 genes, one row each |
| `aoryzae_secretory_clients.csv` | 2,269 proteins predicted to travel *through* that machinery |
| `identifier_crosswalk.csv` | Old gene IDs mapped to current UniProt / KEGG / NCBI IDs |
| `unresolved_identifiers.csv` | IDs needing manual review |
| `transcriptomic_shortlist.csv` | Genes that changed consistently in high-secretion strains |
| `pathway_nodes.csv` | Layout data for the pathway diagram |

All land in `data/processed/`.

**Machinery and clients are kept separate on purpose.** Machinery is the
shipping department; clients are the cargo. Mixing them produces a misleading
table.

## How it works

| Step | Plain English | Main output |
|---|---|---|
| **1. Profile** | Read and clean the 2014 data. Separate machinery from cargo, check IDs, keep all experimental columns. | Cleaned source tables + quality report |
| **2. Update IDs** | Match 2014 gene IDs to current database IDs. | Crosswalk + unresolved list |
| **3. Build** | Add functions, pathway positions, glycosylation roles, evidence tiers. | Final atlas + diagram inputs |

## Project status

**All three notebooks have run. The current build uses UniProt-backed mapping;
KEGG REST was unavailable, so KO and KEGG pathway enrichment remain incomplete.**

- [x] Profile the 2014 supplementary tables
- [x] Confirm subsystem vocabulary (19 labels found)
- [x] Confirm significance criteria — reproduces the published result exactly
- [x] Build the identifier crosswalk
- [x] Produce the machinery and client atlas tables
- [x] Produce pathway-node layout data
- [ ] Add KO and KEGG pathway enrichment when licensed API access is available
- [ ] Render the pathway visualization

### What step 1 found

- 369 machinery genes, 2,269 predicted clients — no duplicate or missing IDs
- Every ID matches the `AO090` + 9 digits format
- 51 machinery genes changed significantly across all three high-secretion
  strains (48 up, 3 down) — matches the published figures
- 6 yeast genes map to more than one *A. oryzae* gene
- **260 of 369 machinery genes have no subsystem label**, which limits how many
  can be placed on a pathway diagram without further assignment
- Clients: this workbook yields 116 significant genes where the paper reports
  111. Recorded as a discrepancy, not silently overridden.

## Quickstart

**Windows (PowerShell)**
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pytest -q
```

**macOS / Linux**
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q
```

Then open the notebooks in VS Code or Jupyter and run them in order:

1. `notebooks/01_profile_liu2014.ipynb`
2. `notebooks/02_fetch_and_crosswalk.ipynb`
3. `notebooks/03_build_component_table.ipynb`

## What a row looks like

Liu's Table S1 has 14 original columns per gene. Notebook 01 preserves all 14
and appends `sig_all_three` and `direction_all_three` only after verifying the
exact adjusted-p-value columns and the `< 0.05` criterion.

This is a real row from Liu's workbook, not placeholder data. It is the ALG2
ortholog `AO090120000461`, one of the 51 machinery genes significant in all
three strains:

| Liu Table S1 column | Actual value |
|---|---|
| `ID` | `AO090120000461` |
| `S. cerevisiae ortholog` | `ALG2` |
| `Subsystems or function` | `Dolichol pathway` |
| `Description` | Dolichol pathway; glycosylation/modification, mannosyltransferase |
| `CF1.1 vs A1560_logFC` | `0.560207` |
| `A16 vs A1560_logFC` | `0.456661` |
| `CF32 vs A1560_logFC` | `0.781801` |
| `CF1.1vs A1560_adj.P.Val` | `0.004760` |
| `A16 vs A1560_adj.P.Val` | `0.006819` |
| `CF32 vs A1560_adj.P.Val` | `0.000714` |
| `1st SOURCE` | `inparanoid+besthit` |
| `2nd SOURCE` | `wang et al 2010-AO` |
| `3rd SOURCE` | blank |
| `4th SOURCE` | blank |
| `sig_all_three` *(derived by Notebook 01)* | `True` |
| `direction_all_three` *(derived by Notebook 01)* | `up` |

To inspect the full 51-gene machinery shortlist:

```python
import pandas as pd

df = pd.read_csv("data/interim/liu_components_raw_cleaned.csv")
df.loc[df["sig_all_three"]].head(10)
```

This shortlist is worth reviewing before adding database mappings or new
biological inference.

### What the final atlas adds

| Liu gives us | The final atlas adds |
|---|---|
| 2014 `ID` | Preserved `liu_ao_locus_tag`, current `ao_locus_tag`, and available UniProt, KEGG-gene, and NCBI identifiers; KO remains blank in this build |
| Yeast ortholog | Current *A. oryzae* gene name when available |
| Subsystem label on 109 of 369 rows | Normalized subsystem plus `pathway_order` for diagram layout; blank source assignments remain visible |
| Free-text description | Normalized function, compartment, and `glycosylation_role` |
| No provenance specific to the glycosylation assignment | `glycosylation_role_source`, distinguishing Liu, KEGG-orthology, and literature assignments |
| Six raw expression values | The six values preserved plus precomputed `sig_all_three` and `direction_all_three` |
| Four free-text `SOURCE` columns | Normalized `evidence_source`, `record_origin`, and citation fields while retaining the original columns upstream |
| No mapping outcome | `mapping_status`, `mapping_method`, and `manual_review_required` |
| Machinery list and client list in separate workbook sheets | Separate machinery/client outputs and an explicit `record_type` on any combined export |
| No unresolved-ID deliverable | An unresolved and ambiguous identifier report |

The same ALG2 row now appears in the generated component table as follows:

| Selected final field | Generated value for the real ALG2 source row |
|---|---|
| `record_id` | `AOR00008` |
| `liu_ao_locus_tag` | `AO090120000461` |
| `ao_locus_tag` | `AO090120000461` |
| `gene_name` | blank; UniProt has no primary gene name for this entry |
| `function` | Dolichol pathway; glycosylation/modification, mannosyltransferase |
| `record_type` | `machinery_component` |
| `subsystem` | `dolichol_pathway` |
| `pathway_order` | `20` |
| `glycosylation_role` | `n_glycan_assembly` |
| `glycosylation_role_source` | `liu2014_table_s1_subsystem` |
| `sig_all_three` | `True` |
| `direction_all_three` | `up` |
| `record_origin` | `liu2014` |
| `mapping_status` | `exact` |
| `mapping_method` | direct AO090 locus-tag match |
| `uniprot_accession` | `Q2U5Y1` |
| `kegg_gene_id` | `aor:AO090120000461` (from UniProt's cross-reference) |
| `kegg_ko` | blank; KEGG enrichment was unavailable |
| `ncbi_gene_id` | `5996426` |

This row is illustrative, not evidence that every atlas field has the same
coverage. See `data/processed/build_audit.json` and
`data/processed/unresolved_identifiers.csv` before downstream analysis.


## How we verify against the source

**Aggregate checks** — row counts, identifier format, and recomputation of
Liu's significance criteria, which reproduces the published result exactly:
51 machinery genes changed significantly across all three high-secretion
strains, 48 up and 3 down.

**Row-level fidelity** — for every gene, fields derived from Liu are checked
against Liu's original cells, which the cleaned tables retain alongside the
derived columns. This catches the failure mode aggregate checks miss: a
correct total with individual genes in the wrong place.

**Composite labels are flagged, not silently split.** Some of Liu's subsystem
cells name more than one process (e.g. "ERADL; Protein folding/UPR"). These are
reported and assigned deliberately rather than by taking the first token.

**Eyeball sample** — `data/processed/source_fidelity_sample.csv` shows 20 genes
with Liu's original values beside ours for manual inspection.

Known discrepancy: for the client list, the workbook yields 116 genes
significant in all three comparisons where the paper reports 111. Documented,
not overridden. The machinery result matches exactly.