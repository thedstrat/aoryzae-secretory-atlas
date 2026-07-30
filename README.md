# A. oryzae Secretory Pathway Atlas

A traceable, current map of the machinery *Aspergillus oryzae* uses to fold, glycosylate, transport, and secrete proteins.

## Start here

| File | Purpose |
|---|---|
| [`pathway_overview.svg`](data/processed/pathway_overview.svg) | Primary pathway placement with unassigned genes shown explicitly. |
| [`FINDINGS.md`](data/processed/FINDINGS.md) | A concise, evidence-bounded interpretation of the current atlas. |
| [`high_secretion_responsive_genes.csv`](data/processed/high_secretion_responsive_genes.csv) | The 51 machinery genes consistently changed in all three secretion strains. |
| [`secretion_machinery_genes.csv`](data/processed/secretion_machinery_genes.csv) | The complete 369-gene machinery atlas, one biological gene per row. |
| [`column_descriptions.csv`](data/processed/column_descriptions.csv) | Definitions, allowed values, and provenance for every atlas column. |

![A. oryzae secretory pathway overview](data/processed/pathway_overview.png)

## Why this exists

*A. oryzae* is an important industrial protein producer, but its standard secretory-machinery reference was published in 2014. This project preserves that source evidence, resolves the genes against current databases, separates machinery from predicted cargo, and makes uncertainty visible rather than guessing.

## Current numbers

| Measure | Current build |
|---|---:|
| Machinery genes | 369 |
| Predicted secretory clients | 2,269 |
| Transcriptomic shortlist | 51 (48 up, 3 down; matches the published result) |
| Primary subsystem coverage | 247/369 (122 unassigned) |
| KEGG KO coverage | 330/369 |

## How it works

| Step | Plain English | Main output |
|---|---|---|
| **1. Profile** | Read Liu's workbook, separate machinery from cargo, preserve source cells, and verify the expression criteria. | Clean source tables and parsing audit |
| **2. Update IDs** | Resolve Liu-era locus tags against current KEGG, UniProt, and NCBI identifiers. | Identifier crosswalk and KO links |
| **3. Build** | Assign documented pathway positions, compartments, glycosylation roles, and evidence provenance. | Final atlas, shortlist, and visualization |

## Project status

All three notebooks run successfully. KEGG access works, KO enrichment is complete for the current release, and both SVG and PNG visualizations are rendered.

Remaining gaps:

- 189 machinery genes still have compartment `unknown`.
- 122 machinery genes remain without a primary subsystem.
- Liu's workbook yields 116 significant clients while the paper reports 111; the discrepancy is recorded, not overridden.

## Outputs

### Deliverables (`data/processed/`)

| File | Grain | Contents |
|---|---|---|
| `secretion_machinery_genes.csv` | One row per machinery gene | The complete 369-gene atlas with identifiers, annotations, expression flags, and provenance. |
| `secreted_proteins_predicted.csv` | One row per predicted client gene | The 2,269 predicted cargo proteins, kept separate from machinery. |
| `high_secretion_responsive_genes.csv` | One row per responsive machinery gene | The 51 genes significant in all three Liu comparisons. |
| `glycosylation_genes.csv` | One row per machinery gene with an assigned glycosylation role | Conservative assembly, transfer, trimming, Golgi, O-glycosylation, and GPI annotations. |
| `gene_pathway_roles.csv` | One row per gene–subsystem relationship | Primary placements plus documented secondary roles and cross-links. |
| `gene_id_mapping.csv` | One row per gene–identifier mapping | Current locus, UniProt, KEGG, and NCBI identifier links, including alternate candidates. |
| `pathway_diagram_data.csv` | One row per machinery gene | Primary visualization placement and secondary-subsystem labels. |
| `pathway_overview.svg` | One diagram per build | Standalone vector overview with definitions, counts, citation, confidence note, and build date. |
| `pathway_overview.png` | One diagram per build | GitHub-friendly raster rendering of the same overview. |
| `column_descriptions.csv` | One row per output field | Data dictionary generated from `atlas/schema.py`, including vocabulary and provenance. |
| `FINDINGS.md` | One interpretation per build | Human-readable response pattern, candidates, complex checks, absent responses, and limitations. |

### Quality assurance (`data/processed/qa/`)

| File | Grain | Contents |
|---|---|---|
| `build_audit.json` | One record per build | Machine-readable counts, coverage, and important limitations. |
| `coverage_summary.md` | One summary per build | Human-readable release coverage, evidence order, and remaining gaps. |
| `genes_needing_review.csv` | One row per reviewed source record | All unresolved identifiers and conflicting annotations requiring review. |
| `subsystem_assignment_audit.csv` | One row per machinery gene | Primary evidence plus rejected subsystem candidates not retained in the component table. |
| `liu_vs_atlas_sample.csv` | One row per sampled machinery gene | Twenty original Liu rows beside their derived atlas fields for inspection. |

`unresolved_identifiers.csv` was consolidated into `genes_needing_review.csv`: its three records were a strict subset of the nine-row review queue.

`subsystem_assignment_audit.csv` is retained because six rows record rejected or conflicting subsystem candidates and evidence summaries that are not present in the component table.

## What a row looks like

ALG2 (`AO090120000461`) is significant in all three Liu comparisons and illustrates the final table without implying equal coverage for every field:

| Field | Value |
|---|---|
| `liu_ao_locus_tag` | `AO090120000461` |
| `liu_table_row` | `S1:10` |
| `yeast_ortholog` | `ALG2` |
| `subsystem` | `dolichol_pathway` |
| `glycosylation_role` | `n_glycan_assembly` |
| `direction_all_three` | `up` |
| `uniprot_accession` | `Q2U5Y1` |
| `kegg_gene_id` | `aor:AO090120000461` |
| `kegg_ko` | `ko:K03843` |
| `ncbi_gene_id` | `5996426` |

## Provenance at a glance

Liu-derived, fetched, and pipeline-generated fields remain distinguishable per row. See the `source` column in [`column_descriptions.csv`](data/processed/column_descriptions.csv) and the table below rather than relying on column-name conventions.

## Data dictionary

<!-- DATA_DICTIONARY_START -->
| Field | Definition | Allowed values | Source |
|---|---|---|---|
| record_id | Unique row ID generated by this project. | free text (AOR/AOC plus five digits) | generated |
| ao_locus_tag | Current A. oryzae RIB40 locus tag. | free text (AO090 plus nine digits) or blank | fetched |
| liu_ao_locus_tag | Locus tag exactly as printed by Liu et al. 2014. | free text | Liu 2014 |
| liu_source_raw | The four Liu SOURCE cells in order, joined with pipes without changing cell text; blank for Table S3, which has no SOURCE columns. | free text or blank | Liu 2014 |
| liu_table_row | Supplementary table and physical workbook row. | free text (S1:n or S3:n) | Liu 2014 |
| gene_name | Current accepted A. oryzae gene symbol where available. | free text or blank | fetched |
| function | Plain-English function selected from the Liu description or current annotation. | free text or blank | generated |
| record_type | Separates secretory machinery from predicted cargo. | machinery_component &#124; secretory_client | generated |
| subsystem | Primary pathway subsystem used for diagram placement. | tc &#124; dolichol_pathway &#124; er_glycosylation &#124; folding &#124; gpi_biosynthesis &#124; erad &#124; copii &#124; copi &#124; golgi_processing &#124; ldsv &#124; hdsv &#124; cpy_pathway &#124; alp_pathway &#124; snare &#124; septin &#124; beta_1_6_glucan_biosynthesis &#124; translation &#124; putative_mitochondria_protein &#124; mitochondrial_m_aaa_protease &#124; blank | generated |
| subsystem_source | Evidence route used for the primary subsystem. | liu2014 &#124; liu2014_description &#124; yeast_scaffold &#124; current_annotation &#124; kegg_pathway &#124; curated_liu_composite &#124; liu2014_description+uniprot+kegg_ko &#124; aspergillus_homolog+uniprot+yeast_ortholog+kegg_ko &#124; unassigned | generated |
| subsystem_confidence | Confidence in the primary subsystem assignment. | high &#124; medium &#124; low | generated |
| subsystem_rationale | Short explanation for the primary subsystem assignment. | free text | generated |
| pathway_order | Numeric diagram position derived from the primary subsystem. | numeric or blank | generated |
| compartment | Controlled primary cellular compartment. | cytosol &#124; ER &#124; Golgi &#124; vesicle &#124; membrane &#124; vacuole &#124; extracellular &#124; unknown | generated |
| compartment_source | Evidence route used for compartment assignment. | uniprot &#124; uniprot_conflict &#124; subsystem_inference &#124; unassigned | generated |
| compartment_confidence | Confidence in compartment assignment. | high &#124; medium &#124; low | generated |
| glycosylation_role | Specific glycosylation process assigned to the gene. | n_glycan_assembly &#124; n_glycan_transfer &#124; n_glycan_trimming &#124; golgi_mannosylation &#124; o_glycosylation &#124; gpi_anchor &#124; none | generated |
| glycosylation_role_source | Evidence route used for glycosylation role. | liu2014 &#124; liu2014_description &#124; yeast_scaffold &#124; current_annotation &#124; kegg_pathway &#124; kegg_ko &#124; liu_description_or_current_annotation &#124; curated_multi_source_evidence &#124; blank | generated |
| glycosylation_role_confidence | Confidence in glycosylation role. | high &#124; medium &#124; low &#124; blank | generated |
| sig_all_three | True when all three Liu adjusted p-values are below 0.05. | true &#124; false | Liu 2014 |
| direction_all_three | Shared log-fold-change direction when significant in all three comparisons. | up &#124; down &#124; blank | Liu 2014 |
| record_origin | Origin of the biological record. | liu2014 &#124; post_2014_literature &#124; new_orthology_inference | generated |
| evidence_source | Normalized strength/type of functional evidence. | a_oryzae_experimental &#124; a_oryzae_transcriptomic &#124; aspergillus_homolog &#124; yeast_inference &#124; database_prediction &#124; unknown | generated |
| yeast_ortholog | S. cerevisiae ortholog exactly as recorded by Liu. | free text or blank | Liu 2014 |
| citation | Publication supporting the record. | free text (DOI or PMID) or blank | generated |
| mapping_status | Outcome of current-identifier resolution. | exact &#124; cross_reference &#124; sequence &#124; orthology &#124; ambiguous &#124; split &#124; merged &#124; unresolved | generated |
| mapping_method | Controlled explanation of identifier resolution. | direct AO090 locus-tag match &#124; direct AO090 match in UniProt cross-reference &#124; AO090 tag has multiple UniProt candidates; no candidate selected &#124; no current locus tag found; try step 2 &#124; no current AO090 cross-reference found | generated |
| manual_review_required | Whether automated evidence requires human review. | true &#124; false | generated |
| uniprot_accession | Current UniProt accession. | free text or blank | fetched |
| kegg_gene_id | KEGG gene identifier. | free text (aor:AO090...) or blank | fetched |
| kegg_ko | KEGG Orthology identifier supplied by KEGG. | free text (ko:Knnnnn) or blank | fetched |
| ncbi_gene_id | NCBI Gene identifier. | free text or blank | fetched |
<!-- DATA_DICTIONARY_END -->

## Verification

Row-level tests join all 369 machinery records back to Liu's retained source cells and verify locus tags, yeast orthologs, direct subsystem normalization, functions, significance, direction, shortlist membership, and composite-label handling. Run:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pytest -q
```

Then execute `notebooks/01_profile_liu2014.ipynb`, `02_fetch_and_crosswalk.ipynb`, and `03_build_component_table.ipynb` in order.

## Licensing

- KEGG data is fetched at runtime and bulk KEGG responses/pathway maps are not redistributed, in accordance with KEGG licensing restrictions. Missing KO values remain blank.
- UniProt-derived annotations are distributed under UniProt's CC BY 4.0 terms.
- NCBI data is public domain.
- Liu and Feizi supplementary source files retain their original publication terms; cite the source publications when reusing them.

## Citation

Primary source: Liu et al. (2014), *Systems-level analysis of the secretion stress response in Aspergillus oryzae*, DOI [`10.1186/1752-0509-8-73`](https://doi.org/10.1186/1752-0509-8-73).

Yeast scaffold: Feizi et al. (2013), *Genome-Scale Modeling of the Protein Secretory Machinery in Yeast*, DOI [`10.1371/journal.pone.0063284`](https://doi.org/10.1371/journal.pone.0063284).
