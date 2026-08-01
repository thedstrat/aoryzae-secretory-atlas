# Independent crosswalk validation

## Result

The audit independently checked all 369 Liu 2014 identifiers against freshly retrieved current records. It did not use gene names or functional similarity to establish identity.

| Audit status | Count |
|---|---:|
| Confirmed | 368 |
| Probable | 0 |
| Ambiguous | 0 |
| Conflicting | 0 |
| Unresolved | 1 |

| Mapping method | Count |
|---|---:|
| Unchanged/exact identifier | 368 |
| Authoritative database cross-reference | 0 |
| Sequence identity | 0 |
| Orthology only | 0 |
| Ambiguous | 0 |
| Split | 0 |
| Merged | 0 |
| Unresolved | 1 |

Of the 368 confirmed loci, 367 were supported by NCBI RefSeq, KEGG, and UniProt. `AO090673000002` was supported by NCBI RefSeq and KEGG but had no matching current UniProt record. No current database identifiers disagreed, no old identifier mapped to multiple current genes, and no multiple old identifiers mapped to one verified current gene.

## Proposed corrections before any existing mapping changes

No proposed AO090 locus needs to be replaced with a different locus. One existing classification should be corrected:

| Original ID | Existing result | Proposed result | Evidence |
|---|---|---|---|
| `AO090011000795` | Ambiguous; current AO090 locus retained but database IDs withheld | Confirmed unchanged locus `AO090011000795`; NCBI GeneID `5998476`; KEGG `aor:AO090011000795`; UniProt `Q2TZM5` | NCBI RefSeq GFF, KEGG, and UniProt agree on the locus. UniProt `Q9Y8E3` uses the separate legacy name `AO090011000795-A`; it is preserved as an alternative record but is not a second current AO090 gene. |

This report and `crosswalk_audit.csv` show the correction first. The existing mapping in `secretion_machinery_genes.csv` has **not** been changed.

## Unresolved identifier

`AO090005001666` is absent from the current NCBI RefSeq annotation, KEGG `aor` gene list, and UniProt organism-wide export. Liu associates it with the yeast ortholog `SEC20`, but orthology or function alone is not enough to assign an exact *A. oryzae* identity. No archived authoritative source sequence was established in this audit, so no sequence comparison or replacement locus is claimed. It remains unresolved and requires manual review.

## Alternative candidates retained

The only additional identifier encountered was UniProt `Q9Y8E3`, named `AO090011000795-A`. It is recorded in `alternative_candidates` for `AO090011000795`. It was not promoted to a current gene candidate because the suffix-bearing name is not a current AO090 locus in NCBI or KEGG.

## Method

For every Liu identifier, the audit compared exact locus identifiers and database cross-references in:

- NCBI RefSeq assembly `GCF_000184455.2` (`ASM18445v3`) genomic GFF;
- the KEGG `aor` organism gene list and its NCBI GeneID and UniProt conversions;
- the UniProt organism-wide export for taxonomy ID `510516`, including KEGG and GeneID cross-references.

The records were retrieved on 2026-08-01 from:

- `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/455/GCF_000184455.2_ASM18445v3/GCF_000184455.2_ASM18445v3_genomic.gff.gz`
- `https://rest.kegg.jp/list/aor`
- `https://rest.kegg.jp/conv/ncbi-geneid/aor`
- `https://rest.kegg.jp/conv/uniprot/aor`
- `https://rest.uniprot.org/uniprotkb/stream` with query `organism_id:510516`

A row is **confirmed** when at least two independent current sources agree on one locus, **probable** when only one current source supports it, **ambiguous** when more than one current candidate remains, **conflicting** when sources disagree, and **unresolved** when no current source supports a candidate. Supporting accessions and GeneIDs are retained in `audit_notes`; alternative candidates are retained separately.

## Limits

This validates identifier continuity and cross-references, not the biological function assigned to each gene. Exact sequence identity was not needed for the 368 unchanged loci because multiple authoritative current sources retained the original identifiers. The unresolved row could not be evaluated by sequence because Liu Table S1 supplies an identifier and ortholog, not the legacy protein sequence.
