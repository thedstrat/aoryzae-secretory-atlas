# Annotation coverage summary

Liu supplied subsystem labels for 109 of 369 machinery genes. The 260 blanks still had yeast orthologs; most were additions from A. niger, A. oryzae, reciprocal-best-hit, or InParanoid source lists rather than rows in the original 16-subsystem yeast scaffold. They were therefore unlabeled, not lost.

This release places **247/369** genes. Assignments preserve their source in this order: Liu's explicit label, an unambiguous Liu description, the Feizi/Liu yeast scaffold, current UniProt annotation text, then KEGG pathway membership. Conflicts remain visible rather than being silently guessed.

KEGG retrieval now validates the organism through `GET https://rest.kegg.jp/list/aor`. **330/369** machinery genes have a supplied KO; unavailable values remain blank. KO family mappings add conservative glycosylation roles without replacing Liu-derived roles.

Four reviewed secondary roles are retained directly in `secretion_machinery_genes.csv`. The pathway figure is built from that table and shows all 122 unassigned genes.
