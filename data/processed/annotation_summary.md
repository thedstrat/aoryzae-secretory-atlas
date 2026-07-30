# Annotation completion summary

Liu supplied subsystem labels for 109 of 369 machinery genes. The 260 blanks still had yeast orthologs; most were additions from A. niger, A. oryzae, reciprocal-best-hit, or InParanoid source lists rather than rows in the original 16-subsystem yeast scaffold. They were therefore unlabeled, not lost.

This release places **247/369** genes. Assignments preserve their source in this order: Liu's explicit label, an unambiguous Liu description, the Feizi/Liu yeast scaffold, current UniProt annotation text, then KEGG pathway membership. Conflicts remain visible rather than being silently guessed.

KEGG retrieval now validates the organism through `GET https://rest.kegg.jp/list/aor`. **330/369** machinery genes have a supplied KO; unavailable values remain blank. KO family mappings add conservative glycosylation roles without replacing Liu-derived roles.

KAR2/BiP is placed primarily in folding, with ERAD retained as a secondary relationship. Seven other genes with explicit cross-subsystem Liu descriptions retain their defensible primary assignment and now expose the additional role in `gene_subsystems.csv`. The pathway figure shows all 122 unassigned genes.
