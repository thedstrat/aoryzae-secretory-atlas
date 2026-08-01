"""Schema and data dictionary for the A. oryzae secretory/glycosylation atlas.

DESIGN DECISIONS (the reasoning, so future-you doesn't relitigate it)
---------------------------------------------------------------------

1. `record_id` is the database primary key; `ao_locus_tag` is the primary
   *biological* identifier. These are different jobs. A surrogate key survives
   the cases that break a natural key: an unresolved source row, one old tag
   splitting into two current records, a yeast component with no A. oryzae
   match yet.

2. AO090... locus tags are alive. NCBI Gene records for RIB40 still use them
   (gene 5999297 is AO090010000120), on RefSeq GCF_000184455.2, and KEGG keys
   on aor:AO090.... AspGD being unmaintained did NOT kill the identifiers.
   Direct locus-tag matching is therefore the primary mapping path.

3. This atlas covers secretion machinery from Liu Table S1. Liu's predicted
   secretome in Table S3 is a different population and remains out of scope.

4. Evidence is a first-class column, not a footnote. Most of this table is
   inference from yeast. Saying so plainly is the difference between a
   resource and a liability.

5. Engineering-hypothesis fields (possible_intervention, engineering_risk,
   ...) are DEFERRED to a later annotation layer. They depend on decisions
   not yet made (target protein, host) and would be empty columns today.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Optional


# --------------------------------------------------------------- vocabularies

class RecordType(str, Enum):
    """Biological population represented by this machinery atlas."""

    MACHINERY = "machinery_component"


class RecordOrigin(str, Enum):
    """Where this row came from. Keeps canonical and new inference distinct."""

    LIU2014 = "liu2014"
    POST_2014_LIT = "post_2014_literature"
    NEW_ORTHOLOGY = "new_orthology_inference"   # our 2026 re-derivation


class MappingStatus(str, Enum):
    """Did the Liu-era identifier resolve, and how cleanly?"""

    EXACT = "exact"                     # direct locus-tag hit
    CROSS_REFERENCE = "cross_reference"  # via UniProt/KEGG/NCBI xref
    SEQUENCE = "sequence"               # resolved by protein sequence match
    ORTHOLOGY = "orthology"             # re-derived from yeast
    AMBIGUOUS = "ambiguous"             # one tag -> several current genes
    SPLIT = "split"                     # gene model was split since 2014
    MERGED = "merged"                   # gene model was merged since 2014
    UNRESOLVED = "unresolved"           # no current gene found - REPORT IT


class EvidenceSource(str, Enum):
    """What kind of evidence backs this gene's assignment?

    Ordered strongest to weakest. For A. oryzae most rows will land at
    YEAST_INFERENCE - which is exactly why the column has to exist.
    """

    AO_EXPERIMENTAL = "a_oryzae_experimental"       # measured in A. oryzae
    AO_TRANSCRIPTOMIC = "a_oryzae_transcriptomic"   # Liu's expression data
    ASPERGILLUS_HOMOLOG = "aspergillus_homolog"     # A. niger / A. nidulans
    YEAST_INFERENCE = "yeast_inference"             # S. cerevisiae ortholog
    DATABASE_PREDICTION = "database_prediction"     # computational only
    UNKNOWN = "unknown"


class GlycosylationRole(str, Enum):
    """Deliberately granular - Vikas's item 3b needs these distinguished."""

    N_GLYCAN_ASSEMBLY = "n_glycan_assembly"     # dolichol/ALG, ER-side
    N_GLYCAN_TRANSFER = "n_glycan_transfer"     # OST complex
    N_GLYCAN_TRIMMING = "n_glycan_trimming"     # ER glucosidases/mannosidases
    GOLGI_MANNOSYLATION = "golgi_mannosylation"  # OCH1/MNN/KTR - the hyperman. step
    O_GLYCOSYLATION = "o_glycosylation"          # PMT family
    GPI_ANCHOR = "gpi_anchor"
    NONE = "none"


# Normalized transcription of every non-null category observed in Liu 2014
# Additional file 2, Table S1. The workbook leaves 260 of 369 rows uncategorized
# and also includes three mitochondria/translation labels beyond the pathway
# categories described in the paper. Keep those source facts visible rather
# than silently forcing rows into the paper's canonical 16-subsystem model.
SUBSYSTEM_ORDER: dict[str, int] = {
    "tc": 10,
    "dolichol_pathway": 20,
    "er_glycosylation": 30,
    "folding": 40,
    "gpi_biosynthesis": 50,
    "erad": 60,
    "copii": 70,
    "copi": 80,
    "golgi_processing": 90,
    "ldsv": 100,
    "hdsv": 110,
    "cpy_pathway": 120,
    "alp_pathway": 130,
    "snare": 140,
    "septin": 150,
    "beta_1_6_glucan_biosynthesis": 160,
    "translation": 170,
    "putative_mitochondria_protein": 180,
    "mitochondrial_m_aaa_protease": 190,
}


# ------------------------------------------------------------------- the record

@dataclass
class GeneRecord:
    """One row of the component table.

    ~20 fields. Every one earns its place by feeding either the pathway
    diagram or a decision Vikas has to make. Anything that does neither is
    deferred, not added "just in case".
    """

    # --- keys --------------------------------------------------------------
    record_id: str                                # surrogate PK, project-generated
    ao_locus_tag: Optional[str] = None            # primary biological ID (current)
    liu_ao_locus_tag: Optional[str] = None        # as printed in Liu 2014
    liu_source_raw: Optional[str] = None           # verbatim SOURCE cells, pipe-delimited
    liu_table_row: Optional[str] = None            # source sheet and physical row

    # --- names and function ------------------------------------------------
    gene_name: Optional[str] = None               # A. oryzae symbol, when one exists
    function: Optional[str] = None                # plain-English description

    # --- crosswalk ---------------------------------------------------------
    uniprot_accession: Optional[str] = None
    kegg_gene_id: Optional[str] = None            # aor:AO090...
    kegg_ko: Optional[str] = None                 # orthology group - finds families
    ncbi_gene_id: Optional[str] = None

    # --- biological classification -----------------------------------------
    record_type: RecordType = RecordType.MACHINERY
    subsystem: Optional[str] = None               # one of SUBSYSTEM_ORDER
    subsystem_source: str = "unassigned"          # controlled assignment provenance
    subsystem_confidence: Optional[str] = None
    subsystem_rationale: Optional[str] = None
    secondary_subsystem: Optional[str] = None
    secondary_subsystem_rationale: Optional[str] = None
    pathway_order: Optional[int] = None           # diagram layout position
    compartment: Optional[str] = None             # ER, Golgi, vesicle, membrane...
    compartment_source: Optional[str] = None
    compartment_confidence: Optional[str] = None
    glycosylation_role: GlycosylationRole = GlycosylationRole.NONE
    glycosylation_role_source: Optional[str] = None  # Liu / KEGG KO / literature provenance
    glycosylation_role_confidence: Optional[str] = None

    # --- evidence and provenance -------------------------------------------
    record_origin: RecordOrigin = RecordOrigin.LIU2014
    evidence_source: EvidenceSource = EvidenceSource.UNKNOWN
    yeast_ortholog: Optional[str] = None          # stable anchor (YMR123W style)
    citation: Optional[str] = None                # PMID or DOI
    mapping_status: MappingStatus = MappingStatus.UNRESOLVED
    mapping_method: Optional[str] = None          # free text: how it resolved
    manual_review_required: bool = False

    # --- Liu's transcriptomic evidence -------------------------------------
    # The single most decision-relevant thing in the source dataset: genes the
    # cell itself ramps up when straining to secrete more alpha-amylase.
    # Measured in the actual host. This is the shortlist to hand Vikas.
    sig_all_three: Optional[bool] = None          # significant in all 3 strains
    direction_all_three: Optional[str] = None     # "up" | "down" | None

    def to_row(self) -> dict:
        d = asdict(self)
        for k, v in d.items():
            if hasattr(v, "value"):
                d[k] = v.value
        return d


COLUMN_ORDER = [
    # keys
    "record_id",
    "ao_locus_tag",
    "liu_ao_locus_tag",
    "liu_source_raw",
    "liu_table_row",
    # names / function
    "gene_name",
    "function",
    # biology
    "record_type",
    "subsystem",
    "subsystem_source",
    "subsystem_confidence",
    "subsystem_rationale",
    "secondary_subsystem",
    "secondary_subsystem_rationale",
    "pathway_order",
    "compartment",
    "compartment_source",
    "compartment_confidence",
    "glycosylation_role",
    "glycosylation_role_source",
    "glycosylation_role_confidence",
    # transcriptomic evidence
    "sig_all_three",
    "direction_all_three",
    # provenance
    "record_origin",
    "evidence_source",
    "yeast_ortholog",
    "citation",
    "mapping_status",
    "mapping_method",
    "manual_review_required",
    # crosswalk
    "uniprot_accession",
    "kegg_gene_id",
    "kegg_ko",
    "ncbi_gene_id",
]


# Fields that must never be null in a shippable table.
REQUIRED_FIELDS = ["record_id", "record_type", "record_origin", "mapping_status"]


def _field(definition: str, allowed_values: str, source: str) -> dict[str, str]:
    """One machine-readable data-dictionary entry."""
    return {"definition": definition, "allowed_values": allowed_values, "source": source}


DATA_DICTIONARY = {
    "record_id": _field("Unique row ID generated by this project.", "free text (AOR/AOC plus five digits)", "generated"),
    "ao_locus_tag": _field("Current A. oryzae RIB40 locus tag.", "free text (AO090 plus nine digits) or blank", "fetched"),
    "liu_ao_locus_tag": _field("Locus tag exactly as printed by Liu et al. 2014.", "free text", "Liu 2014"),
    "liu_source_raw": _field("The four Liu Table S1 SOURCE cells in order, joined with pipes without changing cell text.", "free text", "Liu 2014"),
    "liu_table_row": _field("Table S1 physical workbook row.", "free text (S1:n)", "Liu 2014"),
    "gene_name": _field("Current accepted A. oryzae gene symbol where available.", "free text or blank", "fetched"),
    "function": _field("Plain-English function selected from the Liu description or current annotation.", "free text or blank", "generated"),
    "record_type": _field("Identifies rows as secretion-machinery components.", "machinery_component", "generated"),
    "subsystem": _field("Primary pathway subsystem used for diagram placement.", " | ".join([*SUBSYSTEM_ORDER, "blank"]), "generated"),
    "subsystem_source": _field("Evidence route used for the primary subsystem.", "liu2014 | liu2014_description | yeast_scaffold | current_annotation | kegg_pathway | curated_liu_composite | liu2014_description+uniprot+kegg_ko | aspergillus_homolog+uniprot+yeast_ortholog+kegg_ko | unassigned", "generated"),
    "subsystem_confidence": _field("Confidence in the primary subsystem assignment.", "high | medium | low", "generated"),
    "subsystem_rationale": _field("Short explanation for the primary subsystem assignment.", "free text", "generated"),
    "secondary_subsystem": _field("Reviewed secondary pathway role retained without duplicating the gene row.", "subsystem value or blank", "generated"),
    "secondary_subsystem_rationale": _field("Short explanation for the secondary pathway role.", "free text or blank", "generated"),
    "pathway_order": _field("Numeric diagram position derived from the primary subsystem.", "numeric or blank", "generated"),
    "compartment": _field("Controlled primary cellular compartment.", "cytosol | ER | Golgi | vesicle | membrane | vacuole | extracellular | unknown", "generated"),
    "compartment_source": _field("Evidence route used for compartment assignment.", "uniprot | uniprot_conflict | subsystem_inference | unassigned", "generated"),
    "compartment_confidence": _field("Confidence in compartment assignment.", "high | medium | low", "generated"),
    "glycosylation_role": _field("Specific glycosylation process assigned to the gene.", " | ".join(role.value for role in GlycosylationRole), "generated"),
    "glycosylation_role_source": _field("Evidence route used for glycosylation role.", "liu2014 | liu2014_description | yeast_scaffold | current_annotation | kegg_pathway | kegg_ko | liu_description_or_current_annotation | curated_multi_source_evidence | blank", "generated"),
    "glycosylation_role_confidence": _field("Confidence in glycosylation role.", "high | medium | low | blank", "generated"),
    "sig_all_three": _field("True when all three Liu adjusted p-values are below 0.05.", "true | false", "Liu 2014"),
    "direction_all_three": _field("Shared log-fold-change direction when significant in all three comparisons.", "up | down | blank", "Liu 2014"),
    "record_origin": _field("Origin of the biological record.", " | ".join(value.value for value in RecordOrigin), "generated"),
    "evidence_source": _field("Normalized strength/type of functional evidence.", " | ".join(value.value for value in EvidenceSource), "generated"),
    "yeast_ortholog": _field("S. cerevisiae ortholog exactly as recorded by Liu.", "free text or blank", "Liu 2014"),
    "citation": _field("Publication supporting the record.", "free text (DOI or PMID) or blank", "generated"),
    "mapping_status": _field("Outcome of current-identifier resolution.", " | ".join(value.value for value in MappingStatus), "generated"),
    "mapping_method": _field("Controlled explanation of identifier resolution.", "direct AO090 locus-tag match | direct AO090 match in UniProt cross-reference | AO090 tag has multiple UniProt candidates; no candidate selected | no current locus tag found; try step 2 | no current AO090 cross-reference found", "generated"),
    "manual_review_required": _field("Whether automated evidence requires human review.", "true | false", "generated"),
    "uniprot_accession": _field("Current UniProt accession.", "free text or blank", "fetched"),
    "kegg_gene_id": _field("KEGG gene identifier.", "free text (aor:AO090...) or blank", "fetched"),
    "kegg_ko": _field("KEGG Orthology identifier supplied by KEGG.", "free text (ko:Knnnnn) or blank", "fetched"),
    "ncbi_gene_id": _field("NCBI Gene identifier.", "free text or blank", "fetched"),
}
