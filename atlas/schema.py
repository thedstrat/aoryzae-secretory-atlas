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

3. Machinery components and secretory clients are DIFFERENT POPULATIONS.
   Machinery = the shipping department. Clients = the cargo. Blending them
   into one unlabelled table produces something actively misleading. They get
   separate cleaned files and a `record_type` field on any combined export.

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
    """Machinery vs cargo. Never blend these silently."""

    MACHINERY = "machinery_component"
    CLIENT = "secretory_client"


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
    # names / function
    "gene_name",
    "function",
    # biology
    "record_type",
    "subsystem",
    "subsystem_source",
    "subsystem_confidence",
    "subsystem_rationale",
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


# One-line definition per field. Mirrored into the README so the data
# dictionary and the code cannot drift apart.
DATA_DICTIONARY = {
    "record_id":            "Unique row ID generated by this project. Database primary key.",
    "ao_locus_tag":         "Current A. oryzae RIB40 locus tag (AO090...). Primary biological identifier.",
    "liu_ao_locus_tag":     "Locus tag exactly as printed in Liu et al. 2014, before resolution.",
    "gene_name":            "Accepted A. oryzae gene symbol where one exists. Often absent.",
    "function":             "Plain-English description of what the gene product does.",
    "record_type":          "machinery_component or secretory_client. Never blend these.",
    "subsystem":            "Liu's functional category, e.g. translocation, folding, erad_degradation.",
    "subsystem_source":     "Provenance for subsystem assignment: liu2014, liu2014_description, yeast_scaffold, current_annotation, kegg_pathway, or unassigned.",
    "subsystem_confidence": "Confidence in the primary subsystem assignment: high, medium, or low.",
    "subsystem_rationale":  "Short plain-English explanation for the primary subsystem assignment.",
    "pathway_order":        "Numeric position along the secretory pathway. Drives diagram layout.",
    "compartment":          "Controlled subcellular location: cytosol, ER, Golgi, vesicle, membrane, vacuole, extracellular, or unknown.",
    "compartment_source":   "Evidence used for the controlled compartment assignment, such as UniProt or subsystem inference.",
    "compartment_confidence": "Confidence in the compartment assignment: high, medium, or low.",
    "glycosylation_role":   "Which glycosylation process, if any. Distinguishes N-assembly/transfer/trimming, Golgi mannosylation, O-glycosylation, GPI.",
    "glycosylation_role_source": "Provenance for the glycosylation-role assignment, such as Liu 2014, KEGG orthology, or literature.",
    "glycosylation_role_confidence": "Confidence in the glycosylation-role assignment: high, medium, or low; blank when no role was assigned.",
    "sig_all_three":        "True if significantly changed in all three alpha-amylase overproducing strains (Liu 2014).",
    "direction_all_three":  "Direction of that change: up or down.",
    "record_origin":        "liu2014, post_2014_literature, or new_orthology_inference. Keeps our 2026 inferences distinct from Liu's list.",
    "evidence_source":      "Strength of evidence: A. oryzae experimental / transcriptomic / Aspergillus homolog / yeast inference / database prediction.",
    "yeast_ortholog":       "S. cerevisiae systematic name. Stable anchor - yeast IDs don't rot.",
    "citation":             "PMID or DOI supporting this record.",
    "mapping_status":       "exact, cross_reference, sequence, orthology, ambiguous, split, merged, or unresolved.",
    "mapping_method":       "Free-text note on how the identifier was resolved.",
    "manual_review_required": "True where automated assignment needs a human check.",
    "uniprot_accession":    "UniProt accession. Multiple accessions produce multiple rows, not a delimited cell.",
    "kegg_gene_id":         "KEGG gene ID, aor:AO090... form.",
    "kegg_ko":              "KEGG Orthology group. Use this to find gene families - symbols are mostly absent in A. oryzae.",
    "ncbi_gene_id":         "NCBI Gene identifier.",
}
