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


def _field(
    definition: str,
    allowed_values: str,
    source: str,
    value_meanings: str = "",
    blank_means: str = "This field is not expected to be blank.",
) -> dict[str, str]:
    """One machine-readable data-dictionary entry."""
    return {
        "definition": definition,
        "allowed_values": allowed_values,
        "value_meanings": value_meanings,
        "blank_means": blank_means,
        "source": source,
    }


DATA_DICTIONARY = {
    "record_id": _field("Unique row ID generated by this project.", "free text (AOR/AOC plus five digits)", "generated"),
    "ao_locus_tag": _field("Current A. oryzae RIB40 locus tag selected by the ID-mapping step.", "free text (AO090 plus nine digits) or blank", "fetched", blank_means="No single current locus was independently verified."),
    "liu_ao_locus_tag": _field("Locus tag exactly as printed by Liu et al. 2014.", "free text", "Liu 2014"),
    "liu_source_raw": _field("The four Liu Table S1 SOURCE cells in order, joined with pipes without changing cell text.", "free text", "Liu 2014"),
    "liu_table_row": _field("Table S1 physical workbook row.", "free text (S1:n)", "Liu 2014"),
    "gene_name": _field("Current accepted A. oryzae gene symbol where available.", "free text or blank", "fetched", blank_means="No accepted short gene symbol was found; use ao_locus_tag instead."),
    "function": _field("Plain-English function selected from the Liu description or current annotation.", "free text or blank", "generated", blank_means="Neither Liu nor the current annotation supplied a usable description."),
    "record_type": _field("Row type retained so machinery genes cannot be silently mixed with other biological lists.", "machinery_component", "generated", "machinery_component = one gene from Liu's secretion-machinery list"),
    "subsystem": _field("Primary secretion-pathway stage used for table grouping and diagram placement.", " | ".join([*SUBSYSTEM_ORDER, "blank"]), "generated", "tc = entry into the secretion pathway; dolichol_pathway and er_glycosylation = N-glycan assembly/transfer; folding = protein folding; gpi_biosynthesis = GPI-anchor production; erad = ER quality-control degradation; copii = ER-to-Golgi transport; copi = Golgi-to-ER transport; golgi_processing = Golgi glycan processing; ldsv/hdsv = low/high-density secretory vesicles; cpy_pathway/alp_pathway = vacuolar sorting routes; snare = vesicle fusion; remaining values name specialized processes", "Not enough evidence to choose one primary pathway stage; the gene remains in the atlas."),
    "subsystem_source": _field("Main evidence route used to assign subsystem; see subsystem_rationale for the gene-specific explanation.", "liu2014 | liu2014_description | yeast_scaffold | current_annotation | kegg_pathway | curated_liu_composite | aspergillus_homolog | unassigned", "generated", "liu2014 = Liu supplied the stage directly; liu2014_description = inferred from Liu's description; yeast_scaffold = inferred through the yeast counterpart; current_annotation = inferred from a current database description; kegg_pathway = inferred from broad KEGG pathway membership; curated_liu_composite = manually reviewed agreement among several Liu clues; aspergillus_homolog = inferred from a related Aspergillus gene; unassigned = no route was strong enough"),
    "subsystem_confidence": _field("Strength of support for the primary subsystem assignment, not confidence in the gene's general function or ID mapping.", "high | medium | low", "generated", "high = stated directly by Liu or manually reviewed agreement among multiple specific sources; medium = one specific indirect source such as Liu's description, a unique yeast match, or a current annotation; low = only broad KEGG evidence or no defensible assignment"),
    "subsystem_rationale": _field("Short explanation for the primary subsystem assignment.", "free text", "generated"),
    "secondary_subsystem": _field("Reviewed additional pathway role retained on the same one-gene row.", "subsystem value or blank", "generated", blank_means="No separate secondary role was reviewed and retained; this does not prove the gene has only one function."),
    "secondary_subsystem_rationale": _field("Short explanation for the reviewed secondary pathway role.", "free text or blank", "generated", blank_means="No secondary_subsystem was assigned."),
    "pathway_order": _field("Number used only to place subsystem boxes in the pathway diagram.", "numeric or blank", "generated", "Smaller numbers appear earlier in the diagram; the number is not a confidence score, rank, or biological priority.", "No primary subsystem was assigned, so there is no diagram position."),
    "compartment": _field("Controlled label for the gene product's main cellular location.", "cytosol | ER | Golgi | vesicle | membrane | vacuole | extracellular | unknown", "generated", "unknown = available evidence did not support one location; other values name the assigned location"),
    "compartment_source": _field("Evidence route used for the compartment label.", "uniprot | uniprot_conflict | subsystem_inference | unassigned", "generated", "uniprot = one clear UniProt location; uniprot_conflict = UniProt listed conflicting or multiple locations, so compartment is unknown; subsystem_inference = inferred from the assigned pathway stage; unassigned = no usable location evidence"),
    "compartment_confidence": _field("Strength of support for the compartment label only.", "high | medium | low", "generated", "high = one clear UniProt location; medium = inferred from a pathway stage; low = conflicting or absent location evidence"),
    "glycosylation_role": _field("Specific sugar-modification process assigned to the gene.", " | ".join(role.value for role in GlycosylationRole), "generated", "n_glycan_assembly = build the precursor N-glycan; n_glycan_transfer = attach it to a protein; n_glycan_trimming = remove sugars during processing/quality control; golgi_mannosylation = add mannose in the Golgi; o_glycosylation = attach sugars through oxygen; gpi_anchor = build a GPI membrane anchor; none = no role assigned by this atlas"),
    "glycosylation_role_source": _field("Evidence route used for glycosylation_role.", "liu2014 | liu2014_description | yeast_scaffold | current_annotation | kegg_pathway | kegg_ko | liu_description_or_current_annotation | curated_multi_source_evidence | blank", "generated", "liu2014 = direct Liu pathway label; liu2014_description = Liu description; yeast_scaffold = matching yeast gene; current_annotation = current database description; kegg_pathway = broad KEGG pathway; kegg_ko = KEGG functional group; liu_description_or_current_annotation = one specific matching description; curated_multi_source_evidence = manually reviewed agreement among multiple sources", "Blank because glycosylation_role is none; no source is claimed for an unassigned role."),
    "glycosylation_role_confidence": _field("Strength of support for glycosylation_role only, not general function confidence.", "high | medium | low | blank", "generated", "high = direct Liu assignment or manually reviewed agreement among several specific sources; medium = one specific indirect route such as a description, yeast match, KEGG pathway, or KEGG functional group; low = reserved for weak assignments (none are present in this release)", "Blank because glycosylation_role is none; an unassigned role receives no confidence score."),
    "sig_all_three": _field("Whether all three Liu adjusted p-values are below 0.05.", "true | false", "Liu 2014", "true = statistically significant in all three high-secretion strains; false = at least one comparison was not significant"),
    "direction_all_three": _field("Shared expression direction when the gene is significant in all three comparisons.", "up | down | blank", "Liu 2014", "up = expression increased in all three strains; down = expression decreased in all three strains", "The gene was not significant in all three comparisons, or the three directions did not agree."),
    "record_origin": _field("Publication or process that placed the gene row in the atlas.", " | ".join(value.value for value in RecordOrigin), "generated", "liu2014 = gene was in Liu's 369-gene list; post_2014_literature and new_orthology_inference are reserved for possible future expansion"),
    "evidence_source": _field("Broad type of evidence Liu used to include or support the gene; this is not the subsystem confidence score.", " | ".join(value.value for value in EvidenceSource), "generated", "a_oryzae_experimental = direct experiment in A. oryzae; a_oryzae_transcriptomic = A. oryzae expression evidence; aspergillus_homolog = evidence from a related Aspergillus species; yeast_inference = inferred through a yeast counterpart; database_prediction = computational database prediction; unknown = source text could not be classified"),
    "yeast_ortholog": _field("S. cerevisiae counterpart exactly as recorded by Liu.", "free text or blank", "Liu 2014", blank_means="Liu did not provide a yeast counterpart for this gene."),
    "citation": _field("Liu 2014 DOI supporting every machinery record.", "10.1186/1752-0509-8-73", "Liu 2014"),
    "mapping_status": _field("Result of checking the 2014 gene ID against current database records.", " | ".join(value.value for value in MappingStatus), "generated", "exact = the same AO090 ID is current; cross_reference = linked by an authoritative database cross-reference; sequence = linked by sequence identity; orthology = only an evolutionary counterpart supports the link, not exact identity; ambiguous = more than one current candidate; split = one old gene became several models; merged = several old genes became one model; unresolved = no current gene was verified"),
    "mapping_method": _field("Specific route used to reach mapping_status.", "direct_locus_tag | cross_reference | sequence | orthology | ambiguous | split | merged | unresolved", "generated", "direct_locus_tag = exact AO090 text match; cross_reference = authoritative database link; sequence = sequence identity; orthology = evolutionary relationship only; ambiguous/split/merged/unresolved mirror the unresolved mapping outcome"),
    "manual_review_required": _field("Whether any automated mapping or annotation decision needs a person to inspect it.", "true | false", "generated", "true = at least one mapping or annotation issue needs review; false = no review trigger was raised, not that every biological claim is experimentally proven"),
    "uniprot_accession": _field("Current UniProt protein-record identifier.", "free text or blank", "fetched", blank_means="No single UniProt record was found or selected; alternatives may remain in gene_id_mapping.csv."),
    "kegg_gene_id": _field("Current KEGG A. oryzae gene identifier.", "free text (aor:AO090...) or blank", "fetched", blank_means="No KEGG gene link was supplied or selected for this row."),
    "kegg_ko": _field("KEGG Orthology functional-group identifier supplied by KEGG.", "free text (ko:Knnnnn) or blank", "fetched", "Genes sharing a KO are treated by KEGG as performing an equivalent broad function across species.", "KEGG did not supply a KO link; this is incomplete annotation, not evidence that the gene has no function."),
    "ncbi_gene_id": _field("Current NCBI Gene numeric identifier.", "free text or blank", "fetched", blank_means="No unambiguous NCBI Gene link was available from the retrieved cross-references."),
}
