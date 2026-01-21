from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import List, Optional


def my_dataclass(cls=None, **kwargs):
    if sys.version_info >= (3, 10):
        kwargs["slots"] = True
    return dataclass(cls, **kwargs)


@my_dataclass
class GafStep:
    name: str = ""
    is_reverse: bool = False
    is_stable: bool = False
    is_interval: bool = False
    start: Optional[int] = None
    end: Optional[int] = None


@my_dataclass
class GafRecord:
    query_name: str = ""
    query_length: int = 0
    query_start: int = 0
    query_end: int = 0
    path_length: int = 0
    path_start: int = 0
    path_end: int = 0
    matches: int = 0
    block_length: int = 0
    mapq: int = 0
    strand: str = ""
    path: List[GafStep] = field(default_factory=list)
    opt_fields: dict[str, tuple[str, str]] = field(default_factory=dict)


@my_dataclass
class Edit:
    """
    Edits describe how to generate a new string from elements in the graph. To
    determine the new string, just walk the series of edits, stepping
    from_length distance in the basis node, and to_length in the novel element,
    replacing from_length in the basis node with the sequence. There are
    several types of Edit: - *matches*: from_length == to_length; sequence is
    empty - *snps*: from_length == to_length; sequence = alt - *deletions*:
    to_length == 0 && from_length > to_length; sequence is empty -
    *insertions*: from_length < to_length; sequence = alt
    """

    from_length: int = 0
    to_length: int = 0
    sequence: str = ""


@my_dataclass
class Position:
    node_id: int = 0
    offset: int = 0
    is_reverse: bool = False
    name: str = ""


@my_dataclass
class Mapping:
    """
    A Mapping defines the relationship between a node in system and another
    entity. An empty edit list implies complete match, however it is preferred
    to specify the full edit structure. as it is more complex to handle special
    cases.
    """

    position: Position = Position()
    edit: List[Edit] = field(default_factory=list)
    rank: int = 0


@my_dataclass
class Path:
    """
    Paths are walks through nodes defined by a series of `Edit`s. They can be
    used to represent:    - haplotypes    - mappings of reads, or alignments,
    by including edits    - relationships between nodes    - annotations from
    other data sources, such as:          genes, exons, motifs, transcripts,
    peaks
    """

    name: str = ""
    mapping: List[Mapping] = field(default_factory=list)
    is_circular: bool = False
    length: int = 0


@my_dataclass
class Alignment:
    """
    Alignments link query strings, such as other genomes or reads, to Paths.
    """

    sequence: str = ""
    path: Path = Path()
    name: str = ""
    quality: bytes = b""
    mapping_quality: int = 0
    score: int = 0
    query_position: int = 0
    sample_name: str = ""
    read_group: str = ""
    is_secondary: bool = False
    identity: float = 0.0
    fragment_prev: Optional[Alignment] = None
    fragment_next: Optional[Alignment] = None
    fragment: List[Path] = field(default_factory=list)
    locus: List[Locus] = field(default_factory=list)
    refpos: List[Position] = field(default_factory=list)
    read_paired: bool = False
    """SAMTools-style flags"""

    read_mapped: bool = False
    mate_unmapped: bool = False
    read_on_reverse_strand: bool = False
    mate_on_reverse_strand: bool = False
    soft_clipped: bool = False
    discordant_insert_size: bool = False
    uniqueness: float = 0.0
    correct: float = 0.0
    secondary_score: List[int] = field(default_factory=list)
    fragment_score: float = 0.0
    mate_mapped_to_disjoint_subgraph: bool = False
    fragment_length_distribution: str = ""
    time_used: float = 0.0
    to_correct: Position = Position()
    correctly_mapped: bool = False
    annotation: dict = field(default_factory=dict)


@my_dataclass
class MultipathAlignment:
    """
    A subgraph of the unrolled Graph in which each non-branching path is
    associated with an alignment of part of the read and part of the graph such
    that any path through the MultipathAlignment indicates a valid alignment of
    a read to the graph
    """

    sequence: str = ""
    quality: bytes = b""
    name: str = ""
    sample_name: str = ""
    read_group: str = ""
    subpath: List[Subpath] = field(default_factory=list)
    """
    non-branching paths of the multipath alignment, each containing an
    alignment of part of the sequence to a Graph IMPORTANT: downstream
    applications will assume these are stored in topological order
    """

    mapping_quality: int = 0
    """-10 * log_10(probability of mismapping)"""

    start: List[int] = field(default_factory=list)
    """
    optional: indices of Subpaths that align the beginning of the read (i.e.
    source nodes)
    """

    paired_read_name: str = ""
    annotation: dict = field(default_factory=dict)


@my_dataclass
class Subpath:
    """A non-branching path of a MultipathAlignment"""

    path: Path = Path()
    """describes node sequence and edits to the graph sequences"""

    next: List[int] = field(default_factory=list)
    """
    the indices of subpaths in the multipath alignment that are to the right of
    this path where right is in the direction of the end of the read sequence
    """

    score: int = 0
    """score of this subpath's alignment"""

    connection: List[Connection] = field(default_factory=list)
    """
    connections to other subpaths that are not necessarily contiguous in the
    graph
    """


@my_dataclass
class Connection:
    """
    An edge in a MultipathAlignment between Subpaths that may not be contiguous
    in the graph
    """

    next: int = 0
    """the index of the Subpath that this connection points to"""

    score: int = 0
    """the score of this connection"""


@my_dataclass
class Support:
    """Aggregates information about the reads supporting an allele."""

    quality: float = 0.0
    """
    The overall quality of all the support, as -10 * log10(P(all support is
    wrong))
    """

    forward: float = 0.0
    """
    The number of supporting reads on the forward strand (which may be
    fractional)
    """

    reverse: float = 0.0
    """
    The number of supporting reads on the reverse strand (which may be
    fractional)
    """

    left: float = 0.0
    """TODO: what is this?"""

    right: float = 0.0
    """TODO: What is this?"""


@my_dataclass
class Locus:
    """
    Describes a genetic locus with multiple possible alleles, a genotype, and
    observational support.
    """

    name: str = ""
    """A locus may have an identifying name."""

    allele: List[Path] = field(default_factory=list)
    """
    These are all the alleles at the locus, not just the called ones. Note that
    a primary reference allele may or may not appear.
    """

    support: List[Support] = field(default_factory=list)
    """These supports are per-allele, matching the alleles above"""

    genotype: List[Genotype] = field(default_factory=list)
    """sorted by likelihood or posterior  the first one is the "call"""

    overall_support: Support = Support()
    """
    We also have a Support for the locus overall, because reads may have
    supported multiple alleles and we want to know how many total there were.
    """

    allele_log_likelihood: List[float] = field(default_factory=list)
    """
    We track the likelihood of each allele individually, in addition to
    genotype likelihoods. Stores the likelihood natural logged.
    """


@my_dataclass
class Genotype:
    """Describes a genotype at a particular locus."""

    allele: List[int] = field(default_factory=list)
    """These refer to the offsets of the alleles in the Locus object."""

    is_phased: bool = False
    likelihood: float = 0.0
    log_likelihood: float = 0.0
    log_prior: float = 0.0
    log_posterior: float = 0.0
