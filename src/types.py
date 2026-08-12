from typing import Literal, NamedTuple


class GeneLocation(NamedTuple):
    id: str
    seq_id: (
        str  # the sequence ID (e.g., chromosome or scaffold) where the gene is located
    )
    start: int  # 0-based
    end: int  # 0-based, inclusive
    strand: Literal["+", "-"]  # the strand of the gene
