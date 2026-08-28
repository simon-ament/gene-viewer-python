from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class GeneLocation:
    id: str
    seq_id: (
        str  # the sequence ID (e.g., chromosome or scaffold) where the gene is located
    )
    start: int  # 1-based
    end: int  # 1-based, inclusive
    strand: Literal["+", "-"]  # the strand of the gene
