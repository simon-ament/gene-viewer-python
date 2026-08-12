import hashlib
from abc import abstractmethod
from collections import defaultdict

from src.types import GeneLocation

from .loader import Loader


class ProbesLoader(Loader):
    @property
    @abstractmethod
    def gene_list(self):
        """Return a list of gene IDs for which probes are available."""


class ProbesLoaderManual(ProbesLoader):
    def __init__(self):
        super().__init__()
        self._probes = defaultdict(
            lambda: defaultdict(list)
        )  # {gene_id: {probeset_id: [probe1, probe2, ...]}}
        self._cache_id_hash = hashlib.sha256()

    def add_probe(self, gene_id: str, probeset_id: str, probe):
        self._probes[gene_id][probeset_id].append(probe)
        self._cache_id_hash.update(f"{gene_id}:{probeset_id}:{probe}".encode())

    @property
    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        return self._cache_id_hash.hexdigest()

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        return self._probes.get(gene.id, {})

    @property
    def gene_list(self):
        return list(self._probes.keys())
