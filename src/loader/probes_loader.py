import hashlib
from abc import ABC, abstractmethod
from collections import defaultdict


class ProbesLoader(ABC):
    @abstractmethod
    def load_gene(self, gene_id: str):
        pass

    @abstractmethod
    def gene_list(self):
        pass

    @abstractmethod
    def cache_id(self):
        pass


class ProbesLoaderManual(ProbesLoader):
    def __init__(self):
        self.probes = defaultdict(
            lambda: defaultdict(list)
        )  # {gene_id: {probeset_id: [probe1, probe2, ...]}}
        self.cache_id_hash = hashlib.sha256()

    def add_probe(self, gene_id: str, probeset_id: str, probe):
        self.probes[gene_id][probeset_id].append(probe)
        self.cache_id_hash.update(f"{gene_id}:{probeset_id}:{probe}".encode())

    def load_gene(self, gene_id: str):
        return self.probes.get(gene_id, {})

    def gene_list(self):
        return list(self.probes.keys())

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        return self.cache_id_hash.hexdigest()
