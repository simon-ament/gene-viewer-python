from abc import ABC, abstractmethod

from gene_viewer.types import GeneLocation


class Loader(ABC):
    def __init__(self):
        self._lazy_init_done = False

    @abstractmethod
    def cache_id(self):
        pass

    def _lazy_init(self):
        """Perform any deferred initialization here. Cache_id should be available before this is called."""

    @abstractmethod
    def load_gene(self, gene: GeneLocation):
        if not self._lazy_init_done:
            self._lazy_init()
            self._lazy_init_done = True
