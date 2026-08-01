import shelve
from abc import ABC, abstractmethod

from pybedtools import BedTool

from helpers import _hash_file


class TrackLoader(ABC):
    @abstractmethod
    def load_gene(self, gene_id: str):
        pass

    @abstractmethod
    def gene_list(self):
        pass

    @abstractmethod
    def cache_id(self):
        pass


class TrackLoaderBED(TrackLoader):
    def __init__(self, bed_file_path: str):
        self.bed_file_path = bed_file_path
        self.cache_id_str = None

        # create a persistent index for BED features and keep it open for later access
        self.features_index = shelve.open(
            f"{bed_file_path}.index", flag="c", writeback=True
        )

        # Load BED file
        self.features = BedTool(bed_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash to create a unique cache_id
            self.cache_id_str = f"{_hash_file(self.bed_file_path)}"
        return self.cache_id_str


class TrackLoaderGTF(TrackLoader):
    def __init__(self, gtf_file_path: str):
        self.gtf_file_path = gtf_file_path
        self.cache_id_str = None

        # create a persistent index for GTF features and keep it open for later access
        self.features_index = shelve.open(
            f"{gtf_file_path}.index", flag="c", writeback=True
        )

        # Load GTF file
        self.features = BedTool(gtf_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash to create a unique cache_id
            self.cache_id_str = f"{_hash_file(self.gtf_file_path)}"
        return self.cache_id_str
