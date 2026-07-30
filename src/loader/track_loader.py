from abc import ABC, abstractmethod

from pybedtools import BedTool

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

        # Load BED file
        self.features = BedTool(bed_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        pass

class TrackLoaderGTF(TrackLoader):
    def __init__(self, gtf_file_path: str):
        self.gtf_file_path = gtf_file_path

        # Load GTF file
        self.features = BedTool(gtf_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        pass
