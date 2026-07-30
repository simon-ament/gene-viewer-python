from abc import ABC, abstractmethod

from Bio import SeqIO
from pybedtools import BedTool
from oligo_designer_toolsuite.utils import FastaParser
from collections import defaultdict

class RegionsLoader(ABC):
    @abstractmethod
    def load_gene(self, gene_id: str):
        pass

    @abstractmethod
    def gene_list(self):
        pass

    @abstractmethod
    def cache_id(self):
        pass

class RegionsLoaderGTF(RegionsLoader):
    def __init__(self, gtf_file_path: str, region_types: list = []):
        self.gtf_file_path = gtf_file_path
        self.region_types = region_types

        # Load GTF file and filter features based on region_types
        self.features = BedTool(gtf_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        pass

class RegionsLoaderODTFasta(RegionsLoader):
    def __init__(self, odt_fasta_file_path: str, region_types: list = []):
        self.odt_fasta_file_path = odt_fasta_file_path
        self.region_types = region_types
        self.records_index = SeqIO.index(odt_fasta_file_path, "fasta")
        # create map of Gene IDs to their records to avoid loading all sequences into memory at once
        self.gene_records_map = defaultdict(list)  # {gene_id: [idx, idx, ...]}

        # iterate through the ODTFasta file and parse the headers using FastaParser
        fasta_parser = FastaParser()
        for idx in self.records_index:
            region_name, additional_info, coordinates = fasta_parser.parse_fasta_header(idx)
            gene_id = region_name.lstrip(">")
            # only include sequences that match the specified region_types (if provided)
            if gene_id and (additional_info.get("type") in self.region_types or not self.region_types):
                self.gene_records_map[gene_id].append(idx)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        pass
