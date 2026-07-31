from abc import ABC, abstractmethod
from collections import defaultdict

from Bio import SeqIO
from oligo_designer_toolsuite.utils import FastaParser
from pybedtools import BedTool

from helpers import _hash_file


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
    def __init__(self, gtf_file_path: str, region_types: list | None = None):
        self.gtf_file_path = gtf_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        # Load GTF file and filter features based on region_types
        self.features = BedTool(gtf_file_path)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = (
                f"{_hash_file(self.gtf_file_path)}_{'_'.join(self.region_types)}"
            )
        return self.cache_id_str


class RegionsLoaderODTFasta(RegionsLoader):
    def __init__(
        self,
        odt_fasta_file_path: str,
        region_types: list | None = None,
    ):
        self.odt_fasta_file_path = odt_fasta_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        self.records_index = SeqIO.index(odt_fasta_file_path, "fasta")
        # create map of Gene IDs to their records to avoid loading all sequences into memory at once
        self.gene_records_map = defaultdict(list)  # {gene_id: [idx, idx, ...]}

        # iterate through the ODTFasta file and parse the headers using FastaParser
        fasta_parser = FastaParser()
        for idx in self.records_index:
            region_name, additional_info, _coordinates = fasta_parser.parse_fasta_header(
                idx
            )
            gene_id = region_name.lstrip(">")
            # only include sequences that match the specified region_types (if provided)
            if gene_id and (
                additional_info.get("type") in self.region_types
                or not self.region_types
            ):
                self.gene_records_map[gene_id].append(idx)

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = (
                f"{_hash_file(self.odt_fasta_file_path)}_{'_'.join(self.region_types)}"
            )
        return self.cache_id_str
