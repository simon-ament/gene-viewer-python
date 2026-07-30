from abc import ABC, abstractmethod
import re

from Bio import SeqIO
from pybedtools import BedTool
from oligo_designer_toolsuite.utils import FastaParser
from collections import defaultdict

class SequencesLoader(ABC):
    @abstractmethod
    def load_gene(self, gene_id: str):
        pass

    @abstractmethod
    def gene_list(self):
        pass

    @abstractmethod
    def cache_id(self):
        pass

class SequencesLoaderFastaGTF(SequencesLoader):
    """Loader class for loading sequences from a Fasta file and gene annotations from a GTF file.

    Args:
        fasta_file_path (str): Path to the Fasta file containing the sequences.
        gtf_file_path (str): Path to the GTF file containing gene annotations.
        region_types (list, optional): List of region types to include. Defaults to [].
    """

    def __init__(self, fasta_file_path: str, gtf_file_path: str, region_types: list = []):
        self.fasta_file_path = fasta_file_path
        self.gtf_file_path = gtf_file_path
        self.region_types = region_types

        self.records_index = SeqIO.index(fasta_file_path, "fasta")
        self.gene_features_map = defaultdict(list)  # {gene_id: [int]}

        # Load GTF file and filter features based on region_types
        self.features = BedTool(gtf_file_path)
        for feature_idx, feature in enumerate(self.features):
            if feature[2] in self.region_types or not self.region_types: # feature[2] is the feature type (e.g., "exon", "intron")
                # .attrs fails parsing some GTF files, so we use a regex to extract the gene_id from the attributes string
                # gene_id = feature.attrs.get("gene_id")
                gene_id = self._get_feature_attribute(feature[8], "gene_id")
                if gene_id:
                    self.gene_features_map[gene_id].append(feature_idx)

    @staticmethod
    def _get_feature_attribute(attributes: str, key: str):
        pattern = rf'(?:^|;\s*){re.escape(key)}\s+"([^"]*)"'
        match = re.search(pattern, attributes)
        if match:
            return match.group(1)
        return None

    def load_gene(self, gene_id: str):
        """Loads sequences associated with a specific gene ID from the Fasta file.

        Args:
            gene_id (str): The gene ID for which to load sequences.

        Returns:
            list: A list sequence records (start, sequence) associated with the specified gene ID.
        """

        sequences = []
        for feature_idx in self.gene_features_map.get(gene_id, []):
            feature = self.features[feature_idx]
            seq_record = self.records_index.get(feature[0]) # feature[0] is the sequence ID (e.g., chromosome or scaffold name)
            if seq_record:
                start = int(feature.start) # one-based start position
                end = int(feature.end) # one-based end position
                sequence = str(seq_record.seq[start-1:end]) # extract the sequence from the Fasta record (convert to zero-based indexing)
                sequences.append((start, sequence))
        return sequences


    def gene_list(self):
        """Returns a list of gene IDs available in the GTF file.

        Returns:
            list: A list of gene IDs.
        """
        return list(self.gene_features_map.keys())

    def cache_id(self):
        pass
        

class SequencesLoaderODTFasta(SequencesLoader):
    """Loader class for loading sequences from an ODTFasta file.

    Args:
        odt_fasta_file_path (str): Path to the ODTFasta file containing the sequences.
        region_types (list, optional): List of region types to include. Defaults to [].
    """

    def __init__(self, odt_fasta_file_path: str, region_types: list = []):
        self.odt_fasta_file_path = odt_fasta_file_path
        self.region_types = region_types
        self.records_index = SeqIO.index(odt_fasta_file_path, "fasta")
        # create map of Gene IDs to their records to avoid loading all sequences into memory at once
        self.gene_records_map = defaultdict(list)  # {gene_id: [{idx: str, start: int}, {idx: str, start: int}, ...]}

        # iterate through the ODTFasta file and parse the headers using FastaParser
        fasta_parser = FastaParser()
        for idx in self.records_index:
            region_name, additional_info, coordinates = fasta_parser.parse_fasta_header(idx)
            gene_id = region_name.lstrip(">")
            # only include sequences that match the specified region_types (if provided)
            if gene_id and (additional_info.get("type") in self.region_types or not self.region_types):
                self.gene_records_map[gene_id].append({"idx": idx, "start": coordinates.get("start")})

    def load_gene(self, gene_id: str):
        """Loads sequences associated with a specific gene ID from the ODTFasta file.

        Args:
            gene_id (str): The gene ID for which to load sequences.

        Returns:
            list: A list of sequence records (start, sequence) associated with the specified gene ID.
        """
        return [(record_info["start"], str(self.records_index[record_info["idx"]].seq)) for record_info in self.gene_records_map.get(gene_id, [])]
        

    def gene_list(self):
        """Returns a list of gene IDs available in the ODTFasta file.

        Returns:
            list: A list of gene IDs.
        """
        return list(self.gene_records_map.keys())

    def cache_id(self):
        pass
        
