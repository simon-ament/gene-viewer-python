import shelve
from collections import defaultdict

from Bio import SeqIO
from oligo_designer_toolsuite.utils import FastaParser
from pybedtools import BedTool

from helpers import _get_feature_attribute, _hash_file
from .loader import FileBasedLoader


class SequencesLoader(FileBasedLoader):
    pass


class SequencesLoaderFastaGTF(SequencesLoader):
    """Loader class for loading sequences from a Fasta file and gene annotations from a GTF file.

    Args:
        fasta_file_path (str): Path to the Fasta file containing the sequences.
        gtf_file_path (str): Path to the GTF file containing gene annotations.
        region_types (list, optional): List of region types to include. Defaults to [].
    """

    def __init__(
        self, fasta_file_path: str, gtf_file_path: str, region_types: list | None = None
    ):
        super().__init__()
        self.fasta_file_path = fasta_file_path
        self.gtf_file_path = gtf_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        self.records_index = None  # will be initialized in load_files()
        self.gene_features_map = defaultdict(list)  # {gene_id: [int]}
        self.features_index = None  # will be initialized in load_files()
        self.features = None  # will be initialized in load_files()

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = f"{_hash_file(self.fasta_file_path)}_{_hash_file(self.gtf_file_path)}_{'_'.join(self.region_types)}"
        return self.cache_id_str

    def load_files(self):
        self.records_index = shelve.open(f"{self.fasta_file_path}_sequencesFastaGTF.shelve", flag="c", writeback=True)
        for record in SeqIO.parse(self.fasta_file_path, "fasta"):
            self.records_index[record.id] = record.seq
        # create a persistent index for GTF features and keep it open for later access
        self.features_index = shelve.open(
            f"{self.gtf_file_path}_sequencesFastaGTF.shelve", flag="c", writeback=True
        )

        # Load GTF file and filter features based on region_types
        self.features = BedTool(self.gtf_file_path)
        for feature_idx, feature in enumerate(self.features):
            if (
                feature[2] in self.region_types or not self.region_types
            ):  # feature[2] is the feature type (e.g., "exon", "intron")
                gene_id = _get_feature_attribute(
                    feature[8], "gene_id"
                )  # .attrs fails parsing some GTF files
                # gene_id = feature.attrs.get("gene_id")
                if gene_id:
                    self.gene_features_map[gene_id].append(feature_idx)
                    self.features_index[str(feature_idx)] = {
                        "seq_id": feature[0],
                        "start": feature.start,
                        "end": feature.end,
                    }

    def load_gene(self, gene_id: str):
        """Loads sequences associated with a specific gene ID from the Fasta file.

        Args:
            gene_id (str): The gene ID for which to load sequences.

        Returns:
            list: A list sequence records (start, sequence) associated with the specified gene ID.
        """

        super().load_gene()  # ensure files are loaded
        sequences = []
        for feature_idx in self.gene_features_map.get(gene_id, []):
            feature = self.features_index.get(
                str(feature_idx)
            )  # retrieve the feature from the persistent index
            seq_record = self.records_index.get(
                feature["seq_id"]
            )  # feature["seq_id"] is the sequence ID (e.g., chromosome or scaffold name)
            if seq_record:
                start = int(feature["start"])  # one-based start position
                end = int(feature["end"])  # one-based end position
                sequence = str(
                    seq_record[start - 1 : end]
                )  # extract the sequence from the Fasta record (convert to zero-based indexing)
                sequences.append({"start": start, "sequence": sequence})

        # deduplicate sequences by removing sequences fully contained within other sequences and merging overlapping sequences
        sorted_sequences = sorted(sequences, key=lambda x: (x["start"], len(x["sequence"])))
        deduplicated_sequences = []
        last_end = -1
        for seq in sorted_sequences:
            seq_start = seq["start"]
            seq_end = seq_start + len(seq["sequence"]) - 1
            if seq_start > last_end:
                deduplicated_sequences.append(seq)
                last_end = seq_end
            elif seq_end > last_end:
                # merge overlapping sequences
                merged_sequence = deduplicated_sequences[-1]["sequence"] + seq["sequence"][last_end - seq_start + 1 :]
                deduplicated_sequences[-1]["sequence"] = merged_sequence
                last_end = seq_end
            # else: fully contained sequence, skip it
        
        return deduplicated_sequences

    def gene_list(self):
        """Returns a list of gene IDs available in the GTF file.

        Returns:
            list: A list of gene IDs.
        """
        super().gene_list()  # ensure files are loaded
        return list(self.gene_features_map.keys())

    def delete(self):
        records_index = getattr(self, "records_index", None)
        if records_index is not None:
            try:
                records_index.close()
            except Exception:
                pass
            finally:
                self.records_index = None

        features_index = getattr(self, "features_index", None)
        if features_index is not None:
            try:
                features_index.close()
            except Exception:
                pass
            finally:
                self.features_index = None


class SequencesLoaderODTFasta(SequencesLoader):
    """Loader class for loading sequences from an ODTFasta file.

    Args:
        odt_fasta_file_path (str): Path to the ODTFasta file containing the sequences.
        region_types (list, optional): List of region types to include. Defaults to [].
    """

    def __init__(self, odt_fasta_file_path: str, region_types: list | None = None):
        super().__init__()
        self.odt_fasta_file_path = odt_fasta_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        self.records_index = None  # will be initialized in load_files()
        # create map of Gene IDs to their records to avoid loading all sequences into memory at once
        self.gene_records_map = defaultdict(
            list
        )  # {gene_id: [{idx: str, start: int}, {idx: str, start: int}, ...]}

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = (
                f"{_hash_file(self.odt_fasta_file_path)}_{'_'.join(self.region_types)}"
            )
        return self.cache_id_str

    def load_files(self):
        self.records_index = shelve.open(f"{self.odt_fasta_file_path}_sequencesODT.shelve", flag="c", writeback=True)
        for record in SeqIO.parse(self.odt_fasta_file_path, "fasta"):
            self.records_index[record.id] = record

        # iterate through the ODTFasta file and parse the headers using FastaParser
        fasta_parser = FastaParser()
        for idx in self.records_index:
            region_name, additional_info, coordinates = fasta_parser.parse_fasta_header(
                idx
            )
            gene_id = region_name.lstrip(">")
            # only include sequences that match the specified region_types (if provided)
            if gene_id and (
                additional_info.get("type") in self.region_types
                or not self.region_types
            ):
                self.gene_records_map[gene_id].append(
                    {"idx": idx, "start": coordinates.get("start")}
                )

    def load_gene(self, gene_id: str):
        """Loads sequences associated with a specific gene ID from the ODTFasta file.

        Args:
            gene_id (str): The gene ID for which to load sequences.

        Returns:
            list: A list of sequence records (start, sequence) associated with the specified gene ID.
        """
        super().load_gene()  # ensure files are loaded
        return [
            {
                "start": record_info["start"],
                "sequence": str(self.records_index[record_info["idx"]].seq),
            }
            for record_info in self.gene_records_map.get(gene_id, [])
        ]

    def gene_list(self):
        """Returns a list of gene IDs available in the ODTFasta file.

        Returns:
            list: A list of gene IDs.
        """
        super().gene_list()  # ensure files are loaded
        return list(self.gene_records_map.keys())

    def delete(self):
        records_index = getattr(self, "records_index", None)
        if records_index is not None:
            try:
                records_index.close()
            except Exception:
                pass
            finally:
                self.records_index = None
