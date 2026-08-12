from Bio.Seq import Seq
from oligo_designer_toolsuite.utils import FastaParser

from helpers import _hash_file
from src.index import FastaFileIndex, ODTFastaFileIndex
from src.types import GeneLocation

from .loader import Loader


class SequencesLoader(Loader):
    def __init__(self):
        super().__init__()
        self._gene_locations = None

    def set_gene_locations(self, gene_locations: dict[str, GeneLocation]):
        self._gene_locations = gene_locations


class SequencesLoaderFasta(SequencesLoader):
    def __init__(self, fasta_file_path: str):
        super().__init__()
        self._fasta_file_path = fasta_file_path
        self._fasta_file_index = None  # will be initialized lazily
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = f"{_hash_file(self._fasta_file_path)}"
        return self._cache_id_str

    def _lazy_init(self):
        if not self._gene_locations:
            raise ValueError("Gene locations must be set before loading sequences.")
        self._fasta_file_index = FastaFileIndex(
            self._fasta_file_path, gene_locations=self._gene_locations
        )

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        sequences = []
        sequence = self._fasta_file_index.get(gene.id)
        if gene.strand == "-":
            sequence = str(Seq(sequence).complement())
        if sequence:
            start = gene.start
            sequences.append({"start": start, "sequence": sequence})

        return deduplicate_sequences(sequences)


class SequencesLoaderODTFasta(SequencesLoader):
    def __init__(self, odt_fasta_file_path: str, region_types: list | None = None):
        super().__init__()
        self._odt_fasta_file_path = odt_fasta_file_path
        self._odt_fasta_file_index = None  # will be initialized lazily
        self._region_types = region_types or []
        self._fasta_parser = None  # will be initialized lazily
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = f"{_hash_file(self._odt_fasta_file_path)}_{'_'.join(self._region_types)}"
        return self._cache_id_str

    def _lazy_init(self):
        if not self._gene_locations:
            raise ValueError("Gene locations must be set before loading sequences.")
        gene_list = [gene.id for gene in self._gene_locations]
        self._odt_fasta_file_index = ODTFastaFileIndex(
            self._odt_fasta_file_path, gene_list=gene_list
        )
        self._fasta_parser = FastaParser()

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        sequences = []
        for header, sequence in self._odt_fasta_file_index.get(gene.id):
            _, additional_info, coordinates = self._fasta_parser.parse_fasta_header(
                header
            )
            if gene.strand == "-":
                sequence = str(Seq(sequence).complement())
            if additional_info.get("regiontype") in self._region_types:
                start = coordinates.get(
                    "start"
                )  # TODO: check indexing, should be 0-based
                sequences.append({"start": start, "sequence": sequence})
        return deduplicate_sequences(sequences)


def deduplicate_sequences(sequences):
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
            merged_sequence = (
                deduplicated_sequences[-1]["sequence"]
                + seq["sequence"][last_end - seq_start + 1 :]
            )
            deduplicated_sequences[-1]["sequence"] = merged_sequence
            last_end = seq_end
        # else: fully contained sequence, skip it

    return deduplicated_sequences
