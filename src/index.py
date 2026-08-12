from abc import ABC, abstractmethod
from io import StringIO  # for creating file-like objects from strings

from oligo_designer_toolsuite.utils import FastaParser
from pybedtools import BedTool
from pybedtools import helpers as _pybedtools_helpers

from helpers import _get_GTF_attribute
from src.types import GeneLocation


class FileIndex(ABC):
    def __init__(self, file_path, gene_locations: dict[str, GeneLocation]):
        self._file_path = file_path
        self._file = open(file_path, "r")
        self._index = None
        self._gene_locations = gene_locations

    @abstractmethod
    def _build_index(self):
        # use self._gene_locations to build the index
        pass

    @abstractmethod
    def get(self, key):
        if self._index is None:
            self._index = self._build_index()

    def keys(self):
        if self._index is None:
            self._index = self._build_index()
        return list(self._index.keys())

    def __del__(self):
        self._file.close()


class GeneAwareFileIndex(FileIndex):
    def __init__(
        self,
        file_path,
        gene_list: list[str] | None = None,
        collect_gene_locations: bool = False,
    ):
        super().__init__(file_path, gene_locations={})
        self._gene_list = gene_list
        self._collect_gene_locations = collect_gene_locations

    @abstractmethod
    def _build_index(self):
        # use self._gene_list to build the index
        # if self._collect_gene_locations is True, populate self._gene_locations
        pass

    @property
    def gene_locations(self):
        if self._index is None:
            self._index = self._build_index()
        if not self._collect_gene_locations:
            raise ValueError("Gene location collection is disabled.")
        return self._gene_locations


### Simple file indexes


# key: SeqID, return: sequence
class FastaFileIndex(FileIndex):
    def _build_index(self):
        index = {}
        position = 0
        seq_len_counter = 0
        next_gene_index = 0
        for line in self._file:
            if line.startswith(">"):
                seq_id = (
                    line[1:].strip().split()[0]
                )  # Use the first word after '>' as the seq_id
                seq_len_counter = 0
                next_gene_index = 0
                sorted_genes = self._gene_locations.get(
                    seq_id, []
                )  # already sorted by start position
                next_gene = sorted_genes[next_gene_index] if sorted_genes else None
            else:
                while (
                    next_gene and seq_len_counter + len(line.strip()) > next_gene.start
                ):
                    index[next_gene.id] = (
                        position + next_gene.start - seq_len_counter,
                        next_gene.end - next_gene.start + 1,
                    )  # position and length
                    next_gene_index += 1
                    next_gene = (
                        sorted_genes[next_gene_index]
                        if next_gene_index < len(sorted_genes)
                        else None
                    )
                seq_len_counter += len(line.strip())
            position += len(line)
        return index

    def get(self, key):
        super().get(key)  # Ensure the index is built
        if key not in self._index:
            raise KeyError(f"Key '{key}' not found in index.")
        (position, length) = self._index[key]
        self._file.seek(position)
        # TODO: check if length is correct
        sequence = self._file.read(length).replace("\n", "")
        return sequence


# key: gene ID, return: list of features
class BEDFileIndex(FileIndex):
    def _build_index(self):
        index = {}
        position = 0
        for line in self._file:
            if line.startswith(("#", "track", "browser")):
                position += len(line)
                continue
            string_io = StringIO(line)
            bedtool = BedTool(string_io)
            # clean up the tag to avoid memory leaks in pybedtools
            # TODO: find a better way to clean up
            try:
                _pybedtools_helpers._tags.pop(bedtool._tag, None)
            except Exception:
                pass
            for feature in bedtool:
                start = feature.start
                end = feature.end
                # find genes from self._gene_locations that overlaps with the feature
                # use binary search to find the first gene that overlaps with the feature (self._gene_locations is sorted by start position)
                gene_locations = self._gene_locations.get(feature.chrom, [])
                first_gene_index = binary_search_gene_locations(
                    gene_locations, start, end
                )
                if first_gene_index is not None:
                    for gene_location in gene_locations[first_gene_index:]:
                        if gene_location.start > end:
                            break
                        if gene_location.end >= start:
                            if gene_location.id not in index:
                                index[gene_location.id] = []
                            index[gene_location.id].append(position)

            position += len(line)
        return index

    def get(self, key):
        super().get(key)  # Ensure the index is built
        if key not in self._index:
            raise KeyError(f"Key '{key}' not found in index.")
        positions = self._index[key]
        features = []
        for pos in positions:
            self._file.seek(pos)
            line = self._file.readline().strip()
            string_io = StringIO(line)
            bedtool = BedTool(string_io)
            try:
                _pybedtools_helpers._tags.pop(bedtool._tag, None)
            except Exception:
                pass
            for feature in bedtool:
                features.append(feature)
        return features


### Gene-aware file indexes


# key: gene ID, return: list of features
class GTFFileIndex(GeneAwareFileIndex):
    def _build_index(self):
        index = {}
        position = 0
        for line in self._file:
            if line.startswith(("#", "track", "browser")):
                position += len(line)
                continue
            string_io = StringIO(line)
            bedtool = BedTool(string_io)
            # clean up the tag to avoid memory leaks in pybedtools
            try:
                _pybedtools_helpers._tags.pop(bedtool._tag, None)
            except Exception:
                pass

            for feature in bedtool:
                gene_id = _get_GTF_attribute(feature[8], "gene_id")
                if gene_id:
                    if gene_id not in index:
                        index[gene_id] = []
                    index[gene_id].append(position)
                    type = feature[2]
                    if type == "gene" and self._collect_gene_locations:
                        gene_location = GeneLocation(
                            id=gene_id,
                            seq_id=feature.chrom,
                            start=feature.start,
                            end=feature.end,
                            strand=feature.strand,
                        )
                        if gene_location.seq_id not in self._gene_locations:
                            self._gene_locations[gene_location.seq_id] = []
                        self._gene_locations[gene_location.seq_id].append(gene_location)
            position += len(line)

        if self._collect_gene_locations:
            for seq_id in self._gene_locations:
                self._gene_locations[seq_id] = sorted(
                    self._gene_locations[seq_id], key=lambda g: g.start
                )
        return index

    def get(self, key):
        super().get(key)  # Ensure the index is built
        if key not in self._index:
            raise KeyError(f"Key '{key}' not found in index.")
        positions = self._index[key]
        features = []
        for pos in positions:
            self._file.seek(pos)
            line = self._file.readline().strip()
            string_io = StringIO(line)
            bedtool = BedTool(string_io)
            # clean up the tag to avoid memory leaks in pybedtools
            try:
                _pybedtools_helpers._tags.pop(bedtool._tag, None)
            except Exception:
                pass
            for feature in bedtool:
                features.append(feature)
        return features


# key: gene ID, return: list of (header, sequence)
class ODTFastaFileIndex(GeneAwareFileIndex):
    def _build_index(self):
        index = {}
        position = 0
        fasta_parser = FastaParser()
        for line in self._file:
            if line.startswith(">"):
                region_name, additional_info, coordinates = (
                    fasta_parser.parse_fasta_header(line.strip())
                )
                gene_id = region_name.lstrip(">")
                if gene_id:
                    if gene_id not in index:
                        index[gene_id] = []
                    index[gene_id].append(position)
                    # NOTE: ODT Fasta must include gene entries if used as sequence source
                    if (
                        additional_info.get("regiontype") == "gene"
                        and self._collect_gene_locations
                    ):
                        gene_location = GeneLocation(
                            id=gene_id,
                            seq_id=coordinates.get("chromosome"),
                            start=coordinates.get("start"),
                            end=coordinates.get("end"),
                            strand=additional_info.get("strand"),
                        )
                        if gene_location.seq_id not in self._gene_locations:
                            self._gene_locations[gene_location.seq_id] = []
                        self._gene_locations[gene_location.seq_id].append(gene_location)
                    # TODO: infer gene start and end from other regions if gene entry is not present
            position += len(line)

        if self._collect_gene_locations:
            for seq_id in self._gene_locations:
                self._gene_locations[seq_id] = sorted(
                    self._gene_locations[seq_id], key=lambda g: g.start
                )
        return index

    def get(self, key):
        super().get(key)  # Ensure the index is built
        if key not in self._index:
            raise KeyError(f"Key '{key}' not found in index.")
        positions = self._index[key]
        sequences = []
        for pos in positions:
            self._file.seek(pos)
            header = self._file.readline().strip()
            sequence_lines = []
            while True:
                line = self._file.readline()
                if not line or line.startswith(">"):
                    break
                sequence_lines.append(line.strip())
            sequences.append((header, "".join(sequence_lines)))
        return sequences


def binary_search_gene_locations(gene_locations, start, end):
    """Binary search to find the first gene location that overlaps with the given start and end."""
    low = 0
    high = len(gene_locations) - 1
    while low <= high:
        mid = (low + high) // 2
        gene_location = gene_locations[mid]
        if gene_location.end < start:
            low = mid + 1
        elif gene_location.start > end:
            high = mid - 1
        else:
            # Found an overlapping gene location, now find the first one
            while mid > 0 and gene_locations[mid - 1].end >= start:
                mid -= 1
            return mid
    return None
