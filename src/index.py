from abc import ABC, abstractmethod
from collections import defaultdict
from io import StringIO  # for creating file-like objects from strings

from oligo_designer_toolsuite.utils import FastaParser
from pybedtools import BedTool
from pybedtools import helpers as _pybedtools_helpers

from helpers import _get_GTF_attribute
from src.types import GeneLocation


class FileIndex(ABC):
    def __init__(self, file_path, gene_locations: dict[str, list[GeneLocation]]):
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
        super().__init__(file_path, gene_locations=defaultdict(list))
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
        position = 0  # 0-based position of the first character in the file
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
                    next_gene
                    and seq_len_counter + len(line.strip()) > next_gene.start - 1
                ):
                    index[next_gene.id] = (
                        position + (next_gene.start - 1 - seq_len_counter),
                        next_gene.end - next_gene.start,
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
        seq = ""
        # read the sequence line by line and ignore any newline characters
        while len(seq) < length:
            line = self._file.readline()
            if not line or line.startswith(">"):
                break
            line = line.strip()
            if len(seq) + len(line) > length:
                seq += line[: length - len(seq)]
            else:
                seq += line.strip()
        return seq


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
        inferred_gene_locations = {}
        explicit_gene_locations = {}
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
                if not gene_id:
                    continue
                if gene_id not in index:
                    index[gene_id] = []
                index[gene_id].append(position)

                if not self._collect_gene_locations:
                    continue

                feature_type = feature[2]
                start = int(feature.start)
                end = int(feature.end)
                strand = feature.strand

                if feature_type == "gene":
                    explicit_gene_locations[gene_id] = GeneLocation(
                        id=gene_id,
                        seq_id=feature.chrom,
                        start=start,
                        end=end,
                        strand=strand,
                    )
                    continue

                # non-gene feature, use to infer gene locations
                gene_location = inferred_gene_locations.setdefault(
                    gene_id,
                    {
                        "seq_id": feature.chrom,
                        "start": start,
                        "end": end,
                        "strand": strand
                        if strand != "."
                        else "+",  # could be intron without a strand
                    },
                )
                gene_location["seq_id"] = feature.chrom
                gene_location["start"] = min(gene_location["start"], start)
                gene_location["end"] = max(gene_location["end"], end)
                if strand not in (None, "."):
                    gene_location["strand"] = strand

            position += len(line)

        if self._collect_gene_locations:
            self._gene_locations = gather_gene_locations(
                explicit_gene_locations, inferred_gene_locations
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
        inferred_gene_locations = {}
        explicit_gene_locations = {}
        for line in self._file:
            if line.startswith(">"):
                region_name, additional_info, coordinates = (
                    fasta_parser.parse_fasta_header(line.strip())
                )
                gene_id = region_name.lstrip(">")
                if not gene_id:
                    position += len(line)
                    continue
                if gene_id not in index:
                    index[gene_id] = []
                index[gene_id].append(position)

                if not self._collect_gene_locations:
                    position += len(line)
                    continue

                seq_id = coordinates["chromosome"][0]
                start = int(coordinates["start"][0])
                end = int(coordinates["end"][0])
                strand = additional_info.get("strand", ["+"])[0]
                region_type = additional_info.get("regiontype", ["unknown"])[0]

                if region_type == "gene":
                    explicit_gene_locations[gene_id] = GeneLocation(
                        id=gene_id,
                        seq_id=seq_id,
                        start=start,
                        end=end,
                        strand=strand,
                    )
                    position += len(line)
                    continue

                # non-gene feature, use to infer gene locations
                gene_location = inferred_gene_locations.setdefault(
                    gene_id,
                    {
                        "seq_id": seq_id,
                        "start": start,
                        "end": end,
                        "strand": strand
                        if strand not in (None, ".")
                        else "+",  # could be intron without a strand
                    },
                )
                gene_location["seq_id"] = seq_id
                gene_location["start"] = min(gene_location["start"], start)
                gene_location["end"] = max(gene_location["end"], end)
                if strand not in (None, "."):
                    gene_location["strand"] = strand
            position += len(line)

        if self._collect_gene_locations:
            self._gene_locations = gather_gene_locations(
                explicit_gene_locations, inferred_gene_locations
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


def gather_gene_locations(explicit_gene_locations, inferred_gene_locations):
    """Combine explicit and inferred gene locations into a single dictionary."""
    combined_gene_locations = defaultdict(list)

    # add explicit gene locations
    for gene_id, gene_location in explicit_gene_locations.items():
        combined_gene_locations[gene_location.seq_id].append(gene_location)

    # add inferred gene locations, avoiding duplicates
    for gene_id, gene_location in inferred_gene_locations.items():
        if gene_id not in explicit_gene_locations:
            seq_id = gene_location["seq_id"]
            inferred = GeneLocation(
                id=gene_id,
                seq_id=seq_id,
                start=gene_location["start"],
                end=gene_location["end"],
                strand=gene_location["strand"],
            )
            combined_gene_locations[seq_id].append(inferred)

    # sort the gene locations by start position for each sequence ID
    for seq_id in combined_gene_locations:
        combined_gene_locations[seq_id] = sorted(
            combined_gene_locations[seq_id], key=lambda g: g.start
        )

    return combined_gene_locations
