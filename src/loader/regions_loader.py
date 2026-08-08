import shelve
from collections import defaultdict

from Bio import SeqIO
from oligo_designer_toolsuite.utils import FastaParser
from pybedtools import BedTool

from helpers import _get_feature_attribute, _hash_file
from .loader import FileBasedLoader


class RegionsLoader(FileBasedLoader):
    pass


class RegionsLoaderGTF(RegionsLoader):
    def __init__(self, gtf_file_path: str, region_types: list | None = None):
        super().__init__()
        self.gtf_file_path = gtf_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        self.gene_features_map = defaultdict(list)  # {gene_id: [feature, feature, ...]}
        self.features_index = None  # will be initialized in load_files()
        self.features = None  # will be initialized in load_files()

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = (
                f"{_hash_file(self.gtf_file_path)}_{'_'.join(self.region_types)}"
            )
        return self.cache_id_str

    def load_files(self):
        # create a persistent index for GTF features and keep it open for later access
        self.features_index = shelve.open(
            f"{self.gtf_file_path}_regionsGTF.shelve", flag="c", writeback=True
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
                    # TODO: handle lists
                    transcript_id = _get_feature_attribute(feature[8], "transcript_id")
                    self.features_index[str(feature_idx)] = {
                        "start": feature.start,
                        "end": feature.end,
                        "type": feature[2],
                        "strand": feature.strand,
                        "transcript_id": transcript_id,
                    }

    def load_gene(self, gene_id: str):
        super().load_gene()  # ensure files are loaded
        regions = defaultdict(list)  # {transcript_id: [(start, end, type), ...]}
        for feature_idx in self.gene_features_map.get(gene_id, []):
            feature = self.features_index.get(
                str(feature_idx)
            )  # retrieve the feature from the persistent index
            start = int(feature["start"])  # one-based start position
            end = int(feature["end"])  # one-based end position
            regions[feature["transcript_id"]].append(
                {"start": start, "end": end, "type": feature["type"], "strand": feature["strand"]}
            )
        return regions

    def gene_list(self):
        super().gene_list()  # ensure files are loaded
        return list(self.gene_features_map.keys())

    def delete(self):
        features_index = getattr(self, "features_index", None)
        if features_index is not None:
            try:
                features_index.close()
            except Exception:
                pass
            finally:
                self.features_index = None


class RegionsLoaderODTFasta(RegionsLoader):
    def __init__(
        self,
        odt_fasta_file_path: str,
        region_types: list | None = None,
    ):
        super().__init__()
        self.odt_fasta_file_path = odt_fasta_file_path
        self.region_types = region_types or []
        self.cache_id_str = None

        self.records_index = None  # will be initialized in load_files()
        # create map of Gene IDs to their records to avoid loading all sequences into memory at once
        self.gene_records_map = defaultdict(
            list
        )  # {gene_id: [{idx: str, start: int}, {idx: str, start: int}, ...]}

        # iterate through the ODTFasta file and parse the headers using FastaParser
        fasta_parser = FastaParser()
        for idx in self.records_index:
            region_name, additional_info, _coordinates = (
                fasta_parser.parse_fasta_header(idx)
            )
            gene_id = region_name.lstrip(">")
            # only include sequences that match the specified region_types (if provided)
            if gene_id and (
                additional_info.get("type") in self.region_types
                or not self.region_types
            ):
                self.gene_records_map[gene_id].append(idx)

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # combine file hash and region_types to create a unique cache_id
            self.cache_id_str = (
                f"{_hash_file(self.odt_fasta_file_path)}_{'_'.join(self.region_types)}"
            )
        return self.cache_id_str

    def load_files(self):
        self.records_index = shelve.open(f"{self.odt_fasta_file_path}_regionsODT.shelve", flag="c", writeback=True)
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
                    {
                        # TODO: check if these are actually lists
                        "start": coordinates.get("start"),
                        "end": coordinates.get("end"),
                        "type": additional_info.get("type"),
                        "strand": additional_info.get("strand"),
                        "transcript_id": additional_info.get("transcript_id"),
                    }
                )

    def load_gene(self, gene_id: str):
        super().load_gene()  # ensure files are loaded
        regions = defaultdict(list)  # {transcript_id: [(start, end, type), ...]}
        for record_info in self.gene_records_map.get(gene_id, []):
            regions[record_info["transcript_id"]].append(
                {
                    "start": record_info["start"],
                    "end": record_info["end"],
                    "type": record_info["type"],
                    "strand": record_info["strand"],
                }
            )
        return regions

    def gene_list(self):
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
