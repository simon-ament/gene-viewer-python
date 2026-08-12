from abc import abstractmethod
from collections import defaultdict

from oligo_designer_toolsuite.utils import FastaParser

from helpers import _get_GTF_attribute, _hash_file
from src.index import GTFFileIndex, ODTFastaFileIndex
from src.types import GeneLocation

from .loader import Loader


class RegionsLoader(Loader):
    def __init__(self):
        super().__init__()
        self._gene_list = None

    def limit_to_genes(self, gene_list: list[str]):
        self._gene_list = gene_list

    @abstractmethod
    def gene_locations(self) -> dict[str, list[GeneLocation]]:
        if not self._lazy_init_done:
            self._lazy_init()
            self._lazy_init_done = True


class RegionsLoaderGTF(RegionsLoader):
    def __init__(self, gtf_file_path: str, region_types: list | None = None):
        super().__init__()
        self._gtf_file_path = gtf_file_path
        self._gtf_file_index: GTFFileIndex | None = None  # will be inizialized lazily
        self._region_types = region_types or []
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = (
                f"{_hash_file(self._gtf_file_path)}_{'_'.join(self._region_types)}"
            )
        return self._cache_id_str

    def _lazy_init(self):
        self._gtf_file_index = GTFFileIndex(
            self._gtf_file_path, self._gene_list, collect_gene_locations=True
        )

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        regions = defaultdict(list)  # {transcript_id: [(start, end, type), ...]}
        for feature in self._gtf_file_index.get(gene.id):
            type = feature[2]
            if self._region_types and type not in self._region_types:
                continue
            start = int(feature.start)  # 1-based start position
            end = int(feature.end)  # 1-based end position
            exon_number = _get_GTF_attribute(feature[8], "exon_number")
            region = {
                "start": start,
                "end": end,
                "type": type,
                "strand": feature.strand,
            }
            if exon_number is not None:
                region["exon_number"] = int(exon_number)
            regions[_get_GTF_attribute(feature[8], "transcript_id")].append(region)
        return regions

    @property
    def gene_locations(self):
        super().gene_locations()
        return self._gtf_file_index.gene_locations


class RegionsLoaderODTFasta(RegionsLoader):
    def __init__(
        self,
        odt_fasta_file_path: str,
        region_types: list | None = None,
    ):
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
        self._odt_fasta_file_index = ODTFastaFileIndex(
            self._odt_fasta_file_path, self._gene_list, collect_gene_locations=True
        )
        self._fasta_parser = FastaParser()

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        regions = defaultdict(list)  # {transcript_id: [(start, end, type), ...]}
        for header, sequence in self._odt_fasta_file_index.get(gene.id):
            _, additional_info, coordinates = self._fasta_parser.parse_fasta_header(
                header
            )
            type = additional_info.get("regiontype", ["unknown"])[0]
            if self._region_types and type not in self._region_types:
                continue
            for transcript_id in additional_info.get("transcript_id", ["unknown"]):
                regions[transcript_id].append(
                    {
                        "start": coordinates["start"][0],
                        "end": coordinates["end"][0],
                        "type": type,
                        "strand": additional_info["strand"][0],
                    }
                )
        return regions

    @property
    def gene_locations(self):
        super().gene_locations()
        return self._odt_fasta_file_index.gene_location_list
