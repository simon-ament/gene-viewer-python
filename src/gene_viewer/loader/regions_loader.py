from abc import abstractmethod
from collections import defaultdict

from oligo_designer_toolsuite.utils import FastaParser

from gene_viewer.helpers import _get_GTF_attribute, _hash_file
from gene_viewer.index import GTFFileIndex, ODTFastaFileIndex
from gene_viewer.types import GeneLocation

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
    def __init__(self, gtf_file_path: str, region_types: list[str], gene_id_attribute: str):
        super().__init__()
        self._gtf_file_path = gtf_file_path
        self._gtf_file_index: GTFFileIndex | None = None  # will be inizialized lazily
        self._region_types = region_types
        self._gene_id_attribute = gene_id_attribute
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
            self._gtf_file_path, self._gene_id_attribute, gene_list=self._gene_list, collect_gene_locations=True
        )

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        regions = defaultdict(list)  # {transcript_id: [(start, end, type), ...]}
        for feature in self._gtf_file_index.get(gene.id, default=[]):
            type = feature["feature"]
            if self._region_types and type not in self._region_types:
                continue
            start = int(feature["start"]) + 1  # 0-based -> 1-based start position
            end = int(feature["end"])  # 1-based end position
            exon_number = _get_GTF_attribute(feature["attributes"], "exon_number")
            region = {
                "start": start,
                "end": end,
                "type": type,
                "strand": feature["strand"],
            }
            if exon_number is not None:
                region["exon_number"] = int(exon_number)
            regions[_get_GTF_attribute(feature["attributes"], "transcript_id")].append(region)
        return regions

    @property
    def gene_locations(self):
        super().gene_locations()
        return self._gtf_file_index.gene_locations


class RegionsLoaderODTFasta(RegionsLoader):
    def __init__(
        self,
        odt_fasta_file_path: str,
        region_types: list[str]
    ):
        super().__init__()
        self._odt_fasta_file_path = odt_fasta_file_path
        self._odt_fasta_file_index = None  # will be initialized lazily
        self._region_types = region_types
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
        for header, sequence in self._odt_fasta_file_index.get(gene.id, default=[]):
            _, additional_info, coordinates = self._fasta_parser.parse_fasta_header(
                header
            )
            type = additional_info.get("regiontype", ["unknown"])[0]
            if self._region_types and type not in self._region_types:
                continue
            for transcript_index, transcript_id in enumerate(additional_info.get("transcript_id", ["unknown"])):
                subregions_indices = [0, 1] if type == "exonexonjunction" else [0]
                for idx in subregions_indices:
                    region = {
                        "start": coordinates["start"][idx],
                        "end": coordinates["end"][idx],
                        "type": type,
                        "strand": coordinates["strand"][idx],
                    }
                    exon_number = additional_info.get("exon_number", [None])[transcript_index]
                    if exon_number is not None and type != "exonexonjunction":
                        region["exon_number"] = int(exon_number)
                    if type == "exonexonjunction":
                        region["exon_number"] = int(exon_number.split("__JUNC__")[idx])
                    regions[transcript_id].append(region)
        return regions

    @property
    def gene_locations(self):
        super().gene_locations()
        return self._odt_fasta_file_index.gene_locations
