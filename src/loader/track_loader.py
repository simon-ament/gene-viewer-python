from helpers import _get_GTF_attribute, _hash_file
from src.index import BEDFileIndex, GTFFileIndex
from src.types import GeneLocation

from .loader import Loader


class TrackLoader(Loader):
    def __init__(self):
        super().__init__()
        self._gene_locations = None

    def set_gene_locations(self, gene_locations: dict[str, GeneLocation]):
        self._gene_locations = gene_locations


class TrackLoaderBED(TrackLoader):
    def __init__(
        self,
        bed_file_path: str,
        opacity_from_score: bool = False,
        max_score: float = 1.0,
    ):
        self._bed_file_path = bed_file_path
        self._bed_file_index = None  # will be initialized lazily
        self._opacity_from_score = opacity_from_score
        self._max_score = max_score
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = f"{_hash_file(self._bed_file_path)}_{self._opacity_from_score}_{self._max_score}"
        return self._cache_id_str

    def _lazy_init(self):
        if not self._gene_locations:
            raise ValueError("Gene locations must be set before loading tracks.")
        self._bed_file_index = BEDFileIndex(
            self._bed_file_path, gene_locations=self._gene_locations
        )

    def load_gene(self, gene: GeneLocation):
        track = []
        for feature in self._bed_file_index.get(gene.id):
            start = int(feature["start"])  # one-based start position
            end = int(feature["end"])  # one-based end position
            track.append(
                {
                    "start": start,
                    "end": end,
                    "type": feature[3],
                    "score": feature[4] / self._max_score
                    if self._opacity_from_score
                    else 1.0,
                    "item_rgb": feature[5],
                }
            )
        return track


class TrackLoaderGTF(TrackLoader):
    def __init__(
        self,
        gtf_file_path: str,
        opacity_from_score: bool = False,
        max_score: float = 1.0,
    ):
        self._gtf_file_path = gtf_file_path
        self._gtf_file_index = None  # will be initialized lazily
        self._opacity_from_score = opacity_from_score
        self._max_score = max_score
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = f"{_hash_file(self._gtf_file_path)}_{self._opacity_from_score}_{self._max_score}"
        return self._cache_id_str

    def _lazy_init(self):
        if not self._gene_locations:
            raise ValueError("Gene locations must be set before loading tracks.")
        gene_list = [gene.id for gene in self._gene_locations]
        self._gtf_file_index = GTFFileIndex(self._gtf_file_path, gene_list=gene_list)

    def load_gene(self, gene: GeneLocation):
        track = []
        for feature in self._gtf_file_index.get(gene.id):
            start = int(feature["start"])  # one-based start position
            end = int(feature["end"])  # one-based end position
            track.append(
                {
                    "start": start,
                    "end": end,
                    "type": feature[2],
                    "opacity": feature["score"] / self._max_score
                    if self._opacity_from_score
                    else 1.0,
                    "item_rgb": _get_GTF_attribute(feature[8], "item_rgb"),
                }
            )
        return track
