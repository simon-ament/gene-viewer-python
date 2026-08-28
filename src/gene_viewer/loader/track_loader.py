from gene_viewer.helpers import _get_GTF_attribute, _hash_file
from gene_viewer.index import BEDFileIndex, GTFFileIndex
from gene_viewer.types import GeneLocation

from .loader import Loader


class TrackLoader(Loader):
    def __init__(self):
        super().__init__()
        self._gene_locations: dict[str, list[GeneLocation]] | None = None

    def set_gene_locations(self, gene_locations: dict[str, list[GeneLocation]]):
        self._gene_locations = gene_locations


class TrackLoaderBED(TrackLoader):
    def __init__(
        self,
        bed_file_path: str,
        opacity_from_score: bool,
        max_score: float,
    ):
        super().__init__()
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
        super().load_gene(gene)
        track = []
        for feature in self._bed_file_index.get(gene.id, default=[]):
            start = feature.start + 1  # 0-based -> 1-based start position
            end = feature.end  # 1-based end position
            track.append(
                {
                    "start": start,
                    "end": end,
                    "type": feature[3],
                    "score": feature[4] / self._max_score
                    if self._opacity_from_score
                    else 1.0,
                    "item_rgb": feature[8],
                }
            )
        return track


class TrackLoaderGTF(TrackLoader):
    def __init__(
        self,
        gtf_file_path: str,
        feature_types: list[str] | None,
        opacity_from_score: bool,
        max_score: float,
        gene_id_attribute: str
    ):
        super().__init__()
        self._gtf_file_path = gtf_file_path
        self._gtf_file_index = None  # will be initialized lazily
        self._feature_types = feature_types
        self._opacity_from_score = opacity_from_score
        self._max_score = max_score
        self._gene_id_attribute = gene_id_attribute
        self._cache_id_str = None

    @property
    def cache_id(self):
        if self._cache_id_str is None:
            self._cache_id_str = f"{_hash_file(self._gtf_file_path)}_{self._opacity_from_score}_{self._max_score}"
        return self._cache_id_str

    def _lazy_init(self):
        if not self._gene_locations:
            raise ValueError("Gene locations must be set before loading tracks.")
        gene_list = [gene.id for gene_list in self._gene_locations.values() for gene in gene_list]
        self._gtf_file_index = GTFFileIndex(self._gtf_file_path, self._gene_id_attribute, gene_list=gene_list, collect_gene_locations=True)

    def load_gene(self, gene: GeneLocation):
        super().load_gene(gene)
        track = []
        for feature in self._gtf_file_index.get(gene.id, default=[]):
            type_ = feature["feature"]
            if self._feature_types is not None and type_ not in self._feature_types:
                continue
            start = feature["start"] + 1  # 0-based -> 1-based start position
            end = feature["end"]  # 1-based end position
            track.append(
                {
                    "start": start,
                    "end": end,
                    "type": feature["feature"],
                    "opacity": feature["score"] / self._max_score
                    if self._opacity_from_score
                    else 1.0,
                    "item_rgb": _get_GTF_attribute(feature["attributes"], "item_rgb"),
                }
            )
        return track
