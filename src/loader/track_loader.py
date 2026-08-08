from helpers import _get_feature_attribute
from collections import defaultdict
import shelve

from pybedtools import BedTool

from helpers import _hash_file
from .loader import FileBasedLoader


class TrackLoader(FileBasedLoader):
    pass


class TrackLoaderBED(TrackLoader):
    def __init__(self, bed_file_path: str, track_name: str | None = None):
        super().__init__()
        self.bed_file_path = bed_file_path
        self.track_name = track_name
        self.cache_id_str = None

        self.features_index = None  # will be initialized in load_files()
        self.features = None  # will be initialized in load_files()

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # use file hash to create a unique cache_id
            self.cache_id_str = f"{_hash_file(self.bed_file_path)}"
        return self.cache_id_str

    def load_files(self):
        # create a persistent index for BED features and keep it open for later access
        self.features_index = shelve.open(
            f"{self.bed_file_path}_tracksBED.shelve", flag="c", writeback=True
        )

        # Load BED file
        self.features = BedTool(self.bed_file_path)
        # TODO

    def load_gene(self, gene_id: str):
        pass

    def gene_list(self):
        pass

    def delete(self):
        features_index = getattr(self, "features_index", None)
        if features_index is not None:
            try:
                features_index.close()
            except Exception:
                pass
            finally:
                self.features_index = None


class TrackLoaderGTF(TrackLoader):
    def __init__(self, gtf_file_path: str, track_name: str | None = None):
        super().__init__()
        self.gtf_file_path = gtf_file_path
        self.track_name = track_name
        self.cache_id_str = None

        self.gene_features_map = defaultdict(list)  # {gene_id: [feature, feature, ...]}
        self.features_index = None  # will be initialized in load_files()
        self.features = None  # will be initialized in load_files()

    def cache_id(self):
        # only compute cache_id_str when requested and if not already set
        if self.cache_id_str is None:
            # use file hash to create a unique cache_id
            self.cache_id_str = f"{_hash_file(self.gtf_file_path)}"
        return self.cache_id_str

    def load_files(self):
        # create a persistent index for GTF features and keep it open for later access
        self.features_index = shelve.open(
            f"{self.gtf_file_path}_tracksGTF.shelve", flag="c", writeback=True
        )

        # Load GTF file
        self.features = BedTool(self.gtf_file_path)
        for feature_idx, feature in enumerate(self.features):
            gene_id = _get_feature_attribute(
                feature[8], "gene_id"
            )  # .attrs fails parsing some GTF files
            # gene_id = feature.attrs.get("gene_id")
            if gene_id:
                self.gene_features_map[gene_id].append(feature_idx)
                # TODO: handle lists
                self.features_index[str(feature_idx)] = {
                    "start": feature.start,
                    "end": feature.end,
                    "type": feature[2],
                    "score": feature.score,
                    "item_rgb": _get_feature_attribute(feature[8], "item_rgb"),
                }

    def load_gene(self, gene_id: str):
        super().load_gene()  # ensure files are loaded
        track = {self.track_name: []}  # {track_name: [(start, end, type), ...]}
        for feature_idx in self.gene_features_map.get(gene_id, []):
            feature = self.features_index.get(
                str(feature_idx)
            )  # retrieve the feature from the persistent index
            start = int(feature["start"])  # one-based start position
            end = int(feature["end"])  # one-based end position
            track[self.track_name].append(
                {"start": start, "end": end, "type": feature["type"], "score": feature["score"], "item_rgb": feature["item_rgb"]}
            )
        return track

    def gene_list(self):
        super().load_files()  # ensure files are loaded
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
