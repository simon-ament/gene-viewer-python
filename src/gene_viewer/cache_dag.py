import hashlib
import json
from collections import defaultdict
from pathlib import Path

import zstandard as zstd

from gene_viewer.processor import Processor
from gene_viewer.types import GeneLocation


def hash(string: str) -> str:
    return hashlib.sha256(string.encode()).hexdigest()


class ProcessorNode:
    # can process data from children
    _last_gene_id: str = None
    _last_gene_data: dict[str, any] = None

    def __init__(self, processor: Processor, children: list["DataNode | LoaderNode"]):
        self._processor = processor
        self._children = sorted(children, key=lambda x: x.type)

        self.computation_path = (
            self._processor.id
            + "("
            + ",".join([child.computation_path for child in self._children])
            + ")"
        )
        self.cache_id = hash(self.computation_path)

    def load_gene(self, gene: GeneLocation, loaders: dict[str, list]) -> dict[str, any]:
        if self._last_gene_id == gene.id:
            return self._last_gene_data
        # Load the gene data from the processor
        gene_data = {}
        for child in self._children:
            gene_data[child.type] = child.load_gene(gene, loaders)
        self._last_gene_id = gene.id
        self._last_gene_data = self._processor.process(gene_data)
        return self._last_gene_data


class CachedNode:
    cache_id: str = ""

    _last_index_id: str = None
    _last_index_data: any = None

    def __init__(self, dir_path: Path, type: str):
        self.dir_path = dir_path
        self.type = type

    def _load_cached_gene_ids(self) -> set[str]:
        cache_metadata_path = (
            self.dir_path / f"{self.type}_cache" / self.cache_id / "_metadata.json"
        )
        if cache_metadata_path.exists():
            with open(cache_metadata_path, "r") as f:
                metadata = json.load(f)
                return set(metadata.get("genes_cached", []))
        else:
            return set()

    def _load_cached_gene_data(
        self, gene_id: str, return_ref: bool = False
    ) -> dict[str, any]:
        if return_ref:
            return {"_ref": self.cache_id}

        ref_dir_path = self.dir_path / f"{self.type}_cache" / self.cache_id
        if ref_dir_path.exists():
            if self._last_index_id == self.cache_id:
                index = self._last_index_data
            else:
                with open(ref_dir_path / "_index.json", "r") as index_file:
                    index = json.load(index_file)
                    self._last_index_id = self.cache_id
                    self._last_index_data = index
            gene_info = index.get(gene_id)
            gene_offset = gene_info["offset"]
            gene_length = gene_info["length"]
            with open(ref_dir_path / "data.blob", "rb") as blob_file:
                blob_file.seek(gene_offset)
                cache_data = zstd.ZstdDecompressor().decompress(
                    blob_file.read(gene_length)
                )

            return json.loads(cache_data)
        else:
            return {}


class DataNode(CachedNode):
    # can be cached, otheriwse loads data from a child processor node
    _last_gene_id: str = None
    _last_gene_data: any = None

    def __init__(self, type: str, dir_path: Path, child: ProcessorNode):
        super().__init__(dir_path, type)
        self._child = child

        self.computation_path = self.type + ":" + self._child.computation_path
        self.cache_id = hash(self.computation_path)

        # Load the cached gene IDs from the metadata file if it exists
        self._cached_gene_ids = self._load_cached_gene_ids()

    def load_gene(
        self, gene: GeneLocation, loaders: dict[str, list], return_ref: bool = False
    ) -> any:
        if self._last_gene_id == gene.id:
            return self._last_gene_data

        if gene.id in self._cached_gene_ids:
            # Load the gene data from the cache
            self._last_gene_data = self._load_cached_gene_data(gene.id, return_ref)
        else:
            # Load the gene data from the child processor node
            self._last_gene_data = self._child.load_gene(gene, loaders)[self.type]

        self._last_gene_id = gene.id
        return self._last_gene_data


class LoaderNode(CachedNode):
    # can be cached, otheriwse loads data from loaders
    _last_gene_id: str = None
    _last_gene_data: any = None

    def __init__(
        self,
        type: str,
        dir_path: Path,
        loader_cache_ids: list[str],
        region_loader_cache_ids: list[str],
    ):
        super().__init__(dir_path, type)
        self._loader_cache_ids = sorted(loader_cache_ids)
        self._region_loader_cache_ids = sorted(region_loader_cache_ids)

        self.computation_path = (
            self.type
            + "_loaders("
            + ",".join(self._loader_cache_ids)
            + "|"
            + ",".join(self._region_loader_cache_ids)
            + ")"
        )
        self.cache_id = hash(self.computation_path)

        # Load the cached gene IDs from the metadata file if it exists
        self._cached_gene_ids = self._load_cached_gene_ids()

    def load_gene(
        self, gene: GeneLocation, loaders: dict[str, list], return_ref: bool = False
    ) -> any:
        if self._last_gene_id == gene.id:
            return self._last_gene_data

        if gene.id in self._cached_gene_ids:
            # Load the gene data from the cache
            self._last_gene_id = gene.id
            self._last_gene_data = self._load_cached_gene_data(gene.id, return_ref)
            return self._last_gene_data

        # Load the gene data from the loaders
        if self.type == "sequences":
            gene_data = [
                item
                for loader in loaders["sequences"]
                for item in loader.load_gene(gene)
            ]
        elif self.type == "tracks":
            gene_data = defaultdict(list)
            for track_name, track_loaders in loaders["tracks"].items():
                for loader in track_loaders:
                    loaded_data = loader.load_gene(gene)
                    for item in loaded_data:
                        gene_data[track_name].append(item)
        else:
            gene_data = defaultdict(list)
            for loader in loaders[self.type]:
                loaded_data = loader.load_gene(gene)
                for key, value in loaded_data.items():
                    gene_data[key].extend(value)
        self._last_gene_id = gene.id
        self._last_gene_data = gene_data
        return self._last_gene_data
