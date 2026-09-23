import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Literal

import zstandard as zstd

from gene_viewer.cache_dag import DataNode, LoaderNode, ProcessorNode
from gene_viewer.loader.probes_loader import ProbesLoaderManual
from gene_viewer.loader.regions_loader import RegionsLoader
from gene_viewer.loader.sequences_loader import SequencesLoader
from gene_viewer.loader.track_loader import TrackLoader
from gene_viewer.processor import Processor
from gene_viewer.types import GeneLocation


class GeneViewer:
    def __init__(
        self,
        viewer_id: str,
        dir_path: str,
        genes_without_probes: Literal["ignore", "cache", "visualize"] = "ignore",
        species: str | None = None,
        source: str | None = None,
    ):
        self.viewer_id = viewer_id
        self.dir_path = Path(dir_path)
        self.genes_without_probes = genes_without_probes
        self.species = species
        self.source = source

        self.regions_loaders: list[RegionsLoader] = []
        self.sequences_loaders: list[SequencesLoader] = []
        self.track_loaders: defaultdict[list[TrackLoader]] = defaultdict(list)
        self.probes_loader: ProbesLoaderManual = (
            ProbesLoaderManual()
        )  # default probe loader
        self.processors: list[Processor] = []

        self.loaders: dict[str, list] = {
            "regions": self.regions_loaders,
            "sequences": self.sequences_loaders,
            "tracks": self.track_loaders,
            "probes": [self.probes_loader],
        }

    # Regions Loaders

    def add_regions_loader(self, regions_loader: RegionsLoader):
        """Adds a regions loader to the GeneViewer."""
        self.regions_loaders.append(regions_loader)

    def load_regions_GTF(
        self,
        gtf_file_path: str,
        region_types: list[str] | None = None,
        gene_id_attribute: str = "gene_id",
    ):
        from gene_viewer.loader.regions_loader import RegionsLoaderGTF

        if region_types is None:
            region_types = ["intron", "exon"]

        self.add_regions_loader(
            RegionsLoaderGTF(gtf_file_path, region_types, gene_id_attribute)
        )

    def load_regions_ODTFasta(
        self, odt_fasta_file_path: str, region_types: list[str] | None = None
    ):
        from gene_viewer.loader.regions_loader import RegionsLoaderODTFasta

        if region_types is None:
            region_types = ["intron", "exon"]

        self.add_regions_loader(
            RegionsLoaderODTFasta(odt_fasta_file_path, region_types)
        )

    # Sequences Loaders

    def add_sequences_loader(self, sequences_loader: SequencesLoader):
        """Adds a sequences loader to the GeneViewer."""
        self.sequences_loaders.append(sequences_loader)

    def load_sequences_Fasta(self, fasta_file_path: str):
        from gene_viewer.loader.sequences_loader import SequencesLoaderFasta

        self.add_sequences_loader(SequencesLoaderFasta(fasta_file_path))

    def load_sequences_ODTFasta(
        self, odt_fasta_file_path: str, region_types: list[str] | None = None
    ):
        from gene_viewer.loader.sequences_loader import SequencesLoaderODTFasta

        if region_types is None:
            region_types = ["intron", "exon"]

        self.add_sequences_loader(
            SequencesLoaderODTFasta(odt_fasta_file_path, region_types)
        )

    # Track Loaders

    def add_track_loader(self, track_loader: TrackLoader, track_name: str):
        """Adds a track loader to the GeneViewer."""
        self.track_loaders[track_name].append(track_loader)

    def load_track_BED(
        self,
        bed_file_path: str,
        track_name: str,
        opacity_from_score: bool = False,
        max_score: float = 1000.0,
    ):
        from gene_viewer.loader.track_loader import TrackLoaderBED

        self.add_track_loader(
            TrackLoaderBED(bed_file_path, opacity_from_score, max_score), track_name
        )

    def load_track_GTF(
        self,
        gtf_file_path: str,
        track_name: str,
        feature_types: list[str] | None = None,
        opacity_from_score: bool = False,
        max_score: float = 1000.0,
        gene_id_attribute: str = "gene_id",
    ):
        from gene_viewer.loader.track_loader import TrackLoaderGTF

        self.add_track_loader(
            TrackLoaderGTF(
                gtf_file_path,
                feature_types,
                opacity_from_score,
                max_score,
                gene_id_attribute,
            ),
            track_name,
        )

    # Probes

    def add_probe(self, gene_id: str, probeset_id: str, probe):
        """Adds a probe to the GeneViewer."""
        self.probes_loader.add_probe(gene_id, probeset_id, probe)

    # Processors

    def add_processor(self, processor: Processor):
        """Adds a processor to the GeneViewer."""
        self.processors.append(processor)

    def merge_exon_junctions(self):
        from gene_viewer.processor import ProcessorExonJunctions

        self.add_processor(ProcessorExonJunctions())

    def fill_gaps_with_introns(self):
        from gene_viewer.processor import ProcessorIntronGaps

        self.add_processor(ProcessorIntronGaps())

    def restrict_to_exon_sequences(self):
        from gene_viewer.processor import ProcessorExonSequencesOnly

        self.add_processor(ProcessorExonSequencesOnly())

    # Save

    def save(self):
        ### 1. Build Cache DAG
        dag_roots = {}
        region_loader_cache_ids = [loader.cache_id for loader in self.regions_loaders]
        for output_type, loaders in self.loaders.items():
            if output_type == "tracks":
                cache_ids = [
                    f"{track_name}:{loader.cache_id}"
                    for track_name, loader in loaders.items()
                ]
            else:
                cache_ids = [loader.cache_id for loader in loaders]
            dag_roots[output_type] = LoaderNode(
                output_type, self.dir_path, cache_ids, region_loader_cache_ids
            )

        for processor in self.processors:
            processor_node = ProcessorNode(
                processor, [dag_roots[input_type] for input_type in processor.input]
            )
            for output_type in processor.output:
                dag_roots[output_type] = DataNode(
                    output_type, self.dir_path, processor_node
                )

        ### 2. Determine which genes need to be visualized and cached.
        if self.genes_without_probes != "visualize":
            genes_with_probes: set[str] = set(self.probes_loader.gene_list)

        regions_cache_metadata_path = (
            self.dir_path
            / "regions_cache"
            / f"{dag_roots['regions'].cache_id}"
            / "_metadata.json"
        )
        if regions_cache_metadata_path.exists():
            with open(regions_cache_metadata_path, "r") as f:
                regions_cache_metadata = json.load(f)
            all_gene_locations_raw: dict[str, list[dict]] = regions_cache_metadata.get(
                "gene_locations", {}
            )
            all_gene_locations = {
                seq_id: [
                    GeneLocation(**gene_location) for gene_location in gene_locations
                ]
                for seq_id, gene_locations in all_gene_locations_raw.items()
            }
        else:
            # collect gene locations from all region_loaders
            # from duplicate entries choose the first
            all_gene_locations = defaultdict(list)
            genes_seen = set()
            for loader in self.regions_loaders:
                for gene_locations in loader.gene_locations.values():
                    for gene_location in gene_locations:
                        if gene_location.id not in genes_seen:
                            genes_seen.add(gene_location.id)
                            all_gene_locations[gene_location.seq_id].append(
                                gene_location
                            )

        all_gene_locations_flattened = {}
        for gene_locations in all_gene_locations.values():
            for gene_location in gene_locations:
                all_gene_locations_flattened[gene_location.id] = gene_location
        all_gene_ids: set[str] = set(all_gene_locations_flattened.keys())

        if self.genes_without_probes == "ignore":
            gene_ids_to_process = genes_with_probes
            gene_ids_to_visualize = genes_with_probes
        elif self.genes_without_probes == "cache":
            gene_ids_to_process = all_gene_ids
            gene_ids_to_visualize = genes_with_probes
        elif self.genes_without_probes == "visualize":
            gene_ids_to_process = all_gene_ids
            gene_ids_to_visualize = all_gene_ids

        gene_locations_to_process = {
            seq_id: [
                gene_location
                for gene_location in gene_locations
                if gene_location.id in gene_ids_to_process
            ]
            for seq_id, gene_locations in all_gene_locations.items()
        }

        ### 3. Provide gene locations to loaders that require them.
        for loader in self.regions_loaders:
            loader.limit_to_genes(list(gene_ids_to_visualize))

        for loader in self.sequences_loaders:
            loader.set_gene_locations(gene_locations_to_process)

        for track_loaders in self.track_loaders.values():
            for loader in track_loaders:
                loader.set_gene_locations(gene_locations_to_process)

        ### 4. Evaluate Cache DAG and write data to disk

        blobs = {}
        indices = {}
        offsets = {}
        for output_type in self.loaders:
            blob_location = (
                self.dir_path
                / f"{output_type}_cache"
                / f"{dag_roots[output_type].cache_id}"
                / "data.blob"
            )
            blob_location.parent.mkdir(parents=True, exist_ok=True)
            blobs[output_type] = open(blob_location, "ab")
            indices[output_type] = {}
            offsets[output_type] = 0

        cctx = zstd.ZstdCompressor(level=3)

        for gene_id in gene_ids_to_process:
            if gene_id not in all_gene_locations_flattened:
                continue  # Skip genes that are not found in the gene locations

            gene_location = all_gene_locations_flattened[gene_id]
            gene_data = {}
            for output_type in self.loaders:
                gene_data[output_type] = dag_roots[output_type].load_gene(
                    gene_location, self.loaders, return_ref=True
                )

            if gene_id in gene_ids_to_visualize:
                gene_location = all_gene_locations_flattened[gene_id]
                gene_data_visualization = {
                    "id": gene_id,
                    "seq_id": gene_location.seq_id,
                    "start": gene_location.start,
                    "end": gene_location.end,
                    "strand": gene_location.strand,
                    "species": self.species,
                    "source": self.source,
                }
                for output_type in self.loaders:
                    gene_data_visualization[output_type] = {
                        "_ref": dag_roots[output_type].cache_id,
                    }

                self._write_json(
                    gene_data_visualization,
                    self.dir_path
                    / "visualizations"
                    / f"{self.viewer_id}"
                    / f"{gene_id}.json",
                )

            for output_type in self.loaders:
                cache_data = gene_data[output_type]
                if "_ref" in cache_data:
                    # If the data is a reference, we don't need to store it again
                    continue
                serialized_data = json.dumps(cache_data).encode("utf-8")
                compressed_data = cctx.compress(serialized_data)
                length = len(compressed_data)
                blobs[output_type].write(compressed_data)
                indices[output_type][gene_id] = {
                    "offset": offsets[output_type],
                    "length": length,
                }
                offsets[output_type] += length

        for output_type in self.loaders:
            blobs[output_type].close()
            index_location = (
                self.dir_path
                / f"{output_type}_cache"
                / f"{dag_roots[output_type].cache_id}"
                / "_index.json"
            )
            if index_location.exists():
                with open(index_location, "r") as f:
                    index_content = json.load(f)
                    index_content.update(indices[output_type])
            else:
                index_content = indices[output_type]
            self._write_json(index_content, index_location)

        ### 5. Save the gene list and metadata for the visualizations.

        for output_type in self.loaders:
            metadata_path = (
                self.dir_path
                / f"{output_type}_cache"
                / f"{dag_roots[output_type].cache_id}"
                / "_metadata.json"
            )
            if metadata_path.exists():
                with open(metadata_path, "r") as f:
                    metadata = json.load(f)
                genes_cached = set(metadata.get("genes_cached", []))
            else:
                genes_cached = set()

            metadata = {
                "genes_cached": list(gene_ids_to_process | genes_cached),
            }

            if output_type == "regions":
                metadata["gene_locations"] = {
                    seq_id: [asdict(gene_location) for gene_location in gene_locations]
                    for seq_id, gene_locations in all_gene_locations.items()
                }

            self._write_json(metadata, metadata_path)

        self._write_json(
            {
                "gene_list": list(gene_ids_to_visualize),
            },
            self.dir_path / "visualizations" / f"{self.viewer_id}" / "_metadata.json",
        )

    def save_raw(self):
        ### 1. Determine which genes need to be visualized.

        all_gene_locations = defaultdict(list)
        genes_seen = set()
        for loader in self.regions_loaders:
            for gene_locations in loader.gene_locations.values():
                for gene_location in gene_locations:
                    if gene_location.id not in genes_seen:
                        genes_seen.add(gene_location.id)
                        all_gene_locations[gene_location.seq_id].append(gene_location)

        all_gene_locations_flattened = {}
        for gene_locations in all_gene_locations.values():
            for gene_location in gene_locations:
                all_gene_locations_flattened[gene_location.id] = gene_location
        all_gene_ids: set[str] = set(all_gene_locations_flattened.keys())

        if self.genes_without_probes == "visualize":
            gene_ids_to_visualize: set[str] = all_gene_ids
        else:
            gene_ids_to_visualize: set[str] = set(self.probes_loader.gene_list)

        gene_locations_to_visualize = {
            seq_id: [
                gene_location
                for gene_location in gene_locations
                if gene_location.id in gene_ids_to_visualize
            ]
            for seq_id, gene_locations in all_gene_locations.items()
        }

        ### 2. Provide gene locations to loaders that require them.

        for loader in self.regions_loaders:
            loader.limit_to_genes(list(gene_ids_to_visualize))

        for loader in self.sequences_loaders:
            loader.set_gene_locations(gene_locations_to_visualize)

        for track_loaders in self.track_loaders.values():
            for loader in track_loaders:
                loader.set_gene_locations(gene_locations_to_visualize)

        ### 3. Generate and save raw gene data for all genes in the gene lists.

        for gene_id in gene_ids_to_visualize:
            gene_data = self._get_gene_data(all_gene_locations_flattened[gene_id])

            gene_location = all_gene_locations_flattened[gene_id]
            gene_data_visualization = {
                "id": gene_id,
                "seq_id": gene_location.seq_id,
                "start": gene_location.start,
                "end": gene_location.end,
                "strand": gene_location.strand,
                "species": self.species,
                "source": self.source,
            }
            for output_type in self.loaders:
                gene_data_visualization[output_type] = gene_data[output_type]

            self._write_json(
                gene_data_visualization,
                self.dir_path
                / "visualizations"
                / f"{self.viewer_id}"
                / f"{gene_id}.json",
            )

        self._write_json(
            {
                "gene_list": list(gene_ids_to_visualize),
            },
            self.dir_path / "visualizations" / f"{self.viewer_id}" / "_metadata.json",
        )

    # Helpers

    def _get_gene_data(self, gene: GeneLocation):
        """Retrieves the gene data for a specific gene ID, using only the required inputs and processors"""
        gene_data = {}

        for data_type in self.loaders:
            if data_type == "sequences":
                gene_data[data_type] = [
                    item
                    for loader in self.sequences_loaders
                    for item in loader.load_gene(gene)
                ]
            elif data_type == "tracks":
                gene_data[data_type] = defaultdict(list)
                for track_name, loaders in self.track_loaders.items():
                    for loader in loaders:
                        loaded_data = loader.load_gene(gene)
                        for item in loaded_data:
                            gene_data[data_type][track_name].append(item)
            else:
                gene_data[data_type] = defaultdict(list)
                for loader in self.loaders[data_type]:
                    loaded_data = loader.load_gene(gene)
                    for key, value in loaded_data.items():
                        gene_data[data_type][key].extend(value)

        for processor in self.processors:
            new_gene_data = processor.process(gene_data)
            for data_type in processor.output:
                gene_data[data_type] = new_gene_data[data_type]

        return gene_data

    def _write_json(self, data, file_path: Path):
        """Writes data to a JSON file at the specified file path."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
