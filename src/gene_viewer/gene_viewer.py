import copy
import hashlib
import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Literal

import zstandard as zstd

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

        self.loaders = {
            "regions": self.regions_loaders,
            "sequences": self.sequences_loaders,
            "tracks": [
                loader for loaders in self.track_loaders.values() for loader in loaders
            ],
            "probes": [self.probes_loader],
        }

    # Regions Loaders

    def add_regions_loader(self, regions_loader: RegionsLoader):
        """Adds a regions loader to the GeneViewer."""
        self.regions_loaders.append(regions_loader)

    def load_regions_GTF(self, gtf_file_path: str, region_types: list[str] = ["intron", "exon"], gene_id_attribute: str = "gene_id"):
        from gene_viewer.loader.regions_loader import RegionsLoaderGTF

        self.add_regions_loader(RegionsLoaderGTF(gtf_file_path, region_types, gene_id_attribute))

    def load_regions_ODTFasta(
        self, odt_fasta_file_path: str, region_types: list[str] = ["intron", "exon"]
    ):
        from gene_viewer.loader.regions_loader import RegionsLoaderODTFasta

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
        self, odt_fasta_file_path: str, region_types: list[str] = ["intron", "exon"]
    ):
        from gene_viewer.loader.sequences_loader import SequencesLoaderODTFasta

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
            TrackLoaderGTF(gtf_file_path, feature_types, opacity_from_score, max_score, gene_id_attribute), track_name
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
        ### 1. Determine all dependencies between loaders and processors, and the required inputs for each output type.

        dependencies = {}
        regions_loaders_cache_ids = tuple(
            loader.cache_id for loader in self.regions_loaders
        )

        for data_type, loaders in self.loaders.items():
            dependencies[data_type] = {
                "depends_on": {data_type},
                "processors_used": set(),
                "computation_path": (
                    "base",
                    regions_loaders_cache_ids,
                    tuple(loader.cache_id for loader in loaders),
                ),
            }

        for processor in self.processors:
            # create a snapshot of the current nested dependencies to read from while we update the nested_dependencies dictionary
            deps_snapshot = copy.deepcopy(dependencies)

            for output_type in processor.output:
                # ensure deterministic iteration over inputs
                inputs_sorted = sorted(processor.input)
                dependencies[output_type]["depends_on"].update(
                    dep for dep in inputs_sorted
                )
                dependencies[output_type]["processors_used"].add(processor.id)
                dependencies[output_type]["computation_path"] = (processor.id,) + tuple(
                    deps_snapshot[dep] for dep in inputs_sorted
                )

        ### 2. Determine the cache directories for each data type based on the computation path

        cache_dirs = {
            data_type: Path(
                f"{hashlib.sha256(str(deps['computation_path']).encode()).hexdigest()}"
            )
            for data_type, deps in dependencies.items()
        }

        ### 3. Collect all genes that need to be visualized and cached.

        if self.genes_without_probes != "visualize":
            genes_with_probes: set[str] = set(self.probes_loader.gene_list)

        regions_cache_metadata_path = (
            self.dir_path
            / "regions_cache"
            / f"{cache_dirs['regions']}"
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
                for seq_id, gene_locations in loader.gene_locations.items():
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
            seq_id: [gene_location for gene_location in gene_locations]
            for seq_id, gene_locations in all_gene_locations.items()
            if any(
                gene_location.id in gene_ids_to_process
                for gene_location in gene_locations
            )
        }

        ### 4. Provide gene locations to loaders that require them.

        for loader in self.regions_loaders:
            loader.limit_to_genes(list(gene_ids_to_process))

        for loader in self.sequences_loaders:
            loader.set_gene_locations(gene_locations_to_process)

        for track_loaders in self.track_loaders.values():
            for loader in track_loaders:
                loader.set_gene_locations(gene_locations_to_process)

        ### 5. For every input type, determine which genes are already cached and which need to be processed.

        unprocessed_genes_by_output_type = {}

        for data_type in dependencies:
            cache_metadata_path = (
                self.dir_path
                / f"{data_type}_cache"
                / f"{cache_dirs[data_type]}"
                / "_metadata.json"
            )
            if cache_metadata_path.exists():
                with open(cache_metadata_path, "r") as f:
                    cache_metadata = json.load(f)
                cached_genes: set[str] = set(cache_metadata.get("genes_cached", []))
            else:
                cached_genes: set[str] = set()

            # Update the unprocessed genes for this data type based on the genes to process
            unprocessed_genes = gene_ids_to_process - cached_genes
            unprocessed_genes_by_output_type[data_type] = unprocessed_genes

        ### 6. Build a mapping of output type combinations to the genes, inputs, and processors required to generate them.

        genes_to_process = defaultdict(lambda: (set(), set(), set()))

        for gene_id in gene_ids_to_process:
            required_outputs = set()
            required_inputs = set()
            required_processors = set()

            for (
                data_type,
                unprocessed_genes,
            ) in unprocessed_genes_by_output_type.items():
                if gene_id in unprocessed_genes:
                    required_outputs.add(data_type)
                    required_inputs.update(dependencies[data_type]["depends_on"])
                    required_processors.update(
                        dependencies[data_type]["processors_used"]
                    )

            required_outputs_sorted = tuple(sorted(required_outputs))
            genes_to_process[required_outputs_sorted][0].add(gene_id)
            genes_to_process[required_outputs_sorted][1].update(required_inputs)
            genes_to_process[required_outputs_sorted][2].update(required_processors)

        ### 7. Generate and save gene data for all genes in the gene lists.

        for output_types, (
            gene_set,
            required_inputs,
            required_processors,
        ) in genes_to_process.items():
            for gene_id in gene_set:
                if gene_id not in all_gene_locations_flattened:
                    continue  # Skip genes that are not found in the gene locations
                
                gene_data = self._get_gene_data(
                    all_gene_locations_flattened[gene_id],
                    required_inputs,
                    required_processors,
                )

                for output_type in output_types:
                    cache_data = gene_data[output_type]
                    cache_file_path = (
                        self.dir_path
                        / f"{output_type}_cache"
                        / f"{cache_dirs[output_type]}"
                        / f"{gene_id}.json.zst"
                    )
                    self._write_zstd_json(cache_data, cache_file_path)

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
                    for output_type, cache_dir in cache_dirs.items():
                        gene_data_visualization[output_type] = {
                            "_ref": str(cache_dir / f"{gene_id}.json.zst")
                        }

                    self._write_json(
                        gene_data_visualization,
                        self.dir_path
                        / "visualizations"
                        / f"{self.viewer_id}"
                        / f"{gene_id}.json",
                    )

        ### 8. Save the gene list and metadata for the visualizations.

        for output_type in self.loaders:
            metadata_path = (
                self.dir_path
                / f"{output_type}_cache"
                / f"{cache_dirs[output_type]}"
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

        all_gene_locations = {
            seq_id: [gene_location for gene_location in gene_locations]
            for loader in self.regions_loaders
            for seq_id, gene_locations in loader.gene_locations.items()
        }

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
            seq_id: [gene_location for gene_location in gene_locations]
            for seq_id, gene_locations in all_gene_locations.items()
            if any(
                gene_location.id in gene_ids_to_visualize
                for gene_location in gene_locations
            )
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

    def _get_gene_data(
        self,
        gene: GeneLocation,
        required_inputs: set[str] | None = None,
        required_processors: set[str] | None = None,
    ):
        """Retrieves the gene data for a specific gene ID, using only the required inputs and processors"""
        gene_data = {}
        filtered_inputs = [
            data_type
            for data_type in self.loaders
            if required_inputs is None or data_type in required_inputs
        ]

        for data_type in filtered_inputs:
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
            if required_processors is None or processor.id in required_processors:
                new_gene_data = processor.process(gene_data)
                for data_type in processor.output:
                    gene_data[data_type] = new_gene_data[data_type]

        return gene_data

    def _write_json(self, data, file_path: Path):
        """Writes data to a JSON file at the specified file path."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)

    def _write_zstd_json(self, data, file_path: Path):
        """Writes data as a zstandard-compressed JSON file."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "wb") as f:
            cctx = zstd.ZstdCompressor()
            compressed = cctx.compress(json.dumps(data).encode("utf-8"))
            f.write(compressed)
        