import copy
import hashlib
import json
from pathlib import Path

from src.loader.probes_loader import ProbesLoaderManual


class GeneViewer:
    def __init__(self, viewer_id: str, dir_path: str, genes_from: str = "probes"):
        self.viewer_id = viewer_id
        self.dir_path = Path(dir_path)
        self.genes_from = (
            genes_from  # "probes", "regions", "sequences", "tracks", "all"
        )
        self.gene_list = set()  # unique gene IDs from all sources
        self.regions_loaders = []
        self.sequences_loaders = []
        self.track_loaders = []
        self.probes_loader = ProbesLoaderManual()  # default probe loader
        self.processors = []

    # Regions Loaders

    def add_regions_loader(self, regions_loader):
        """Adds a regions loader to the GeneViewer."""
        self.regions_loaders.append(regions_loader)
        if self.genes_from == "regions" or self.genes_from == "all":
            self._append_gene_list(regions_loader.gene_list())

    def add_regions_GTF(self, gtf_file_path: str, region_types: list | None = None):
        from src.loader.regions_loader import RegionsLoaderGTF

        self.add_regions_loader(RegionsLoaderGTF(gtf_file_path, region_types or []))

    def add_regions_ODTFasta(self, odt_fasta_file_path: str, region_types: list | None = None):
        from src.loader.regions_loader import RegionsLoaderODTFasta

        self.add_regions_loader(
            RegionsLoaderODTFasta(odt_fasta_file_path, region_types or [])
        )

    # Sequences Loaders

    def add_sequences_loader(self, sequences_loader):
        """Adds a sequences loader to the GeneViewer."""
        self.sequences_loaders.append(sequences_loader)
        if self.genes_from == "sequences" or self.genes_from == "all":
            self._append_gene_list(sequences_loader.gene_list())

    def add_sequences_FastaGTF(
        self, fasta_file_path: str, gtf_file_path: str, region_types: list | None = None
    ):
        from src.loader.sequences_loader import SequencesLoaderFastaGTF

        self.add_sequences_loader(
            SequencesLoaderFastaGTF(fasta_file_path, gtf_file_path, region_types or [])
        )

    def add_sequences_ODTFasta(self, odt_fasta_file_path: str, region_types: list | None = None):
        from src.loader.sequences_loader import SequencesLoaderODTFasta

        self.add_sequences_loader(
            SequencesLoaderODTFasta(odt_fasta_file_path, region_types or [])
        )

    # Track Loaders

    def add_track_loader(self, track_loader):
        """Adds a track loader to the GeneViewer."""
        self.track_loaders.append(track_loader)
        if self.genes_from == "tracks" or self.genes_from == "all":
            self._append_gene_list(track_loader.gene_list())

    def add_track_BED(self, bed_file_path: str):
        from src.loader.track_loader import TrackLoaderBED

        self.add_track_loader(TrackLoaderBED(bed_file_path))

    def add_track_GTF(self, gtf_file_path: str):
        from src.loader.track_loader import TrackLoaderGTF

        self.add_track_loader(TrackLoaderGTF(gtf_file_path))

    # Probes

    def add_probe(self, gene_id: str, probeset_id: str, probe):
        """Adds a probe to the GeneViewer."""
        self.probes_loader.add_probe(gene_id, probeset_id, probe)
        if self.genes_from == "probes" or self.genes_from == "all":
            self._append_gene_list([gene_id])

    # Processors

    def add_processor(self, processor):
        """Adds a processor to the GeneViewer."""
        self.processors.append(processor)

    def merge_exon_junctions(self):
        from src.processor import ProcessorExonJunctions

        self.add_processor(ProcessorExonJunctions())

    def fill_gaps_with_introns(self):
        from src.processor import ProcessorIntronGaps

        self.add_processor(ProcessorIntronGaps())

    def restrict_to_exon_sequences(self):
        from src.processor import ProcessorExonSequencesOnly

        self.add_processor(ProcessorExonSequencesOnly())

    # Save

    def save(self):
        """Saves the gene data for all genes in the gene list to JSON files in the specified directory.

        Cached files are used if available, otherwise new files are generated.
        """

        cache_dirs, required_outputs, required_inputs, required_processors = (
            self._get_cache_dirs()
        )

        for gene_id in self.gene_list:
            gene_data = self._get_gene_data_cached(
                gene_id, required_inputs, required_processors
            )

            for output_type in required_outputs:
                cache_data = gene_data[output_type]
                cache_file_path = cache_dirs[output_type] / f"{gene_id}.json"
                self._write_json(cache_data, cache_file_path)

            # TODO: collision risk: "ref"
            gene_data_refs = {
                output_type: {"ref": str(dir / f"{gene_id}.json")}
                for output_type, dir in cache_dirs.items()
            }
            self._write_json(
                gene_data_refs,
                self.dir_path
                / "visualizations"
                / f"{self.viewer_id}"
                / f"{gene_id}.json",
            )
            self._write_json(
                list(self.gene_list),
                self.dir_path / "visualizations" / f"{self.viewer_id}_gene_list.json",
            )

    def save_raw(self):
        """Saves the raw gene data for all genes in the gene list to JSON files in the specified directory."""

        for gene_id in self.gene_list:
            gene_data = self._get_gene_data_raw(gene_id)

            self._write_json(
                gene_data,
                self.dir_path
                / "visualizations"
                / f"{self.viewer_id}"
                / f"{gene_id}.json",
            )
            self._write_json(
                list(self.gene_list),
                self.dir_path / "visualizations" / f"{self.viewer_id}_gene_list.json",
            )

    # Helpers

    def _append_gene_list(self, gene_ids):
        """Appends gene IDs to the gene list, ensuring uniqueness."""
        self.gene_list.update(gene_ids)

    def _get_cache_dirs(self):
        """Determines the cache directories for each data type based on the loaders and processors used."""

        nested_dependencies = {
            "regions": {
                "depends_on": {"regions"},
                "processors_used": set(),  # TODO: same processor twice?
                "computation_path": tuple(
                    loader.cache_id() for loader in self.regions_loaders
                ),
                "is_cached": False,
            },
            "sequences": {
                "depends_on": {"sequences"},
                "processors_used": set(),
                "computation_path": tuple(
                    loader.cache_id() for loader in self.sequences_loaders
                ),
                "is_cached": False,
            },
            "tracks": {
                "depends_on": {"tracks"},
                "processors_used": set(),
                "computation_path": tuple(
                    loader.cache_id() for loader in self.track_loaders
                ),
                "is_cached": False,
            },
            "probes": {
                "depends_on": {"probes"},
                "processors_used": set(),
                "computation_path": (self.probes_loader.cache_id(),),
                "is_cached": False,
            },
        }

        for processor in self.processors:
            # create a snapshot of the current nested dependencies to read from while we update the nested_dependencies dictionary
            deps_snapshot = copy.deepcopy(nested_dependencies)

            for output_type in processor.output:
                # ensure deterministic iteration over inputs
                inputs_sorted = sorted(processor.input)

                nested_dependencies[output_type]["depends_on"].update(
                    dep for dep in inputs_sorted
                )
                nested_dependencies[output_type]["processors_used"].add(processor.id)
                nested_dependencies[output_type]["computation_path"] = (
                    processor.id,
                ) + tuple(deps_snapshot[dep] for dep in inputs_sorted)

        # check for existing cached files and use them if available, otherwise generate new files
        cache_dirs = {
            data_type: self.dir_path
            / f"{data_type}_cache"
            / f"{hashlib.sha256(str(deps['computation_path']).encode()).hexdigest()}"
            for data_type, deps in nested_dependencies.items()
        }
        for data_type, deps in nested_dependencies.items():
            if cache_dirs[data_type].exists():
                deps["is_cached"] = True

        required_outputs = set()
        required_inputs = set()
        required_processors = set()
        for data_type, deps in nested_dependencies.items():
            if not deps["is_cached"]:
                required_outputs.add(data_type)
                required_inputs.update(deps["depends_on"])
                required_processors.update(deps["processors_used"])

        return cache_dirs, required_outputs, required_inputs, required_processors

    def _get_gene_data_cached(self, gene_id: str, required_inputs, required_processors):
        """Retrieves the gene data for a specific gene ID, skipping cached data if available, otherwise generating new data."""
        gene_data = {}

        for data_type in required_inputs:
            if data_type == "probes":
                gene_data[data_type] = self.probes_loader.load_gene(gene_id)
            else:
                loaders = {
                    "regions": self.regions_loaders,
                    "sequences": self.sequences_loaders,
                    "tracks": self.track_loaders,
                }
                gene_data[data_type] = [
                    item
                    for loader in loaders[data_type]
                    for item in loader.load_gene(gene_id)
                ]

        for processor in self.processors:
            if processor.id in required_processors:
                gene_data = processor.process(gene_data)

        return gene_data

    def _get_gene_data_raw(self, gene_id: str):
        """Retrieves the raw gene data for a specific gene ID without caching."""
        gene_data = {
            "regions": [
                region
                for regions_loader in self.regions_loaders
                for region in regions_loader.load_gene(gene_id)
            ],
            "sequences": [
                sequence
                for sequences_loader in self.sequences_loaders
                for sequence in sequences_loader.load_gene(gene_id)
            ],
            "tracks": [
                track
                for track_loader in self.track_loaders
                for track in track_loader.load_gene(gene_id)
            ],
            "probes": self.probes_loader.load_gene(gene_id),
        }

        for processor in self.processors:
            gene_data = processor.process(gene_data)

        return gene_data

    def _write_json(self, data, file_path: Path):
        """Writes data to a JSON file at the specified file path."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
