import copy
import hashlib
import json
from pathlib import Path
from collections import defaultdict

from src.loader.probes_loader import ProbesLoaderManual


class GeneViewer:
    def __init__(self, viewer_id: str, dir_path: str, genes_from: str = "probes"):
        self.viewer_id = viewer_id
        self.dir_path = Path(dir_path)
        self.genes_from = (
            genes_from  # "probes", "regions", "sequences", "tracks", "all"
        )
        self.regions_loaders = []
        self.sequences_loaders = []
        self.track_loaders = []
        self.probes_loader = ProbesLoaderManual()  # default probe loader
        self.processors = []

        self.loaders = {
            "regions": self.regions_loaders,
            "sequences": self.sequences_loaders,
            "tracks": self.track_loaders,
            "probes": [self.probes_loader],
        }

    # Regions Loaders

    def add_regions_loader(self, regions_loader):
        """Adds a regions loader to the GeneViewer."""
        self.regions_loaders.append(regions_loader)

    def add_regions_GTF(self, gtf_file_path: str, region_types: list | None = None):
        from src.loader.regions_loader import RegionsLoaderGTF

        self.add_regions_loader(RegionsLoaderGTF(gtf_file_path, region_types or []))

    def add_regions_ODTFasta(
        self, odt_fasta_file_path: str, region_types: list | None = None
    ):
        from src.loader.regions_loader import RegionsLoaderODTFasta

        self.add_regions_loader(
            RegionsLoaderODTFasta(odt_fasta_file_path, region_types or [])
        )

    # Sequences Loaders

    def add_sequences_loader(self, sequences_loader):
        """Adds a sequences loader to the GeneViewer."""
        self.sequences_loaders.append(sequences_loader)

    def add_sequences_FastaGTF(
        self, fasta_file_path: str, gtf_file_path: str, region_types: list | None = None
    ):
        from src.loader.sequences_loader import SequencesLoaderFastaGTF

        self.add_sequences_loader(
            SequencesLoaderFastaGTF(fasta_file_path, gtf_file_path, region_types or [])
        )

    def add_sequences_ODTFasta(
        self, odt_fasta_file_path: str, region_types: list | None = None
    ):
        from src.loader.sequences_loader import SequencesLoaderODTFasta

        self.add_sequences_loader(
            SequencesLoaderODTFasta(odt_fasta_file_path, region_types or [])
        )

    # Track Loaders

    def add_track_loader(self, track_loader):
        """Adds a track loader to the GeneViewer."""
        self.track_loaders.append(track_loader)

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

    def save(self, store_genes_eagerly: bool = True):
        """Saves the gene data for all genes in the gene list to JSON files in the specified directory.

        Args:
            store_genes_eagerly (bool): If True, gene data is generated and stored for all genes in a loader's gene list.
            If False, gene data is generated and stored only for genes that will be visualized (as specified by the `genes_from` parameter).

        Cached files are used if available, otherwise new files are generated.
        """

        # TODO: switch all gene_lists to sets, implement the store_genes_eagerly option, and only generate json files when data exists

        genes_to_process, cache_dirs, gene_set_visualized, gene_set_by_output = (
            self._get_genes_to_process(store_genes_eagerly)
        )

        # generate and save gene data for all genes in the gene lists
        for output_types, (gene_set, required_inputs, required_processors) in genes_to_process.items():
            for gene_id in gene_set:
                gene_data = self._get_gene_data_cached(
                    gene_id, required_inputs, required_processors
                )

                for output_type in output_types:
                    cache_data = gene_data[output_type]
                    cache_file_path = cache_dirs[output_type] / f"{gene_id}.json"
                    self._write_json(cache_data, cache_file_path)

                # TODO: collision risk: "_ref"
                if gene_id in gene_set_visualized:
                    gene_data_refs = {
                        output_type: {"_ref": str(dir / f"{gene_id}.json")}
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
            {
                "gene_list": list(gene_set_visualized),
            },
            self.dir_path / "visualizations" / f"{self.viewer_id}" / "_metadata.json",
        )

        # update metadata for each output type
        for output_type in self.loaders.keys():
            # if store_genes_eagerly, store all genes in the loader's gene list, otherwise add genes that were newly processed
            if store_genes_eagerly:
                metadata = {
                    "genes_cached": list(gene_set_by_output[output_type]),
                    "genes_uncached": [],
                }
            else:
                # load existing metadata if it exists
                metadata_file_path = cache_dirs[output_type] / "_metadata.json"
                if metadata_file_path.exists():
                    with open(metadata_file_path, "r") as f:
                        metadata = json.load(f)
                else:
                    metadata = {"genes_cached": [], "genes_uncached": []}

                # update cached and uncached gene lists
                metadata["genes_cached"].extend(
                    list(gene_set_by_output[output_type])
                )
                metadata["genes_uncached"] = list(
                    set(metadata["genes_uncached"])
                    - set(gene_set_by_output[output_type])
                )
            
            self._write_json(
                metadata,
                cache_dirs[output_type] / "_metadata.json",
            )

    def save_raw(self, store_genes_eagerly: bool = True):
        """Saves the raw gene data for all genes in the gene list to JSON files in the specified directory."""

        gene_set = set()
        for loader in self.loaders[self.genes_from]:
            gene_set.update(loader.gene_list())
        gene_list = list(gene_set)

        for gene_id in gene_list:
            gene_data = self._get_gene_data_raw(gene_id)

            self._write_json(
                gene_data,
                self.dir_path
                / "visualizations"
                / f"{self.viewer_id}"
                / f"{gene_id}.json",
            )
            self._write_json(
                gene_list,
                self.dir_path / "visualizations" / f"{self.viewer_id}_gene_list.json",
            )

    # Helpers

    def _get_genes_to_process(self, store_genes_eagerly: bool = True):
        """Determines which genes need to be processed based on the loaders, processors, and caching status.
        
        Also returns other helpful information collected during the process, such as cache directories and the set of genes that will be visualized.
        """
        nested_dependencies = {}

        for data_type, loaders in self.loaders.items():
            nested_dependencies[data_type] = {
                "depends_on": {data_type},
                "processors_used": set(),
                "computation_path": tuple(loader.cache_id() for loader in loaders),
                "unprocessed_genes": set(),
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

        # determine the cache directories for each data type based on the computation path
        cache_dirs = {
            data_type: self.dir_path
            / f"{data_type}_cache"
            / f"{hashlib.sha256(str(deps['computation_path']).encode()).hexdigest()}"
            for data_type, deps in nested_dependencies.items()
        }

        gene_set_full: set[str] = set()
        gene_set_visualized: set[str] = set()
        gene_set_by_output = {data_type: set() for data_type in self.loaders.keys()}

        # load gene sets
        for output_type in self.loaders.keys():
            if cache_dirs[output_type].exists():
                # load cached gene set
                metadata = json.load(open(cache_dirs[output_type] / "_metadata.json"))
                gene_set_for_output = set(metadata["genes_cached"] + metadata["genes_uncached"])
                gene_set_full.update(gene_set_for_output)
                gene_set_by_output[output_type] = gene_set_for_output
                nested_dependencies[output_type]["unprocessed_genes"] = set(metadata["genes_uncached"])
                if self.genes_from == output_type or self.genes_from == "all":
                    gene_set_visualized.update(gene_set_for_output)
            else:
                # load uncached gene set
                gene_set_for_output = set()
                for loader in self.loaders[output_type]:
                    gene_set_for_output.update(loader.gene_list())
                gene_set_full.update(gene_set_for_output)
                gene_set_by_output[output_type] = gene_set_for_output
                nested_dependencies[output_type]["unprocessed_genes"] = gene_set_for_output
                if self.genes_from == output_type or self.genes_from == "all":
                    gene_set_visualized.update(gene_set_for_output)

        # if genes are not stored eagerly, only process genes that will be visualized
        if not store_genes_eagerly:
            for output_type in self.loaders.keys():
                nested_dependencies[output_type]["unprocessed_genes"] &= gene_set_visualized

        # Genes to Process
        # generate genes_to_process, which is a dictionary
        # its keys are tuples of output_types
        # its values are tuples containing a gene set, the required inputs and the required processors

        genes_to_process = defaultdict(lambda: (set(), set(), set()))

        for gene_id in gene_set_full:
            required_outputs = set()
            required_inputs = set()
            required_processors = set()

            for data_type, deps in nested_dependencies.items():
                if gene_id in deps["unprocessed_genes"]:
                    required_outputs.add(data_type)
                    required_inputs.update(deps["depends_on"])
                    required_processors.update(deps["processors_used"])

            required_outputs_sorted = tuple(sorted(required_outputs))
            genes_to_process[required_outputs_sorted][0].add(gene_id)
            genes_to_process[required_outputs_sorted][1].update(required_inputs)
            genes_to_process[required_outputs_sorted][2].update(
                required_processors
            )

        return genes_to_process, cache_dirs, gene_set_visualized, gene_set_by_output

    def _get_gene_data_cached(self, gene_id: str, required_inputs, required_processors):
        """Retrieves the gene data for a specific gene ID, skipping cached data if available, otherwise generating new data."""
        gene_data = {}

        for data_type in required_inputs:
            if data_type == "probes":
                gene_data[data_type] = self.probes_loader.load_gene(gene_id)
            else:
                gene_data[data_type] = [
                    item
                    for loader in self.loaders[data_type]
                    for item in loader.load_gene(gene_id)
                ]

        for processor in self.processors:
            if processor.id in required_processors:
                gene_data = processor.process(gene_data)

        return gene_data

    def _get_gene_data_raw(self, gene_id: str):
        """Retrieves the raw gene data for a specific gene ID without caching."""
        gene_data = {}

        for data_type in self.loaders.keys():
            if data_type == "probes":
                gene_data[data_type] = self.probes_loader.load_gene(gene_id)
            else:
                gene_data[data_type] = [
                    item
                    for loader in self.loaders[data_type]
                    for item in loader.load_gene(gene_id)
                ]

        for processor in self.processors:
            gene_data = processor.process(gene_data)

        return gene_data

    def _write_json(self, data, file_path: Path):
        """Writes data to a JSON file at the specified file path."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
