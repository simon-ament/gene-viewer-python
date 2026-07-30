from pathlib import Path
from collections import defaultdict
import json

class GeneViewer:
    def __init__(self, viewer_id: str, dir_path: str, genes_from: str = "probes"):
        self.viewer_id = viewer_id
        self.dir_path = Path(dir_path)
        self.genes_from = genes_from  # "probes", "regions", "sequences", "tracks", "all"
        self.gene_list = set()  # unique gene IDs from all sources
        self.regions_loaders = []
        self.sequences_loaders = []
        self.track_loaders = []
        self.probes = defaultdict(lambda: defaultdict(list)) # {gene_id: {probeset_id: [probe1, probe2, ...]}}
        self.processors = []

    # Regions Loaders

    def add_regions_loader(self, regions_loader):
        """Adds a regions loader to the GeneViewer."""
        self.regions_loaders.append(regions_loader)
        if (self.genes_from == "regions" or self.genes_from == "all"):
            self._append_gene_list(regions_loader.gene_list())

    def add_regions_GTF(self, gtf_file_path: str, region_types: list = []):
        from src.loader.regions_loader import RegionsLoaderGTF
        self.add_regions_loader(RegionsLoaderGTF(gtf_file_path, region_types))

    def add_regions_ODTFasta(self, odt_fasta_file_path: str, region_types: list = []):
        from src.loader.regions_loader import RegionsLoaderODTFasta
        self.add_regions_loader(RegionsLoaderODTFasta(odt_fasta_file_path, region_types))

    # Sequences Loaders

    def add_sequences_loader(self, sequences_loader):
        """Adds a sequences loader to the GeneViewer."""
        self.sequences_loaders.append(sequences_loader)
        if (self.genes_from == "sequences" or self.genes_from == "all"):
            self._append_gene_list(sequences_loader.gene_list())

    def add_sequences_FastaGTF(self, fasta_file_path: str, gtf_file_path: str, region_types: list = []):
        from src.loader.sequences_loader import SequencesLoaderFastaGTF
        self.add_sequences_loader(SequencesLoaderFastaGTF(fasta_file_path, gtf_file_path, region_types))

    def add_sequences_ODTFasta(self, odt_fasta_file_path: str, region_types: list = []):
        from src.loader.sequences_loader import SequencesLoaderODTFasta
        self.add_sequences_loader(SequencesLoaderODTFasta(odt_fasta_file_path, region_types))

    # Track Loaders

    def add_track_loader(self, track_loader):
        """Adds a track loader to the GeneViewer."""
        self.track_loaders.append(track_loader)
        if (self.genes_from == "tracks" or self.genes_from == "all"):
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
        self.probes[gene_id][probeset_id].append(probe)

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

        for gene_id in self.gene_list:
            gene_data, cached = self._get_gene_data_cached(gene_id)

            for dir, cached_file in cached.items():
                self._write_json(cached_file.data, self.dir_path / dir / f"{cached_file.filename}.json")

            self._write_json(gene_data, self.dir_path / "visualizations" / f"{self.viewer_id}" / f"{gene_id}.json")

    def save_raw(self):
        """Saves the raw gene data for all genes in the gene list to JSON files in the specified directory."""

        for gene_id in self.gene_list:
            gene_data = self._get_gene_data_raw(gene_id)

            self._write_json(gene_data, self.dir_path / "visualizations" / f"{self.viewer_id}" / f"{gene_id}.json")

    # Helpers

    def _append_gene_list(self, gene_ids):
        """Appends gene IDs to the gene list, ensuring uniqueness."""
        self.gene_list.update(gene_ids)

    def _get_gene_data_cached(self, gene_id: str):
        """Retrieves the gene data for a specific gene ID, using cached files if available."""

        # TODO: cached file names are bases on cache ids, generator ids and dependencies (recursively)
        # TODO: check for existing cached files and use them if available, otherwise generate new files
        pass

    def _get_gene_data_raw(self, gene_id: str):
        """Retrieves the raw gene data for a specific gene ID without caching."""
        gene_data = {
            "gene_id": gene_id,
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
            "probes": self.probes.get(gene_id, {}),
        }

        for processor in self.processors:
            gene_data = processor.process(gene_data)

        return gene_data

    def _write_json(self, data, file_path: Path):
        """Writes data to a JSON file at the specified file path."""
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
