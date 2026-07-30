from pathlib import Path

class GeneViewerServer:
    def __init__(self, dir_path: Path):
        self.dir_path = dir_path

    def serve(self, viewer_id: str, gene_id: str):
        # load json file
        # resolve references
        # special case: empty gene_id (provide list of available genes together with first gene)
        pass
