import json
from pathlib import Path


class GeneViewerServer:
    def __init__(self, dir_path: Path):
        self.dir_path = dir_path

    def serve(self, viewer_id: str, gene_id: str):
        resolved_gene_id = gene_id

        # special case: empty gene_id (provide list of available genes together with first gene)
        if not resolved_gene_id:
            gene_list_file_path = (
                self.dir_path / "visualizations" / f"{viewer_id}_gene_list.json"
            )
            with open(gene_list_file_path, "r") as f:
                gene_list = json.load(f)
            resolved_gene_id = gene_list[0] if gene_list else None

        # load json file (self.dir_path / "visualizations" / f"{self.viewer_id}" / f"{gene_id}.json")
        gene_data_file_path = (
            self.dir_path
            / "visualizations"
            / f"{viewer_id}"
            / f"{resolved_gene_id}.json"
        )
        with open(gene_data_file_path, "r") as f:
            gene_data = json.load(f)

        # resolve references
        for key in ["regions", "sequences", "tracks"]:
            for item in gene_data[key]:
                if isinstance(item, dict) and "_ref" in item:
                    ref_file_path = (
                        self.dir_path
                        / "visualizations"
                        / f"{key}"
                        / f"{item['_ref']}.json"
                    )
                    with open(ref_file_path, "r") as ref_f:
                        ref_data = json.load(ref_f)
                    item.update(ref_data)
                    del item["_ref"]

        if not gene_id:
            return {"gene_list": gene_list, "data": gene_data}

        return gene_data
