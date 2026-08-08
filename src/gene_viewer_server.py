import json
from pathlib import Path


class GeneViewerServer:
    def __init__(self, dir_path: str):
        self.dir_path = Path(dir_path)

    def serve(self, viewer_id: str, gene_id: str):
        resolved_gene_id = gene_id

        # special case: empty gene_id (provide list of available genes together with first gene)
        if not resolved_gene_id:
            gene_list_file_path = (
                self.dir_path / "visualizations" / f"{viewer_id}" / "_metadata.json"
            )
            if not gene_list_file_path.exists():
                return { "error": f"Gene list metadata file not found for viewer_id: {viewer_id}" }
            with open(gene_list_file_path, "r") as f:
                metadata = json.load(f)
            resolved_gene_id = metadata["gene_list"][0] if metadata["gene_list"] else None

        # load json file (self.dir_path / "visualizations" / f"{self.viewer_id}" / f"{gene_id}.json")
        gene_data_file_path = (
            self.dir_path
            / "visualizations"
            / f"{viewer_id}"
            / f"{resolved_gene_id}.json"
        )

        if not gene_data_file_path.exists():
            return { "error": f"Gene data file not found for gene_id: {resolved_gene_id}" }

        with open(gene_data_file_path, "r") as f:
            gene_data = json.load(f)

        # resolve references
        for data_type in gene_data:
            item = gene_data[data_type]
            if isinstance(item, dict) and "_ref" in item:
                # check if the reference file exists
                ref_file_path = self.dir_path / f"{data_type}_cache" / item["_ref"]
                if ref_file_path.exists():
                    with open(ref_file_path, "r") as ref_f:
                        ref_data = json.load(ref_f)
                    gene_data[data_type] = ref_data

        if not gene_id:
            return { "geneList": metadata["gene_list"], "geneId": resolved_gene_id, "gene": gene_data }

        return { "geneId": resolved_gene_id, "gene": gene_data}
