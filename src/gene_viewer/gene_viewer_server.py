import json
from pathlib import Path

import zstandard as zstd


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
                return {"error": f"Viewer with ID {viewer_id} not found"}
            with open(gene_list_file_path, "r") as f:
                metadata = json.load(f)
            resolved_gene_id = (
                metadata["gene_list"][0] if metadata["gene_list"] else None
            )

        # load json file (self.dir_path / "visualizations" / f"{self.viewer_id}" / f"{gene_id}.json")
        gene_data_file_path = (
            self.dir_path
            / "visualizations"
            / f"{viewer_id}"
            / f"{resolved_gene_id}.json"
        )

        if not gene_data_file_path.exists():
            return {"error": f"Data for gene {resolved_gene_id} not found"}

        with open(gene_data_file_path, "r") as f:
            gene_data = json.load(f)

        # resolve references
        for data_type in gene_data:
            item = gene_data[data_type]
            if isinstance(item, dict) and "_ref" in item:
                # check if the reference file exists
                ref_dir_path = self.dir_path / f"{data_type}_cache" / item["_ref"]
                if ref_dir_path.exists():
                    with open(ref_dir_path / "_index.json", "r") as index_file:
                        index = json.load(index_file)
                        gene_info = index.get(resolved_gene_id)
                        gene_offset = gene_info["offset"]
                        gene_length = gene_info["length"]
                    with open(ref_dir_path / "data.blob", "rb") as blob_file:
                        blob_file.seek(gene_offset)
                        ref_data = zstd.ZstdDecompressor().decompress(
                            blob_file.read(gene_length)
                        )
                    ref_data = json.loads(ref_data)
                    gene_data[data_type] = ref_data

        if not gene_id:
            return {
                "geneList": metadata["gene_list"],
                "geneId": resolved_gene_id,
                "gene": gene_data,
            }

        return {"geneId": resolved_gene_id, "gene": gene_data}
