from src.gene_viewer import GeneViewer

viewer = GeneViewer("test_viewer", "output", genes_from="all")
viewer.add_sequences_FastaGTF("data/GCF_009729015.1_ASM972901v1_genomic.fna", "data/GCF_009729015.1_ASM972901v1_genomic.gtf")
viewer.save_raw()
