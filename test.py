from src.gene_viewer import GeneViewer

viewer = GeneViewer("test_viewer-odt", "output", genes_from="all")
viewer.add_sequences_FastaGTF(
    "data/GCF_009729015.1_ASM972901v1_genomic.fna",
    "data/GCF_009729015.1_ASM972901v1_genomic.gtf",
)
# viewer.add_sequences_ODTFasta("data/exon_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38.fna")
viewer.save()
