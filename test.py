from gene_viewer import GeneViewer

# homo_sapiens: GCF_000001405.40_GRCh38.p14_genomic
# drosophila: GCF_000001215.4_Release_6_plus_ISO1_MT_genomic
# odt: exon_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38
# odt2: exon_exon_junction_annotation_source-NCBI_species-Homo_sapiens_annotation_release-110_genome_assemly-GRCh38

if __name__ == "__main__":
    viewer = GeneViewer("odt", "visualizations", genes_without_probes="visualize")
    viewer.load_regions_ODTFasta(
        "data/exon_annotation_source-NCBI.fna",
    )
    viewer.load_sequences_ODTFasta(
        "data/exon_annotation_source-NCBI.fna",
    )
    viewer.load_track_BED("data/hgTables.bed", track_name="tables")
    # viewer.load_track_GTF(
    #     "data/GCF_000001405.40_GRCh38.p14_genomic.gtf",
    #     "start_stop_codons",
    #     feature_types=["start_codon", "stop_codon"]
    # )
    # viewer.merge_exon_junctions()
    # viewer.fill_gaps_with_introns()
    # viewer.restrict_to_exon_sequences()
    viewer.save()
