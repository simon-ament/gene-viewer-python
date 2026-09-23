import random

from src.gene_viewer.index import FastaFileIndex
from src.gene_viewer.types import GeneLocation

if __name__ == "__main__":
    gene_location = GeneLocation(
        id="Dmel_CG3082",
        seq_id="NT_033778.4",
        start=23071552,
        end=23071555,
        strand="+",
    )
    gene_locations = {gene_location.id: gene_location}
    fastaIndex = FastaFileIndex(
        "data/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        gene_locations=gene_locations,
    )

    fasta_keys = fastaIndex.keys()
    random.shuffle(fasta_keys)
    for key in fasta_keys:
        print(f"Fetching sequence for gene {key}...")
        print(fastaIndex.get(key))
        fastaIndex.get(key)

    # gtfIndex = GTFFileIndex("data/GCF_009729015.1_ASM972901v1_genomic.gtf")
    # gtf_keys = gtfIndex.keys()
    # random.shuffle(gtf_keys)
    # for key in gtf_keys:
    #     gtfIndex.get(key)
