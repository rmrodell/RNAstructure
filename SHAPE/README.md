# Preparing SHAPE Data

Start from a bam file labeled with a barcode from Nanopore basecalling outputs,
ideally with super-high-accuracy basecalling.

The pipeline (process_single_file.sh) goes through the following steps:
1. Convert bam to fastq, renaming from barcode in the process.
2. Trim Nanopore adapters from the fastq in sense and antisense directions. 
3. Reverse complement the antisense to unify orientation of reads.
4. Extract UMIs.
5. Trim the Pool adapters in the sense direction.
6. Map reads with bowtie2 & pool1 bowtie2 index.
7. Convert to sorted bam & index.
8. Deduplicate reads.
9. Convert to fastq for input into shapemapper

The following variables are hardcoded into process_single_file.sh:
- BOWTIE2_INDEX="/home/groups/nicolemm/rodell/pool1/pool1_bowtie_index/pool1_bowtie2"
- SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
- SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
- ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
- ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
- POOL_SENSE_ADAPTER_5PRIME="GACGCTCTTCCGATCT"
- POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
- UMI_PATTERN="NNNNNNNNNN"

To use process_single_file.sh, you need to submit it as a slurm array using submit_slurm_array.sh:

Example command:
```bash
bash $HOME/RNAstructure/SHAPE/submit_slurm_array.sh \
    --mail-user "rodell@stanford.edu" \
    --script-path $HOME/RNAstructure/SHAPE/process_single_file.sh \
    --sample-map /scratch/groups/nicolemm/rodell/SHAPE_InCell/barcodes.txt \
    --bam-dir /scratch/groups/nicolemm/rodell/SHAPE_InCell/merged_bams \
    --output-dir /scratch/groups/nicolemm/rodell/SHAPE_InCell/combined
```

The barcodes file should be formatted as follows:

HEK293T_SHAPE_Rep1:barcode75 

HEK293T_DMSO_Rep1:barcode76 

This will submit an individual job for each bam file in the bam directory.
The final deduplicated fastq will be located in a directory inside the output directory.
