# RNAstructure
Code for analysis of RNA structure

## RNAfold
To gather information of pairing probabilities and folding energies, Vienna RNA 2.5.1 is used. The following scripts are optimized to use on Sherlock, Stanford's HPC.

Use as follows:

```
bash run_RNAfold.sh <fasta_file> <output_dir>
```

This runs RNAfold on every sequence in the given fasta and outputs all files to the specified directory.

To concatenate these to a single csv file more suitable for downstream analysis, use `extract_pairing_prob.R` and `fold_summary.R` as described below.

```
Rscript extract_pairing_prob.R \
    -i <input directory> \
    -o <output csv> \
    -l <sequence length>
```

The above script extracts the pairing probabilites from the "ubox" lines in all dp.ps files output by RNAfold. This represents the total of all pairing probabilities for a given sequence, NOT the pairing probabilites of the MFE structures. The sequence length argument dictates how many columns will be created and should be equivalent to the length of the sequence. To extract the pairing probabilities of the MFE structure, use the script `RNAFold_MFE.R`.

```
Rscript fold_summary.R \
    -i <input directory> \
    -o <output csv>
```

The above script extracts information from all of the .fold files produced by RNAfold. This includes the sequence and dot-bracket notation for the MFE, MEA, and centroid structures, as well as the free energy (delta G) for those respective structure.

## Nano-SHAPE-Amp Pipeline (SHAPE subfolder)

### Mapping Reads: Nanopore Pipeline
This pipeline is optimized to run on Stanford's HPC Sherlock system. Make sure you adjust have the following installed and adjust package loading appropriately:

- samtools 1.16.1
- cutadapt 1.18 (requires Python 3.6)
- java 11
- bowtie2 2.3.4.1
- UMICollapse: https://github.com/Daniel-Liu-c0deb0t/UMICollapse

Start from a bam file labeled with a barcode from Nanopore basecalling outputs, ideally with super-high-accuracy basecalling since this is SHAPE and expected mutation rates are very low.

The pipeline `SHAPE/SHAPE_trim_map_dedup_mpra.sh` goes through the following steps:

1. Convert bam to fastq, renaming from barcode in the process.
2. Trim Nanopore adapters from the fastq in sense and antisense directions.
3. Reverse complement the antisense to unify orientation of reads.
4. Extract UMIs.
5. Trim the MPRA adapters in the sense direction.
6. Map reads with bowtie2.
7. Convert to sorted bam & index.
8. Deduplicate reads.
9. Convert to fastq for input into shapemapper

The following variables are hardcoded into `SHAPE/process_single_file.sh`:

- SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
- SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
- ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
- ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
- POOL_SENSE_ADAPTER_5PRIME="GACGCTCTTCCGATCT"
- POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
- UMI_PATTERN="NNNNNNNNNN"

To use, you need to submit it as a slurm array using `SHAPE/submit_slurm_array.sh`.

Example command:

```
bash SHAPE/submit_slurm_array.sh \
    --mail-user you@institution.edu \
    --script-path SHAPE/SHAPE_trim_map_dedup_mpra.sh \
    --map-file <barcodes.txt> \
    --input-dir <basecalling output directory> \
    --output-dir <output dir> \
    --bowtie-index <bowtie index>
```

The barcodes file should be formatted as follows:

```
HEK293T_SHAPE_Rep1:barcode75 
HEK293T_DMSO_Rep1:barcode76 
```

This will submit an individual job for each bam file in the bam directory. The final deduplicated fastq will be located in a directory inside the output directory.

**Note**: This mapping pipeline is NOT optimized for a mutagenesis pool with closely related sequences, even with barcodes. That is an outstanding problem that I am still looking to solve.

### Mapping Reads: Illumina Pipeline

The Illumina data starts in a slightly different place, namely pair-ended reads, than the Nanopore data, so it requires a slightly different pipeline.

The same packages are still required:

- samtools 1.16.1
- cutadapt 1.18 (requires Python 3.6)
- java 11
- bowtie2 2.3.4.1
- UMICollapse: https://github.com/Daniel-Liu-c0deb0t/UMICollapse

The pipeline `SHAPE/process_single_file_illumina_downsample.sh` goes through the following steps:

1. Starts with pair-ended .fastq.gz files.
2. Extract UMIs from the start of R2 and appends it to R1 headers.
3. Trims 3' adapters from R1.
4. Downsamples to ~20M reads to approximately match the read counts in the Nanopore sample.
5. Map reads with bowtie2.
6. Convert to sorted bam & index.
7. Deduplicate reads.
8. Convert to fastq for input into shapemapper

The following variables are hardcoded into process_single_file.sh:

- SENSE_ADAPTER_5PRIME="TTTCTGTTGGTGCTGATATTGCG"
- SENSE_ADAPTER_3PRIME="GAAGATAGAGCGACAGGCAAGT"
- ANTISENSE_ADAPTER_5PRIME="ACTTGCCTGTCGCTCTATCTTC"
- ANTISENSE_ADAPTER_3PRIME="CGCAATATCAGCACCAACAGAAA"
- POOL_SENSE_ADAPTER_5PRIME="GACGCTCTTCCGATCT"
- POOL_SENSE_ADAPTER_3PRIME="CACTCGGGCACCAAGGAC"
- UMI_PATTERN="NNNNNNNNNN"

To use, you need to submit it as a slurm array using `SHAPE/submit_slurm_array_illumina.sh`

Example command:

```
bash /SHAPE/submit_slurm_array_illumina.sh \
    --mail-user you@institution.edu \
    --script-path SHAPE/process_single_file_illumina_downsample.sh" \
    --sample-map <map_file.tsv> \
    --output-dir <output dir> \
    --ref-index <bowtie index> 
``` 

The sample map file should be formatted as follows:

```
No_PUS_DMSO No_PUS_DMSO_S4_L001_R1_001.fastq.gz No_PUS_DMSO_S4_L001_R2_001.fastq.gz
No_PUS_SHAPE    No_PUS_SHAPE_S3_L001_R1_001.fastq.gz    No_PUS_SHAPE_S3_L001_R2_001.fastq.gz
```

This will submit an individual job for each pair of fastq files in the map file. The final deduplicated fastq will be located in a directory inside the output directory.

### SHAPE-Mapper
This pipeline `SHAPE/shapemapper_pipeline.sbatch` goes through the following steps:

1. Runs shapemapper2-2.3 paired fastq files (DMSO and SHAPE)
2. Annotates RNAs with poor quality scores
3. Produces correlation plots of the replicates
4. Averages SHAPE reactivities across replicates
5. Performs SHAPE-informed RNA fold
6. Extracts pairing probabilities and structure feature summaries

This pipeline works the same for reads of Nanopore or Illumina origin, once they have been processed through the above pipelines.

One limitations of shapemapper is the number of sequences it can process at a given time. It will crash if your reference FASTA has more than ~100 sequences (represented by > header lines). To get around this, you should "chunk" your fasta into smaller lists. This script will run the chunks separately and then concatenate the results.

Shapemapper log output includes a list of "poor quality" RNAs, indicating an issue identified during processing. This typically results from inadequate read coverage for a given sequence. Low coverage can skew reactivities, so these should be excluded from downstream analysis, which is what the annotation section does.

In the presence of multiple replicates, inverse variance weighting averaging is used to create one set of shape reactivites to go into RNAfold.

After running SHAPE-informed RNAfold, pairing probabilities and folding information is extracted using the same scripts as in the RNAfold folder. Structures are also visualized with SHAPE reactivities as described below.

This script is optimized to work on Sherlock, Stanford's HPC. The following packages are required for this to run:

- Vienna RNA 2.5.1
- shapemapper2-2.3 https://github.com/Weeks-UNC/shapemapper2

Example command to run a single condition with two replicates:

```
sbatch --array=1-2 \
    --mail-user=you@institution.edu \
    ~/PUS7regulation2026/Figure4/SHAPE_mapper/shapemapper_pipeline.sbatch \
        --sample-map <sample_map.tsv> \
        --output-dir <sample output dir> \
        --ref-fasta-dir <reference_fasta_chunks> \
        --ref-rnafold-fasta <reference_fasta> \
        --num-samples 1
```


The sample map should look as follows:

```
GroupName       SampleName      Untreated       Treated
PUS7    PUS7_Rep1       PUS7_1_DMSO_deduplicated_for_shapemapper.fastq    PUS7_1_SHAPE_deduplicated_for_shapemapper.fastq
PUS7    PUS7_Rep2       PUS7_2_DMSO_deduplicated_for_shapemapper.fastq    PUS7_2_SHAPE_deduplicated_for_shapemapper.fastq
```

## Structure Visualization
These script automates the batch generation of 2D RNA MEA structure vector-based SVG diagrams using VARNA, without or without mapping SHAPE reactivity data to nucleotide colors (low, medium, high). It can also draw a custom bounding box and apply special formatting to highlight a specific nucleotide position of interest, such as the target uridine.

These scripts require python 3.9.

Activate the virtual environment and use as follows:

```
source varna/varna-env/bin/activate

python3 varna/run_varna_highlight_SHAPE.py \
  --input-csv <fold summary csv> \
  --shape-dir <shape directory> \
  --highlight-csv <highlight csv> \
  --output-dir <output dir> 

python3 varna/run_varna_highlight.py \
  --input-csv <fold summary csv> \
  --highlight-pos <int>> \
  --output-dir <output dir> 
```

The fold summary csv should contain the sequence and dot-bracket notation of the MEA structure, and which can be produced by `fold_summary.R`.

The highlight csv should contain the site name and highlight position (1-indexed) (optional arguments to specify) for the residue that should be highlighted in pink for the SHAPE script. For the non-SHAPE script, the highlight is done at a single consistent position.


## Mutagenesis Generation

### Accessibility (psi_pairing)

The goal of all scripts in this section is to alter the pairing status of target nucleotides - either 1, 2, 3, or 5 at a time. This is done by either breaking an existing pair to create unpaired residues, or by making a pair out of existing unpaired residues. 

An input directory of .fold files, generated by Vienna RNA 2.5.1 ```RNAfold -p -o --MEA <fasta_file>```, is required for these script to run. The script will attempt mutations on all sequences in the directory.

**Break pair** is straight forward: you take the target nucleotide, identify it's pairing partner based on dot-bracket notation, and then mutate the pairing partner to become the same nucleotide as the target nucleotide. So if you start with a target C nucleotide that pairs with a G (C-G pair), you will end up with a C-C pair.

**Make pair** takes one of two strategies. 
1. The first strategy is to "extend stem", wherein the script looks for a stem within 3 nucleotides of the target position. If one is found, the script then makes mutations on the opposing strand from the target, turning residues on the opposing strand into potential pairing partners for the target position(s).
2. The second strategy is to "seed hairpin", taken only if no nearby pairing region is found. In this scenario, mutations are made to the soon-to-be opposing strand to form a 2 base pair stem of a hairpin loop, locking up the target residue.

Structure mutations do not always work as intended. To check on the impact of the mutation, RNAfold is run on the mutated sequence. The script then assesses if the target positions were impacted as anticipated: break pair are now unpaired, make pair are now paired. This only examines pairing status and does not account for potentially larger structure rearrangements.

#### Target Uridine: alter_psi_pairing.R

The script alters the pairing status of a single nucleotide at a constant position across all sequences within the fold directory. 

During make pair - extend stem, the script looks for a maximum of 3 unpaired nucleotide between the target position and start of the stem, non-inclusive.

Example usage:
```
Rscript psi_pairing/alter_psi_pairing.R \
  --fold_dir <fold_dir> \
  --output_dir <output_dir> \
  --target_pos 65
```

#### Target Uridine and Plus 1: alter_psi_pairing_01.R

The script alters the pairing status of two nucleotides at a constant position across all sequences within the fold directory. The nucleotides are defined by a target position, and the +1 base immediately downstream of that position.

Break pair is applied to any targets that are both paired or of mixed status (one paired, one unpaired). Make pair is applied to any targets that are both unpaired or of mixed status.

During make pair - extend stem, the script looks for a maximum of 3 unpaired nucleotide between the target position and start of the stem, non-inclusive.

During make pair - seed hairpin, the script creates a 1 nucleotide hairpin with a two nucleotide stem: ((.)) This should probably be updated before the script is used again.

Example usage, targeting positions 65 and 66:
```
Rscript psi_pairing/alter_psi_pairing_01.R \
  --fold_dir <fold_dir> \
  --output_dir <output_dir> \
  --target_pos 65
```

#### Unpaired2: alter_3mer_pairing.R

The script alters the pairing status of three nucleotides at a constant position across all sequences within the fold directory. The nucleotides are defined by a target position, and the +1 base immediately downstream of that position, and the -1 base immediately upstream of the that position. The target position should be the center of the block.
This was applied to modulate the pairing status of the unpaired2 region downstream of the target uridine.

There is an additional ```--protect_pos``` argument that prevents any mutations to a 5mer region center on the positions. This is applied to prevent the 5mer sequence motif centered on the pseudouridine from being mutated while trying to alter the target 3mer pairing status.

Break pair is applied to any targets that are both paired or of mixed status (paired and unpaired nucleotides). Make pair is applied to any targets that are both unpaired or of mixed status.

During make pair - extend stem, the script looks for a maximum of 3 unpaired nucleotide between the target position and start of the stem, non-inclusive.

During make pair - seed hairpin, the script creates a 4 nucleotide hairpin with a 3 nucleotide stem: (((....)))

Example usage, targeting positions 72, 73, 74 and protecting positions 63 - 67:
```
Rscript psi_pairing/alter_3mer_pairing.R \
  --fold_dir <fold_dir> \
  --output_dir <output_dir> \
  --target_pos 73 \
  --protect_pos 65
```

#### 5mer motif: alter_5mer_pairing.R

The script alters the pairing status of five nucleotides at a constant position across all sequences within the fold directory. The nucleotides are defined by a target position, and the +1 and +2 bases immediately downstream of that position, and the -1 and -2 bases immediately upstream of the that position. The target position should be the center of the block. 
This was applied to modulate the pairing status of the 5mer containing the target uridine.

Break pair is applied to any targets that are both paired or of mixed status (paired and unpaired nucleotides). Make pair is applied to any targets that are both unpaired or of mixed status.

During make pair - extend stem, the script looks for a stem up to 10 nucleotides away from the end of the 5mer.

During make pair - seed hairpin, the script creates a 4 nucleotide hairpin with a 5 nucleotide stem: (((((....)))))

Example usage, targeting positions 63 - 67:
```
Rscript psi_pairing/alter_3mer_pairing.R \
  --fold_dir <fold_dir> \
  --output_dir <output_dir> \
  --target_pos 65
```

### Structure

The goal of `structure_mutation.R` is to disrupt all base pairs within a defined window (the Pupup signature), and then to restore those base pairs with compensatory mutations.

An input directory of .fold files, generated by Vienna RNA 2.5.1 ```RNAfold -p -o --MEA <fasta_file>```, is required for these script to run. The script will attempt mutations on all sequences in the directory.

The script is given a target position, the target uridine. It protects a 5 nucleotide block center on this position from mutation. It identifies the upstream window as -10 to -3, inclusive, from the target position. It identifies the downstream window as +3 to +15, inclusive, from the target positions.

First, the script generates disruption mutations. It identifies all paired residues within the target window and, similar to the accessibility scripts above, mutates the partner to be identical to the window base. So an A-U pairs, where A is in the window, becomes and A-A mismatch.

Next, the script generates compensation mutations. From the paired residues identified above within the target window, it now mutates the window base to complement its mutated partner. So the A-A mismatch above now becomes an U-A match. This is the opposite of the original A-U pair.

Finally, the script uses RNAfold to predict the dot-bracket notation of the mutated sequences and saves those into a single CSV file `all_mutagenesis_results.csv` for future analysis.

Example Usage:

```bash
Rscript structure_mutation/structure_mutation.R \
  --input_dir <fold_dir> \
  --output_dir <output_dir> \
  --position_zero 65
```


## Parameter Sweep / motif matcher

Code to identify structure motifs in RNA structures where .fold files are available. None of this code was actually used in mutagenesis for Pool2, so proceed with caution and an informed perspective.

### original version
Created by Ronit Jain.

motif_matcher.R

#### Description

Takes as input RNA sequences IDs and corresponding RNAfold .fold files. Searches for motifs that match a paired-unpaired-paired pattern based on dot notation from RNAfold. Specify the position in the sequence to start the search from, the offset, and the min and max number of base pairs for each region. Returns a list of sequences that match the pattern.

Search is performed by querying each sequence for the specified regions, starting with paired1 at the location specified by input_position. Start with the maximum length given for paired1 (default = 6 nts), queries the sequence for that number of sequentially paired bases. If not found, decreases range by 1, tries again, and repeats until the minimum is reached. Takes the maximum value where all sequential bases are paired as part of the motif, starts search for unpaired region immediately downstream of that. Again, starts with the max range given, this time for unpaired bases, and tests until a full range of sequential unpaired bases is found. Takes that maximum range of unpaired bases, starts search for paired2 region downstream of that. Again, starts with max range given, searches down until only sequential paired bases are found. Returns all sequences that are found with this, though individual sequences are likley to have different lengths of paired-unpaired-paired.

If sequence downstream of the input_position is not longer than the combined maximums for each range, the sequence will not be considered.

#### Parameters

Flag            Description

-i, --input     Path to the input CSV file containing RNA sequence IDs and associated data. This file is passed into motif_matcher.R. The csv file must have a column called "id" which contains sequence names that match a substring in the name of the associated .fold file. This is how the script pulls out the relevant .fold files which correspond to your sequences of interest. 

-f, --fold_dir	Directory containing .fold files from running RNAFold on the input sequences (output from RNAfold). Used by motif_matcher.R to identify candidate motifs. Note that this needs to be done before running the pipeline. 

-o, --output	Path to the final output CSV file that will include the original and mutant data, along with motif annotations after re-folding.

--input_position	Pseudouridine position in your sequences. Default is 65 (based on the pool1 sequences where psi is centered at position 65). 

--offset_min	Minimum offset from the input_position to start searching for the upstream paired region (paired1). Default is 1.

--offset_max	Maximum offset from the input_position to search for the start of the paired1 region. Default is 3.

--min_paired1	Minimum number of base pairs in the upstream (paired1) region. Default is 3.

--max_paired1	Maximum number of base pairs in the paired1 region. Default is 6.

--min_unpaired	Minimum length of the unpaired region connecting paired1 and paired2. Default is 3.

--max_unpaired	Maximum length of the unpaired stretch between the two paired regions. Default is 7.

--min_paired2	Minimum number of base pairs in the downstream (paired2) region. Default is 3.

--max_paired2	Maximum number of base pairs in the paired2 region. Default is 10. Note that the paired region could be greater than this maximum. For instance, the paired region could be actually 13 nucleotides long, but it would still be extracted since there is still a 10 nt paired region within that 13 nt region. 

### version two
Edited by Rebecca Rodell. 

motif_matcher_v2.R

#### Description

Same purpose as motif_matcher.R, original version, but with some key changes:
1. includes another unpaired region before the original paired-unpaired-paired pattern
2. includes optionality for paired2 with --include_paried2

Takes as input RNA sequences IDs and corresponding RNAfold .fold files. Searches for motifs that match a unpaired-paired-unpaired*-paired) pattern based on dot notation from RNAfold. Specify the position in the sequence to start the search from, the offset, and the min and max number of base pairs for each region. Returns a list of sequences that match the pattern.

Code searches for motifs with same method as described in original motif_matcher.R, just with extra unpaired region.

#### Parameters

Flag                Description

-i, --input	        Path to the input CSV file containing RNA sequence IDs and associated data. This file is passed into motif_matcher.R. The csv file must have a column called "id" which contains sequence names that match a substring in the name of the associated .fold file. This is how the script pulls out the relevant .fold files which correspond to your sequences of interest. 

-f, --fold_dir	    Directory containing .fold files from running RNAFold on the input sequences (output from RNAfold). Used by motif_matcher.R to identify candidate motifs. Note that this needs to be done before running the pipeline. 

-o, --output	    Path to the final output CSV file that will include the original and mutant data, along with motif annotations after re-folding.

--input_position	Pseudouridine position in your sequences. 1-indexed position to start searching. **Default is 59** (based on the pool1 sequences where psi is centered at position 59, a key difference from original version). 

--offset_min	    Minimum offset from the input_position to start searching for the upstream unpaired region (unpaired1). Default is 1.

--offset_max	    Maximum offset from the input_position to search for the start of the unpaired1 region. Default is 3.

--min_unpaired1     Minimum number of base pairs in the first unpaired region. Default is 0.

--max_unpaired1     Maximum number of base pairs in the first unpaired region. Default is 1.

--min_paired1	    Minimum number of base pairs in the first paired region. Default is 2.

--max_paired1	    Maximum number of base pairs in the first paired region. Default is 6.

--min_unpaired2	    Minimum number of base pairs in the second unpaired region. Default is 2.

--max_unpaired2	    Maximum number of base pairs in the second unpaired region. Default is 6.

--include_paired2   Optionality to include paired2 in the search. Default is TRUE, where paired2 would be included in the motif.

--min_paired2	    Minimum number of base pairs in the second paired region. Default is 3.

--max_paired2	    Maximum number of base pairs in the second paired region. Default is 10. Note that the paired region could be greater than this maximum. For instance, the paired region could be actually 13 nucleotides long, but it would still be extracted since there is still a 10 nt paired region within that 13 nt region. 

### parameter sweep
Created by Rebecca Rodell.

parametersweep_v2.R

*Note: This is hard-coded for my purposes, including the parameters to screen and how the f1 score is calculated. Please review and manually alter for your purposes before using.*

#### Description

Performs a grid search to test a range of parameters, with computational parallelization. For each parameters range, calculates an f1 score for the true positive rate (recall) and positive predictor rate (precision) of a given input dataframe with binary values. Returns various f1 scores and corresponding parameters. Outputs individual motif files. This script works, but it could use some improvement in how results are output.

Most important output: all_analysis_results.csv

Perform analysis on this output dataframe with **parametersweep_analysis.R** in R Studio. Extract the range of parameters that give you the greatest f1 score for structure alone. This is the range of values representative of the most predicitive and precise motif boundaries for the given data.


parametersweep_v2_test.R: *This is hard-coded for my purposes, using a more limited version of the parameters to ensure the script is functional. Please review and manually alter for your purposes before using.*

#### Parameters

Flag                Description

-i, --input         Input CSV file with RNA sequences IDs. What will be used by motif_matcher_v2.R. default="input_pool1.csv"

-f, --fold_dir      Directory containing .fold files. What will be used by motif_matcher_v2. R. default="/scratch/users/rodell/RNAfold_psipos"

-o, --output_dir    Output directory for results.

-d, --dataset       Dataset file with binary values for evaluation by f1 score. default="pool1_psipos_info.csv"