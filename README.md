# RNAstructure
Code for analysis of RNA structure

## motif matcher

Code to identify structure motifs in RNA structures where .fold files are available.

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