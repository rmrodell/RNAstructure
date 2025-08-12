Process Summary:

1. Identify a window +/- from psipos
2. Identify AU / GC pairs where at least ONE nucleotide is within that range. Both do not have to be in the range.
      - Exclude anything within +/- 2 of psipos (the UNUAG 5mer): these bases or what pair to them will not be mutated
3. Mutate pairs to AU, GU, or GC. Only one mutation type is made in a given run (example: all AU mutated to GC, or all AU mutated to GU)
4. Refold mutated sequences with RNAfold.
5. Compare dot-bracket notation between original and mutated sequences to see if the same pairing status (paired or unpaired) is preserved.
6. This does NOT check if a base is pairing to the same thing, just if its still paired or unpaired.
7. Return metrics:
    a. the total number of sequences with increased / decreased delta G (> |5|)
    b. if mutations are upstream, downstream, or both from psipos
    c. structure conservation based on pairing status
        i. fraction different between original and mutant
        ii. number without changes in psipos 5mer (UNUAG)
        iii. number without changes at mutated sites

Sample usage: 
Rscript RNAstructure/deltaG_mutations/alter_deltaG.R --fold_dir <directory with .fold files> --output_dir <dir>/stabilizeGC_window15 --mode stabilize --window 15 --mutation_strategy GC

Required:
--fold_dir: location of folder with .fold files
--output_dir: where all output files should be placed

Defaults:
--psipos = 59
--window = 5
--deltaG_threshold = 5.0
--mode = destabilize
--mutation_strategy = AU
