Do mapping with dmel r5.57

Sawtooth

 - Stick with current implementation now, improve if needed later
 - Look for overlaps between samples and time points
 - Need to determine exactly which set of genes to run this on!
   - Noticed higher false positive rates on 

SJR

 - Do all alignments using HISAT2
 - Plot motif score versus number of reads
 - Potentially plan some sort of cutoff?
 - figure out why chrX:21082827-21152072 aligns improperly
 - seems that there can sometimes be slight mapping issues
    ... get around by looking a few bases on each side for motif?

BPs

 - Write software, ignoring inverted reads for now
 - Later cut things out based on intronic location or motif information


Combining Data

 - Lean towards high recovery rate in all individual analysis
 - Later filter based on agreement between data types