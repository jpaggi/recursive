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
 ... seems to be an indel close to the sites? Or a micro exon???
 - seems that there can sometimes be slight mapping issues
    ... get around by looking a few bases on each side for motif?

 - have read pair straddling these sites, indicative of other avoidable misses
 ... 2R	15461832	15488386	15476248	-
 ... 3R	25256963	25297372	25272648	+
 ... 2R	6878630	6987036	6940011	- (seems to be mappibility issues)
 ... 2R	2391194	2501609	2433805	+ (btw igv won't show very long pairs)
 ... 3L	14310438	14359847	14339217	+
 ... 2L	18320731	18383401	18333385	+
 ... 2L	4032770	4074513	4056254	+
 ... 3R	5701946	5756085	5718362	+  
 ... 3R	5701946	5756085	5727133	+  check 
 ... 3R	5701946	5756085	5745466	+  check read ID  1327:35050


 - 3R	4906745	4982839	4918040 has no support but in expressed intron

 - has support in total rna-seq
 ... 3R	4270410	4303140	4289597	
 ... 2L	17197087	17251508	17226473	- 
 ... X	21082827	21152072	21109455	- and straddling read in 5 minute sample

- 2R	13163569	13238598	13193491	+ is an annotated 3'ss (ostensibly not in use though, very much detected)

- putative non-ratchet rs intron in middle of 3R:747,829-769,025

BPs

 - Write software, ignoring inverted reads for now
 - Later cut things out based on intronic location or motif information


Combining Data

 - Lean towards high recovery rate in all individual analysis
 - Later filter based on agreement between data types