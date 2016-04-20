Do mapping with dmel r5.57

Overall
- Helpful to rule out rs as well as confirming it
- Need to come up with set of genes to look in
... based on expression cut offs
... number of reads supporting intron being spliced
... length of introns
... presence in graveley study
- trend towards global statistics instead of individual examples


Sequence Motifs
- Figure out what a bit score for an individual sequence means
... try summing up information content of all characters of equal or greater probability???
- Might want to cut off low scoring annotated motifs before including in consensus motif
- Develop background model of rs site motifs
... look around all AG dinucleotides in long introns
... 


Sawtooth
 - go forward here despite being unsure
 - increase temperature
 - Stick with current implementation now, improve if needed later
 - Look for overlaps between samples and time points
 - Need to determine exactly which set of genes to run this on!
   - Noticed higher false positive rates on 

Sequence

SJR

- incorporate PE read information to numbers
- seems that should be able to count excessively long insert length PE
... reads on equal footing as split reads
- eventualy need to incorporate graveley data

 - figure out why chrX:2108827-21152072 aligns improperly
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
 ... 11 total (counting one from below)


 - 3R	4906745	4982839	4918040 has no support but in expressed intron

 - has support in total rna-seq
 ... 3R	4270410	4303140	4289597	
 ... 2L	17197087	17251508	17226473	- 
 ... X	21082827	21152072	21109455	- and straddling read in 5 minute sample

- 2R	13163569	13238598	13193491	+ is an annotated 3'ss (ostensibly not in use though, easily detected)

- putative non-ratchet rs intron in middle of 3R:747,829-769,025

threshold = 0.8
novel [157, 68, 31, 34, 22, 18, 17, 16, 73] 436
gravely [17, 4, 7, 5, 6, 7, 11, 9, 92] 158
total [5366, 1732, 937, 561, 424, 291, 216, 150, 300] 9977

threshold = 0.85
novel [44, 25, 14, 15, 8, 6, 10, 10, 65] 197
gravely [17, 4, 7, 5, 6, 7, 11, 9, 92] 158
total [5366, 1732, 937, 561, 424, 291, 216, 150, 300] 9977

threshold = 0.87  # min score of a graveley site
novel [24, 15, 12, 9, 5, 6, 9, 8, 62] 150
gravely [17, 4, 7, 5, 6, 7, 11, 9, 92] 158
total [5366, 1732, 937, 561, 424, 291, 216, 150, 300] 9977

threshold = 0.9
novel   [16, 9, 9, 6, 3, 3, 6, 5, 52] 109
gravely [16, 3, 6, 4, 5, 5, 10, 8, 82] 139
total [5366, 1732, 937, 561, 424, 291, 216, 150, 300] 9977


BPs

 - Write software, ignoring inverted reads for now
 - Later cut things out based on intronic location or motif information


Combining Data

 - Lean towards high recovery rate in all individual analysis
 - Later filter based on agreement between data types
 