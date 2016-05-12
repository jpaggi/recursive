Do mapping with dmel r5.57

something very strange is happening at 2L 20967746 21051772

Do we get all here: X 21082826 21152071


Could look for splicing accuracy using 4sU data
Look for non-recursive splicing events that are present in 5 
minute sample but not total (need to find background in total)

Do not necessarily expect to see 'abherrant' spliced transcripts removed from
total RNA-seq

Overall
- Helpful to rule out rs as well as confirming it
- Need to come up with set of genes to look in
... based on expression cut offs
... number of reads supporting intron being spliced
... length of introns
... presence in graveley study
- trend towards global statistics instead of individual examples
- chr2L   2890915 2924913 + is clearly being used as a ratchet site!, but marked as annotated

Need to be able to:
1) compute overlap between sjr and mcmc reads
2) look for nearest good motif / find best motif in a given range
3) compare found sites with annotated ss and such


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


BPs

 - Write software, ignoring inverted reads for now
 - Later cut things out based on intronic location or motif information


Combining Data

 - Lean towards high recovery rate in all individual analysis
 - Later filter based on agreement between data types
