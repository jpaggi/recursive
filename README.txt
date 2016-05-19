
Defining set of expressed genes
... the fact that there are "unexpressed genes" that get detected indicates that
you are potentially not doing something right...

But it could indicate that there is simply a large subset of genes with RS
that are at low coverage levels. In fact, this seems like the most sensible choice?

44 -> 14 not detected and 153 -> 143 detected. This seems like a pretty good enrichment?
Certainly goes along with the "there is a lot of RS hypothesis"


something very strange is happening at 2L 20967746 21051772

Do we get all here: X 21082826 21152071


Could look for splicing accuracy using 4sU data
Look for non-recursive splicing events that are present in 5 
minute sample but not total (need to find background in total)

Do not necessarily expect to see 'abherrant' spliced transcripts removed from
total RNA-seq

Overall
- Helpful to rule out rs as well as confirming it
- trend towards global statistics instead of individual examples
- chr2L   2890915 2924913 + is clearly being used as a ratchet site!, but marked as annotated


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
