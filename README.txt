
Unsolved Mysteries.
	New pipeline misses [('2R', '-', 6940011), ('3R', '-', 4918040)]
	... first has an unmappable region, second is likely not actually used
	... both only supported by pe r

	something very strange is happening at 2L 20967746 21051772

Defining Intron Set.
	Scripts to make introns from GTF file and get coverage data are in sequence/coverage.

	get_introns.py parses the GTF file and outputs a bed file of long introns

	get_intron_expression.py takes the long intron file and a BAM file
	... gives pileup of coverage at each postion + number of sjr spanning intron

	merge.py sums coverage data from multiple BAM files

	expressed.py takes a sjr, min size and max size cutoff and prints the corresponding set of introns

Preparing Coverage for MCMC.
	Still in sequence/coverage.

	remove_exons.sh does everything... workflow is:
	  1) Extract reads in different size ranges.
	  2) Use bedtools subtract to remove exonic regions
	  3) Merge together intronic regions, summing the sjr counts
	  4) Recompute coverage of each intronic segment for all replicates
	  5) Merged all replicate expression levels together
	  6) Mask repeats

Running MCMC.
	All scripts are in sawtooth/mcmc_core.

	run_mcmc.py feeds coverage data into mcmc.py
	... tunable parameters are set in run_mcmc.py

	call_sites.py transforms the MCMC probabilities into RS predicitons
	... setting here to instead produce random peaks for FDR analysis

	merge_sites.py merges individual sites implicated by seperate peaks

Running SJR pipeline.

Motif Analysis.

Rates.

Combining Data.




