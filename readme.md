# Identification of recursive splice sites using metabolic labeling data

## General Notes.
	load_genome.py and get_motifs.py are used throughout
	... must generally add sjr/core_pipeline to path or copy these files into working directory

## Defining Intron Set.
	Scripts to make introns from GTF file and get coverage data are in sequence/coverage.

	get_introns.py parses the GTF file and outputs a bed file of long introns

	get_intron_expression.py takes the long intron file and a BAM file
	... gives pileup of coverage at each postion + number of sjr spanning intron

	merge.py sums coverage data from multiple BAM files

	expressed.py takes a sjr, min size and max size cutoff and prints the corresponding set of introns

## Preparing Coverage for MCMC.
	Still in sequence/coverage.

	remove_exons.sh does everything... workflow is:
	  1) Extract reads in different size ranges.
	  2) Use bedtools subtract to remove exonic regions
	  3) Merge together intronic regions, summing the sjr counts
	  4) Recompute coverage of each intronic segment for all replicates
	  5) Merged all replicate expression levels together
	  6) Mask repeats

## Running MCMC.
	All scripts are in sawtooth/mcmc_core.

	run_mcmc.py feeds coverage data into mcmc.py
	... tunable parameters are set in run_mcmc.py

	call_sites.py transforms the MCMC probabilities into RS predicitons
	... setting here to instead produce random peaks for FDR analysis

	merge_sites.py merges individual sites implicated by seperate peaks

## Running SJR pipeline.
	All necessary scripts are in sjr/core_pipeline/
	Read run_all.sh and change hard coded inputs as appropriate.
	This script will automatically run both the sjr and read pair pipeline
	I honestly don't see why you would need to rerun this though.
	There is the raw output of this on UTR/jpaggi/reads/sjr

## Motif Analysis.
	To prepare intron annotations follow directions above,
	then use bedtools merge (this stops intervals from being double counted).
	Then use sjr/get_background.py to get random motifs

	Look in code for how to load in the 

## Splicing Rates.
	First use combine/standard_table_reader.py to extract the set of recursive sites that you want to use
	Then, run sequence/make_introns.py with the return statement uncommented. This merges recursive sites
	together into groups by intron.
	Finally, call rates/run_all.py with this file as an argument. ( you also need a gtf file of exons.. just use grep exons)

	You can use plot.py with the summary table as an argument to plot the mean psi values
	Use plot_transformed.py to plot coverage in transformed space (for QC)

## Combining Data.
	Use standard_table.py to create a combined table of graveley events, sawtooth events, and sjr detected sites
	... read through to replace hardcoded file paths

	Call combine/define_introns.py followed by combine/get_sjr.py to fill in exon detection data
