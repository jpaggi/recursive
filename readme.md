# Identification of recursive splice sites using metabolic labeling data

Three independent methods for identification of recursive splice sites from RNA-seq data are developed:

- RatchetJunction, a previously described method using splice junction reads,
- RatchetPair, which uses paired end reads straddling a splice junction,
- RatchetScan, which infers recursive splice site locations from patterns in the read coverage of introns.

We applied these methods to identify recursive sites in Drosophila, but they could be applied to study recursive splicing in any organism. The results of our study are described [here](https://www.biorxiv.org/content/early/2017/02/13/107995).

The script *run_all_methods.sh* was used to run all of these methods in our study.

## General notes.
- The scripts load_genome.py and get_motifs.py are used throughout, add RatchetJunction_RatchetPair/core_pipeline to your PYTHONPATH or copy these files into working directory
- Code for aligning reads is not included in this repository. We used hisat2 with default settings in our study, but feel free to use your favorite *spliced*-read aligner.
- Many of the high-level scripts are hardcoded with file names particular to our study. If you want to apply our method to a new dataset, you will need to change these.

## Defining the set of introns.
Scripts to make introns from GTF file and get coverage data are in intron_expression_processing/coverage.

1. Download a gtf file for the organism you are studying.
2. Extract just the 'exon' entries and sort them by transcript name.
3. Run "cat <file>.gtf | get_introns.py > introns.bed" to parse the GTF file into a bed file of 1000+ base pair introns
4. Run "python get_intron_expression.py reads.bam introns.bed > intron_expression.bed" to obtain a pileup of coverage at each postion + number of sjr spanning intron
5. Merge coverage data from mulitple BAM files using merge.py.

## Preparing coverage for RatchetScan.
Still in intron_expression_processing/coverage.

remove_exons.sh does everything, note that this script contains hardcoded paths, so you will have to edit it... workflow is:
1. Extract reads in different size ranges.
2. Use bedtools subtract to remove exonic regions
3. Merge together intronic regions, summing the sjr counts
4. Recompute coverage of each intronic segment for all replicates
5. Merged all replicate expression levels together
6. Replace expression values in repeat regions with average of neighboring regions.

## Running RatchetScan.
All scripts are in RatchetScan/mcmc_core.

run.sh will run all steps of the RatchetScan pipeline. Substeps are

1. run_mcmc.py feeds coverage data into mcmc.py. tunable parameters are set in run_mcmc.py, they are currently set as we used in our study.
2. call_sites.py transforms the MCMC probabilities into RS predicitons. setting here to instead produce random peaks for FDR analysis
3. merge_sites.py merges individual sites implicated by seperate peaks

## Running RatchetJunction and RatchetPair.
RatchetJunction identifies recursive sites using reads that span a 5'SS - recursve site splice junction. These reads directly implicate a recursive splice site in the same way that a standard splice junction read gives away the position of normal splice sites.

RatchetPair uses paired end reads that have one read pair aligning upstream of a 5'SS and the other end aligning far into the intron to implicate recursive sites. This method relies on our nowldge of the distribution of distances between read pairs (commonly called insert lengths) and the expected sequence motif at ratchet sites. Using this information, RatchetPair uses a variant of the GEM algorithm to find a sparse set of recursive sites that can explain the recursive splice junction reads and putative recursive splice site straddling reads.

All necessary scripts are in RatchetJunction_RatchetPair/core_pipeline/.

run_all.sh will run all steps of both pipelines. Substeps are

1. Extract reads that potentially straddle recursive splice junctions.
2. Extract putative recursive splice junction reads
3. Group putative recursive splice junction reads that have shared a 5'ss
4. Merge together data from all timepoints and replicates
5. Run RatchetPair algorithm (See straddle_gem.py).
6. Assign groups of reads to an intron based on intron expression levels

## Combining output of methods
Use combine/standard_table.py to create a combined table of graveley events, sawtooth events, and sjr detected sites
... read through to replace hardcoded file paths

Call combine/define_introns.py followed by combine/get_sjr.py to fill in exon detection data.

The combine directory contains several tools for visualizing expression levels in long introns.