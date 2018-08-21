# Identification of recursive splice sites using metabolic labeling data

Three independent methods for identification of recursive splice sites from RNA-seq data are developed:

- RatchetJunction, a previously described method using splice junction reads,
- RatchetPair, which uses paired end reads straddling a splice junction,
- RatchetScan, which infers recursive splice site locations from patterns in the read coverage of introns.

We applied these methods to identify recursive sites in Drosophila, but they could be applied to study recursive splicing in any organism. The results of our study are described [here](https://www.biorxiv.org/content/early/2017/02/13/107995). For a snapshot of our code at the time of this paper, please check out the branch 'drosophila', but I expect that the current code base will be much more useful for someone seaking to use our methods on a new dataset.

Below, I go through an example application of our methods, which I hope will be straightfowardly adaptable to a new dataset. The example folder contains
example input and output files (excluding the input reads files and results for all but chomosome 2L to meet githubs space restraints).

(In practice, it may be best to run some of the analysis steps in parallel.
In particular, anything that is performed for all input datasets can be run
on a separate core (i.e. any expression wrapped in "for LIB in $LIBRARIES; do ... done").
Additionally, RatchetScan can be parallelized by splitting the input file 
($EXPRESS_DIR/all.masked.bed) into multiple files and running each separately.)

##  Procure Input Files

For this example, we will use nascent RNA-seq data from Drosophila S2 cells. You can get the reads
using tools from the SRA tool kit. Here we use hisat2 to align the reads (any spliced read aligner
would do), then samtools to produce a sorted and indexed bam file.

```
prefetch -v SRX2500203
fastq-dump --outdir reads/ --split-files /Users/jpaggi/ncbi/public/sra/SRX2500203.sra

hisat2 -x dm6/genome -1 SRX2500203_1.fastq -2 SRX2500203_2.fastq -p 2 -S SRX2500203.sam
samtools view -bS SRX2500203.sam > SRX2500203.bam
samtools sort SRX2500203.bam -o SRX2500203_sorted.bam
samtools index SRX2500203_sorted.bam SRX2500203_sorted.bai
```

You will need 3 annotation files
- A GTF annotation file (i.e. dmel-all-r6.21.gtf from ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gtf/)
- A FASTA file of the organism of study's genome (i.e. dmel-all-chromosome-r6.21.fasta from ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/)
- A RepeatMasker annotation file available for download at repeatmasker.org (i.e. dm6.fa.out (renamed dm6.repeats.txt) from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html)


## Environment set-up
In the example below we will use the following environment variables:

```
CODE=/path/to/here/
ANNO=anno
GTF=dmel-all-r6.21.gtf
FASTA=dmel-all-chromosome-r6.21.fasta
REPEATS=dm6.repeats.txt
LIBRARIES='SRX2500203_sorted'

READS_DIR=reads
EXPRESS_DIR=expression
SJR_DIR=sjr
SAWTOOTH_DIR=sawtooth
```

Additionally, you will need to add the utils directory to your PYTHONPATH, which you can do by
changing to the utils directory and executing "export PYTHONPATH=$PYTHONPATH:\`pwd\`".

## Process annotations

Convert GTF format annotations into bed files of introns and exons.

```
python $CODE/expression/get_introns.py $ANNO/$GTF | sort -k1,1 -k2,3n -u > $ANNO/introns.bed
awk ' $3=="exon" {print $1"\t"$4"\t"$5"\t.\t.\t"$7}' $ANNO/$GTF | sort -k1,1 -k2,3n -u > $ANNO/exons.bed
```
## RatchetJunction

Identifies reads that span a putative 5'ss-recursive site junction.
Specifically, finds reads where the upstream half aligns just upstream of an
annotated 5'ss and the downstream end aligns inside an intron at a site
with sequence AG|GT (juxtaposed 3' and 5' ss motif).

Each library (bam file) is processed independently, then all detected sites are merged together.

```
for LIB in $LIBRARIES; do
    python $CODE/RatchetJunction/putative_ratchet_sjr.py \
    	   $READS_DIR/$LIB.bam $ANNO/introns.bed $ANNO/$FASTA > $SJR_DIR/$LIB.sjr.bed
    sort -k1,1 -k2,3n $SJR_DIR/$LIB.sjr.bed \ 
    	   | python $CODE/RatchetJunction/group_introns.py sjr_$LIB > $SJR_DIR/$LIB.sjr_groups.bed
done

cat $SJR_DIR/*.sjr_groups.bed | sort -k1,1 -k2n,3n \ 
    | python $CODE/RatchetJunction/merge_sjr_reps.py  > $SJR_DIR/all.sjr_groups.bed
```
## RatchetPair

Identifies reads where one mate is close upstream of a 5'ss and the other mate is far into the intron and
applies a probabalistic model to assign these to a sparse set of putative recursive sites. This is accomplished
 using an expectation maximization framework.

Note that this method takes as input the output of RatchetJunction. Paired end reads are preferentially assigned
to putative recursive sites that have support from junction reads.

```
# Extract putative 5'ss - recursive site straddling reads.
for LIB in $LIBRARIES; do
    python $CODE/RatchetPair/straddle_jxn.py \
    	   $READS_DIR/$LIB.bam $ANNO/introns.bed $ANNO/$FASTA pe_$LIB > $SJR_DIR/$LIB.straddle.bed
done
cat $SJR_DIR/*.straddle.bed > $SJR_DIR/all.straddle.bed

# Expectation - maximization algorithm to assign reads to recursive sites.
python $CODE/RatchetPair/straddle_gem.py \
       $SJR_DIR/all.straddle.bed $SJR_DIR/all.sjr_groups.bed $ANNO/introns.bed $ANNO/$FASTA \
       | sort -k1,1 -k2,4n > $SJR_DIR/all.straddle_gem.bed
```
## RatchetScan

This algorithm recognizes the "sawtooth pattern" in RNA-seq read counts present in introns that
are recursively spliced. Generally, there is an approximately linear decay in read coverage between
the 5'ss and the first recursive site, a jump in coverage at the recursive site, followed by another
linear decay in coverage to the next recursive site, and so on. This pattern arises because the upstream
end of introns exist longer before the intron is spliced and degraded than the downstream end.

Our method detects this pattern by using a Markov Chain Monte Carlo (MCMC) algorithm to sample
potential combinations of recursive sites and assess how well they explain the rna-seq coverage.
This procedure is used to compute a probability that each position in an intron is a recursive site.
Recursive sites are then called by finding a sites near peaks in probability with strong juxtaposed 3'ss - 5'ss motifs.

### Compute coverage in introns and number of junction spanning reads

Compute rna-seq read coverage in all introns and find the subset of introns that
are expressed i.e. have at least one read spanning the intron.
```
for LIB in $LIBRARIES;
    do python $CODE/expression/get_intron_expression.py \
       $READS_DIR/$LIB.bam $ANNO/introns.bed > $EXPRESS_DIR/$LIB.full.bed
done
sort -k1,1 -k2,3n $EXPRESS_DIR/*.full.bed | python $CODE/expression/merge.py > $EXPRESS_DIR/all.full.bed
awk '$5 > 0 && $3 - $2 > 8000' $EXPRESS_DIR/all.full.bed | \  # $5 is junction read count
    cut -f 1-6 > $EXPRESS_DIR/expressed.full.bed
```
### Compute coverage without exons

Introns with internal exons present a major challenge for this approach. Here we get around this problem by
only considering subsets of introns that don't overlap with an exon. RatchetScan should be applied to this
cleaner subset of the data.

```
bedtools subtract -s -a $EXPRESS_DIR/expressed.full.bed -b $ANNO/exons.bed | sort -k1,1 -k2n,3n | \
	 python $CODE/expression/merge_subtracted.py > $ANNO/expressed_introns_no_exons.bed
for LIB in $LIBRARIES; do
    python $CODE/expression/get_subtract_expression.py \
    	   $READS_DIR/$LIB.bam $ANNO/expressed_introns_no_exons.bed > $EXPRESS_DIR/$LIB.noexon.bed
done
cat $EXPRESS_DIR/*.noexon.bed | sort -k1,1 -k2n,3n | \
    python $CODE/expression/merge.py > $EXPRESS_DIR/all.noexon.bed
```
###  Mask repeats

Repeat regions represent another challenge for this approach because they often have erratic rna-seq coverage
as it is impossible to uniquely map reads to them. Here we mask out these regions by replacing the coverage values
in repeat regions by the coverage values in the regions flanking them.

```
python $CODE/expression/mask_repeats.py \
       $EXPRESS_DIR/all.noexon.bed $ANNO/$REPEATS > $EXPRESS_DIR/all.masked.bed`
```

### Core RatchetScan Algorithm

Note that this first script will have the longest runtime of the whole pipeline. It could trivially be parallelized by splitting the all.masked.bed file into smaller files and then merging the results.
```
python $CODE/RatchetScan/run_mcmc.py $EXPRESS_DIR/all.masked.bed > $SAWTOOTH_DIR/mcmc.bed
python $CODE/RatchetScan/call_sites.py \
       $SAWTOOTH_DIR/all.mcmc.bed $ANNO/$FASTA $ANNO/introns.bed > $SAWTOOTH_DIR/sites.bed
sort -k1,1 -k2,3n $SAWTOOTH_DIR/sites.bed \
       | python $CODE/RatchetScan/merge_peak_calls.py > $SAWTOOTH_DIR/merged_sites.bed
```
##  Post Processing / Visualization

The goal of these scripts is to tabulate the results of the above 3 methods and visualize the results.

### Combine the output of all the methods into a single file
```
python $CODE/utils/combine_results.py $SJR_DIR/all.sjr_groups.bed $SJR_DIR/all.straddle_gem.bed \
       $SAWTOOTH_DIR/sites.bed $ANNO/$FASTA $ANNO/introns.bed | sort -k1,1 -k2,3n > results.bed
```

### Plot the results of RatchetScan, annotated with locations of strong motifs and splice junction reads.
```
python $CODE/RatchetScan/plot_mcmc.py \
       $ANNO/$FASTA $ANNO/introns.bed $SAWTOOTH_DIR/mcmc.bed
```