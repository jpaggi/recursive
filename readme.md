
##############################################################################
#################    Download Files, align reads, etc.   #####################
##############################################################################

# Required software
# HISAT2
# SRA Toolkit
# samtools
# bedtools

# Get Raw Data (SRA Toolkit scripts)
prefetch -v SRX2500203
fastq-dump --outdir reads/ --split-files /Users/jpaggi/ncbi/public/sra/SRX2500203.sra
# dmel-all-r6.21.gtf              from ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gtf/
# dmel-all-chromosome-r6.21.fasta from ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/
# dm6 (hisat2 index file)         from https://ccb.jhu.edu/software/hisat2/manual.shtml
# dm6.fa.out -> dm6.repeats.txt   from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html

# Align reads
hisat2 -x dm6/genome -1 SRX2500203_1.fastq -2 SRX2500203_2.fastq -p 2 -S SRX2500203.sam
samtools view -bS SRX2500203.sam > SRX2500203.bam
samtools sort SRX2500203.bam -o SRX2500203_sorted.bam
samtools index SRX2500203_sorted.bam SRX2500203_sorted.bai

#############################################################################
##################          Process annotations          ####################
#############################################################################

CODE=~/Documents/burge/recursive_example/recursive
ANNO=anno
GTF=dmel-all-r6.21.gtf
FASTA=dmel-all-chromosome-r6.21.fasta
REPEATS=dm6.repeats.txt
LIBRARIES='SRX2500203_sorted'
READS_DIR=reads
EXPRESS_DIR=expression
SJR_DIR=sjr
SAWTOOTH_DIR=sawtooth

python $CODE/expression/get_introns.py $ANNO/$GTF | sort -k1,1 -k2,3n -u > $ANNO/introns.bed
awk ' $3=="exon" {print $1"\t"$4"\t"$5"\t.\t.\t"$7}' $ANNO/$GTF | sort -k1,1 -k2,3n -u > $ANNO/exons.bed

#############################################################################
##################    Compute Expression in Introns      ####################
#############################################################################

# Compute coverage in introns and number of junction spanning reads
for LIB in $LIBRARIES; do python $CODE/expression/get_intron_expression.py $READS_DIR/$LIB.bam $ANNO/introns.bed > $EXPRESS_DIR/$LIB.full.bed; done;
sort -k1,1 -k2,3n $EXPRESS_DIR/*.full.bed | python $CODE/expression/merge.py > $EXPRESS_DIR/all.full.bed
awk '$5 > 0 && $3 - $2 > 8000' $EXPRESS_DIR/all.full.bed | cut -f 1-6 > $EXPRESS_DIR/expressed.full.bed  # $5 is junction read count

# Get expression without exons
bedtools subtract -s -a $EXPRESS_DIR/expressed.full.bed -b $ANNO/exons.bed | sort -k1,1 -k2n,3n | python $CODE/expression/merge_subtracted.py > $ANNO/expressed_introns_no_exons.bed
for LIB in $LIBRARIES; do
    python $CODE/expression/get_subtract_expression.py $READS_DIR/$LIB.bam $ANNO/expressed_introns_no_exons.bed > $EXPRESS_DIR/$LIB.noexon.bed
done
cat $EXPRESS_DIR/*.noexon.bed | sort -k1,1 -k2n,3n | python $CODE/expression/merge.py > $EXPRESS_DIR/all.noexon.bed

# Mask repeats
python $CODE/expression/mask_repeats.py $EXPRESS_DIR/all.noexon.bed $ANNO/$REPEATS > $EXPRESS_DIR/all.masked.bed

##############################################################################
#################             RatchetJunction            #####################
##############################################################################

for LIB in $LIBRARIES; do
    python $CODE/RatchetJunction/putative_ratchet_sjr.py $READS_DIR/$LIB.bam $ANNO/introns.bed $ANNO/$FASTA > $SJR_DIR/$LIB.sjr.bed
    sort -k1,1 -k2,3n $SJR_DIR/$LIB.sjr.bed | python $CODE/RatchetJunction/group_introns.py sjr_$LIB > $SJR_DIR/$LIB.sjr_groups.bed
done

cat $SJR_DIR/*.sjr_groups.bed | sort -k1,1 -k2n,3n | python $CODE/RatchetJunction/merge_sjr_reps.py  > $SJR_DIR/all.sjr_groups.bed

##############################################################################
#################               RatchetPair              #####################
##############################################################################

for LIB in $LIBRARIES; do
    python $CODE/RatchetPair/straddle_jxn.py $READS_DIR/$LIB.bam $ANNO/introns.bed $ANNO/$FASTA pe_$LIB > $SJR_DIR/$LIB.straddle.bed
done
cat $SJR_DIR/*.straddle.bed > $SJR_DIR/all.straddle.bed

python $CODE/RatchetPair/straddle_gem.py $SJR_DIR/all.straddle.bed $SJR_DIR/all.sjr_groups.bed $ANNO/introns.bed $ANNO/$FASTA | sort -k1,1 -k2,4n > $SJR_DIR/all.straddle_gem.bed

##############################################################################
#################               RatchetScan                   ################
##############################################################################
# Note that this first script will have the longest runtime of the whole pipeline
# It could trivially be parallelized by splitting the all.masked.bed file into
# smaller files and then merging the results.

python $CODE/RatchetScan/run_mcmc.py $EXPRESS_DIR/all.masked.bed > $SAWTOOTH_DIR/mcmc.bed
python $CODE/RatchetScan/call_sites.py $SAWTOOTH_DIR/all.mcmc.bed $ANNO/$FASTA $ANNO/introns.bed > $SAWTOOTH_DIR/sites.bed
sort -k1,1 -k2,3n $SAWTOOTH_DIR/sites.bed | python $CODE/RatchetScan/merge_peak_calls.py > $SAWTOOTH_DIR/merged_sites.bed
