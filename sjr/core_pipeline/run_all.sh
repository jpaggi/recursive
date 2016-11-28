ss=$1 # 1 long intron file #ss="../genomes/dmel/long_introns.bed"
fa=$2 # 2 genome sequence fa="../genomes/dmel/downloaded/dmel-all-chromosome-r5.57.fasta"
exp=$3 # 3 expression data #exp="../reads/coverage/all_merged.bed"
short_ss=$4 # 4 splice site file #short_ss="../genomes/dmel/dmel-splice-sites-r5.57.ss"
DIR=$5 # 5 output directory # DIR="../reads/junctions"
IN_DIR=$6 # 6 input directory
SUFFIX=$7


# Adelman_4sU_RNA-seq_10min_rep1.sub0.2.bam

samples='5min_rep1
5min_rep2
5min_rep3
10min_rep1
10min_rep2
10min_rep3
20min_rep1
20min_rep2
20min_rep3
total_rep1
total_rep2'

for BAM in $samples
do
    sh   run.sh $IN_DIR/Adelman_4sU_RNA-seq_$BAM'.'$SUFFIX $DIR/$BAM $BAM $ss $fa $short_ss
done
wait

cat $DIR/*/sjr_groups.bed | sort -k1,1 -k2,3n | python merge_sjr_reps.py  > $DIR'/all_sjr_groups.bed'

python assign_introns.py $DIR'/all_sjr_groups.bed' $exp > $DIR'/all_sjr_introns.bed'

cat $DIR/*/straddle.bed | python straddle_gem.py $DIR'/all_sjr_introns.bed' $short_ss $fa | sort -k1,1 -k2,4n > $DIR'/temp_straddle.bed'

python assign_introns.py $DIR'/temp_straddle.bed' $exp > $DIR'/all_straddle.bed'

rm $DIR'/temp_straddle.bed' $DIR'/all_sjr_groups.bed'

cat $DIR/'all_sjr_introns.bed' $DIR/'all_straddle.bed' | sort -k1,1 -k2,4n | python merge_introns.py > $DIR'/all.bed'
