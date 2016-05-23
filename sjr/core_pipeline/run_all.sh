ss="../genomes/dmel/long_introns.bed"
fa="../genomes/dmel/downloaded/dmel-all-chromosome-r5.57.fasta"
exp="../reads/coverage/all_merged.bed"
short_ss="../genomes/dmel/dmel-splice-sites-r5.57.ss"
DIR="../reads/rs_jxns"

samples='5_repA
5_repB
5_repC
10_repA
10_repB
10_repC
20_repA
20_repB
20_repC
total_repA
total_repB'

for BAM in $samples
do
    echo run.sh ../reads/timecourse/$BAM.bam $DIR/$BAM $BAM $ss $fa $short_ss
    sh   run.sh ../reads/timecourse/$BAM.bam $DIR/$BAM $BAM $ss $fa $short_ss &
done
wait

cat $DIR/*/sjr_groups.bed | sort -k1,1 -k2,3n | python merge_sjr_reps.py  > $DIR'/all_sjr_groups.bed'

python assign_introns.py $DIR'/all_sjr_groups.bed' $exp > $DIR'/all_sjr_introns.bed'

cat $DIR/*/straddle.bed | python straddle_gem.py $DIR'/all_sjr_introns.bed' $short_ss $fa | sort -k1,1 -k2,4n > $DIR'/temp_straddle.bed'

python assign_introns.py $DIR'/temp_straddle.bed' $exp > $DIR'/all_straddle.bed'

rm $DIR'/temp_straddle.bed' $DIR'/all_sjr_groups.bed'

cat $DIR/'all_sjr_introns.bed' $DIR/'all_straddle.bed' | sort -k1,1 -k2,4n | python merge_introns.py > $DIR'/all.bed'
