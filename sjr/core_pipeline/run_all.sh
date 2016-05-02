ss="../genomes/dmel/dmel-splice-sites-r5.57.ss"
fa="../genomes/dmel/downloaded/dmel-all-chromosome-r5.57.fasta"
exp="../reads/coverage/all_merged.bed"

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
    echo run.sh ../reads/timecourse/$BAM.bam ../reads/jxns/$BAM $BAM $ss $fa $exp
    sh run.sh ../reads/timecourse/$BAM.bam ../reads/jxns/$BAM $BAM $ss $fa $exp
done

