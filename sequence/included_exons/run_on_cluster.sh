BEGIN="sh sequence/included_exons/pipeline.sh genomes/dmel/downloaded/dmel-allElementsCDS-r5.57.gtf rs_exons/"
MIDDLE=" genomes/dmel/downloaded/dmel-all-chromosome-r5.57.fasta reads/timecourse/"

SAMPLES="5_repA
5_repB
5_repC
10_repA
10_repB
10_repC
20_repA
20_repB
20_repC
total_repA
total_repB"

for sample in $SAMPLES
do
	$BEGIN$sample$MIDDLE$sample'.bam'
done