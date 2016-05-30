# 1) introns bed
# 2) exons gtf
# 3) output base
# 4) bam base


samples="5_repA
5_repB
5_repC
10_repA
10_repB
10_repC
20_repA
20_repB
20_repC"

for sample in $samples
do
	sh run.sh $1 $2 $4/$sample'.bam' $3/$sample 220 80 51
done


samples="total_repA
total_repB"

for sample in $samples
do
	sh run.sh $1 $2 $4/$sample'.bam' $3/$sample 320 80 101
done

python extract_psi.py $3 $1 > $3/table.tsv


