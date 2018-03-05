# 1 introns bed file
# 2 output directory
# 3 bam directory
# 4 bam suffix


LIBRARIES='5min_rep1
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

for lib in $LIBRARIES
do
	python get_intron_expression.py $3/Adelman_4sU_RNA-seq_$lib'.'$4 $1 > $2'/'$lib'.bed'
done
wait

cat $2/* | sort -k1,1 -k2,2n -k3,3n | python merge.py > $2'/merged.bed'
