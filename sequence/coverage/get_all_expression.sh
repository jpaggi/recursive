# 1 introns bed file
# 2 output directory
# 3 bam directory


LIBRARIES='5_repA
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

for lib in $LIBRARIES
do
	python get_intron_expression.py $3$lib'.bam' $1 > $2'/'$bin'_'$lib'.bed'
done


cat $2/* | sort -k1,1 -k2,2n -k3,3n | merge.py > $2'/merged.bed'
