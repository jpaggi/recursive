# 1 intron file
# 2 BAM directory
# 3 exon file
# 4 repeat file
# 5 full expression output
# 6 final output 
# 7 bam suffix

# create expression file
sh get_all_expression.sh $1 $5 $2 $7

# compute exon-free expression file
sh remove_exons.sh $5'/merged.bed' $6 $3 $2 $4 $7