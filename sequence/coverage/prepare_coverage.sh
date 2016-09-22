# 1 intron file
# 2 BAM directory
# 3 exon file
# 4 repeat file
# 5 full expression output
# 6 final output 

# create expression file

sh get_all_expression.sh $1 $2 $5

# compute exon-free expression file

sh remove_exons.sh $5'/merged.bed' $6 $3 $2 $4