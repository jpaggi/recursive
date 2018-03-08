# 1 expression file
# 2 output directory
# 3 exon file
# 4 bam directory
# 5 repeat file
# 6 bam suffix

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

# Get expressed introns.
cat $1 | python expressed.py 0 8000  15000   | cut -f 1-6 > $2'/8_15.bed'
cat $1 | python expressed.py 0 15000 30000   | cut -f 1-6 > $2'/15_30.bed'
cat $1 | python expressed.py 0 30000 9999999 | cut -f 1-6 > $2'/30_plus.bed'

for bin in 8_15 15_30 30_plus
do
    # Remove exonic regions from introns (otherwise this could cause false positives in RatchetScan)
    bedtools subtract -s -a $2'/'$bin'.bed' -b $3 | sort -k1 -k2n,3n | python merge_subtracted.py > $2'/'$bin'_subtract.bed'
    
    # Recompute expression values with the exonic portions removed
    for lib in $LIBRARIES
    do
	python get_subtract_expression.py $4'/Adelman_4sU_RNA-seq_'$lib'.'$6 $2'/'$bin'_subtract.bed' > $2'/'$bin'_'$lib'.bed'
    done

    # Process each timepoint seperately
    for time in 5 10 20 total
    do
	f=$2'/'$bin'_'$time'*'
	# Merge replicates for each timepoint.
	cat $f | sort -k1,1 -k2n,3n | python merge.py > $2'/'$bin'_'$time'_merged.bed'
	# Replace expression levels in repeat regions by average of neighboring regions.
	python mask_repeats.py $2'/'$bin'_'$time'_merged.bed' $5 > $2'/'$bin'_'$time'_merged_masked.bed'
    done
    
    # Process merged timepoints. Steps same as above
    f=$2'/'$bin'*merged.bed'
    cat $f | sort -k1,1 -k2n,3n | python merge.py > $2'/'$bin'_all_merged.bed'
    python mask_repeats.py $2'/'$bin'_all_merged.bed' $5 > $2'/'$bin'_all_merged_masked.bed'
done