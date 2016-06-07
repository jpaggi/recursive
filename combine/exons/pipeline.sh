
# Goal is to arrive at a set of exons that are either expressed or do not 
# contain a recursive splice site at thier 5' end.

# These will be used to filter out sjr otherwise fitting a rs pattern from consideration.

# Information you plan to use includes:
#  1) The relative expression levels in 5 minute replicates vs. total replicates
#  2) The presence of sjr leaving from the 3' end of introns
#  3) The presence of a RS motif
#  4) Whether the exon is a first, middle or last exon


# $1 annotation file in gtf format
# $2 output directory
# $3 fasta file
# $4 alignment file

mkdir $2

echo $1
echo $2
echo $3
echo $4

# # Extract introns over 1000 bps
sort -k10,10 $1 | python sequence/get_introns.py | sort -u -k1,1 -k2,3n > $2'/long_introns.bed'

# # Extract exons for consideration, intersect with long introns
sort -k10,10 $1 | python sequence/included_exons/exons_to_bed.py | sort -u -k1,1 -k2,3n -k4,4 > $2'/exons.bed'

bedtools intersect -u -s -a $2'/exons.bed' -b $2'/long_introns.bed' > $2'/overlapping_exons.bed'

# # Get sjr from exon
python sequence/included_exons/get_expression.py $2/'overlapping_exons.bed' $4  > $2'/sjr.bed'

# # # Look for motifs.
python sequence/included_exons/motif_checker.py $2'/sjr.bed' $3 $2'/motif.bed' $2'/no_motif.bed'

# # Compare expression between 5 and total samples
# python sequence/included_exons/compare_expression.py motif_expression motif_expressed.bed motif_not_expressed.bed
