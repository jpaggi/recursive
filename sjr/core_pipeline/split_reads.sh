# 1 sorted, indexed bam file of reads
# 2 output directory name
# 3 sample name
# 4 splice site file
# 5 fasta file of genome
# 6 intron expression data

# Gives directory of this file...
# copied from http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python $DIR/putative_ratchet_sjr.py $1 $2'/sjr.bed' $4

echo 'Finished extracting splice junction reads'

cat $2'/sjr.bed' | sort -k1,1 -k2n,3n | python $DIR/group_introns.py - $2'/sjr_groups.bed' $3

echo 'made groups'

python $DIR/seq_for_groups.py $2'/sjr_groups.bed' $5 $2'/sjr_seq.bed' $4

echo 'got recursive site sequences'

python $DIR/assign_introns.py $2'/sjr_seq.bed' $6 > $2'/sjr_introns.bed'

echo 'assigned to introns'