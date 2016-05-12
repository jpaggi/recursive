# 1 sorted, indexed bam file of reads
# 2 output directory name
# 3 sample name
# 4 splice site file
# 5 fasta file of genome
# 6 intron expression data
# 7 sjr file

# Gives directory of this file...
# copied from http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#python straddle_jxn.py $1 $4 $3 > $2'/straddle.bed'

echo 'Finished extracting splice junction reads'

python straddle_assign.py $2'/straddle.bed' $7 $4 $5 | sort -k1 -k2n,3n > $2'/straddle_groups.bed'

echo 'made groups'

python seq_for_groups.py $2'/straddle_groups.bed' $5 $2'/straddle_seq.bed' $4

echo 'got recursive site sequences'

python assign_introns.py $2'/straddle_seq.bed' $6 > $2'/straddle_introns.bed'

echo 'assigned to introns'