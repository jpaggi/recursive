# 1 sorted, indexed bam file of reads
# 2 output directory name
# 3 sample name
# 4 splice site file
# 5 fasta file of genome
# 6 intron expression data

# Gives directory of this file...
# copied from http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir $2

#sh split_reads.sh $1 $2 'sjr_'$3 $4 $5 $6

sh straddle_reads.sh $1 $2 'pe_'$3 $4 $5 $6 ../reads/jxns/all_sjr_seq.bed