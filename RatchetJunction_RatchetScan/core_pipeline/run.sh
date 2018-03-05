# 1 sorted, indexed bam file of reads
# 2 output directory name
# 3 sample name
# 4 long intron splice site file
# 5 fasta file of genome
# 6 all intron splice site file

mkdir $2

python straddle_jxn.py $1 $6 'pe_'$3 $5 > $2'/straddle.bed' &

python putative_ratchet_sjr.py $1 $4 $5 > $2'/sjr.bed'
cat $2'/sjr.bed' | sort -k1,1 -k2n,3n | python group_introns.py - 'sjr_'$3 > $2'/sjr_groups.bed'

wait
