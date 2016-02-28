mkdir $2

python putative_ratchet_sjr.py $1 $2'/sjr.sam' $2'/sjr.bed'

echo 'Finished extracting splice junction reads'

cat $2'/sjr.bed' | sort -k1,1 -k2n,3n | python group_introns.py - $2'/groups.bed' $3

echo 'made groups'

python seq_for_groups.py $2'/groups.bed' ~/Downloads/dmel-all-chromosome-r5.57.fasta $2'/seq.bed'

echo 'got recursive site sequences'