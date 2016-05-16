
cat $1 | python sequence/insert/get_long_exons.py $2 | sort -k1,1 -k4,4n > temp_long_exons.gtf

grep $'\texon\t' $1 | sort -k1,1 -k4,4n > temp_exons.gtf

bedtools intersect -sorted -c -s -a temp_long_exons.gtf -b temp_exons.gtf | grep $'\t1$' | cut -f1-9

rm temp_exons.gtf
rm temp_long_exons.gtf
