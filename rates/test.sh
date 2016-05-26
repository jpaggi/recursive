
rm ../data/test_grav.* ../data/test_annotations.gff
rm -rf ../data/indexed

python rates/paired_end.py ../data/test_introns.bed ../data/downloaded/dmel-Exons.gtf

echo 'wrote sam plus annotations'

samtools view -bS ../data/test_grav.sam | samtools sort - -o ../data/test_grav.bam

echo 'sorted bam'

samtools index ../data/test_grav.bam ../data/test_grav.bam.bai

echo 'indexed bam'

index_gff --index ../data/test_annotations.gff ../data/indexed/

echo 'indexed annotations'

rm -rf ../data/test_miso ../data/summarize_output

miso --run ../data/indexed/ ../data/test_grav.bam --paired-end 220 80 --read-len 51 --output-dir ../data/test_miso

echo 'ran miso'

summarize_miso --summarize-samples ../data/test_miso/ ../data/summary_output/

echo 'summarized miso'

