sh sequence/included_exons/pipeline.sh ../data/downloaded/dmel-allElements-r5.57.gtf ../data/five_included ../data/downloaded/dmel-all-chromosome-r5.57.fasta ../data/Adelman_4sU_RNA-seq_5min_merged.bam

echo 'five'

sh sequence/included_exons/pipeline.sh ../data/downloaded/dmel-allElements-r5.57.gtf ../data/ten_included ../data/downloaded/dmel-all-chromosome-r5.57.fasta ../data/Adelman_4sU_RNA-seq_10min_merged.bam

echo 'ten'

sh sequence/included_exons/pipeline.sh ../data/downloaded/dmel-allElements-r5.57.gtf ../data/twenty_included ../data/downloaded/dmel-all-chromosome-r5.57.fasta ../data/Adelman_4sU_RNA-seq_20min_merged.bam

echo 'twenty'

sh sequence/included_exons/pipeline.sh ../data/downloaded/dmel-allElements-r5.57.gtf ../data/total_included ../data/downloaded/dmel-all-chromosome-r5.57.fasta ../data/Adelman_4sU_RNA-seq_total_merged.bam