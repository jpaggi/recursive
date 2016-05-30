# source of example unaligning reads
#hisat2 -p 12 -u 1000000 -x indexes/dmel5.57/dmel5.57 -1 ../athma/Adelman/timecourse/recursive/TH_5_A_L1.D702_501_1.fastq.gz -2 ../athma/Adelman/timecourse/recursive/TH_5_A_L1.D702_501_2.fastq.gz -S temp/aligning.sam --un-conc temp/unaligning%.fq


# # align both ends individually
# hisat2 --sp 1000,999 -x ~/Downloads/dmel-r5.57 -U $1 -S data/branchpoints/end1.sam
# hisat2 --sp 1000,999 -x ~/Downloads/dmel-r5.57 -U $2 -S data/branchpoints/end2.sam

# # sort each by ID
# sort -k1,1 data/branchpoints/single_end1.sam > data/branchpoints/sorted_single_end1.sam
# sort -k1,1 data/branchpoints/single_end2.sam > data/branchpoints/sorted_single_end2.sam

# # merge each file together using merge_alignments.py
# python code/branchpoints/merge_alignments.py data/branchpoints/sorted_single_end1.sam > data/branchpoints/single_end_merged1.joe
# python code/branchpoints/merge_alignments.py data/branchpoints/sorted_single_end2.sam > data/branchpoints/single_end_merged2.joe

# separate out discordant and reads for which only one end aligns
python code/branchpoints/merge_ends.py data/branchpoints/single_end_merged1.joe data/branchpoints/single_end_merged2.joe data/branchpoints/test_discordant.joe data/branchpoints/test1.fq data/branchpoints/test2.fq

# align broken reads
# hisat2 --rf --no-discordant  --no-mixed -x ~/Downloads/dmel-r5.57 -1 data/branchpoints/test1.fq -2 data/branchpoints/test2.fq -S data/branchpoints/pe_findbps_aligning.sam

# # sort by ID (could use --reorder instead)
# sort -k1,1 data/branchpoints/pe_findbps_aligning.sam >  data/branchpoints/pe_sorted_findbps.sam 

# # extend alignments only working for broken read2 on + strand
# python code/branchpoints/parse_initial_mapping.py data/branchpoints/pe_sorted_findbps.sam data/downloaded/dmel-all-chromosome-r5.57.fasta


# hisat2 very aggressively soft trims reads by default!!! should disable this behaviour for both the inital 
# mapping and the intermediate mapping

# also seems to be missing some alignments???? 

# maybe should switch to bowtie for later alignments???

# 1 create test cases from reads found with findbps
# 2 figure out best (efficiency/correctness) way to get unmapped reads
# 3 think about what software to use for alignment of each end
# 4 software / settings to use for PE alignment
# 5 alignment requirements


# hisat2 --no-discordant --no-mixed -p 12 -u 500000 --sp 1000,999 -x indexes/dmel5.57/dmel5.57 -1 ../athma/Adelman/timecourse/recursive/TH_5_A_L1.D702_501_1.fastq.gz -2 ../athma/Adelman/timecourse/recursive/TH_5_A_L1.D702_501_2.fastq.gz -S test_align.sam --un-conc unaligning%.fq




