
#normal alignment!
# might want to increase the max insert size here?
# how does this work with transcriptome mapping???
hisat2 -p 12 <index> -1 <in1.fq> -2 <in2.fq> --un-conc unaligning%.fq



# get all discordantly aligning reads
# set max insert length to be longest expected intron length
hisat2 -p 12 --no-mixed --no-discordant -X 100000 <index> -1 unaligning2.fq -2 unaligning1.fq --un-conc nodiscordant%.fq

# split reads such that LLLLLLLXXXXXXXXXXXXX--------------YYYYYYYYYYYYYYYYRRRRRRRRRRRRRRRRRRR
# become:
# XXXXXXXXXXXXXXX      YYYYYYYYYYYYRRRRRRRRRRRRRRR
# YYYYYYYYYYYYYYRRRRRRRRRRRRRR     LLLLLLLLLLLLLLLLLL
#
# LLLLLLLLLLLLLLLLLLLLLXXXXXXXXXXXXXXXX     YYYYYYYYYYYYYYYYYYYY
# RRRRRRRRRRRRR   LLLLLLLLLLLLLXXXXXXXXXXXXXXXXX
python seperate.py no_discordant split_reads

hisat2 -p 12 --no-discordant <index> -1 split_reads1.fq -2 split_reads2.fq --un nosplit%.fq


# run findbps on nosplit%.fq reads

# process split reads
# extend using linear algorithm
# see if GT or GC close to where 5' end extends to...



#align normally
#hisat2 -u 1000000 -x Downloads/dmel-r5.57 -1 Downloads/TH_5_B_L1_1.fastq.gz -2 Downloads/TH_5_B_L1_2.fastq.gz -S Documents/recursive_splicing/data/five_minute_repB_aligning.sam --un-conc Documents/recursive_splicing/data/five_minute_repB_unaligning%.fq

#align discordantly

hisat2 -X 1500000 --rf --no-discordant --no-mixed -x Downloads/dmel-r5.57 -1 Documents/recursive_splicing/data/five_minute_repB_unaligning2.fq -2 Documents/recursive_splicing/data/five_minute_repB_unaligning1.fq -S Documents/recursive_splicing/data/five_minute_repB_discordant.sam --un-conc Documents/recursive_splicing/data/five_minute_repB_no_discordant%.fq
