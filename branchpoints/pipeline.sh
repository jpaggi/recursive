
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
