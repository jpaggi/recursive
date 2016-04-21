import sys
import subprocess

# ran with command                                                                                                                                                                 
# python splicing_rates/process_annotations/get_unique_non_overlapping_introns.py /net/utr/data/atf/jpaggi/genomes/dmel/dmel-sorted-introns-r5.57.bed
#/net/utr/data/atf/jpaggi/gen
#omes/dmel/dmel-unique-introns-r5.57.bed /net/utr/data/atf/jpaggi/genomes/dmel/dmel-nonoverlapping-introns-r5.57.bed                                                                

# used command                                                                                                                                                                     
# bedtools intersect -v -wa -a  dmel-nonoverlapping-introns-r5.57.bed -b dmel-allElementsCDS-r5.57.gtf > dmel-introns-constituative-r5.57.bed                                      
# to check for overlaps with exons                                                                                                                                                 


introns = open(sys.argv[1], 'r')
unique_introns = open(sys.argv[2], 'w')
non_overlapping_introns = open(sys.argv[3], 'w')

unique = subprocess.Popen("sort -u -k1,1 -k2,2n -k3,3n", stdin = introns, stdout = unique_introns, shell = True)

unique.wait()
unique_introns.close()
unique_int = sys.argv[2]

counts = subprocess.Popen("bedtools intersect -c -s -a " + unique_int + " -b " + unique_int, stdout = subprocess.PIPE, shell = True)

non_overlapping = subprocess.Popen("grep -P '\t1$'", stdin = counts.stdout, stdout = subprocess.PIPE, shell = True)

final = subprocess.Popen("cut --complement -f7", stdin = non_overlapping.stdout, stdout = non_overlapping_introns, shell = True)

final.wait()
non_overlapping_introns.close()