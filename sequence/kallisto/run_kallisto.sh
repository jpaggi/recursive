#!/bin/bash'

BASE="quant -i genomes/dmel/downloaded/kallisto_index -o run_kallisto/"
PATH="../athma/Adelman/timecourse/recursive/"

echo kallisto $BASE'total_repA' $PATH'LZ1.1.fastq' $PATH'LZ1.2.fastq&'
echo kallisto $BASE'total_repB' $PATH'LZ2.1.fastq' $PATH'LZ2.2.fastq&'
 
echo kallisto $BASE'5_repA' $PATH'TH_5_A_L1.D702_501_1.fastq.gz' $PATH'TH_5_A_L1.D702_501_2.fastq.gz&'
echo kallisto $BASE'5_repB' $PATH'TH_5_B_L1_1.fastq.gz' $PATH'TH_5_B_L1_2.fastq.gz&'
echo kallisto $BASE'5_repC' $PATH'TH_5_C_L1.D704_501_1.fastq.gz' $PATH'TH_5_C_L1.D704_501_2.fastq.gz&'

echo kallisto $BASE'10_repA' $PATH'TH_10_A_L1.D705_501_1.fastq.gz' $PATH'TH_10_A_L1.D705_501_2.fastq.gz&'
echo kallisto $BASE'10_repB' $PATH'TH_10_B_L1.D706_501_1.fastq.gz' $PATH'TH_10_B_L1.D706_501_2.fastq.gz&'
echo kallisto $BASE'10_repC' $PATH'TH_10_C_L1.D707_501_1.fastq.gz' $PATH'TH_10_C_L1.D707_501_2.fastq.gz&'

echo kallisto $BASE'20_repA' $PATH'TH_20_A_L1.D708_1.fastq.gz' $PATH'TH_20_A_L1.D708_2.fastq.gz&'
echo kallisto $BASE'20_repB' $PATH'TH_20_B_L1.D709_1.fastq.gz' $PATH'TH_20_B_L1.D709_2.fastq.gz&'
echo kallisto $BASE'20_repC' $PATH'TH_20_C_L1.D711_1.fastq.gz' $PATH'TH_20_C_L1.D711_2.fastq.gz&'
