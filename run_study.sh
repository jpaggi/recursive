LONG_INTRONS='../annotations/dmel/long-introns-r5.57.bed'
BAM_DIR='../../athma/Adelman/timecourse/subsampled/'$1
EXONS='../annotations/dmel/dmel-Exons.gtf'
REPEATS='../annotations/dmel/repeats.txt'
EXPRESSION_DIR='../subsamples/expression/'$1
DENSITY_DIR='../subsamples/read_density/'$1
BAM_SUFFIX=$1'.bam'
GENOME='../annotations/dmel/downloaded/dmel-all-chromosome-r5.57.fasta'
SHORT_SS='../annotations/dmel/dmel-anno-splice-sites.ss'
SJR_DIR='../subsamples/sjr/'$1
MCMC_DIR='../subsamples/mcmc/'$1
GRAV='../annotations/dmel/graveley.bed'
OUT='../subsamples/output/'$1

cd coverage
# < 20 minutes
mkdir $EXPRESSION_DIR
mkdir $DENSITY_DIR
rm $EXPRESSION_DIR/*
rm $DENSITY_DIR/*

sh prepare_coverage.sh $LONG_INTRONS $BAM_DIR $EXONS $REPEATS $EXPRESSION_DIR $DENSITY_DIR $BAM_SUFFIX

cd ../sjr_pipeline
# < 10 minutes
mkdir $SJR_DIR
rm $SJR_DIR/*
sh run_all.sh $LONG_INTRONS $GENOME $EXPRESSION_DIR'/merged.bed' $SHORT_SS $SJR_DIR $BAM_DIR $BAM_SUFFIX

cd mcmc_core
mkdir $MCMC_DIR
rm $MCMC_DIR/*
sh run.sh $DENSITY_DIR $MCMC_DIR $GENOME $SHORT_SS

cd ../combine
python standard_table.py $SJR_DIR/'all.bed' $MCMC_DIR/'sites_all.bed' $GRAV $EXONS $GENOME $SHORT_SS $EXPRESSION_DIR'/merged.bed' > $OUT
