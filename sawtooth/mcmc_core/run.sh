# 1 expression directory name
# 2 output directory name

for bin in 8_15 15_30 30_plus
do
	python run_mcmc.py $1/$bin'_all_merged_masked.bed' > $2'/mcmc_'$bin'.bed'

	python call_sites.py $2'/mcmc_'$bin'.bed' > $2'/sites_'$bin'.bed'
done

cat $2'/sites_'* | sort -k1,1 -k2,3n | python merge_peak_calls.py > $2'/sites_all.bed'