for SAMPLE in 1 2 3 4 5
do
	sh run_study.sh 'sub0.'$SAMPLE &
done
wait

# for SAMPLE in 6 7 8 9
# do
# 	sh run_study.sh 'sub0.'$SAMPLE &
# done
# wait