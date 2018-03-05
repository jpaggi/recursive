for SAMPLE in 1 2 3 4 5
do
	sh run_study.sh 'sub0.'$SAMPLE &
done
wait
