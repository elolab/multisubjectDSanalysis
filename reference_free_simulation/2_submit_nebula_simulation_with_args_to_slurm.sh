cd ..

RESULT_DIR="simulated_analysis/nebula_data"

i=1
for sig2 in 0.05 0.10 0.20 0.50
do
	for sig3 in 0.1 1 10 100
	do
		sbatch -J R --mem 50000 -p normal -t 720 --wrap="module add R/4.1.0 ; Rscript --vanilla simulated_analysis/simulate_nebula_with_args.R $sig2 $sig3 $i $RESULT_DIR"
		i=$((i + 1))
		#exit 0
	done
done

cd simulated_analysis
