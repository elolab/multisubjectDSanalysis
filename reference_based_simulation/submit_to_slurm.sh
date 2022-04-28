cd ..

for DATASET in $(find simulated_analysis/data_main_comparison/ -maxdepth 1 -type f -printf "%f\n")
do
	for METHOD in "pseudobulk_ROTS_sum" "pseudobulk_ROTS_mean" "pseudobulk_Limma_sum" "pseudobulk_Limma_mean" "pseudobulk_edgeR_sum" "pseudobulk_DESeq2_sum" "Seurat_wilcoxon" "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR"  "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent" "muscat_MM" "NEBULA-LN" "MAST_RE"
	do
		OUTPUT_FILE=simulated_analysis/results_main_comparison/${METHOD}_${DATASET}
		if [ ! -f "$OUTPUT_FILE" ]; then
    			echo "$OUTPUT_FILE does not exist."
		
			sbatch -J R --mem 50000 -p normal -t 720 --wrap="module add R/4.1.0 ; Rscript --vanilla simulated_analysis/run_DE_method_with_args_simulated_analysis.R $METHOD simulated_analysis/data_main_comparison/$DATASET simulated_analysis/results_main_comparison/${METHOD}_${DATASET}"
			#exit 0
		fi
		#sbatch -J R --mem 75000 -p long -t 720 --wrap="module add R/4.1.0 ; Rscript --vanilla simulated_analysis/run_DE_method_with_args_simulated_analysis.R $METHOD simulated_analysis/data_main_comparison/$DATASET simulated_analysis/results_main_comparison/${METHOD}_${DATASET}"
		#exit 0
	done
done

cd simulated_analysis
