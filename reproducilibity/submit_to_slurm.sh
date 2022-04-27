cd ..

for DATASET in $(find simulated_analysis/data_Liu/ -maxdepth 1 -type f -printf "%f\n")
do
	for METHOD in "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR"  "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent" "muscat_MM" "MAST_RE"
	#for METHOD in "pseudobulk_ROTS_sum" "pseudobulk_ROTS_mean" "pseudobulk_Limma_sum" "pseudobulk_Limma_mean" "pseudobulk_edgeR_sum" "pseudobulk_DESeq2_sum" "Seurat_wilcoxon" "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR"  "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent" "muscat_MM" "NEBULA-LN" "MAST_RE"
	do
		OUTPUT_FILE=simulated_analysis/reproduciblity_Liu_results/${METHOD}_${DATASET}
		if [ ! -f "$OUTPUT_FILE" ]; then
    			echo "$OUTPUT_FILE does not exist."
		
			sbatch -J R --mem 100000 -p long --wrap="module add R/4.1.0 ; Rscript --vanilla simulated_analysis/run_DE_method_with_args_simulated_reproducibility_analysis.R $METHOD simulated_analysis/data_Liu/$DATASET simulated_analysis/reproduciblity_Liu_results/${METHOD}_${DATASET}"
			#exit 0
		fi
		#exit 0
	done
done

cd simulated_analysis
