cd ..

for DATASET in $(find simulated_analysis/Liu/ -maxdepth 1 -type f -printf "%f\n")
do
	#for METHOD in "pseudobulk_ROTS_sum" "pseudobulk_ROTS_mean" "pseudobulk_Limma_sum" "pseudobulk_Limma_mean" "pseudobulk_edgeR_sum" "pseudobulk_DESeq2_sum" "Seurat_wilcoxon" "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR"  "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent" "NEBULA-LN" "MAST_RE" "muscat_MM"
	do
		sbatch -J R --mem 100000 -p normal -t 720 --wrap="module add R/4.1.0 ; Rscript --vanilla simulated_analysis/run_DE_method_with_args_simulated_analysis.R $METHOD simulated_analysis/Liu/$DATASET simulated_analysis/results_Liu/${METHOD}_${DATASET}"
		#exit 0
	done
done

cd simulated_analysis
