# for METHOD in "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR" "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent"
# for DATASET in "Combes_B_cells.RData" "Lee_B_cells.RData" "Liu_B_cells.RData" "Combes_CD14_Mono.RData" "Lee_CD14_Mono.RData" "Liu_CD14_Mono.RData" "Combes_NK.RData" "Lee_NK.RData" "Liu_NK.RData"

cd ..

for DATASET in "Liu_B_cells.RData"
do
	for METHOD in "Seurat_wilcoxon" "Seurat_MAST" "Seurat_MAST_latent" "Seurat_LR"  "Seurat_LR_latent" "Seurat_negbinom" "Seurat_negbinom_latent" "Seurat_poisson" "Seurat_poisson_latent" "pseudobulk_ROTS_sum" "pseudobulk_ROTS_mean" "pseudobulk_Limma_sum" "pseudobulk_Limma_mean" "pseudobulk_edgeR_sum" "pseudobulk_DESeq2_sum" "MAST_RE" "MAST_RE_no_ngeneson" "muscat_MM"
	do
		for SEED in `seq 1 30`
		do
			sbatch -J R --mem 200000 -p normal --wrap="module add R/4.1.0 ; Rscript --vanilla covid_mock_analysis/run_DE_method_with_args_covid_mock_analysis.R $METHOD covid_mock_analysis/$DATASET covid_mock_analysis/results/${METHOD}_${SEED}_${DATASET} $SEED"
			#exit 0
		done
	done
done

cd covid_mock_analysis
