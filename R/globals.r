# ============================================================================
# R/globals.r
# Declares data.table and ggplot2 NSE column names plus data.table
# special symbols, so R CMD check does not flag them as unbound globals.
# ============================================================================

utils::globalVariables(c(
  # ---- data.table special symbols ----
  ".", ":=", ".SD", ".N", ".I", ".GRP", ".BY", ".EACHI",

  # ---- parent-scope variable references (double-dot prefix) ----
  "..cols", "..cols_chunk", "..dedup_cols", "..keep", "..keep_cols",
  "..pub_cols", "..geno_cols", "..sample_cols",

  # ---- genotype / variant table ----
  "SNP", "CHR", "POS", "REF", "ALT",
  "CHR_original", "REF_original", "ALT_original",
  "CHR_normalized", "is_biallelic_snp", "variant_key", "CHR_N",
  "Chr", "Start", "End", "Strand", "Width",
  "CDS_Pos_Start", "CDS_Pos_End", "Transcript_ID",
  "Pos_in_CDS", "Codon_Num", "Codon_Pos",
  "AA_Change", "AA_REF", "AA_ALT", "Codon_REF", "Codon_ALT",
  "Is_Nonsynonymous", "AA_Validated", "REF_Matches_FASTA",
  "Variant_Label", "Genomic_Order", "Region",
  "NS_Transcripts", ".rank", "liftover_status",
  "A1", "A2", "BP", "overlap", "dist", "poskey",

  # ---- haplotype table ----
  "Haplotype", "NT_Sequence", "AA_Sequence", "Count", "Frequency",
  "NT_Seq", "AA_Seq", "Mean", "Median", "SE", "N", "Letter",

  # ---- LD / regions ----
  "Original_GWAS_Row", ".SNP_ID_INTERNAL", ".CHR_INTERNAL", ".POS_INTERNAL",
  "Region_Start", "Region_End", "Region_Length", "GWAS_SNP_Count",
  "GWAS_SNPs", "GWAS_Positions", "Original_Window_IDs", "Original_GWAS_Rows",
  "Original_Window_ID", "LD_Region_ID",
  "LD_Region_Start", "LD_Region_End", "LD_Region_Length",
  "LD_Decay_bp", "LD_Decay_kb",
  "Cum_End", "Previous_Cum_End", "New_Region", "Region_Group",
  "Distance", "R2", "Distance_Bin",
  "Mean_R2", "Median_R2", "N_Pairs", "Min_Distance", "Max_Distance",
  "CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B",
  ".Region_Chr", ".Region_Start", ".Region_End",
  ".Gene_Chr", ".Gene_Start", ".Gene_End",
  ".join_start", ".join_end",
  "N_Unique_GWAS_SNPs", ".snps_usable", ".rowid",
  "n_positions", "n_chr", "n_pos", "bin", "GenomicPos", "ChrPos_Start",

  # ---- association test ----
  "Duncan_Group", "Tukey_Group", "Significance_Group",
  "Mean_Duncan", "Duncan_Group.new", "Tukey_Group.new",
  "hap", "mean_y", "n",

  # ---- superior-haplotype / pipeline summary ----
  "Is_Superior", "Is_Superior_Nominal", "Superior_Haplotype",
  "Superior_NT_Seq", "Superior_AA_Seq",
  "P_value", "Permutation_P", "Adjusted_P", "Cohens_d",
  "Improvement_Pct", "Superior_Mean", "Superior_Freq_Analysed_Pct",
  "permutation_n_valid", "test_res",
  "Trait", "Direction", "Mode", "Genes_Tested", "Superior_Haplotypes",
  "Output_Dir", "Schema_Version",
  "Significant_FDR", "Significant_Nominal", "Significant_Permutation",
  "Annotation_Source", "Exit_Stage", "i.Position", "ID",

  # ---- genetic map / bins ----
  "chr", "cm", "bp", "source",
  "length_cM", "n_markers", "n_chrom", "total_cM",
  "bin_of", "center", "snp_col",
  "gene_bin", "gene_chr", "gene_cm", "gene_first",

  # ---- cross design / block methods ----
  "units", "block", "carry", "ucarry",
  "n_target_haplotypes", "n_genes", "genes", "traits", "panel_freq",
  "start_cM", "end_cM",
  "n_blocks_fixed", "n_blocks_segregating", "n_blocks_split",
  "n_blocks_incomplete", "n_blocks_covered", "n_blocks_target",
  "n_blocks_full", "n_blocks_partial",
  "p_pyramid_DH", "p_pyramid_F2", "p_pyramid_F2_hom",
  "n_DH_95", "n_F2_95", "n_F2hom_95",
  "mid_index", "sd_index", "uc", "relatedness",
  "parent_a", "parent_b", "rank", "pa", "pb",     # pa/pb: hbr_plot_cross_landscape
  "unique_blocks", "set_blocks_covered", "n_blocks", "coverage",
  "truncation_coverage", "uncovered_blocks",
  "hap_info", "hap_eff", "w_eff", "bins", "pairs", "sel_frac", "h_sel",
  "feasible", "mean_relatedness",

  # ---- CV / comparison / prediction quality ----
  "method", "baseline", "metric", "n_genes",
  "boot_lo", "boot_hi", "mean_diff", "sig", "row_label",
  "comparison", "cmp_order",
  "naive", "honest", "optimism", "kind", "gain", "top_frac",
  "pred", "obs", "rel", "carrier", "value",
  "M0", "M2", "M3",
  "ntrain", "ntest", "n_features", "n_obs",
  "cross_id", "rank_mid", "rank_uc", "shift",     # hbr_plot_cross_ranking

  # ---- breeding decision ----
  "id", "index", "P_top15", "phenos", "index_SE",
  "g_M0", "g_M2", "rel_M0", "rel_M2", "PEV_M2",
  "sigma2_g", "sigma2_e", "w_hat",
  "p1", "p2", "complement", "n_regions", "score",
  "trait", "direction", "mode",
  "tested", "superior", "action", "mean_rel", "haplo_cov", "n_crosses",
  "pred_base", "pred_m", "is_superior",
  "Donor", "Complete_Region",
  "N_Superior_In_Region", "N_Superior_Genes_In_Region",
  "N_Total_Genes_In_Region", "n_traits",
  "parent", "panel", "selected",                  # hbr_plot_parent_profile / selected_set_coverage

  # ---- cross-annotation / publication table ----
  "Gene", "Gene_ID", "Gene_Function", "Gene_Chromosome",
  "Gene_Start", "Gene_End", "Gene_Strand", "Gene_Name",
  "Gene_Relationship", "Annotation_Status",
  "MSU7_Gene_ID", "RAPDB_Gene_ID", "MSU_Gene_ID",
  "Cross_Gene_ID", "Cross_GO_IDs", "Cross_GO_BP", "Cross_GO_MF",
  "Cross_GO_CC",
  "GWAS_SNP", "GWAS_Position", "Candidate_Gene_ID", "GWAS_SNP_Position",
  "N_GWAS_SNPs", "Chromosome", "i.SNP", "i.CHR", "i.POS", "i.REF",
  "REF_reference", "x.transcript_id", "x.cds_start",
  ".Gene_ID_Source", "annotation", "model",
  "Mean_TPM", "Max_TPM", "N_Tissues_Expressed", "Expression_Breadth",
  "Coexpression_Module",
  "Original", "Cleaned",

  # ---- GO enrichment ----
  "GO_ID", "GO_Term", "GO_Aspect", "OR", "N_Superior", "N_Universe",
  "pct", "aspect", "GO", "ONTOLOGY", "Aspect",

  # ---- subpopulation ----
  "Total", "Freq_pct", "Subpopulation",

  # ---- gene model / GFF3 ----
  "Type", "Attributes", "Parent", "Length_Source", "UTR_Type",
  "type", "strand", "gene_id", "gene_function", "transcript_id",
  "start", "end", "width", "CDS_bp", "grp",
  "has_cds_fasta", "has_pep_fasta", "cds_length", "cds_length_ok",
  "translated_length", "pep_length", "translation_matches",
  "cds_start", "cds_end",
  "Exon_CDS_Start", "Exon_Start", "Exon_End",

  # ---- ggplot2 NSE / plot-specific ----
  ".data", ".gene", "x", "y", "xend", "yend", "label",
  "trait_val", "fill", "fill_col", "is_signal", "lwd",
  "i", "j", "vshape", "vlabel",
  "neglogP", "P", "r2", "in_highlight",
  "Expected", "Observed", "density",
  "cap_out", "species",
  
  # ---- plot-layer  ----
  "pa", "pb",                     # hbr_plot_cross_landscape
  "rank_uc", "rank_mid", "shift", "cross_id",  # hbr_plot_cross_ranking
  "parent",                       # hbr_plot_parent_profile
  "panel", "selected",            # hbr_plot_selected_set_coverage
  "top_frac",                     # .hgsi_compare (after fix 6b)

  # ---- haplotypes / workbook ----
  "REF_ref", "concordant", "n_unique_hap",
  "gene", "haplotype", "nt_sequence", "region"
))