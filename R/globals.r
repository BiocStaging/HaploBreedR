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
  "A1", "A2", "BP", "overlap", "dist",

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
  "N_Unique_GWAS_SNPs", ".snps_usable", ".rowid", "poskey",
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

  # ---- breeding decision ----
  "id", "index", "rank", "P_top15", "phenos", "index_SE",
  "g_M0", "g_M2", "rel_M0", "rel_M2", "PEV_M2",
  "sigma2_g", "sigma2_e", "w_hat", "n_features", "n_obs",
  "p1", "p2", "mid_index", "relatedness", "coverage", "complement",
  "n_regions", "score", "trait", "direction", "mode",
  "tested", "superior", "action", "mean_rel", "haplo_cov", "n_crosses",
  "pred", "pred_base", "pred_m", "method", "is_superior",
  "Donor", "Complete_Region",
  "N_Superior_In_Region", "N_Superior_Genes_In_Region",
  "N_Total_Genes_In_Region", "n_traits",

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
  "value", "Total", "Freq_pct", "Subpopulation",

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

  # ---- haplotypes / workbook ----
  "REF_ref", "concordant", "n_unique_hap",
  "gene", "haplotype", "nt_sequence", "region"
))