# ============================================================================
# R/publication-tables.r
# Manuscript-ready multi-sheet Excel workbook from a pipeline run.
# ============================================================================

#' Build a manuscript-ready multi-sheet Excel workbook from a pipeline run
#'
#' Reads `pipeline_result_for_breeding.rds` and the CSVs written by
#' [haplobreed_analysis()] in the same directory, and writes one Excel
#' workbook containing up to 18 sheets covering study summary, gene
#' results, donor tables, breeding decisions, and annotation.
#'
#' @param rds_file Path to `pipeline_result_for_breeding.rds`. Required.
#' @param output_dir Directory holding the pipeline CSVs. `NULL` uses
#'   `dirname(rds_file)`.
#' @param output_file Path to the Excel workbook. `NULL` auto-names to
#'   `<trait>_manuscript_tables.xlsx` in `output_dir`.
#' @param phenotype_data Optional phenotype table.
#' @param id_col ID column in `phenotype_data`.
#' @param profile One of `"manuscript"`, `"supplementary"`, `"full"`.
#' @param include_sheets Character vector of sheet names to include.
#' @param exclude_sheets Character vector of sheet names to exclude.
#' @param overwrite Overwrite `output_file` if it exists.
#' @param verbose Print progress.
#' @return Invisibly, a list with `file`, `sheets`, `skipped`.
#' @examples
#' res <- structure(
#'   list(
#'     traits = list(GY = list(
#'       name      = "GY",
#'       source    = "toy",
#'       trait_col = "GY",
#'       direction = "higher",
#'       info      = list(
#'         source         = "toy",
#'         rds_path       = NULL,
#'         trait_col      = "GY",
#'         direction      = "higher",
#'         mode           = "GWAS",
#'         genes_tested   = 10L,
#'         superior_count = 2L,
#'         n_perm         = 1000L,
#'         ld_decay_bp    = NA_real_))),
#'     parents  = data.frame(rank = 1L, id = "IND01", index = 0.5,
#'                           SE = 0.1, P_top15 = 0.9, haplo_cov = 0.8,
#'                           mean_rel = 0.1, n_crosses = 2L,
#'                           action = "advance"),
#'     crosses  = data.frame(rank = 1L, p1 = "IND01", p2 = "IND02",
#'                           mid_index = 0.4, relatedness = 0.1,
#'                           complement = 0.5, score = 1.4),
#'     evidence = data.frame(trait = "GY", direction = "higher",
#'                           mode = "GWAS", tested = 10L,
#'                           superior = 2L, n_obs = 100L,
#'                           sigma2_g = 0.5, sigma2_e = 0.5,
#'                           h2 = 0.5, w_hat = 0.3, n_features = 5L),
#'     G            = list(G = diag(1), n_total = 5L, n_used = 5L),
#'     fits         = list(GY = list(ids = "IND01",
#'                                   g_M0 = 0.1, g_M2 = 0.1,
#'                                   rel_M2 = 0.8, sigma2_g = 0.5)),
#'     weights      = c(GY = 1),
#'     directions   = c(GY = "higher"),
#'     output_dir   = tempdir(),
#'     plink_prefix = "toy",
#'     maf_min      = 0.05,
#'     parameters   = list(relatedness_max = 0.3,
#'                         cv = list(folds = 3L,
#'                                   repeats = 1L,
#'                                   n_perm = 50L)),
#'     sup_sets     = list(),
#'     index        = data.frame(id = "IND01", index = 0.5,
#'                               SE = 0.1, rank = 1L,
#'                               P_top15 = 0.9),
#'     pipeline_summary = data.frame(
#'       Trait               = "GY",
#'       Direction           = "higher",
#'       Mode                = "GWAS",
#'       Genes_Tested        = 10L,
#'       Superior_Haplotypes = 2L,
#'       LD_Decay_bp         = NA_real_,
#'       Output_Dir          = tempdir(),
#'       stringsAsFactors    = FALSE)),
#'   class = c("hgsi_result", "list"))
#'
#' rds <- file.path(tempdir(), "fake_result.rds")
#' saveRDS(res, rds)
#' suppressWarnings(
#'   hbr_build_manuscript_tables(
#'     rds_file    = rds,
#'     output_file = file.path(tempdir(), "tables.xlsx")))
#' @export
hbr_build_manuscript_tables <- function(
    rds_file       = NULL,
    output_dir     = NULL,
    output_file    = NULL,
    phenotype_data = NULL,
    id_col         = "ID",
    profile        = c("manuscript", "supplementary", "full"),
    include_sheets = NULL,
    exclude_sheets = NULL,
    overwrite      = TRUE,
    verbose        = TRUE) {

  profile <- match.arg(profile)

  if (missing(rds_file) || is.null(rds_file) || !file.exists(rds_file))
    stop("`rds_file` must exist. Expected ",
         "pipeline_result_for_breeding.rds.", call. = FALSE)

  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package `openxlsx` is required. Install with ",
         "install.packages('openxlsx').", call. = FALSE)

  if (is.null(output_dir)) output_dir <- dirname(rds_file)
  if (!dir.exists(output_dir))
    stop("Output directory does not exist: ", output_dir, call. = FALSE)

  rds <- tryCatch(readRDS(rds_file), error = function(e)
    stop("Could not read RDS: ", conditionMessage(e), call. = FALSE))

  trait_col <- as.character(rds$pipeline_summary$Trait[1])
  if (is.na(trait_col) || !nzchar(trait_col)) trait_col <- "trait"

  if (is.null(output_file))
    output_file <- file.path(output_dir,
                             paste0(trait_col, "_manuscript_tables.xlsx"))

  if (file.exists(output_file) && !isTRUE(overwrite))
    stop("Output file exists and `overwrite = FALSE`: ", output_file,
         call. = FALSE)

  .t0 <- Sys.time()
  if (verbose) {
    cat("\n", strrep("=", 68), "\n", sep = "")
    cat("MANUSCRIPT TABLES\n")
    cat(strrep("=", 68), "\n", sep = "")
    cat("  RDS   :", rds_file, "\n")
    cat("  Dir   :", output_dir, "\n")
    cat("  Out   :", output_file, "\n")
    cat("  Trait :", trait_col, "\n")
    cat("  Profile:", profile, "\n\n")
  }

  sheet_specs <- .hbr_mt_sheet_specs()

  if (!is.null(include_sheets))
    sheet_specs <- sheet_specs[names(sheet_specs) %in% include_sheets]
  if (!is.null(exclude_sheets))
    sheet_specs <- sheet_specs[!names(sheet_specs) %in% exclude_sheets]

  sheets  <- list()
  skipped <- list()

  wb <- openxlsx::createWorkbook()

  for (sheet_name in names(sheet_specs)) {
    spec <- sheet_specs[[sheet_name]]
    if (verbose) cat(sprintf("  [%-24s] ", sheet_name))

    built <- tryCatch(
      spec$builder(output_dir, rds, trait_col, phenotype_data, id_col),
      error = function(e) {
        structure(list(msg = conditionMessage(e)), class = "mt_build_fail")
      })

    if (inherits(built, "mt_build_fail")) {
      if (verbose) cat("skipped: ", built$msg, "\n", sep = "")
      skipped[[sheet_name]] <- built$msg
      next
    }
    if (is.null(built) || !nrow(built)) {
      if (verbose) cat("skipped: no data\n")
      skipped[[sheet_name]] <- "no data"
      next
    }

    curated <- tryCatch(
      .hbr_mt_curate_columns(built, sheet_name, profile),
      error = function(e) built)

    if (!nrow(curated)) {
      if (verbose) cat("skipped: no columns survive profile\n")
      skipped[[sheet_name]] <- "no columns survive profile"
      next
    }

    excel_df <- as.data.frame(curated, check.names = FALSE)
    names(excel_df) <- .hbr_mt_excel_colnames(names(excel_df))

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, excel_df,
                        headerStyle = .hbr_mt_header_style())

    n_cols <- ncol(excel_df)
    if (n_cols > 0) {
      col_widths <- vapply(excel_df, function(col) {
        w <- max(nchar(as.character(col)), na.rm = TRUE)
        if (!is.finite(w) || w < 8) w <- 8
        min(w + 2, 60)
      }, numeric(1))
      openxlsx::setColWidths(wb, sheet_name,
                             cols = seq_len(n_cols),
                             widths = col_widths)
    }
    openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)

    sheets[[sheet_name]] <- curated
    if (verbose) cat(sprintf("%d rows x %d cols\n",
                             nrow(curated), ncol(curated)))
  }

  if (length(sheets) == 0) {
    warning("No sheets could be built; workbook not written.",
            call. = FALSE)
    return(invisible(list(file = NULL, sheets = list(),
                          skipped = skipped)))
  }

  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

  if (verbose) {
    cat("\n", strrep("=", 68), "\n", sep = "")
    cat(sprintf("Written: %s\n", output_file))
    cat(sprintf("  Sheets written : %d\n", length(sheets)))
    cat(sprintf("  Sheets skipped : %d\n", length(skipped)))
    if (length(skipped)) {
      cat("\n  Skipped:\n")
      for (nm in names(skipped))
        cat(sprintf("    %-24s %s\n", nm, skipped[[nm]]))
    }
    cat(strrep("=", 68), "\n", sep = "")
    .hbr_log_time(.t0, "manuscript tables", verbose)
  }

  invisible(list(file = output_file, sheets = sheets, skipped = skipped))
}


# ============================================================================
# Sheet specifications
# ============================================================================

#' @noRd
.hbr_mt_sheet_specs <- function() {
  list(
    Study_Summary           = list(builder = .hbr_mt_build_study_summary),
    Superior_Genes          = list(builder = .hbr_mt_build_superior_genes),
    Superior_Haplotypes     = list(builder = .hbr_mt_build_superior_haplotypes),
    All_Tested_Genes        = list(builder = .hbr_mt_build_all_genes),
    Gene_Attrition          = list(builder = .hbr_mt_build_attrition),
    Parents_Ranked          = list(builder = .hbr_mt_build_parents_ranked),
    Cross_Predictions       = list(builder = .hbr_mt_build_cross_predictions),
    Haplotype_Blocks        = list(builder = .hbr_mt_build_blocks_table),
    Selected_Parents_Blocks = list(builder = .hbr_mt_build_selected_blocks),
    Cross_Validation        = list(builder = .hbr_mt_build_cv_summary),
    GO_Enrichment           = list(builder = .hbr_mt_build_go_enrichment),
    Subpop_Haplotype_Freq   = list(builder = .hbr_mt_build_subpop_freq)
  )
}


# ============================================================================
# Sheet builders
# ============================================================================

#' @noRd
.hbr_mt_build_study_summary <- function(output_dir, rds, trait_col, ...) {
  ps_path <- file.path(output_dir, "pipeline_summary.csv")
  if (!file.exists(ps_path)) return(NULL)
  ps <- data.table::fread(ps_path, na.strings = c("", "NA"))
  if (!nrow(ps)) return(NULL)

  long <- data.table::data.table(
    Parameter = names(ps),
    Value     = vapply(ps[1], function(x) as.character(x), character(1)))
  long[]
}

#' @noRd
.hbr_mt_build_superior_genes <- function(output_dir, rds, trait_col, ...) {
  p <- .hbr_mt_read_pub(output_dir)
  if (is.null(p)) return(NULL)
  if (!"Is_Superior" %in% names(p)) return(NULL)
  out <- p[Is_Superior == TRUE]
  if (!nrow(out)) return(NULL)
  out[]
}

#' @noRd
.hbr_mt_build_all_genes <- function(output_dir, rds, trait_col, ...) {
  p <- .hbr_mt_read_pub(output_dir)
  if (is.null(p)) return(NULL)
  p[]
}

#' @noRd
.hbr_mt_build_attrition <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "gene_attrition_log.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_superior_haplotypes <- function(output_dir, rds, trait_col, ...) {
  superior_ids <- .hbr_mt_superior_ids(output_dir, rds)
  if (!length(superior_ids)) return(NULL)

  rows <- list()
  for (g in superior_ids) {
    hm <- .hbr_mt_find_per_gene(output_dir, g, "_hap_map.csv")
    if (is.null(hm)) next
    hap <- data.table::fread(hm, na.strings = c("", "NA"))
    if (!nrow(hap)) next
    hap[, Gene := g]
    rows[[g]] <- hap
  }
  if (!length(rows)) return(NULL)
  data.table::rbindlist(rows, fill = TRUE)
}

#' @noRd
.hbr_mt_build_parents_ranked <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "hgsi", "advanced_parents_ranked.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_cross_predictions <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "hgsi", "advanced_crosses_predicted.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_blocks_table <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "hgsi", "advanced_blocks_table.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_selected_blocks <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "hgsi", "advanced_selected_set_blocks.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_cv_summary <- function(output_dir, rds, trait_col, ...) {
  p <- file.path(output_dir, "hgsi", "advanced_cv_summary.csv")
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, na.strings = c("", "NA"))
}

#' @noRd
.hbr_mt_build_go_enrichment <- function(output_dir, rds, trait_col, ...) {
  rows <- list()
  for (asp in c("BP", "MF", "CC")) {
    p <- file.path(output_dir, paste0("GO_enrichment_", asp, ".csv"))
    if (file.exists(p)) {
      dt <- tryCatch(data.table::fread(p, na.strings = c("", "NA")),
                     error = function(e) NULL)
      if (!is.null(dt) && nrow(dt)) {
        if (!"Aspect" %in% names(dt)) dt[, Aspect := asp]
        rows[[asp]] <- dt
      }
    }
  }
  if (!length(rows)) return(NULL)
  data.table::rbindlist(rows, fill = TRUE)
}

#' @noRd
.hbr_mt_build_subpop_freq <- function(output_dir, rds, trait_col, ...) {
  superior_ids <- .hbr_mt_superior_ids(output_dir, rds)
  if (!length(superior_ids)) return(NULL)

  rows <- list()
  for (g in superior_ids) {
    f <- .hbr_mt_find_per_gene(output_dir, g, "_subpop_freq.csv")
    if (is.null(f)) next
    dt <- tryCatch(data.table::fread(f, na.strings = c("", "NA")),
                   error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) next
    dt[, Gene := g]
    rows[[g]] <- dt
  }
  if (!length(rows)) return(NULL)
  data.table::rbindlist(rows, fill = TRUE)
}


# ============================================================================
# Column profiles
# ============================================================================

#' @noRd
.hbr_mt_col_spec <- function() {
  list(
    Superior_Genes = list(
      manuscript = c(
        "Gene", "Gene_Function", "Chr", "Gene_Start", "Gene_End",
        "LD_Region_ID", "LD_Region_Start", "LD_Region_End",
        "GWAS_SNP", "GWAS_Position",
        "N_Functional_SNPs_Final", "N_Haplotypes", "PIC",
        "Test", "P_value", "Adjusted_P", "Permutation_P",
        "Cohens_d", "Superior_Mean", "Superior_SE",
        "Superior_Haplotype", "Superior_NT_Seq", "Superior_AA_Seq",
        "Superior_Freq_Analysed_Pct", "Improvement_Pct",
        "Functional_Changes", "Cross_Gene_ID"),
      supplementary = c(
        "GFF3_Gene_ID", "Transcript_ID", "N_Transcripts_Considered",
        "CDS_bp", "UTR5_bp", "UTR3_bp",
        "N_NS_CDS_SNPs", "N_5UTR_SNPs", "N_3UTR_SNPs",
        "N_REF_Mismatch", "N_Functional_SNPs_Raw", "N_Functional_SNPs_QC",
        "N_Individuals_Genotyped",
        "Nei_Hd", "Neff",
        "PC_Corrected", "PCs_Used",
        "Significant_Nominal", "Significant_Permutation",
        "Significant_FDR", "Omega_Squared", "CV_Percent",
        "Second_Haplotype", "Descriptive_Letter",
        "Gate_A_Sig", "Gate_B_Perm_Sig", "Gate_C_Effect",
        "Gate_D_Threshold",
        "Diff_Mean", "Diff_CI_Lower", "Diff_CI_Upper",
        "Superior_CI_Lower", "Superior_CI_Upper", "Superior_Freq_Panel_Pct",
        "LD_Decay_bp", "LD_Decay_kb", "Annotation_Source",
        "Phenotype_Type")),

    All_Tested_Genes = list(
      manuscript = c(
        "Gene", "Gene_Function", "Chr", "Gene_Start", "Gene_End",
        "LD_Region_ID",
        "N_Haplotypes", "PIC",
        "Test", "P_value", "Adjusted_P", "Permutation_P",
        "Cohens_d", "Superior_Haplotype", "Is_Superior",
        "Improvement_Pct"),
      supplementary = c(
        "GFF3_Gene_ID", "Transcript_ID", "CDS_bp",
        "N_NS_CDS_SNPs", "N_5UTR_SNPs", "N_3UTR_SNPs",
        "N_Functional_SNPs_Final", "N_Individuals_Genotyped",
        "Nei_Hd", "Neff", "Omega_Squared", "CV_Percent",
        "GWAS_SNP", "GWAS_Position",
        "Significant_Nominal", "Significant_Permutation",
        "Significant_FDR", "Annotation_Source")),

    Superior_Haplotypes = list(
      manuscript = c("Gene", "Haplotype", "Count", "Frequency",
                     "NT_Sequence", "AA_Sequence"),
      supplementary = c()),

    Gene_Attrition = list(
      manuscript = c("Gene", "Exit_Stage", "Exit_Reason",
                     "N_Functional_SNPs", "N_Genotyped_SNPs",
                     "N_After_QC", "N_Individuals", "N_Haplotypes"),
      supplementary = c("Index")),

    Donor_Summary = list(
      manuscript = c(
        "Donor", "N_Superior_Genes", "N_Unique_Superior_Haplotypes",
        "N_LD_Regions", "N_Complete_LD_Regions",
        "Genes", "Superior_Haplotypes", "Amino_Acid_Changes"),
      supplementary = c(
        "N_Functional_Variants_Total",
        "Complete_Regions", "Regions_Carried",
        "N_Total_Superior_In_Regions", "N_Total_Region_Genes",
        "Rank_Effect_Context")),

    Shortlist = list(
      manuscript = c("Donor", "Shortlist_Rank", "Composite_Rank",
                     "Rank_Count", "Rank_Region", "Rank_GEBV",
                     "N_Superior_Genes", "N_Complete_LD_Regions"),
      supplementary = c("Score_Count", "Score_Region", "Score_GEBV",
                        "N_Ranks", "N_LD_Regions",
                        "Genes", "Superior_Haplotypes")),

    Cross_Design = list(
      manuscript = c("Cross", "Trait", "Composite_Rank", "Composite_Score",
                     "F2_Mean", "F2_P95",
                     "N_Linkage_Aware", "Usefulness",
                     "N_F2_for_Usefulness_95", "N_Genes_Segregating",
                     "Segregating_Genes"),
      supplementary = c("Donor_A", "Donor_B",
                        "Norm_Pyramid", "Norm_Usefulness",
                        "Norm_Parent", "Norm_P95",
                        "F2_SD_from_Hap", "Parent_Mean_A",
                        "Parent_Mean_B", "Segregating_Blocks",
                        "P_F2_Full_Carrier", "P_F2_Full_Homozygous",
                        "N_F2_for_Full_Homozygous",
                        "N_F2_for_Carrier_95")),

    Cross_Predictions = list(
      manuscript = c("Cross", "Trait",
                     "F2_Mean", "F2_SD_from_Hap", "F2_P95",
                     "N_Genes_Segregating", "Segregating_Blocks",
                     "P_F2_Full_Homozygous", "N_Linkage_Aware",
                     "Usefulness", "N_F2_for_Usefulness_95",
                     "Segregating_Genes"),
      supplementary = c("Donor_A", "Donor_B",
                        "Parent_Mean_A", "Parent_Mean_B",
                        "P_F2_Full_Carrier", "N_F2_for_Full_Homozygous",
                        "N_F2_for_Carrier_95")),

    MAS_Panel = list(
      manuscript = c("SNP", "Gene", "Chr", "POS", "REF", "ALT",
                     "Region", "AA_Change", "Effect_d",
                     "Assay_Priority", "Amplicon_Size"),
      supplementary = c("Flank_Left", "Flank_Right", "Context_Ref",
                        "F_REF_Primer", "F_ALT_Primer", "R_Primer",
                        "F_Tm", "R_Tm",
                        "Multiplex_Conflict", "In_Repeat_Region",
                        "Paralog_Risk", "Design_Status")),

    MAS_Plate_Map = list(
      manuscript = c("Well", "Plate", "SNP", "Gene", "Chr", "POS",
                     "REF", "ALT", "Effect_d", "Assay_Priority"),
      supplementary = c("F_REF_Primer", "F_ALT_Primer", "R_Primer",
                        "Amplicon_Size",
                        "Multiplex_Conflict", "In_Repeat_Region",
                        "Paralog_Risk")),

    MABB_Plan = list(
      manuscript = c("Year", "Generation", "Cross", "Action",
                     "Population_Size", "Expected_RPG_Pct",
                     "Expected_Target_Loci",
                     "Foreground_Markers", "Background_Markers"),
      supplementary = c("N_Independent", "N_LinkageAware", "Notes")),

    MABB_Pop_Size = list(
      manuscript = c("Per_Locus_Recomb_Fraction",
                     "N_Independent_Loci", "N_LD_Blocks",
                     "N_Independent", "N_LinkageAware",
                     "Success_Prob_Target", "Flanking_cM"),
      supplementary = c("Per_Locus_Double_Recomb_Prob",
                        "P_All_Loci_Independent",
                        "P_All_Blocks_LinkageAware")),

    GO_Enrichment = list(
      manuscript = c("Aspect", "GO_ID", "GO_Term", "GO_Aspect",
                     "N_Superior", "N_Universe", "OR",
                     "P", "Adjusted_P"),
      supplementary = c()),

    Subpop_Haplotype_Freq = list(
      manuscript = c("Gene", "Haplotype", "Subpopulation",
                     "N", "Total", "Freq_pct"),
      supplementary = c())
  )
}

#' @noRd
.hbr_mt_curate_columns <- function(dt, sheet_name, profile) {
  if (identical(profile, "full")) return(dt)

  spec <- .hbr_mt_col_spec()[[sheet_name]]
  if (is.null(spec)) return(dt)

  if (identical(profile, "manuscript")) {
    keep <- intersect(spec$manuscript, names(dt))
  } else {
    keep <- intersect(c(spec$manuscript, spec$supplementary), names(dt))
  }

  if (!length(keep)) return(dt)
  extra <- setdiff(names(dt), c(spec$manuscript, spec$supplementary))
  keep  <- c(keep, intersect(extra, names(dt)))
  dt[, ..keep]
}


# ============================================================================
# Helpers
# ============================================================================

#' @noRd
.hbr_mt_read_pub <- function(output_dir) {
  for (cand in c(
      file.path(output_dir, "publication_table_functional.csv"),
      file.path(output_dir, "superior_genes",
                "publication_table_functional.csv"))) {
    if (file.exists(cand))
      return(tryCatch(data.table::fread(cand, na.strings = c("", "NA")),
                      error = function(e) NULL))
  }
  NULL
}

#' @noRd
.hbr_mt_superior_ids <- function(output_dir, rds) {
  if (!is.null(rds$superior_results) && length(rds$superior_results))
    return(names(rds$superior_results))

  p <- .hbr_mt_read_pub(output_dir)
  if (is.null(p) || !"Is_Superior" %in% names(p)) return(character(0))
  gc <- intersect(c("RAPDB_Gene_ID", "MSU7_Gene_ID", "MSU_Gene_ID",
                    "Gene_ID", "Gene"), names(p))[1]
  if (is.na(gc)) return(character(0))
  unique(as.character(p[Is_Superior == TRUE, get(gc)]))
}

#' @noRd
.hbr_mt_find_per_gene <- function(output_dir, gene_id, suffix) {
  candidates <- c(
    file.path(output_dir, "superior_genes", paste0(gene_id, suffix)),
    file.path(output_dir, paste0(gene_id, suffix)))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) NULL else hit
}



#' @noRd
.hbr_mt_excel_colnames <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

#' @noRd
.hbr_mt_header_style <- function() {
  openxlsx::createStyle(
    textDecoration = "bold",
    halign  = "left",
    valign  = "center",
    fgFill  = "#2C3E50",
    fontColour = "#FFFFFF",
    border  = "bottom",
    borderStyle = "thin",
    borderColour = "#1A252F")
}
