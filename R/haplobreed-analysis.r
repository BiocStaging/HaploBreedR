# ============================================================================
# R/haplobreed-analysis.r
# Single-trait identification of genes carrying superior haplotypes.
#
# Memory model
#   The genotype is held as an `hbr_genotype` whose `$genotypes` slot is a
#   `BEDMatrix::BEDMatrix` streaming handle when the package is installed,
#   or a plain numeric matrix otherwise. Every downstream read goes through
#   `.hbr_pull_genotypes()`, which reads only the requested columns and
#   applies the row subset (`$row_keep`) and strand swap (`$swap_cols`) at
#   read time. Canonicalization stores a `canonical_map` (canonical SNP
#   name -> bed column index) instead of copying the matrix. The RDS
#   packer streams the matrix chunk by chunk.
#
# ============================================================================

#' Single-trait identification of genes carrying superior haplotypes
#'
#' Runs reference and genotype preparation, candidate-gene definition
#' (GWAS LD windows or a gene list), haplotype construction, association
#' testing, the superior-haplotype gates (with BH-FDR), functional
#' annotation of the publication table and, optionally, GO enrichment,
#' for ONE trait. Results are written to `output_dir`. To analyse
#' several traits call this function once per trait with distinct
#' `trait_col` / `output_dir`.
#'
#' @param plink_prefix PLINK file prefix (no extension). `.bed`, `.bim`,
#'   `.fam`, and a `.raw` file are required.
#' @param plink_path Path to the PLINK executable.
#' @param ref_alt Optional `data.frame` with columns `SNP, CHR, POS, REF, ALT`.
#' @param harmonize_strand Correct REF/ALT against the reference genome.
#' @param drop_ambiguous Drop A/T and C/G variants after harmonization.
#' @param gff3,cds,pep,genome_fasta Local reference file paths. Any left
#'   `NULL` is fetched from `reference_source`.
#' @param reference_source One of `"local"`, `"msu7"`, `"rapdb"`, `"ensembl"`.
#' @param crop_species,crop_assembly Ensembl Plants species and assembly.
#' @param ensembl_release Ensembl Plants release. `NULL` queries the
#'   current release.
#' @param reference_cache Cache directory for downloaded reference files.
#' @param pheno_id_col ID column in `phenotype_data`.
#' @param on_fail One of `"stop"` or `"warn"` on compatibility failure.
#' @param qc_report_dir Directory for the input-QC report.
#' @param skip_validation Skip compatibility and ID-overlap checks.
#' @param work_dir Directory for intermediate files.
#' @param gwas_snps GWAS hit table. Mutually exclusive with `gene_list`.
#' @param gene_list Character vector of gene IDs. Direct-gene mode.
#' @param phenotype_data Phenotype table with an ID and the trait column.
#' @param trait_col Trait column name. Exactly one trait per run.
#' @param trait_label Human-readable trait label for figures.
#' @param direction One of `"higher"` or `"lower"`.
#' @param phenotype_type One of `"single_env"`, `"blup"`, `"blue"`.
#' @param pc_cols Principal-component columns. `NULL` skips correction.
#' @param subpop_col Subpopulation column for frequency summaries.
#' @param posthoc Post-hoc tests for three-or-more haplotypes.
#' @param n_perm,n_boot Permutation and bootstrap iterations.
#' @param min_cohens_d Minimum Cohen's d for the superior-haplotype gate.
#' @param alpha Significance threshold.
#' @param trait_threshold Absolute trait value required for the superior
#'   haplotype.
#' @param use_all_transcripts Use all transcripts or only the canonical one.
#' @param regions_to_include Functional regions for haplotype construction.
#' @param max_miss_snp,maf_min,max_het_ind,max_miss_ind QC thresholds.
#' @param min_individuals_qc Minimum individuals after QC.
#' @param min_ind_per_hap Minimum individuals per haplotype group.
#' @param min_hap_freq Minimum haplotype frequency.
#' @param ld_decay_distance Optional LD-decay distance in bp.
#' @param ld_decay_file Path to a saved LD-decay RDS from a previous run.
#' @param output_dir Output directory.
#' @param save_only_superior Write per-gene files only for superior genes.
#' @param keep_intergenic Keep LD regions with no overlapping gene.
#' @param figure_formats Figure output formats.
#' @param aa_full_names Use three-letter amino-acid codes.
#' @param gwas_snp_col,gwas_pos_col,gwas_chr_col Column names in `gwas_snps`.
#' @param require_fdr If `TRUE` (default) a gene is `Is_Superior` only if
#'   it also passes BH-FDR (< `alpha`) across the genes tested in this run.
#' @param plot_style One of `"gene_list"`, `"per_gene"`, `"both"`, `"none"`.
#' @param plot_palette One of [hbr_palette_names()].
#' @param plot_theme One of [hbr_theme_names()].
#' @param plot_table_style `"color_coded"` or `"row_highlight"`.
#' @param plot_panels Panels to draw.
#' @param plot_width,plot_height Figure dimensions in inches.
#' @param plot_combined_pdf Write a combined multi-page PDF.
#' @param go_enrichment Run GO enrichment and annotation.
#' @param go_enrichment_source One of `"msu7"` or `"rapdb"`.
#' @param go_enrichment_min_term,go_enrichment_min_superior Enrichment
#'   size gates.
#' @param go_enrichment_cache Cache directory for annotation downloads.
#' @param save_breeding_rds Save `pipeline_result_for_breeding.rds`.
#' @param rds_compress Compression of the RDS file(s). One of `"fast"`,
#'   `"none"`, `"xz"`.
#' @param rds_reference Where the reference sequences are saved: 
#'   `"embedded"` (inside the breeding RDS), `"separate"` (a sibling
#'   `pipeline_reference.rds`), or `"none"` (drop the sequences).
#' @param save_plot If `FALSE`, no figures are written.
#' @param store_all_results Retain every per-gene result. `TRUE` keeps the
#'   full per-gene object for every tested gene; memory grows linearly with
#'   the number of genes. Prefer the CSVs on disk unless the caller really
#'   needs every gene in memory.
#' @param verbose Print progress. `TRUE` prints one concise status line
#'   per pipeline phase and per gene; `FALSE` prints nothing except
#'   errors and warnings.
#'
#' @return Invisibly, a list with the pipeline result. Notable elements:
#'   \describe{
#'     \item{`reference`}{the `hbr_reference` object used, minus the
#'       genome FASTA if it was not needed for plotting}
#'     \item{`geno_mat`}{the `hbr_genotype` source (streaming or in-memory)}
#'     \item{`geno_ids`}{character vector of the individuals that survived
#'       ID matching}
#'     \item{`superior_results`}{list per superior gene: `ns_result`,
#'       `geno_data`, `hap_mat`, `hap_data`, `test_res`, `superior`,
#'       `sp_freq`, `gene_annotation`}
#'     \item{`summary`, `superior_summary`, `attrition`}{per-gene tables
#'       also written to CSV}
#'     \item{`phenotype_data`}{the cleaned ID + trait table actually used}
#'     \item{`functional_snps_all`}{every functional SNP considered
#'       across all tested genes}
#'     \item{`gate_params`}{the gate parameter list, so downstream tools
#'       can replay the gate without re-parsing the summary}
#'   }
#'
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' res <- try(haplobreed_analysis(
#'   plink_prefix     = file.path(toy, "toy"),
#'   reference_source = "local",
#'   gff3             = file.path(toy, "toy.gff3"),
#'   cds              = file.path(toy, "toy.fa"),
#'   gwas_snps        = read.csv(file.path(toy, "toy_gwas.csv")),
#'   phenotype_data   = read.csv(file.path(toy, "toy_pheno.csv")),
#'   trait_col        = "trait",
#'   direction        = "higher",
#'   output_dir       = tempdir(),
#'   save_plot        = FALSE,
#'   n_perm           = 100L,
#'   n_boot           = 100L,
#'   verbose          = FALSE))
#'
#' @export
haplobreed_analysis <- function(
    plink_prefix        = NULL,
    plink_path          = "plink",
    ref_alt             = NULL,
    harmonize_strand    = TRUE,
    drop_ambiguous      = FALSE,
    gff3                = NULL,
    cds                 = NULL,
    pep                 = NULL,
    genome_fasta        = NULL,
    reference_source    = c("local", "msu7", "rapdb", "ensembl"),
    crop_species        = NULL,
    crop_assembly       = NULL,
    ensembl_release     = NULL,
    reference_cache     = NULL,
    pheno_id_col        = "ID",
    on_fail             = c("stop", "warn"),
    qc_report_dir       = NULL,
    skip_validation     = FALSE,
    work_dir            = tempdir(),

    gwas_snps           = NULL,
    gene_list           = NULL,
    phenotype_data,
    trait_col,
    trait_label         = NULL,
    direction           = c("higher", "lower"),
    phenotype_type      = c("single_env", "blup", "blue"),
    pc_cols             = c("PC1", "PC2", "PC3", "PC4", "PC5"),
    subpop_col          = "Subpopulation",
    posthoc             = c("tukey", "duncan"),
    n_perm              = 1000,
    n_boot              = 1000,
    min_cohens_d        = 0.3,
    alpha               = 0.05,
    trait_threshold     = NULL,
    use_all_transcripts = TRUE,
    regions_to_include  = c("CDS_nonsynonymous", "5UTR", "3UTR"),
    max_miss_snp        = 0.20,
    maf_min             = 0.05,
    max_het_ind         = 0.10,
    max_miss_ind        = 0.10,
    min_individuals_qc  = 10,
    min_ind_per_hap     = 3,
    min_hap_freq        = 0,
    ld_decay_distance   = NULL,
    ld_decay_file       = NULL,
    output_dir          = "HaploBreedR_results",
    save_only_superior  = TRUE,
    keep_intergenic     = TRUE,
    figure_formats      = c("pdf", "tiff", "png"),
    aa_full_names       = TRUE,
    gwas_snp_col        = "SNP",
    gwas_pos_col        = "Position",
    gwas_chr_col        = "Chromosome",
    require_fdr         = TRUE,

    plot_style          = c("gene_list", "per_gene", "both", "none"),
    plot_palette        = "haplobreed",
    plot_theme          = "haplobreed",
    plot_table_style    = c("color_coded", "row_highlight"),
    plot_panels         = c("gene_model", "haplotype_network",
                            "haplotype_table", "trait_boxplot"),
    plot_width          = 11,
    plot_height         = 13,
    plot_combined_pdf   = TRUE,

    go_enrichment              = FALSE,
    go_enrichment_source       = c("msu7", "rapdb"),
    go_enrichment_min_term     = 5L,
    go_enrichment_min_superior = 2L,
    go_enrichment_cache        = NULL,

    save_breeding_rds   = TRUE,
    rds_compress        = c("fast", "none", "xz"),
    rds_reference       = c("embedded", "separate", "none"),
    save_plot           = TRUE,

    store_all_results   = FALSE,
    verbose             = TRUE) {

  # --------------------------------------------------------------------------
  # Argument validation
  # --------------------------------------------------------------------------
  direction         <- match.arg(direction)
  phenotype_type    <- match.arg(phenotype_type)
  on_fail           <- match.arg(on_fail)
  plot_style        <- match.arg(plot_style)
  rds_compress      <- match.arg(rds_compress)
  rds_reference     <- match.arg(rds_reference)
  if (!isTRUE(save_plot)) plot_style <- "none"
  plot_table_style  <- match.arg(plot_table_style)
  reference_source  <- match.arg(reference_source)
  go_enrichment_source <- match.arg(go_enrichment_source)
  .hbr_check_palette_theme(plot_palette, plot_theme)

  if (is.null(plink_prefix))
    stop("`plink_prefix` is required.", call. = FALSE)
  plink_prefix <- sub("\\.(bed|bim|fam|raw)$", "", plink_prefix,
                      ignore.case = TRUE)
  for (ext in c(".bed", ".bim", ".fam")) {
    f <- paste0(plink_prefix, ext)
    if (!file.exists(f))
      stop("PLINK file not found: ", f, call. = FALSE)
  }

  mode_gwas <- !is.null(gwas_snps)
  mode_gene <- !is.null(gene_list)
  if (mode_gwas && mode_gene)
    stop("Supply `gwas_snps` OR `gene_list`, not both.", call. = FALSE)
  if (!mode_gwas && !mode_gene)
    stop("Supply `gwas_snps` or `gene_list`.", call. = FALSE)

  if (is.null(trait_label) || !nzchar(trait_label))
    trait_label <- trait_col

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(work_dir))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  superior_dir <- file.path(output_dir, "superior_genes")
  if (save_only_superior && !dir.exists(superior_dir))
    dir.create(superior_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(qc_report_dir))
    qc_report_dir <- file.path(work_dir, "qc_report")
  if (!dir.exists(qc_report_dir))
    dir.create(qc_report_dir, recursive = TRUE, showWarnings = FALSE)

  .t_pipeline <- Sys.time()

  # --------------------------------------------------------------------------
  # Header banner (one block)
  # --------------------------------------------------------------------------
  if (verbose) {
    cat("\n")
    cat(strrep("=", 75), "\n", sep = "")
    cat("HaploBreedR :: SUPERIOR-HAPLOTYPE ANALYSIS\n")
    cat(strrep("=", 75), "\n", sep = "")
    cat(sprintf("Mode=%s | Trait=%s | Direction=%s | Ref=%s\n",
                if (mode_gwas) "GWAS" else "DirectGene",
                trait_col, direction, reference_source))
    cat(strrep("-", 75), "\n")
  }

  # ==========================================================================
  # PHASE A -- INPUT PREPARATION
  # ==========================================================================

  # ----- A.1 Reference -----
  .t <- Sys.time()
  if (verbose) cat("[A.1] Reference ..................... ")

  ref_paths <- .hbr_resolve_reference_paths(
    gff3            = gff3,
    cds             = cds,
    pep             = pep,
    genome_fasta    = genome_fasta,
    source          = reference_source,
    crop_species    = crop_species,
    crop_assembly   = crop_assembly,
    ensembl_release = ensembl_release,
    cache_dir       = reference_cache,
    verbose         = verbose)

  gff3         <- ref_paths$gff3
  cds          <- ref_paths$cds
  pep          <- ref_paths$pep
  genome_fasta <- ref_paths$genome_fasta

  if (!is.null(ref_paths$cache_dir))
    gff3 <- .hbr_maybe_untar_gff3(gff3, ref_paths$cache_dir)

  if (is.null(gff3) || !file.exists(gff3))
    stop("GFF3 not found. Supply `gff3` or set reference_source.",
         call. = FALSE)
  if (is.null(cds) || !file.exists(cds))
    stop("CDS FASTA not found. Supply `cds` or set reference_source.",
         call. = FALSE)

  annotation_source <- switch(
    reference_source,
    msu7    = "MSU7",
    rapdb   = "RAPDB",
    local   = "CUSTOM",
    ensembl = "CUSTOM")

  reference <- hbr_load_reference(
    genome       = basename(gff3),
    annotation   = annotation_source,
    gff3         = gff3,
    cds          = cds,
    pep          = pep,
    genome_fasta = genome_fasta,
    verbose      = verbose)

  if (is.null(reference$cds) || is.null(reference$utr))
    stop("Reference missing CDS or UTR features.", call. = FALSE)
  if (is.null(reference$cds_seqs) || is.null(reference$pep_seqs))
    stop("Reference missing CDS or peptide sequences.", call. = FALSE)

  if (verbose)
    cat(sprintf("done (%d genes, %d tx, %.1fs)\n",
                nrow(reference$genes), nrow(reference$transcripts),
                as.numeric(difftime(Sys.time(), .t, units = "secs"))))

  # ----- A.2 Genotype (streaming when possible) -----
  .t <- Sys.time()
  if (verbose) cat("[A.2] Genotype ...................... ")

  genotype <- .hbr_open_genotype_plink(plink_prefix, plink_path,
                                        verbose = verbose)

  if (is.null(genotype$report))
    genotype <- hbr_normalize_genotype(genotype, reference,
                                        verbose = verbose)

  if (!skip_validation) {
    compat <- hbr_validate_genotype_reference(genotype, reference,
                                               verbose = verbose)
    if (compat$status == "FAIL" && !isTRUE(harmonize_strand)) {
      msg <- "Compatibility FAILED and harmonize_strand = FALSE."
      if (on_fail == "stop") stop(msg, call. = FALSE)
      else warning(msg, call. = FALSE)
    }
  }

  sub_out   <- .hbr_subset_genotype(
    genotype       = genotype,
    phenotype_data = phenotype_data,
    pheno_id_col   = pheno_id_col,
    verbose        = verbose)
  genotype  <- sub_out$genotype
  id_report <- sub_out$id_report
  geno_obj  <- genotype

  clean_fun <- function(x)
    gsub("\\s+", "_", gsub("-", "_", trimws(as.character(x))))
  matched_raw    <- id_report$matched_ids
  phenotype_data <- data.table::as.data.table(
    data.table::copy(phenotype_data))
  pheno_clean    <- clean_fun(phenotype_data[[pheno_id_col]])
  phenotype_data <- phenotype_data[pheno_clean %in% matched_raw]

  if (verbose)
    cat(sprintf("done (%d ind, %d snps, %.1fs)\n",
                genotype$n_individuals, genotype$n_variants,
                as.numeric(difftime(Sys.time(), .t, units = "secs"))))

  # ----- A.2b Strand harmonization (records swap_cols, no matrix copy) -----
  if (isTRUE(harmonize_strand)) {
    if (verbose) cat("[A.2b] Strand harmonization ......... ")
    if (is.null(reference$genome_seqs)) {
      if (verbose) cat("skipped (no FASTA)\n")
      warning("No genome FASTA loaded; skipping strand harmonization.",
              call. = FALSE)
    } else {
      genotype <- hbr_harmonize_strand(
        genotype       = genotype,
        reference      = reference,
        drop_ambiguous = drop_ambiguous,
        verbose        = verbose)

      data.table::fwrite(
        genotype$strand_report,
        file.path(output_dir, "strand_harmonization_report.csv"))

      # Keep the geno_obj handle consistent with the updated swap/drop
      geno_obj$variants  <- genotype$variants
      geno_obj$swap_cols <- genotype$swap_cols
      geno_obj$drop_cols <- genotype$drop_cols %||% integer(0)

      if (verbose) cat("done\n")
    }
  }

  # ----- A.3 ref_alt -----
  .t <- Sys.time()
  if (verbose) cat("[A.3] REF/ALT table ................. ")

  if (!is.null(ref_alt)) {
    ref_alt <- data.table::as.data.table(data.table::copy(ref_alt))
    ref_alt <- ref_alt[
      !is.na(CHR) & !is.na(POS) & !is.na(REF) & !is.na(ALT) &
        as.character(REF) != "" & as.character(ALT) != "",
      .(SNP = as.character(SNP),
        CHR = as.character(CHR),
        POS = as.numeric(POS),
        REF = toupper(as.character(REF)),
        ALT = toupper(as.character(ALT)))]
    if (!nrow(ref_alt))
      stop("`ref_alt` empty after filtering.", call. = FALSE)

  } else if (isTRUE(harmonize_strand) &&
             !is.null(genotype$strand_report)) {
    ref_alt <- data.table::as.data.table(genotype$variants)[
      !is.na(CHR) & !is.na(POS) & !is.na(REF) & !is.na(ALT) &
        as.character(REF) != "" & as.character(ALT) != "",
      .(SNP = as.character(SNP),
        CHR = as.character(CHR),
        POS = as.numeric(POS),
        REF = toupper(as.character(REF)),
        ALT = toupper(as.character(ALT)))]
    if (!nrow(ref_alt))
      stop("Harmonized genotype has no valid variants.", call. = FALSE)

  } else {
    ref_alt <- .hbr_ref_alt_from_plink_prefix(plink_prefix)
    if (is.null(ref_alt) || !nrow(ref_alt))
      stop("Failed to derive REF/ALT from PLINK prefix.", call. = FALSE)
  }

  geno_snps <- as.character(genotype$variants$SNP)
  overlap   <- sum(ref_alt$SNP %in% geno_snps)
  if (overlap == 0)
    stop("`ref_alt` shares zero SNP IDs with the genotype.", call. = FALSE)

  tryCatch(
    saveRDS(ref_alt, file.path(output_dir, "snps.rds")),
    error = function(e)
      warning("Could not write snps.rds: ",
              conditionMessage(e), call. = FALSE))

  if (verbose) cat(sprintf("done (%d snps)\n", nrow(ref_alt)))

  # ----- A.4 Input-QC report -----
  if (!skip_validation && !is.null(qc_report_dir)) {
    if (verbose) cat("[A.4] Input-QC report ................ ")
    tryCatch({
      hbr_generate_input_qc_report(
        reference      = reference,
        genotype       = geno_obj,
        phenotype_data = phenotype_data,
        pheno_id_col   = pheno_id_col,
        out_dir        = qc_report_dir)
      if (verbose) cat("done\n")
    }, error = function(e) {
      if (verbose) cat("skipped\n")
    })
  }

  # ==========================================================================
  # PHASE B -- SCIENTIFIC ENGINE
  # ==========================================================================

  cds_dt    <- reference$cds
  utr_dt    <- reference$utr
  cds_seqs  <- reference$cds_seqs
  pep_seqs  <- reference$pep_seqs
  all_genes <- data.table::as.data.table(reference$genes)

  # ----- B.1 LD decay -----
  if (verbose) cat("[B.1] LD decay ...................... ")
  if (mode_gwas) {
    if (!is.null(ld_decay_distance)) {
      ld_decay_distance <- as.numeric(ld_decay_distance)[1]
      ld_decay_summary <- data.table::data.table(
        LD_Decay_bp = ld_decay_distance,
        LD_Decay_kb = ld_decay_distance / 1000,
        Method      = "User_provided")
      if (verbose) cat(sprintf("user %.1f kb\n", ld_decay_distance / 1000))

    } else if (!is.null(ld_decay_file)) {
      if (!file.exists(ld_decay_file))
        stop("LD-decay file not found: ", ld_decay_file, call. = FALSE)
      ldr <- readRDS(ld_decay_file)
      ld_decay_distance <- as.numeric(ldr$decay_distance_bp)
      ld_decay_summary <- data.table::data.table(
        LD_Decay_bp = ld_decay_distance,
        LD_Decay_kb = ld_decay_distance / 1000,
        Method      = ldr$method %||% NA_character_)
      if (verbose) cat(sprintf("file %.1f kb\n", ld_decay_distance / 1000))

    } else {
      hap_prefix <- file.path(work_dir, "hbr_geno_hap")
      hbr_write_plink(genotype, hap_prefix, plink_path, verbose = FALSE)

      ldr <- hbr_calc_ld_decay(
        plink_prefix    = hap_prefix,
        plink           = plink_path,
        max_distance    = 500000,
        bin_size        = 10000,
        decay_threshold = 0.20,
        use_median      = FALSE,
        verbose         = verbose)
      ld_decay_distance <- ldr$decay_distance_bp
      ld_decay_summary <- data.table::data.table(
        LD_Decay_bp = ld_decay_distance,
        LD_Decay_kb = ld_decay_distance / 1000,
        Method      = ldr$method)
      saveRDS(ldr, file.path(output_dir, "ld_decay_result.rds"))
      if (verbose)
        cat(sprintf("computed %.1f kb (%s)\n",
                    ld_decay_distance / 1000, ldr$method))
    }
    data.table::fwrite(ld_decay_summary,
                       file.path(output_dir, "ld_decay_summary.csv"))
  } else {
    ld_decay_distance <- NA_real_
    ld_decay_summary  <- data.table::data.table(
      LD_Decay_bp = NA_real_, LD_Decay_kb = NA_real_,
      Method      = "Not_calculated (direct gene mode)")
    if (verbose) cat("n/a (direct gene mode)\n")
  }

  # ----- B.2 Candidate gene table -----
  if (verbose) cat("[B.2] Candidate gene table ........... ")

  if (mode_gwas) {
    gwas_snps_dt <- hbr_validate_gwas_snps(
      gwas_snps = gwas_snps, snp_col = gwas_snp_col,
      chr_col = gwas_chr_col, pos_col = gwas_pos_col,
      conflict_action = "first", verbose = FALSE)
    gwas_snps <- data.table::as.data.table(
      data.table::copy(gwas_snps_dt))

    gwas_snps[, SNP        := as.character(get(gwas_snp_col))]
    gwas_snps[, Chromosome := as.character(get(gwas_chr_col))]
    gwas_snps[, Position   := as.numeric(get(gwas_pos_col))]

    gwas_snps_original <- data.table::copy(gwas_snps)
    gwas_snps <- gwas_snps[!is.na(SNP) & SNP != "" &
                             !is.na(Chromosome) & Chromosome != "" &
                             !is.na(Position)]

    regions <- hbr_create_ld_windows(
      gwas_snps         = gwas_snps,
      ld_decay_distance = ld_decay_distance,
      snp_col           = "SNP",
      chr_col           = "Chromosome",
      pos_col           = "Position",
      merge_overlapping = TRUE,
      verbose           = FALSE)

    regions[, N_Unique_GWAS_SNPs := lengths(GWAS_SNPs)]
    regions[, LD_Decay_bp := ld_decay_distance]
    regions[, LD_Decay_kb := ld_decay_distance / 1000]

    region_snp_map <- regions[
      , .(GWAS_SNP = unlist(GWAS_SNPs, use.names = FALSE)),
      by = .(LD_Region_ID, Chromosome, Region_Start, Region_End)]
    region_snp_map[, GWAS_SNP := as.character(GWAS_SNP)]
    # GWAS_SNP_to_LD_region_map.csv not written: identical region<->SNP
    # pairs are already in candidate_genes_from_GWAS_LD_regions.csv.

    annotation <- hbr_annotate_regions_with_genes(
      regions = regions, genes = all_genes,
      chr_col = "Chromosome", start_col = "Region_Start",
      end_col = "Region_End", keep_intergenic = keep_intergenic,
      verbose = FALSE)

    annotation <- data.table::as.data.table(annotation)
    if (!"Gene_ID" %in% names(annotation)) {
      if ("gene_id" %in% names(annotation))
        annotation[, Gene_ID := as.character(gene_id)]
      else stop("Gene annotation lacks Gene_ID.", call. = FALSE)
    }
    annotation[, Gene_ID := .hbr_clean_gene_id(Gene_ID)]
    annotation <- annotation[!is.na(Gene_ID)]

    if (!"GWAS_SNPs" %in% names(annotation))
      annotation <- merge(annotation,
                          regions[, .(LD_Region_ID, Chromosome,
                                      Region_Start, Region_End, GWAS_SNPs)],
                          by = c("LD_Region_ID", "Chromosome",
                                 "Region_Start", "Region_End"),
                          all.x = TRUE)

    for (nm in c("Gene_Chromosome", "Gene_Start", "Gene_End"))
      if (!nm %in% names(annotation))
        annotation[, (nm) := if (nm == "Gene_Chromosome")
          NA_character_ else NA_real_]

    annotation[, .snps_usable := vapply(
      GWAS_SNPs, function(s) length(s) > 0 && !anyNA(s), logical(1))]
    usable <- annotation[.snps_usable == TRUE]
    annotation[, .snps_usable := NULL]
    usable[, .snps_usable := NULL]
    usable[, .rowid := .I]

    candidate_genes <- usable[
      , .(GWAS_SNP          = as.character(unlist(GWAS_SNPs,
                                                  use.names = FALSE)),
          Gene_ID           = as.character(Gene_ID),
          LD_Region_ID      = as.character(LD_Region_ID),
          Chromosome        = as.character(Chromosome),
          Region_Start      = as.numeric(Region_Start),
          Region_End        = as.numeric(Region_End),
          Gene_Chromosome   = as.character(Gene_Chromosome),
          Gene_Start        = as.numeric(Gene_Start),
          Gene_End          = as.numeric(Gene_End),
          Gene_Relationship = as.character(Gene_Relationship),
          Annotation_Status = as.character(Annotation_Status)),
      by = .rowid]
    candidate_genes[, .rowid := NULL]
    candidate_genes <- candidate_genes[!is.na(GWAS_SNP) & GWAS_SNP != ""]
    candidate_genes[, Candidate_Gene_ID := Gene_ID]
    candidate_genes[, Gene := Gene_ID]
    candidate_genes[gwas_snps_original,
                    on = c(GWAS_SNP = "SNP"),
                    GWAS_SNP_Position := i.Position]
    candidate_genes <- unique(candidate_genes)
    data.table::fwrite(candidate_genes,
                       file.path(output_dir,
                                 "candidate_genes_from_GWAS_LD_regions.csv"))

    analysis_gene_map <- unique(candidate_genes[
      !is.na(Candidate_Gene_ID) & Candidate_Gene_ID != "" &
        !is.na(GWAS_SNP) & GWAS_SNP != "",
      .(Candidate_Gene_ID = .hbr_clean_gene_id(Candidate_Gene_ID),
        GWAS_SNP          = as.character(GWAS_SNP),
        GWAS_SNP_Position,
        LD_Region_ID      = as.character(LD_Region_ID),
        Chromosome        = as.character(Chromosome))])

  } else {
    gene_list <- unique(as.character(gene_list))
    ref_ids <- as.character(all_genes$gene_id)
    idx <- match(.hbr_clean_gene_id(gene_list), .hbr_clean_gene_id(ref_ids))
    idx <- idx[!is.na(idx)]
    if (!length(idx))
      stop("None of `gene_list` found in reference.", call. = FALSE)

    matched <- all_genes[idx]
    candidate_genes <- data.table::data.table(
      Gene_ID           = as.character(matched$gene_id),
      Candidate_Gene_ID = as.character(matched$gene_id),
      Gene              = as.character(matched$gene_id),
      Chromosome        = as.character(matched$Chr),
      Region_Start      = as.numeric(matched$Start),
      Region_End        = as.numeric(matched$End),
      Gene_Chromosome   = as.character(matched$Chr),
      Gene_Start        = as.numeric(matched$Start),
      Gene_End          = as.numeric(matched$End),
      Gene_Relationship = "GENE_SELF",
      Annotation_Status = "Direct gene lookup",
      LD_Region_ID      = NA_character_,
      GWAS_SNP          = NA_character_,
      GWAS_SNP_Position = NA_real_)

    analysis_gene_map <- data.table::data.table(
      Candidate_Gene_ID = as.character(matched$gene_id),
      GWAS_SNP          = NA_character_,
      GWAS_SNP_Position = NA_real_,
      LD_Region_ID      = NA_character_,
      Chromosome        = as.character(matched$Chr))

    data.table::fwrite(candidate_genes,
                       file.path(output_dir,
                                 "candidate_genes_direct_list.csv"))
  }

  # analysis_gene_snp_map.csv no longer written: its long-form content
  # is already in candidate_genes_from_GWAS_LD_regions.csv / direct list.

  analysis_genes <- unique(analysis_gene_map$Candidate_Gene_ID)
  analysis_genes <- analysis_genes[!is.na(analysis_genes) &
                                     analysis_genes != ""]
  n_genes <- length(analysis_genes)
  if (!n_genes) stop("No valid candidate genes.", call. = FALSE)

  if (verbose) {
    if (mode_gwas)
      cat(sprintf("done (%d genes, %d LD regions)\n",
                  data.table::uniqueN(candidate_genes$Candidate_Gene_ID),
                  nrow(regions)))
    else
      cat(sprintf("done (%d direct genes)\n", nrow(candidate_genes)))
  }

  # ----- B.3 Canonicalize genotype (lazy) -----
  if (verbose) cat("[B.3] Genotype canonicalization ...... ")

  var_dt <- data.table::as.data.table(genotype$variants)
  var_dt[, CHR := as.character(CHR)]
  var_dt[, POS := as.numeric(POS)]
  var_dt[, poskey := paste0(CHR, ":", POS)]

  ref_alt[, CHR := as.character(CHR)]
  ref_alt[, POS := as.numeric(POS)]
  ref_alt[, poskey := paste0(CHR, ":", POS)]

  ref_unique <- ref_alt[!duplicated(poskey)]
  pos2canon  <- stats::setNames(ref_unique$SNP, ref_unique$poskey)

  # Canonical name per bed column (NA when no match)
  canon <- unname(pos2canon[var_dt$poskey])

  # Honour drop_ambiguous that was recorded during harmonization
  drop_cols <- genotype$drop_cols %||% integer(0)
  if (length(drop_cols)) canon[drop_cols] <- NA_character_

  keep  <- !is.na(canon) & canon != ""

  # First occurrence wins for duplicated canonical names
  canon_kept <- canon[keep]
  bed_idx    <- which(keep)
  dup        <- duplicated(canon_kept)
  canon_kept <- canon_kept[!dup]
  bed_idx    <- bed_idx[!dup]

  if (!length(bed_idx))
    stop("No genotype SNPs matched canonical SNP IDs.", call. = FALSE)

  # Store canonical_map = canonical name -> bed column index
  genotype$canonical_map       <- stats::setNames(bed_idx, canon_kept)
  genotype$canonical_snp_names <- canon_kept

  # Re-align row_keep so rows correspond to the trait table
  trait_data <- data.table::as.data.table(
    phenotype_data[, c(pheno_id_col, trait_col), with = FALSE])
  data.table::setnames(trait_data, pheno_id_col, "ID")
  trait_data[, ID := as.character(ID)]

  row_keep <- genotype$row_keep
  row_ids  <- genotype$row_names[row_keep]
  common   <- intersect(row_ids, trait_data$ID)

  keep_row    <- row_ids %in% common
  genotype$row_keep      <- row_keep[keep_row]
  genotype$individuals   <- row_ids[keep_row]
  genotype$n_individuals <- length(genotype$row_keep)

  trait_data <- trait_data[match(genotype$individuals, trait_data$ID), ]

  if (verbose)
    cat(sprintf("done (%d snps x %d ind, lazy)\n",
                length(genotype$canonical_snp_names),
                genotype$n_individuals))

  # Free the whole-genome FASTA when the plotting layer won't need it.
  if (!plot_style %in% c("gene_list", "both"))
    reference$genome_seqs <- NULL

  # ----- B.4 Per-gene loop -----
  outcome_log <- vector("list", as.integer(max(n_genes, 1L)))

  record_outcome <- function(i, gene_id, stage, reason = NA_character_,
                             n_functional = NA_integer_,
                             n_genotyped = NA_integer_,
                             n_after_qc = NA_integer_,
                             n_individuals = NA_integer_,
                             n_haplotypes = NA_integer_) {
    list(Index = i, Gene = gene_id,
         Exit_Stage = stage, Exit_Reason = reason,
         N_Functional_SNPs = n_functional,
         N_Genotyped_SNPs  = n_genotyped,
         N_After_QC        = n_after_qc,
         N_Individuals     = n_individuals,
         N_Haplotypes      = n_haplotypes)
  }

  summary_slots <- vector("list", n_genes)
  summary_i     <- 0L
  all_results           <- list()
  superior_gene_results <- list()
  functional_snps_log   <- vector("list", n_genes)

  if (verbose) cat(sprintf("[B.4] Processing %d genes ...\n", n_genes))

  .gline <- function(i, id, status)
    if (isTRUE(verbose))
      cat(sprintf("[%d/%d] %-22s ... %s\n", i, n_genes, id, status))

  for (i in seq_along(analysis_genes)) {
    gene_id <- analysis_genes[i]

    gene_rows <- candidate_genes[Candidate_Gene_ID == gene_id]
    if (!nrow(gene_rows)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "gene_lookup", "Not found in candidate table")
      .gline(i, gene_id, "not in candidate table")
      next
    }

    r_start <- suppressWarnings(min(gene_rows$Region_Start, na.rm = TRUE))
    r_end   <- suppressWarnings(max(gene_rows$Region_End,   na.rm = TRUE))
    if (!is.finite(r_start)) r_start <- NA_real_
    if (!is.finite(r_end))   r_end   <- NA_real_

    cross_id <- if ("RAPDB_Gene_ID" %in% names(gene_rows))
      as.character(gene_rows$RAPDB_Gene_ID[1]) else
        if ("MSU_Gene_ID" %in% names(gene_rows))
          as.character(gene_rows$MSU_Gene_ID[1]) else NA_character_

    gene_row_for_summary <- data.table::data.table(
      RAPDB_Gene_ID    = as.character(gene_rows$Gene_ID[1]),
      Chromosomes      = as.character(gene_rows$Chromosome[1]),
      Gene_Start       = as.numeric(gene_rows$Gene_Start[1]),
      Gene_End         = as.numeric(gene_rows$Gene_End[1]),
      LD_Region_Start  = r_start,
      LD_Region_End    = r_end,
      LD_Region_Length = if (!is.na(r_start) && !is.na(r_end))
        r_end - r_start + 1 else NA_real_,
      Cross_Gene_ID    = cross_id)

    ns_res <- tryCatch(
      hbr_get_functional_variants(
        rapdb_id = gene_id, cds_dt = cds_dt, utr_dt = utr_dt,
        cds_seqs = cds_seqs, pep_seqs = pep_seqs,
        all_genes = all_genes, ref_alt = ref_alt,
        use_all_transcripts = use_all_transcripts,
        regions_to_include  = regions_to_include,
        verbose = FALSE),
      error = function(e) NULL)

    if (is.null(ns_res)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "functional_variant_identification", "failed")
      .gline(i, gene_id, "no transcript / no CDS overlap")
      next
    }

    n_functional <- nrow(ns_res$functional_snps)
    if (!n_functional) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "no_functional_variants",
        "No requested functional variants", n_functional = 0L)
      .gline(i, gene_id, "no functional SNPs")
      next
    }

    fs <- data.table::as.data.table(ns_res$functional_snps)
    fs[, gene := gene_id]
    fs[, Chr  := as.character(Chr)]
    fs[, SNP  := as.character(SNP)]
    fs[, POS  := as.numeric(POS)]
    fs[, REF  := as.character(REF)]
    fs[, ALT  := as.character(ALT)]
    fs[, Region := as.character(Region)]
    if (!"AA_Change" %in% names(fs))          fs[, AA_Change := NA_character_]
    if (!"Variant_Label" %in% names(fs))      fs[, Variant_Label := NA_character_]
    if (!"Codon_Num" %in% names(fs))          fs[, Codon_Num := NA_real_]
    if (!"Is_Nonsynonymous" %in% names(fs))   fs[, Is_Nonsynonymous := NA]
    if (!"Strand" %in% names(fs))             fs[, Strand := NA_character_]
    if (!"Genomic_Order" %in% names(fs))
      fs[, Genomic_Order := ifelse(Strand %in% "-", -POS, POS)]

    keep_cols <- intersect(
      c("gene", "SNP", "POS", "REF", "ALT", "Chr", "Strand",
        "Region", "AA_Change", "Variant_Label",
        "Codon_Num", "Is_Nonsynonymous", "Genomic_Order"),
      names(fs))
    functional_snps_log[[i]] <- fs[, ..keep_cols]

    # Streaming extract: only this gene's columns are read from disk
    gd <- tryCatch(
      hbr_extract_genotypes(
        ns_res$functional_snps, genotype,
        max_miss_snp, maf_min, max_het_ind, max_miss_ind,
        min_individuals_qc, gene_label = NULL, verbose = FALSE),
      error = function(e) NULL)

    if (is.null(gd)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "genotype_qc", "No SNPs/individuals passed QC",
        n_functional = n_functional)
      .gline(i, gene_id, "QC failed")
      next
    }

    hap_mat <- tryCatch(hbr_build_haplotype_matrix(gd),
                        error = function(e) NULL)
    if (is.null(hap_mat)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "haplotype_matrix", "failed",
        n_functional, gd$n_snps, gd$n_snps, gd$n_individuals)
      .gline(i, gene_id, "haplotype matrix failed")
      next
    }

    hap_data <- tryCatch(
      hbr_define_haplotypes(hap_mat, min_ind_per_hap, min_hap_freq),
      error = function(e) NULL)
    if (is.null(hap_data)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "haplotype_definition", "too few groups",
        n_functional, gd$n_snps, gd$n_snps, gd$n_individuals)
      .gline(i, gene_id, "haplotypes < 2")
      next
    }

    sp_freq <- if (subpop_col %in% names(phenotype_data))
      tryCatch(hbr_subpop_freq(hap_data, phenotype_data, subpop_col),
               error = function(e) NULL) else NULL

    test_res <- tryCatch(
      hbr_test_association(
        hap_data          = hap_data,
        pheno_data        = phenotype_data,
        trait_col         = trait_col,
        pc_cols           = pc_cols,
        posthoc           = posthoc,
        n_perm            = n_perm,
        n_boot            = n_boot),
      error = function(e) NULL)
    if (is.null(test_res)) {
      outcome_log[[i]] <- record_outcome(
        i, gene_id, "association_test", "failed",
        n_functional, gd$n_snps, gd$n_snps,
        hap_data$n_individuals, hap_data$n_haplotypes)
      .gline(i, gene_id, "association test failed")
      next
    }

    sup <- tryCatch(
      hbr_identify_superior(
        test_res, direction, trait_col, min_cohens_d, alpha,
        trait_threshold, n_panel = hap_data$n_individuals),
      error = function(e) NULL)
    if (is.null(sup)) {
      sup <- list(
        haplotype = NA_character_, nt_sequence = NA_character_,
        aa_sequence = NA_character_, mean = NA_real_, se = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_, n = NA_integer_,
        freq_analysed_pct = NA_real_, freq_panel_pct = NA_real_,
        letter = NA_character_, second_haplotype = NA_character_,
        is_superior = FALSE,
        gate_a_sig = NA, gate_b_perm_sig = NA,
        gate_c_effect = NA, gate_d_threshold = NA,
        superior_label = "No superior", improvement_pct = NA_real_,
        cohens_d = NA_real_, diff_mean = NA_real_,
        diff_ci_lower = NA_real_, diff_ci_upper = NA_real_)
    }

    is_superior <- isTRUE(sup$is_superior)

    .gline(i, gene_id, if (is_superior)
      sprintf("SUPERIOR  %s  d=%.2f  P=%s",
              sup$haplotype,
              if (is.finite(sup$cohens_d)) sup$cohens_d else NA_real_,
              hbr_format_p(test_res$p_value))
      else
      sprintf("tested    P=%s  %d hap%s",
              hbr_format_p(test_res$p_value),
              hap_data$n_haplotypes,
              if (hap_data$n_haplotypes == 1L) "" else "s"))

    if (is_superior) {
      superior_gene_results[[gene_id]] <- list(
        ns_result = ns_res, geno_data = gd, hap_mat = hap_mat,
        hap_data = hap_data, test_res = test_res, superior = sup,
        sp_freq = sp_freq, gene_annotation = gene_rows)
    }

    gene_output_dir <- if (save_only_superior && is_superior)
      superior_dir else output_dir

    if (is_superior) {
      ind_out <- merge(data.table::copy(hap_data$haplotypes),
                       trait_data, by = "ID", all.x = TRUE)
      data.table::fwrite(ind_out,
                         file.path(gene_output_dir,
                                   paste0(gene_id, "_individuals.csv")))
      data.table::fwrite(hap_data$haplotype_map,
                         file.path(gene_output_dir,
                                   paste0(gene_id, "_hap_map.csv")))
      data.table::fwrite(test_res$groups,
                         file.path(gene_output_dir,
                                   paste0(gene_id, "_groups.csv")))
      data.table::fwrite(ns_res$functional_snps,
                         file.path(gene_output_dir,
                                   paste0(gene_id, "_functional_snps.csv")))
      # Per-gene _cds_snps.csv and _utr_snps.csv are not written:
      # _functional_snps.csv already carries what the plot layer needs,
      # and the additional variants are always excluded from haplotypes.
      if (!is.null(sp_freq))
        data.table::fwrite(sp_freq$long,
                           file.path(gene_output_dir,
                                     paste0(gene_id, "_subpop_freq.csv")))

      if (plot_style %in% c("per_gene", "both")) {
        tryCatch(
          hbr_plot_haplotype(
            gene_id, hap_data, test_res, trait_col,
            gene_output_dir, sup, hap_mat$functional_snps_used,
            adjusted_p = NA_real_,
            aa_full_names = aa_full_names,
            formats = figure_formats),
          error = function(e) NULL)
      }
    }

    sr <- tryCatch(
      .hbr_build_summary_row(
        gene_row = gene_row_for_summary, ns_res = ns_res,
        gd = gd, hap_mat = hap_mat, hap_data = hap_data,
        test_res = test_res, sup = sup,
        phenotype_type = phenotype_type,
        avail_pcs = intersect(pc_cols, names(phenotype_data)),
        all_genes = all_genes, all_annotation = NULL,
        gwas_snp_col = gwas_snp_col, gwas_pos_col = gwas_pos_col),
      error = function(e) NULL)

    if (!is.null(sr)) {
      sr <- data.table::as.data.table(sr)
      if (!"Gene" %in% names(sr)) sr[, Gene := gene_id]
      if (!"Is_Superior" %in% names(sr))
        sr[, Is_Superior := is_superior]
      sr[, LD_Decay_bp := ld_decay_distance]
      sr[, LD_Decay_kb := ld_decay_distance / 1000]
      sr[, Annotation_Source := annotation_source]

      gene_map_rows <- analysis_gene_map[Candidate_Gene_ID == gene_id]
      ld_ids <- sort(unique(as.character(gene_map_rows$LD_Region_ID)))
      ld_ids <- ld_ids[!is.na(ld_ids)]
      sr[, LD_Region_ID := if (length(ld_ids))
        paste(ld_ids, collapse = ";") else NA_character_]

      snp_tbl <- gene_map_rows[
        !is.na(GWAS_SNP) & GWAS_SNP != "",
        .(GWAS_SNP = as.character(GWAS_SNP), GWAS_SNP_Position)]
      snp_tbl <- snp_tbl[!duplicated(GWAS_SNP)]
      data.table::setorder(snp_tbl, GWAS_SNP)

      sr[, GWAS_SNP := if (nrow(snp_tbl))
        paste(snp_tbl$GWAS_SNP, collapse = ";") else NA_character_]
      sr[, GWAS_Position := if (nrow(snp_tbl))
        paste(snp_tbl$GWAS_SNP_Position, collapse = ";")
        else NA_character_]
      sr[, N_GWAS_SNPs := nrow(snp_tbl)]

      summary_i <- summary_i + 1L
      summary_slots[[summary_i]] <- sr
    }

    if (store_all_results)
      all_results[[gene_id]] <- list(
        ns_result = ns_res, geno_data = gd, hap_mat = hap_mat,
        hap_data = hap_data, test_res = test_res, superior = sup,
        sp_freq = sp_freq, gene_annotation = gene_rows,
        is_superior = is_superior)

    outcome_log[[i]] <- record_outcome(
      i, gene_id,
      if (is_superior) "superior" else "tested_not_superior",
      if (is_superior) "Passed gates" else "Not superior",
      n_functional = n_functional, n_genotyped = gd$n_snps,
      n_after_qc = gd$n_snps,
      n_individuals = hap_data$n_individuals,
      n_haplotypes = hap_data$n_haplotypes)
  }

  functional_snps_all <- data.table::rbindlist(
    functional_snps_log[!vapply(functional_snps_log, is.null, logical(1))],
    fill = TRUE)
  if (nrow(functional_snps_all)) {
    fs_path <- file.path(output_dir,
                         "functional_snps_all_tested_genes.csv")
    data.table::fwrite(functional_snps_all, fs_path)
  }

  # ----- B.5 Summary + FDR + publication table + figures -----
  summary_table <- if (summary_i > 0L)
    data.table::rbindlist(summary_slots[seq_len(summary_i)], fill = TRUE)
  else data.table::data.table()

  superior_summary <- NULL
  summary_table_clean <- NULL

  if (nrow(summary_table) && "P_value" %in% names(summary_table)) {
    prefix  <- if (annotation_source == "CUSTOM") "" else
      paste0(annotation_source, "_")
    gid_col <- paste0(prefix, "Gene_ID")
    gfn_col <- paste0(prefix, "Gene_Function")

    if ("Gene" %in% names(summary_table))
      data.table::setnames(summary_table, "Gene", gid_col)
    if ("Gene_Function" %in% names(summary_table))
      data.table::setnames(summary_table, "Gene_Function", gfn_col)

    summary_table[, Adjusted_P :=
                    hbr_apply_fdr_correction(P_value, method = "BH")]
    summary_table[, Significant_FDR :=
                    !is.na(Adjusted_P) & Adjusted_P < alpha]
    summary_table[, Significant_Nominal :=
                    !is.na(P_value) & P_value < alpha]
    summary_table[, Significant_Permutation :=
                    !is.na(Permutation_P) & Permutation_P < alpha]

    if ("Is_Superior" %in% names(summary_table)) {
      summary_table[, Is_Superior_Nominal := as.logical(Is_Superior)]
      if (isTRUE(require_fdr)) {
        summary_table[, Is_Superior := Is_Superior_Nominal %in% TRUE &
                        Significant_FDR %in% TRUE]
        keep_sup <- as.character(summary_table[Is_Superior == TRUE][[gid_col]])
        superior_gene_results <- superior_gene_results[
          names(superior_gene_results) %in% keep_sup]
      }
    }

    keep_always <- unique(c(
      gid_col, gfn_col, "Is_Superior",
      "P_value", "Adjusted_P", "Permutation_P",
      "Superior_Haplotype", "Annotation_Source",
      "Gene", "Gene_Function"))

    summary_table_clean <- .hbr_drop_blank_cols(summary_table,
                                                 keep_always = keep_always)

    data.table::fwrite(summary_table_clean,
                       file.path(output_dir,
                                 "haplotype_summary_functional.csv"))

    pub_cols <- c(
      gid_col, gfn_col, "Cross_Gene_ID",
      "Cross_GO_IDs", "Cross_GO_BP", "Cross_GO_MF", "Cross_GO_CC",
      "GWAS_SNP", "GWAS_Position",
      "LD_Region_ID", "LD_Region_Start", "LD_Region_End",
      "LD_Region_Length",
      "GFF3_Gene_ID", "Transcript_ID", "N_Transcripts_Considered", "Chr",
      "CDS_bp", "UTR5_bp", "UTR3_bp",
      "N_NS_CDS_SNPs", "N_5UTR_SNPs", "N_3UTR_SNPs", "N_REF_Mismatch",
      "N_Functional_SNPs_Final", "Functional_Changes",
      "N_Haplotypes", "Nei_Hd", "PIC",
      "Test", "PC_Corrected", "PCs_Used",
      "P_value", "Permutation_P", "Adjusted_P",
      "Significant_Nominal", "Significant_Permutation", "Significant_FDR",
      "Omega_Squared", "CV_Percent",
      "Superior_Haplotype", "Superior_NT_Seq", "Superior_AA_Seq",
      "Superior_Mean", "Superior_SE",
      "Superior_CI_Lower", "Superior_CI_Upper",
      "Superior_Freq_Analysed_Pct", "Superior_Freq_Panel_Pct",
      "Second_Haplotype", "Descriptive_Letter",
      "Is_Superior", "Is_Superior_Nominal", "Gate_A_Sig", "Gate_B_Perm_Sig",
      "Gate_C_Effect", "Gate_D_Threshold",
      "Improvement_Pct", "Cohens_d",
      "Diff_Mean", "Diff_CI_Lower", "Diff_CI_Upper",
      "LD_Decay_bp", "LD_Decay_kb", "Annotation_Source")
    extra <- grep("^Candidate_", names(summary_table_clean), value = TRUE)
    pub_cols <- unique(c(pub_cols, extra))
    pub_cols <- pub_cols[pub_cols %in% names(summary_table_clean)]

    pub_table <- summary_table_clean[, ..pub_cols]
    pub_table <- .hbr_drop_blank_cols(pub_table, keep_always = keep_always)

    data.table::fwrite(pub_table,
                       file.path(output_dir,
                                 "publication_table_functional.csv"))

    superior_summary <- summary_table_clean[Is_Superior == TRUE]
    if (nrow(superior_summary)) {
      data.table::fwrite(superior_summary,
                         file.path(output_dir,
                                   "superior_genes_summary.csv"))
      if (save_only_superior) {
        data.table::fwrite(superior_summary,
                           file.path(superior_dir,
                                     "superior_genes_summary.csv"))
        data.table::fwrite(pub_table,
                           file.path(superior_dir,
                                     "publication_table_functional.csv"))
      }

      # superior_genes_list.csv is not written: it was a strict subset of
      # superior_genes_summary.csv. It is still constructed in memory for
      # the console top-5 print below.
      superior_genes_list <- data.table::data.table(
        Gene                       = superior_summary[[gid_col]],
        Superior_Haplotype         = superior_summary$Superior_Haplotype,
        P_value                    = superior_summary$P_value,
        Adjusted_P                 = superior_summary$Adjusted_P,
        Improvement_Pct            = superior_summary$Improvement_Pct,
        Cohens_d                   = superior_summary$Cohens_d,
        Superior_Mean              = superior_summary$Superior_Mean,
        Superior_Freq_Analysed_Pct =
          superior_summary$Superior_Freq_Analysed_Pct)

      if (verbose) {
        cat(sprintf("      %d superior gene(s):\n", nrow(superior_summary)))
        top_hits <- utils::head(superior_genes_list, 5L)
        for (k in seq_len(nrow(top_hits)))
          cat(sprintf("        %-22s %-8s d=%.2f P=%s\n",
                      top_hits$Gene[k],
                      top_hits$Superior_Haplotype[k],
                      top_hits$Cohens_d[k],
                      hbr_format_p(top_hits$P_value[k])))
        if (nrow(superior_genes_list) > 5L)
          cat(sprintf("        ... (+%d more)\n",
                      nrow(superior_genes_list) - 5L))
      }

      if (plot_style %in% c("gene_list", "both")) {
        if (verbose) cat("[B.5b] Plotting superior genes ........ ")
        superior_ids <- as.character(superior_summary[[gid_col]])
        superior_ids <- superior_ids[!is.na(superior_ids) &
                                       nzchar(superior_ids)]
        per_gene_for_plot <- if (save_only_superior)
          superior_dir else output_dir

        if (length(superior_ids)) {
          tryCatch({
            plot_out <- hbr_plot_gene_list(
              gene_list     = superior_ids,
              output_dir    = output_dir,
              per_gene_dir  = per_gene_for_plot,
              reference     = reference,
              trait_col     = trait_col,
              trait_label   = trait_label,
              panels        = plot_panels,
              palette       = plot_palette,
              theme         = plot_theme,
              table_style   = plot_table_style,
              aa_full_names = aa_full_names,
              width         = plot_width,
              height        = plot_height,
              formats       = figure_formats,
              combined_pdf  = plot_combined_pdf,
              figures_dir   = "figures",
              verbose       = FALSE)
            if (verbose)
              cat(sprintf("done (%d plotted, %d missing)\n",
                          length(plot_out$plots),
                          length(plot_out$missing_genes)))
          }, error = function(e)
            if (verbose) cat("skipped\n"))
        } else if (verbose) {
          cat("skipped (no IDs)\n")
        }
      }

      if (isTRUE(go_enrichment)) {
        if (verbose) cat("[B.5c] Publication-table annotation ... ")
        tryCatch({
          ann_out <- hbr_superior_candidate_genes_functional_characterization(
            candidate_genes            = summary_table_clean,
            functional_annotation_file = NULL,
            goslim_file                = NULL,
            interpro_file              = NULL,
            pfam_file                  = NULL,
            source                     = go_enrichment_source,
            cache_dir                  = go_enrichment_cache,
            output_dir                 = output_dir,
            output_prefix              = "publication_table_annotated",
            verbose                    = FALSE)

          if (!is.null(ann_out$table) && nrow(ann_out$table)) {
            data.table::fwrite(
              ann_out$table,
              file.path(output_dir, "publication_table_functional.csv"))
            if (save_only_superior)
              data.table::fwrite(
                ann_out$table,
                file.path(superior_dir,
                          "publication_table_functional.csv"))
            if (verbose) cat("done\n")
          } else if (verbose) cat("done (no rows)\n")
        }, error = function(e)
          if (verbose) cat("skipped\n"))

        if (verbose) cat("[B.5d] GO enrichment ................. ")
        go_out <- tryCatch({
          ann <- hbr_read_osa1_r7_annotations(
            source    = go_enrichment_source,
            cache_dir = go_enrichment_cache,
            verbose   = FALSE)

          if (is.null(ann$goslim) || !nrow(ann$goslim))
            stop("No GOSlim annotation available from ",
                 go_enrichment_source, call. = FALSE)

          sup_genes <- unique(as.character(superior_summary[[gid_col]]))
          sup_genes <- sup_genes[!is.na(sup_genes) & nzchar(sup_genes)]

          if (!is.null(summary_table_clean) &&
              gid_col %in% names(summary_table_clean)) {
            uni_genes <- unique(as.character(
              summary_table_clean[[gid_col]]))
          } else {
            uni_genes <- unique(as.character(summary_table[[gid_col]]))
          }
          uni_genes <- uni_genes[!is.na(uni_genes) & nzchar(uni_genes)]

          if (length(uni_genes) <= length(sup_genes)) {
            NULL
          } else {
            sup_genes <- hbr_normalize_osa1_gene_id(sup_genes)
            uni_genes <- hbr_normalize_osa1_gene_id(uni_genes)

            gs        <- ann$goslim
            go_ids    <- unique(as.character(gs$Gene_ID))
            n_uni_go  <- length(intersect(uni_genes, go_ids))
            if (n_uni_go == 0L) stop("no universe gene matches GOSlim.")

            enr <- list()
            for (asp in c("BP", "MF", "CC")) {
              enr[[asp]] <- tryCatch(
                hbr_go_enrichment(
                  genes_superior = sup_genes,
                  genes_all      = uni_genes,
                  goslim_dt      = ann$goslim,
                  aspect         = asp,
                  min_term_size  = go_enrichment_min_term,
                  min_superior   = go_enrichment_min_superior,
                  verbose        = FALSE),
                error = function(e) data.table::data.table())
            }

            combined <- data.table::rbindlist(enr, idcol = "Aspect")
            if (nrow(combined))
              data.table::fwrite(
                combined,
                file.path(output_dir,
                          "GO_enrichment_superior_genes.csv"))

            for (asp in names(enr)) {
              if (nrow(enr[[asp]]))
                data.table::fwrite(
                  enr[[asp]],
                  file.path(output_dir,
                            paste0("GO_enrichment_", asp, ".csv")))
            }
            combined
          }
        }, error = function(e) NULL)

        if (!is.null(go_out) && nrow(go_out)) {
          n_sig <- sum(go_out$Adjusted_P < 0.05, na.rm = TRUE)
          if (verbose)
            cat(sprintf("done (%d terms, %d FDR-significant)\n",
                        nrow(go_out), n_sig))
        } else if (verbose) {
          cat("done (no terms)\n")
        }
      }

    } else if (verbose) {
      cat("      No superior genes identified.\n")
    }
  } else {
    if (verbose) cat("      No genes reached the association-testing stage.\n")
    data.table::fwrite(data.table::data.table(),
                       file.path(output_dir,
                                 "haplotype_summary_functional.csv"))
    data.table::fwrite(data.table::data.table(),
                       file.path(output_dir,
                                 "publication_table_functional.csv"))
  }

  # ----- B.6 Attrition log -----
  valid_outcomes <- outcome_log[!vapply(outcome_log, is.null, logical(1))]
  outcome_table <- if (length(valid_outcomes))
    data.table::rbindlist(valid_outcomes, fill = TRUE)
  else data.table::data.table(
    Index = integer(), Gene = character(),
    Exit_Stage = character(), Exit_Reason = character(),
    N_Functional_SNPs = integer(), N_Genotyped_SNPs = integer(),
    N_After_QC = integer(), N_Individuals = integer(),
    N_Haplotypes = integer())

  if (isTRUE(require_fdr) && nrow(outcome_table) && nrow(summary_table) &&
      all(c("Is_Superior", "Is_Superior_Nominal") %in% names(summary_table)) &&
      gid_col %in% names(summary_table)) {
    fdr_fail <- as.character(summary_table[
      Is_Superior_Nominal %in% TRUE & !(Is_Superior %in% TRUE)][[gid_col]])
    ix <- outcome_table$Exit_Stage == "superior" &
      outcome_table$Gene %in% fdr_fail
    if (any(ix, na.rm = TRUE))
      outcome_table[which(ix), `:=`(
        Exit_Stage  = "superior_nominal_not_fdr",
        Exit_Reason = "passed the nominal superior gates but not BH-FDR")]
  }

  data.table::fwrite(outcome_table,
                     file.path(output_dir, "gene_attrition_log.csv"))

  if (verbose) {
    att <- outcome_table[, .N, by = Exit_Stage][order(-N)]
    cat("[B.6] Attrition:\n")
    for (k in seq_len(nrow(att)))
      cat(sprintf("        %-32s %d\n", att$Exit_Stage[k], att$N[k]))
  }

  # ----- B.7 Final summary + RDS -----
  n_tested   <- if (nrow(summary_table)) nrow(summary_table) else 0L
  n_fdr      <- if (nrow(summary_table) &&
                    "Significant_FDR" %in% names(summary_table))
    sum(summary_table$Significant_FDR, na.rm = TRUE) else 0L
  n_superior <- if (nrow(summary_table) &&
                    "Is_Superior" %in% names(summary_table))
    sum(summary_table$Is_Superior, na.rm = TRUE) else 0L

  pipeline_summary <- data.table::data.table(
    Mode                 = if (mode_gwas) "GWAS" else "DirectGene",
    Annotation_Source    = annotation_source,
    Reference_Source     = reference_source,
    Crop_Species         = ref_paths$crop_species %||% NA_character_,
    Crop_Assembly        = ref_paths$crop_assembly %||% NA_character_,
    Ensembl_Release      = ref_paths$ensembl_release %||% NA_integer_,
    Trait                = trait_col,
    Direction            = direction,
    PLINK_Prefix         = plink_prefix,
    Ref_Alt_Source       = if (!is.null(ref_alt)) "user"
    else if (isTRUE(harmonize_strand)) "harmonized"
    else "plink_bim",
    Harmonize_Strand     = isTRUE(harmonize_strand),
    Plot_Style           = plot_style,
    Plot_Palette         = plot_palette,
    Plot_Theme           = plot_theme,
    Genes_Requested      = if (mode_gene) length(gene_list) else NA_integer_,
    Genes_Resolved       = nrow(candidate_genes),
    Genes_Analysed       = n_genes,
    Genes_Tested         = n_tested,
    Significant_FDR      = n_fdr,
    Superior_Haplotypes  = n_superior,
    LD_Decay_bp          = ld_decay_distance,
    LD_Decay_kb          = if (is.finite(ld_decay_distance))
      ld_decay_distance / 1000 else NA_real_,
    GO_Enrichment        = isTRUE(go_enrichment),
    Save_Only_Superior   = save_only_superior,
    Validation_Skipped   = isTRUE(skip_validation),
    Plink_Path           = plink_path,
    QC_Report_Dir        = qc_report_dir %||% NA_character_,
    Save_Plot            = isTRUE(save_plot),
    RDS_Compress         = if (isTRUE(save_breeding_rds)) rds_compress else NA_character_,
    RDS_Reference        = if (isTRUE(save_breeding_rds)) rds_reference else NA_character_,
    Stage                = "gene_identification",
    Schema_Version       = "2",
    FDR_Required         = isTRUE(require_fdr),
    Alpha_Gate           = alpha,
    Min_Cohens_D_Gate    = min_cohens_d,
    Min_Ind_Per_Hap_Gate = min_ind_per_hap,
    Min_Hap_Freq_Gate    = min_hap_freq,
    Max_Miss_SNP_Gate    = max_miss_snp,
    MAF_Min_Gate         = maf_min,
    Max_Het_Ind_Gate     = max_het_ind,
    Max_Miss_Ind_Gate    = max_miss_ind,
    Min_Individuals_QC   = min_individuals_qc,
    N_Perm_Gate          = n_perm,
    N_Boot_Gate          = n_boot,
    Posthoc_Gate         = paste(posthoc, collapse = "+"),
    PC_Cols_Gate         = paste(pc_cols, collapse = ","),
    Trait_Threshold_Gate = if (is.null(trait_threshold)) NA_real_ else
                             as.numeric(trait_threshold),
    Require_FDR_Gate     = isTRUE(require_fdr),
    Output_Dir           = output_dir)

  data.table::fwrite(pipeline_summary,
                     file.path(output_dir, "pipeline_summary.csv"))

  # Phenotype table actually used (ID-matched, QC-passed). Saved so that
  # hbr_advanced_decision() can recover the exact per-individual values.
  tryCatch(
    data.table::fwrite(
      phenotype_data[, c(pheno_id_col, trait_col), with = FALSE],
      file.path(output_dir, "phenotype_data_used.csv")),
    error = function(e)
      warning("Could not write phenotype_data_used.csv: ",
              conditionMessage(e), call. = FALSE))

  if (isTRUE(save_breeding_rds)) {
    tryCatch(
      .hbr_write_breeding_rds(
        superior_results = superior_gene_results,
        geno_mat         = genotype,        # hbr_genotype source; streamed pack
        reference        = reference,
        ref_paths        = ref_paths,
        plink_prefix     = plink_prefix,
        pipeline_summary = pipeline_summary,
        output_dir       = output_dir,
        rds_compress     = rds_compress,
        rds_reference    = rds_reference,
        verbose          = verbose),
      error = function(e)
        warning("Could not save breeding RDS: ",
                conditionMessage(e), call. = FALSE))
  }

  if (verbose) {
    cat(strrep("-", 75), "\n")
    cat(sprintf("[DONE] %s | %d analysed, %d tested, %d FDR, %d superior | %.1f min\n",
                pipeline_summary$Mode, n_genes, n_tested, n_fdr, n_superior,
                as.numeric(difftime(Sys.time(), .t_pipeline, units = "mins"))))
    cat("       output: ", output_dir, "\n", sep = "")
    cat(strrep("=", 75), "\n")
  }

  invisible(list(
    mode              = pipeline_summary$Mode,
    reference         = reference,
    reference_paths   = ref_paths,
    plink_prefix      = plink_prefix,
    geno_mat          = genotype,             # hbr_genotype (streaming or in-memory)
    geno_ids          = genotype$individuals,
    strand_report     = genotype$strand_report %||% NULL,
    ld_decay          = ld_decay_summary,
    genomic_regions   = if (mode_gwas) regions else NULL,
    region_snp_map    = if (mode_gwas) region_snp_map else NULL,
    gene_annotation   = if (mode_gwas) annotation else NULL,
    candidate_genes   = candidate_genes,
    analysis_gene_map = analysis_gene_map,
    ref_alt           = ref_alt,
    results           = all_results,
    superior_results  = superior_gene_results,
    summary           = summary_table,
    superior_summary  = superior_summary,
    attrition         = outcome_table,
    pipeline_summary  = pipeline_summary,
    id_report         = id_report %||% NULL,
    qc_report_dir     = qc_report_dir,
    phenotype_data      = phenotype_data,
    functional_snps_all = functional_snps_all,
    gate_params         = list(
      alpha              = alpha,
      min_cohens_d       = min_cohens_d,
      min_ind_per_hap    = min_ind_per_hap,
      min_hap_freq       = min_hap_freq,
      max_miss_snp       = max_miss_snp,
      maf_min            = maf_min,
      max_het_ind        = max_het_ind,
      max_miss_ind       = max_miss_ind,
      min_individuals_qc = min_individuals_qc,
      n_perm             = n_perm,
      n_boot             = n_boot,
      posthoc            = posthoc,
      pc_cols            = pc_cols,
      trait_threshold    = trait_threshold,
      require_fdr        = require_fdr)))
}