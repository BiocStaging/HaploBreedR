# ============================================================================
# R/plot-functions.r
# All user-facing plotting functions.
#
# ----------------------------------------------------------------------------
# Per-gene boxplot + table (hbr_plot_haplotype)
# ----------------------------------------------------------------------------

#' Plot haplotype boxplot with summary table
#'
#' Draws a two-panel figure: a haplotype boxplot (with jittered points,
#' mean +/- SE overlay, and post-hoc letters) and a haplotype summary
#' table (N, frequency, per-SNP alleles, mean/SE, superior marker). The
#' function has two modes: **pipeline mode** (supply `gene_id` +
#' `hap_data` + `test_res`) and **CSV mode** (supply `gene_csv` or
#' `gene_data` + optional individuals/genotype files).
#'
#' @param gene_id Gene ID (pipeline mode).
#' @param hap_data Test objects from the pipeline (pipeline mode).
#' @param test_res Test results from the pipeline (pipeline mode).
#' @param trait_col Trait column name.
#' @param output_dir Output directory. `NULL` writes to the working
#'   directory.
#' @param superior_info Superior-haplotype info (pipeline mode).
#' @param functional_snps_used SNPs used (pipeline mode).
#' @param adjusted_p Optional FDR-adjusted p-value.
#' @param aa_full_names Use three-letter amino-acid codes.
#' @param formats Output formats. Any of `"pdf"`, `"png"`, `"tiff"`,
#'   `"svg"`.
#' @param gene_list Alternative: gene ID (CSV mode).
#' @param gene_data,gene_csv Alternative: in-memory or CSV input (CSV
#'   mode).
#' @param individuals_csv,individuals_data Alternative individuals source
#'   (CSV mode).
#' @param geno_file,geno_matrix,geno_format Alternative genotype source
#'   (CSV mode).
#' @param plink_path Path to PLINK.
#' @return Invisibly, a list with `plot` (the composed grob) and `files`
#'   (named list of paths written).
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' hbr_plot_haplotype(
#'   gene_list       = "gene1",
#'   gene_csv        = file.path(toy, "gene1_viz.csv"),
#'   individuals_csv = file.path(toy, "gene1_individuals.csv"),
#'   trait_col       = "trait_value",
#'   output_dir      = tempdir(),
#'   formats         = "pdf")
#'
#' @export
hbr_plot_haplotype <- function(
    gene_id              = NULL,
    hap_data             = NULL,
    test_res             = NULL,
    trait_col            = "trait",
    output_dir           = NULL,
    superior_info        = NULL,
    functional_snps_used = NULL,
    adjusted_p           = NA_real_,
    aa_full_names        = TRUE,
    formats              = c("pdf", "tiff", "png"),
    gene_list            = NULL,
    gene_data            = NULL,
    individuals_csv      = NULL,
    individuals_data     = NULL,
    geno_file            = NULL,
    geno_matrix          = NULL,
    geno_format          = NULL,
    plink_path           = "plink",
    gene_csv             = NULL) {

  palette <- "haplobreed"
  theme   <- "haplobreed"

  csv_mode <- !is.null(gene_csv) || !is.null(gene_data) || !is.null(gene_list)

  if (csv_mode) {
    g <- if (!is.null(gene_list)) gene_list[1]
    else if (!is.null(gene_data)) as.character(gene_data$gene_id[1])
    else sub("_viz\\.csv$", "", basename(gene_csv))

    res <- hbr_build_gene_object(
      gene_csv         = gene_csv,
      gene_data        = gene_data,
      individuals_csv  = individuals_csv,
      individuals_data = individuals_data,
      geno_file        = geno_file,
      geno_matrix      = geno_matrix,
      geno_format      = geno_format,
      trait_col        = trait_col,
      plink_path       = plink_path)

    hap_data             <- list(haplotypes    = res$individuals_df,
                                 haplotype_map = res$hap_map)
    test_res             <- list(data          = res$individuals_df,
                                 groups        = res$groups,
                                 p_value       = res$p_value,
                                 permutation_p = res$permutation_p,
                                 test          = "CSV")
    superior_info        <- res$superior
    functional_snps_used <- res$functional_snps
    gene_id              <- g
  }

  if (is.null(gene_id))
    stop("Supply `gene_id` (pipeline mode) OR `gene_csv` / `gene_data` ",
         "(manual mode).", call. = FALSE)

  if (is.null(test_res) || !nrow(test_res$data))
    stop("No test data available for '", gene_id, "'.", call. = FALSE)

  pal_anchors <- .hbr_palette_anchors(palette)
  srf         <- .hbr_theme_surface(theme)
  th          <- hbr_theme(theme, base_size = 11)
  accent      <- pal_anchors$superior

  dat     <- data.table::copy(data.table::as.data.table(test_res$data))
  groups  <- data.table::copy(data.table::as.data.table(test_res$groups))
  hap_map <- if (!is.null(hap_data$haplotype_map))
    data.table::as.data.table(hap_data$haplotype_map)
  else data.table::data.table(Haplotype = unique(dat$Haplotype))

  superior_hap <- if (!is.null(superior_info)) superior_info$haplotype else NULL
  if (is.null(superior_hap) || is.na(superior_hap))
    superior_hap <- as.character(groups$Haplotype[which.max(groups$Mean)])

  dat[, Haplotype := as.character(Haplotype)]
  dat[[trait_col]] <- as.numeric(dat[[trait_col]])
  dat <- dat[!is.na(dat[[trait_col]])]

  ord <- groups$Haplotype[order(-as.numeric(groups$Mean))]
  dat[, Haplotype := factor(Haplotype, levels = ord)]
  groups[, Haplotype := factor(as.character(Haplotype), levels = ord)]

  pal_cols <- stats::setNames(rep(srf$text_muted, length(ord)), ord)
  if (superior_hap %in% names(pal_cols)) pal_cols[superior_hap] <- accent

  y_max <- max(dat[[trait_col]], na.rm = TRUE)
  y_min <- min(dat[[trait_col]], na.rm = TRUE)
  y_rng <- y_max - y_min
  if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1

  letter_col <- intersect(
    c("Letter", "letter", "Significance_Group",
      "Duncan_Group", "Tukey_Group"),
    names(groups))[1]

  lbl_df <- NULL
  if (!is.na(letter_col)) {
    grp_max <- dat[, .(y = max(get(trait_col), na.rm = TRUE) +
                         y_rng * 0.05),
                   by = Haplotype]
    ldf <- merge(
      data.table::data.table(
        Haplotype = as.character(groups$Haplotype),
        Letter    = as.character(groups[[letter_col]])),
      grp_max, by = "Haplotype", all.x = TRUE)
    ldf[, Haplotype := factor(Haplotype, levels = ord)]
    lbl_df <- ldf[!is.na(y)]
  }

  n_by_hap <- dat[, .N, by = Haplotype]
  x_labels <- stats::setNames(
    sprintf("%s\n(N = %d)", n_by_hap$Haplotype, n_by_hap$N),
    as.character(n_by_hap$Haplotype))

  p_box <- ggplot2::ggplot(dat,
                           ggplot2::aes(x = Haplotype,
                                        y = .data[[trait_col]],
                                        fill = Haplotype)) +
    ggplot2::geom_boxplot(width = 0.58, linewidth = 0.55, alpha = 0.85,
                          outlier.shape = NA, colour = srf$text_primary)

  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    p_box <- p_box + ggbeeswarm::geom_quasirandom(
      width = 0.22, size = 1.35, alpha = 0.55,
      colour = srf$text_primary, shape = 16)
  } else {
    p_box <- p_box + ggplot2::geom_jitter(
      width = 0.12, size = 1.35, alpha = 0.45,
      colour = srf$text_primary, seed = 42)
  }

  p_box <- p_box +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23,
                          fill = "white", colour = srf$text_primary,
                          size = 2.6, stroke = 0.55) +
    ggplot2::geom_errorbar(
      data = groups,
      ggplot2::aes(x = Haplotype, y = Mean,
                   ymin = Mean - SE, ymax = Mean + SE),
      inherit.aes = FALSE,
      colour = accent, linewidth = 0.6, width = 0.15) +
    ggplot2::scale_fill_manual(values = pal_cols) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::coord_cartesian(ylim = c(y_min - y_rng * 0.05,
                                      y_max + y_rng * 0.20),
                             clip = "off")

  if (!is.null(lbl_df) && nrow(lbl_df)) {
    p_box <- p_box + ggplot2::geom_text(
      data = lbl_df,
      ggplot2::aes(x = Haplotype, y = y, label = Letter),
      inherit.aes = FALSE, fontface = "bold",
      size = 4.5, colour = srf$text_primary, na.rm = TRUE)
  }

  sub_parts <- character(0)
  if (!is.null(test_res$test) && !is.na(test_res$test))
    sub_parts <- c(sub_parts, as.character(test_res$test))
  if (!is.null(test_res$p_value) && is.finite(test_res$p_value))
    sub_parts <- c(sub_parts, paste0("P = ", hbr_format_p(test_res$p_value)))
  if (!is.null(test_res$permutation_p) && is.finite(test_res$permutation_p))
    sub_parts <- c(sub_parts,
                   paste0("P_perm = ",
                          hbr_format_p(test_res$permutation_p)))
  if (!is.na(adjusted_p) && is.finite(adjusted_p))
    sub_parts <- c(sub_parts, paste0("FDR = ", hbr_format_p(adjusted_p)))
  if (!is.null(superior_info) && !is.null(superior_info$cohens_d) &&
      is.finite(superior_info$cohens_d))
    sub_parts <- c(sub_parts, sprintf("d = %.2f", superior_info$cohens_d))

  if (length(sub_parts)) {
    p_box <- p_box + ggplot2::annotate(
      "label", x = Inf, y = Inf,
      label = paste(sub_parts, collapse = "\n"),
      hjust = 1.05, vjust = 1.10, size = 2.7,
      linewidth = 0.25, label.padding = grid::unit(0.3, "lines"),
      fill = srf$table_bg, colour = srf$text_primary)
  }

  p_box <- p_box +
    ggplot2::labs(x = NULL, y = trait_col) +
    th +
    ggplot2::theme(
      legend.position    = "none",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(6, 10, 6, 10),
      axis.text.x        = ggplot2::element_text(size = 10,
                                                 colour = srf$text_primary),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 11))

  has_cds <- !is.null(functional_snps_used) &&
    nrow(functional_snps_used) &&
    "Region" %in% names(functional_snps_used) &&
    any(functional_snps_used$Region == "CDS_nonsynonymous", na.rm = TRUE)

  tbl_df <- data.table::data.table(
    Haplotype  = as.character(hap_map$Haplotype),
    N          = as.integer(hap_map$Count),
    `Freq (%)` = round(100 * as.numeric(hap_map$Frequency), 1),
    NT         = as.character(hap_map$NT_Sequence),
    AA         = if (has_cds) as.character(hap_map$AA_Sequence) else "")

  if (!is.null(functional_snps_used) && nrow(functional_snps_used)) {
    alleles_list <- strsplit(as.character(hap_map$NT_Sequence), "-",
                             fixed = TRUE)
    for (j in seq_len(nrow(functional_snps_used))) {
      snp <- functional_snps_used[j]
      hdr <- .hbr_short_snp_header(snp, j)
      tbl_df[[hdr]] <- vapply(alleles_list,
                              function(a) if (length(a) >= j) a[j] else "",
                              character(1))
    }
  }

  tbl_df[["Superior"]] <- ifelse(
    as.character(tbl_df$Haplotype) == superior_hap, "*", "")

  gd_tbl <- data.table::data.table(
    Haplotype = as.character(groups$Haplotype),
    Mean      = round(as.numeric(groups$Mean), 2),
    SE        = round(as.numeric(groups$SE),   2))
  tbl_df <- merge(tbl_df, gd_tbl, by = "Haplotype", all.x = TRUE)
  tbl_df <- tbl_df[order(-Mean)]

  col_names <- names(tbl_df)
  seq_cols  <- intersect(c("NT", "AA"), col_names)
  meta_cols <- c("Haplotype", "N", "Freq (%)", "Superior", "Mean", "SE")
  snp_cols  <- setdiff(col_names, c(meta_cols, seq_cols))

  n_col <- ncol(tbl_df)
  fnt_core <- if      (n_col <=  5) 10.0
  else if (n_col <=  8) 9.5
  else if (n_col <= 12) 8.8
  else if (n_col <= 16) 7.8
  else if (n_col <= 20) 6.8
  else 6.0

  widths <- vapply(col_names, function(cn) {
    if (cn == "Haplotype")               1.8
    else if (cn %in% seq_cols)           2.5
    else if (cn %in% c("N", "Superior")) 0.7
    else if (grepl("Freq", cn))          0.9
    else if (cn %in% c("Mean", "SE"))    0.9
    else                                 0.85
  }, numeric(1))
  widths <- widths / sum(widths)

  fill_mat <- .hbr_hap_table_fill(
    tbl_df, snp_cols,
    superior_haplotype = superior_hap,
    cell_colours       = hbr_cell_colours(theme, palette),
    neutral            = srf$table_bg,
    superior_row       = srf$table_superior_row)

  sup_row <- which(tbl_df$Haplotype == superior_hap)
  fg_mat  <- matrix("plain", nrow(tbl_df), n_col,
                    dimnames = list(NULL, col_names))
  if (length(sup_row)) fg_mat[sup_row, ] <- "bold"

  tg <- gridExtra::tableGrob(
    tbl_df, rows = NULL,
    widths = grid::unit(widths, "npc"),
    theme = gridExtra::ttheme_minimal(
      core = list(
        fg_params = list(fontsize = fnt_core, fontface = fg_mat,
                         col = srf$text_primary),
        bg_params = list(fill = fill_mat, col = NA)),
      colhead = list(
        fg_params = list(fontface = "bold", fontsize = fnt_core + 0.6,
                         col = srf$text_primary),
        bg_params = list(fill = srf$table_header, col = NA))))

  tg <- gtable::gtable_add_grob(
    tg, grid::segmentsGrob(gp = grid::gpar(lwd = 1.1,
                                            col = srf$rule_strong)),
    1, 1, 1, ncol(tg), name = "rule_top")
  tg <- gtable::gtable_add_grob(
    tg, grid::segmentsGrob(gp = grid::gpar(lwd = 0.8,
                                            col = srf$rule_strong)),
    2, 1, 2, ncol(tg), name = "rule_header")
  tg <- gtable::gtable_add_grob(
    tg, grid::segmentsGrob(gp = grid::gpar(lwd = 1.1,
                                            col = srf$rule_strong)),
    nrow(tg), 1, nrow(tg), ncol(tg), name = "rule_bottom")
  tg <- gtable::gtable_add_padding(tg, grid::unit(c(6, 10, 6, 10), "pt"))

  p_box_t <- .hbr_add_panel_tag(p_box, "A", colour = srf$text_primary)
  tg_t    <- .hbr_add_panel_tag(tg,    "B", colour = srf$text_primary)
  aligned <- .hbr_align_panels(list(p_box_t, tg_t))

  master <- grid::textGrob(
    paste0(gene_id, " \u2014 haplotype summary"),
    gp = grid::gpar(fontface = "bold", fontsize = 13,
                    col = srf$text_primary),
    x = grid::unit(10, "pt"), hjust = 0, vjust = 1)

  fig <- gridExtra::arrangeGrob(
    aligned[[1]], aligned[[2]],
    ncol    = 1,
    heights = c(0.62, 0.38 + nrow(tbl_df) * 0.02),
    top     = master,
    padding = grid::unit(2, "lines"))

  figures_out <- if (!is.null(output_dir)) output_dir else getwd()
  if (!dir.exists(figures_out))
    dir.create(figures_out, recursive = TRUE, showWarnings = FALSE)

  files <- .hbr_save_fig(fig, gene_id, "haplotype",
                         figures_out, width = 12, height = 9, formats)

  invisible(list(plot = fig, files = files))
}


# ----------------------------------------------------------------------------
# Batch 4-panel figure per gene (hbr_plot_gene_list)
# ----------------------------------------------------------------------------

#' Plot publication-quality haplotype figures for a list of genes
#'
#' Builds a four-panel figure per gene (gene model, haplotype network,
#' haplotype summary table, trait distribution) and saves each figure in
#' the requested formats. `gene_list` may be a character vector of gene
#' IDs, a path to a summary CSV, a directory containing the pipeline
#' outputs, or an in-memory `data.table`.
#'
#' @param gene_list IDs, CSV path, directory path, or data.frame.
#' @param output_dir Directory holding per-gene files.
#' @param summary_csv Optional path to the summary CSV.
#' @param per_gene_dir Override for the directory holding per-gene files.
#' @param result Optional in-memory [haplobreed_analysis()] output.
#' @param reference Optional `hbr_reference` object.
#' @param gff3_path Optional GFF3 path; parsed once and cached for the
#'   session.
#' @param trait_col Trait column name.
#' @param trait_label Human-readable trait label.
#' @param n_gene Keep at most this many genes.
#' @param sort_by Character vector of columns to rank by.
#' @param panels Panels to draw.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param table_style `"color_coded"` or `"row_highlight"`.
#' @param include_nt,include_aa Include NT / AA sequence columns.
#' @param aa_full_names Use three-letter amino-acid codes.
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param combined_pdf Also write a multi-page PDF.
#' @param figures_dir Subdirectory of `output_dir` for figures.
#' @param verbose Print progress.
#' @return Invisibly, a list with `plots`, `files`, `missing_genes`,
#'   `combined_pdf`.
#' @keywords internal
hbr_plot_gene_list <- function(
    gene_list,
    output_dir   = NULL,
    summary_csv  = NULL,
    per_gene_dir = NULL,
    result       = NULL,
    reference    = NULL,
    gff3_path    = NULL,
    trait_col    = "NT",
    trait_label  = NULL,
    n_gene       = NULL,
    sort_by      = c("Cohens_d", "Superior_Mean"),
    panels       = c("gene_model", "haplotype_network",
                     "haplotype_table", "trait_boxplot"),
    palette      = "haplobreed",
    theme        = "haplobreed",
    table_style  = c("color_coded", "row_highlight"),
    include_nt   = FALSE,
    include_aa   = FALSE,
    aa_full_names = TRUE,
    width        = 12,
    height       = 14,
    formats      = c("pdf", "png"),
    combined_pdf = TRUE,
    figures_dir  = "figures",
    verbose      = TRUE) {

  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)
  table_style <- match.arg(table_style)

  reference <- .hbr_resolve_reference(
    reference = reference, result = result,
    gff3_path = gff3_path, verbose = verbose)

  if (is.null(output_dir) && !is.null(summary_csv))
    output_dir <- if (dir.exists(summary_csv)) summary_csv
                  else dirname(summary_csv)
  if (is.null(per_gene_dir)) per_gene_dir <- output_dir
  if (is.null(per_gene_dir) || !dir.exists(per_gene_dir))
    stop("Cannot locate per-gene files. Supply `output_dir`.", call. = FALSE)

  if (is.null(summary_csv)) {
    for (fn in c("publication_table_functional.csv",
                 "haplotype_summary_functional.csv",
                 "superior_genes_summary.csv")) {
      cand <- file.path(output_dir, fn)
      if (file.exists(cand)) { summary_csv <- cand; break }
    }
  }

  genes <- .hbr_resolve_genes(
    gene_list   = if (missing(gene_list)) NULL else gene_list,
    summary_csv = summary_csv, output_dir = output_dir,
    n_gene      = n_gene, sort_by = sort_by, verbose = verbose)

  if (!length(genes))
    stop("No genes to plot. Supply `gene_list` or a valid output dir.",
         call. = FALSE)

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (verbose)
    message("Plotting ", length(genes), " gene(s) | palette = ", palette,
            " | theme = ", theme, " | table_style = ", table_style)

  fig_dir <- if (!is.null(figures_dir) && nzchar(figures_dir))
    file.path(output_dir, figures_dir) else output_dir
  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  plots     <- list()
  files_out <- list()
  missing   <- character(0)

  for (g in genes) {
    if (verbose) message("\n[", g, "]")

    pieces <- tryCatch(
      .hbr_load_gene_pieces(gene_id = g, per_gene_dir = per_gene_dir,
                            summary_csv = summary_csv,
                            trait_col = trait_col, verbose = verbose),
      error = function(e) {
        if (verbose) message("  load failed: ", conditionMessage(e))
        NULL
      })
    if (is.null(pieces)) { missing <- c(missing, g); next }

    fig <- tryCatch(
      .hbr_build_4panel(
        gene_id       = g,
        pieces        = pieces,
        panels        = panels,
        trait_col     = trait_col,
        trait_label   = trait_label,
        aa_full_names = aa_full_names,
        palette       = palette,
        theme         = theme,
        include_nt    = include_nt,
        include_aa    = include_aa,
        table_style   = table_style,
        reference     = reference),
      error = function(e) {
        if (verbose) message("  panel build failed: ", conditionMessage(e))
        NULL
      })
    if (is.null(fig)) { missing <- c(missing, g); next }

    files <- .hbr_save_fig(fig, g, "haplotype_figure", fig_dir,
                           width, height, formats)
    plots[[g]]     <- fig
    files_out[[g]] <- files
  }

  combined_path <- .hbr_write_combined(plots, fig_dir, combined_pdf,
                                       width, height, verbose)

  if (verbose) {
    message("\n", strrep("=", 60))
    message("Genes requested : ", length(genes))
    message("Genes plotted   : ", length(plots))
    message("Genes missing   : ", length(missing))
    if (length(missing))
      message("Missing IDs     : ",
              paste(utils::head(missing, 5), collapse = ", "),
              if (length(missing) > 5) " ..." else "")
    message("Figures dir     : ", fig_dir)
    if (!is.null(combined_path)) message("Combined PDF    : ", combined_path)
    message(strrep("=", 60))
  }

  invisible(list(plots         = plots,
                 files         = files_out,
                 missing_genes = missing,
                 combined_pdf  = combined_path))
}


# ----------------------------------------------------------------------------
# Gene structure (hbr_plot_gene_structure)
# ----------------------------------------------------------------------------

#' Plot gene structure for one or more genes
#'
#' Draws UTRs, exons, intron chevrons (strand-aware), and optional
#' functional-variant markers for one or more genes. When
#' `facet_rows = FALSE` and more than one gene is requested, a combined-
#' axis facet figure is also written.
#'
#' @param gene_list Gene IDs, path to a summary CSV, directory, or
#'   data.frame.
#' @param reference A `hbr_reference` object.
#' @param gff3_path Path to a GFF3 file. Parsed once and cached.
#' @param annotation_source `"MSU7"`, `"RAPDB"`, or `"CUSTOM"`.
#' @param cdna_path,cds_path,pep_path,genome_fasta_path Optional FASTA
#'   paths.
#' @param output_dir Directory for figures.
#' @param figures_dir Subdirectory of `output_dir`.
#' @param n_gene Keep at most this many genes.
#' @param sort_by Character vector of columns to rank by.
#' @param fun_variants Optional named list or data.frame.
#' @param variant_panel_dir Directory with per-gene functional SNP files.
#' @param show_variants Draw functional-variant markers.
#' @param variant_shape_by_region Shape by region.
#' @param aa_full_names Expand `R65G` to `Arg65Gly`.
#' @param facet_rows One gene per row or all on a single axis.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param combined_pdf Also write a multi-page PDF.
#' @param verbose Print progress.
#' @return Invisibly, a list with `plots`, `files`, `combined_pdf`,
#'   `missing_genes`.
#'   
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' hbr_plot_gene_structure(
#'   gene_list         = "gene1",
#'   gff3_path         = file.path(toy, "toy.gff3"),
#'   genome_fasta_path = file.path(toy, "toy.fa"),
#'   output_dir        = tempdir(),
#'   formats           = "pdf",
#'   verbose           = FALSE)
#'
#' @export
hbr_plot_gene_structure <- function(
    gene_list,
    reference            = NULL,
    gff3_path            = NULL,
    annotation_source    = c("MSU7", "RAPDB", "CUSTOM"),
    cdna_path            = NULL,
    cds_path             = NULL,
    pep_path             = NULL,
    genome_fasta_path    = NULL,
    output_dir           = NULL,
    figures_dir          = "figures",
    n_gene               = NULL,
    sort_by              = c("Cohens_d", "Superior_Mean"),
    fun_variants         = NULL,
    variant_panel_dir    = NULL,
    show_variants        = TRUE,
    variant_shape_by_region = TRUE,
    aa_full_names        = TRUE,
    facet_rows           = TRUE,
    palette              = "haplobreed",
    theme                = "scientific",
    width                = 10,
    height               = 4,
    formats              = c("pdf", "png"),
    combined_pdf         = TRUE,
    verbose              = TRUE) {

  annotation_source <- match.arg(annotation_source)
  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)

  reference <- .hbr_resolve_reference(
    reference = reference, gff3_path = gff3_path, verbose = verbose)

  if (is.null(reference)) {
    if (!is.null(gff3_path) && file.exists(gff3_path)) {
      reference <- tryCatch(
        hbr_load_reference(
          genome       = basename(gff3_path),
          annotation   = annotation_source,
          gff3         = gff3_path,
          cds          = cds_path %||% cdna_path,
          pep          = pep_path,
          genome_fasta = genome_fasta_path),
        error = function(e) NULL)
      if (!is.null(reference))
        .hbr_set_cached_reference(gff3_path, reference)
    }
  }

  if (is.null(reference))
    stop("Supply `reference` or a valid `gff3_path` with the FASTA paths.",
         call. = FALSE)

  genes <- .hbr_resolve_genes(
    gene_list   = if (missing(gene_list)) NULL else gene_list,
    summary_csv = NULL, output_dir = NULL,
    n_gene      = n_gene, sort_by = sort_by, verbose = verbose)

  if (!length(genes)) stop("No genes to plot.", call. = FALSE)

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (verbose)
    message("Plotting gene structure for ", length(genes), " gene(s).")

  variant_source <- .hbr_gs_resolve_variants(
    genes, fun_variants, variant_panel_dir)

  fig_dir <- if (!is.null(output_dir) && nzchar(figures_dir))
    file.path(output_dir, figures_dir)
  else output_dir
  if (is.null(fig_dir) || !nzchar(fig_dir)) fig_dir <- getwd()
  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  plots  <- list()
  files  <- list()
  missing <- character(0)

  for (g in genes) {
    gene_gff_id <- .hbr_resolve_gene_gff_id(
      g, data.table::data.table(Transcript_ID = character(0)), reference)
    if (is.null(gene_gff_id)) gene_gff_id <- g

    fig <- tryCatch(
      .hbr_gs_build_one(
        gene_id          = g,
        gene_gff_id      = gene_gff_id,
        reference        = reference,
        variants         = variant_source[[g]],
        show_variants    = show_variants,
        variant_shape_by_region = variant_shape_by_region,
        aa_full_names    = aa_full_names,
        palette          = palette,
        theme            = theme),
      error = function(e) {
        if (verbose) message("  [", g, "] failed: ", conditionMessage(e))
        NULL
      })
    if (is.null(fig)) { missing <- c(missing, g); next }

    plots[[g]] <- fig
    out <- .hbr_save_fig(fig, g, "gene_structure",
                         fig_dir, width, height, formats)
    files[[g]] <- out
  }

  combined_path <- if (isTRUE(combined_pdf) && length(plots) >= 2)
    .hbr_write_combined(plots, fig_dir, TRUE, width, height, verbose)
  else NULL

  if (!isTRUE(facet_rows) && length(genes) > 1) {
    if (verbose) message("Building combined-axis facet figure...")
    facet_fig <- tryCatch(
      .hbr_gs_build_facet(genes, reference, variant_source,
                          palette, theme, aa_full_names),
      error = function(e) NULL)
    if (!is.null(facet_fig)) {
      out <- .hbr_save_fig(facet_fig, "multi_gene", "gene_structure",
                           fig_dir, width * 1.4,
                           height * 0.5 * length(genes),
                           formats)
      files[["multi_gene"]] <- out
    }
  }

  if (verbose) {
    message("\n", strrep("=", 60))
    message("Genes plotted : ", length(plots))
    message("Genes missing : ", length(missing))
    message("Figures dir   : ", fig_dir)
    message(strrep("=", 60))
  }

  invisible(list(plots = plots, files = files,
                 combined_pdf = combined_path, missing_genes = missing))
}


# ----------------------------------------------------------------------------
# Haplotype network (hbr_plot_haplotype_network)
# ----------------------------------------------------------------------------

#' Plot a haplotype network for one or more genes
#'
#' Renders a minimal spanning graph of a gene's haplotypes: nodes sized
#' by frequency, edges drawn when two haplotypes differ by exactly one
#' SNP. Optionally highlights and labels the superior haplotype.
#'
#' @inheritParams hbr_plot_gene_list
#' @param colour_by_superior Highlight the superior node.
#' @param label_superior Label the superior node.
#' @param layout `"circle"` or `"grid"`.
#' @param node_scale Multiplier for maximum node radius.
#' @return Invisibly, a list with `plots`, `files`, `missing_genes`,
#'   `combined_pdf`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' hbr_plot_haplotype_network(
#'   gene_list    = "gene1",
#'   output_dir   = tempdir(),
#'   per_gene_dir = toy,
#'   trait_col    = "trait_value",
#'   formats      = "pdf",
#'   combined_pdf = FALSE,
#'   verbose      = FALSE)
#'
#' @export
hbr_plot_haplotype_network <- function(
    gene_list,
    output_dir   = NULL,
    summary_csv  = NULL,
    per_gene_dir = NULL,
    result       = NULL,
    reference    = NULL,
    gff3_path    = NULL,
    trait_col    = "NT",
    n_gene       = NULL,
    sort_by      = c("Cohens_d", "Superior_Mean"),
    colour_by_superior = TRUE,
    label_superior     = TRUE,
    layout             = c("circle", "grid"),
    node_scale         = 1.0,
    palette            = "haplobreed",
    theme              = "haplobreed",
    width              = 8,
    height             = 7,
    formats            = c("pdf", "png"),
    combined_pdf       = TRUE,
    figures_dir        = "figures",
    verbose            = TRUE) {

  layout <- match.arg(layout)
  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)

  .hbr_resolve_reference(reference = reference, result = result,
                          gff3_path = gff3_path, verbose = verbose)

  if (is.null(output_dir) && !is.null(summary_csv))
    output_dir <- if (dir.exists(summary_csv)) summary_csv
                  else dirname(summary_csv)
  if (is.null(per_gene_dir)) per_gene_dir <- output_dir
  if (is.null(per_gene_dir) || !dir.exists(per_gene_dir))
    stop("Cannot locate per-gene files. Supply `output_dir`.", call. = FALSE)

  genes <- .hbr_resolve_genes(
    gene_list = if (missing(gene_list)) NULL else gene_list,
    summary_csv = summary_csv, output_dir = output_dir,
    n_gene = n_gene, sort_by = sort_by, verbose = verbose)
  if (!length(genes)) stop("No genes to plot.", call. = FALSE)
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (verbose)
    message("Plotting haplotype networks for ", length(genes), " gene(s).")

  fig_dir <- if (nzchar(figures_dir %||% ""))
    file.path(output_dir, figures_dir) else output_dir
  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- list(); files_out <- list(); missing <- character(0)

  for (g in genes) {
    if (verbose) message("\n[", g, "]")
    pieces <- tryCatch(
      .hbr_load_gene_pieces(g, per_gene_dir, trait_col = trait_col,
                            verbose = verbose),
      error = function(e) NULL)
    if (is.null(pieces) || nrow(pieces$hap_map) < 2) {
      missing <- c(missing, g); next
    }
    fig <- tryCatch(
      .hbr_build_network_figure(g, pieces$hap_map, pieces$groups,
                                colour_by_superior, label_superior,
                                layout, node_scale,
                                palette = palette, theme = theme),
      error = function(e) NULL)
    if (is.null(fig)) { missing <- c(missing, g); next }
    files <- .hbr_save_fig(fig, g, "haplotype_network", fig_dir,
                           width, height, formats)
    plots[[g]] <- fig; files_out[[g]] <- files
  }

  combined_path <- .hbr_write_combined(plots, fig_dir, combined_pdf,
                                       width, height, verbose)

  if (verbose) {
    message("\n", strrep("=", 60))
    message("Genes plotted   : ", length(plots))
    message("Genes missing   : ", length(missing))
    message("Figures dir     : ", fig_dir)
    if (!is.null(combined_path))
      message("Combined PDF    : ", combined_path)
    message(strrep("=", 60))
  }

  invisible(list(plots = plots, files = files_out,
                 missing_genes = missing, combined_pdf = combined_path))
}


# ----------------------------------------------------------------------------
# Haplotype table (internal)
# ----------------------------------------------------------------------------

#' Plot the haplotype summary table as a figure
#'
#' Renders a single-panel table figure (Haplotype, N, Freq %, per-SNP
#' alleles, mean/SE) with the superior row highlighted. Internal helper.
#'
#' @inheritParams hbr_plot_gene_list
#' @param include_nt,include_aa Include NT / AA sequence columns.
#' @param include_per_snp Include one column per SNP.
#' @param highlight_superior Highlight the superior row.
#' @param height Figure height; `NULL` auto-sizes.
#' @return Invisibly, a list with `plots`, `files`, `missing_genes`,
#'   `combined_pdf`.
#' @examples
#' \dontrun{
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' hbr_plot_haplotype_table(
#'   gene_list  = "gene1",
#'   output_dir = tempdir())
#' }
#'
#' @noRd
hbr_plot_haplotype_table <- function(
    gene_list,
    output_dir   = NULL,
    summary_csv  = NULL,
    per_gene_dir = NULL,
    result       = NULL,
    reference    = NULL,
    gff3_path    = NULL,
    trait_col    = "trait",
    n_gene       = NULL,
    sort_by      = c("Cohens_d", "Superior_Mean"),
    include_nt         = FALSE,
    include_aa         = FALSE,
    include_per_snp    = TRUE,
    highlight_superior = TRUE,
    palette            = "haplobreed",
    theme              = "haplobreed",
    width              = 12,
    height             = NULL,
    formats            = c("pdf", "png"),
    combined_pdf       = TRUE,
    figures_dir        = "figures",
    verbose            = TRUE) {

  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)

  .hbr_resolve_reference(reference = reference, result = result,
                          gff3_path = gff3_path, verbose = verbose)

  if (is.null(output_dir) && !is.null(summary_csv))
    output_dir <- if (dir.exists(summary_csv)) summary_csv
                  else dirname(summary_csv)
  if (is.null(per_gene_dir)) per_gene_dir <- output_dir
  if (is.null(per_gene_dir) || !dir.exists(per_gene_dir))
    stop("Cannot locate per-gene files. Supply `output_dir`.", call. = FALSE)

  genes <- .hbr_resolve_genes(
    gene_list = if (missing(gene_list)) NULL else gene_list,
    summary_csv = summary_csv, output_dir = output_dir,
    n_gene = n_gene, sort_by = sort_by, verbose = verbose)
  if (!length(genes)) stop("No genes to plot.", call. = FALSE)
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (verbose)
    message("Plotting haplotype tables for ", length(genes), " gene(s).")

  fig_dir <- if (nzchar(figures_dir %||% ""))
    file.path(output_dir, figures_dir) else output_dir
  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- list(); files_out <- list(); missing <- character(0)

  for (g in genes) {
    if (verbose) message("\n[", g, "]")
    pieces <- tryCatch(
      .hbr_load_gene_pieces(g, per_gene_dir, trait_col = trait_col,
                            verbose = verbose),
      error = function(e) NULL)
    if (is.null(pieces)) { missing <- c(missing, g); next }

    tbl <- tryCatch(
      .hbr_build_hap_table(pieces, highlight_superior,
                           include_nt, include_aa, include_per_snp,
                           palette),
      error = function(e) NULL)
    if (is.null(tbl)) { missing <- c(missing, g); next }

    fig <- .hbr_render_hap_table_grob(g, tbl, highlight_superior,
                                       pieces$groups,
                                       palette = palette, theme = theme)

    h <- if (is.null(height)) max(2.5, 0.9 + nrow(tbl) * 0.32) else height
    files <- .hbr_save_fig(fig, g, "haplotype_table", fig_dir,
                           width, h, formats)
    plots[[g]] <- fig; files_out[[g]] <- files
  }

  combined_path <- .hbr_write_combined(plots, fig_dir, combined_pdf,
                                       width, height %||% 6, verbose)

  if (verbose) {
    message("\n", strrep("=", 60))
    message("Genes plotted   : ", length(plots))
    message("Genes missing   : ", length(missing))
    message("Figures dir     : ", fig_dir)
    if (!is.null(combined_path))
      message("Combined PDF    : ", combined_path)
    message(strrep("=", 60))
  }

  invisible(list(plots = plots, files = files_out,
                 missing_genes = missing, combined_pdf = combined_path))
}


# ----------------------------------------------------------------------------
# LD heatmap (hbr_plot_ld_heatmap)
# ----------------------------------------------------------------------------

#' Plot an LD heatmap for a gene's functional SNPs
#'
#' Renders a lower-triangular LD heatmap (r-squared between pairs of the
#' gene's functional SNPs) with the gene model drawn above it. Requires
#' genotype data for the gene's SNPs, supplied either in memory or via
#' a file.
#'
#' @inheritParams hbr_plot_gene_list
#' @param geno_file Genotype file path.
#' @param geno_matrix In-memory genotype matrix.
#' @param geno_format Genotype file format override.
#' @param plink_path Path to the PLINK executable.
#' @param stat `"r2"` (currently only r-squared is supported).
#' @return Invisibly, a list with `plots`, `files`, `missing_genes`,
#'   `combined_pdf`, `ld_matrix`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' geno <- matrix(sample(c(0L, 1L, 2L), 30, replace = TRUE),
#'                nrow = 10, ncol = 3,
#'                dimnames = list(sprintf("IND%02d", 1:10),
#'                                c("snp005", "snp020", "snp045")))
#' hbr_plot_ld_heatmap(
#'   gene_list    = "gene1",
#'   output_dir   = tempdir(),
#'   per_gene_dir = toy,
#'   geno_matrix  = geno,
#'   formats      = "pdf",
#'   combined_pdf = FALSE,
#'   verbose      = FALSE)
#'
#' @noRd
hbr_plot_ld_heatmap <- function(
    gene_list,
    output_dir   = NULL,
    summary_csv  = NULL,
    per_gene_dir = NULL,
    result       = NULL,
    reference    = NULL,
    gff3_path    = NULL,
    geno_file    = NULL,
    geno_matrix  = NULL,
    geno_format  = NULL,
    plink_path   = "plink",
    trait_col    = "NT",
    n_gene       = NULL,
    sort_by      = c("Cohens_d", "Superior_Mean"),
    stat         = c("r2", "dprime"),
    palette      = "haplobreed",
    theme        = "haplobreed",
    width        = 8,
    height       = 6,
    formats      = c("pdf", "png"),
    combined_pdf = TRUE,
    figures_dir  = "figures",
    verbose      = TRUE) {

  stat <- match.arg(stat)
  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)

  reference <- .hbr_resolve_reference(reference = reference, result = result,
                                       gff3_path = gff3_path,
                                       verbose = verbose)

  if (is.null(output_dir) && !is.null(summary_csv))
    output_dir <- if (dir.exists(summary_csv)) summary_csv
                  else dirname(summary_csv)
  if (is.null(per_gene_dir)) per_gene_dir <- output_dir
  if (is.null(per_gene_dir) || !dir.exists(per_gene_dir))
    stop("Cannot locate per-gene files. Supply `output_dir`.", call. = FALSE)

  genes <- .hbr_resolve_genes(
    gene_list = if (missing(gene_list)) NULL else gene_list,
    summary_csv = summary_csv, output_dir = output_dir,
    n_gene = n_gene, sort_by = sort_by, verbose = verbose)
  if (!length(genes)) stop("No genes to plot.", call. = FALSE)
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (verbose)
    message("Plotting LD heatmaps for ", length(genes), " gene(s).")

  fig_dir <- if (nzchar(figures_dir %||% ""))
    file.path(output_dir, figures_dir) else output_dir
  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- list(); files_out <- list(); missing <- character(0)
  ld_list <- list()

  for (g in genes) {
    if (verbose) message("\n[", g, "]")
    pieces <- tryCatch(
      .hbr_load_gene_pieces(g, per_gene_dir,
                            genotype_m  = geno_matrix,
                            geno_file   = geno_file,
                            geno_format = geno_format,
                            plink_path  = plink_path,
                            trait_col   = trait_col,
                            verbose     = verbose),
      error = function(e) NULL)
    if (is.null(pieces) || is.null(pieces$geno) ||
        ncol(pieces$geno) < 2) {
      missing <- c(missing, g); next
    }
    ld_mat <- tryCatch(.hbr_compute_ld_matrix(pieces$geno, stat = stat),
                       error = function(e) NULL)
    if (is.null(ld_mat)) { missing <- c(missing, g); next }
    ld_list[[g]] <- ld_mat

    fig <- tryCatch(
      .hbr_build_ld_figure(g, ld_mat, pieces$func_snps, stat,
                           palette   = palette,
                           theme     = theme,
                           reference = reference),
      error = function(e) NULL)
    if (is.null(fig)) { missing <- c(missing, g); next }

    files <- .hbr_save_fig(fig, g, "ld_heatmap", fig_dir,
                           width, height, formats)
    plots[[g]] <- fig; files_out[[g]] <- files
  }

  combined_path <- .hbr_write_combined(plots, fig_dir, combined_pdf,
                                       width, height, verbose)

  if (verbose) {
    message("\n", strrep("=", 60))
    message("Genes plotted   : ", length(plots))
    message("Genes missing   : ", length(missing))
    message("Figures dir     : ", fig_dir)
    if (!is.null(combined_path))
      message("Combined PDF    : ", combined_path)
    message(strrep("=", 60))
  }

  invisible(list(plots = plots, files = files_out,
                 missing_genes = missing, combined_pdf = combined_path,
                 ld_matrix = ld_list))
}


# ----------------------------------------------------------------------------
# Superior-haplotype narrative (hbr_plot_superior_haplotype)
# ----------------------------------------------------------------------------

#' Plot a superior GWAS signal as a narrative multi-panel figure
#'
#' Assembles a full-page figure combining: a Manhattan plot with Q-Q
#' inset, three progressive zoom tiers to the candidate locus, a gene
#' model, and an LD heatmap.
#'
#' @param gene_id Top superior gene to highlight.
#' @param output_dir Pipeline output directory.
#' @param gwas_data Full GWAS SNP table.
#' @param result Optional in-memory [haplobreed_analysis()] output.
#' @param reference Optional `hbr_reference` object.
#' @param gff3_path Optional GFF3 path.
#' @param geno_file Genotype file for the LD heatmap.
#' @param geno_format Format override for `geno_file`.
#' @param plink_path PLINK executable path.
#' @param trait_col Trait column name.
#' @param suggestive_threshold Suggestive significance line.
#' @param window_bp Regional window for the tier-2 zoom.
#' @param zoom_gap_label Optional text on the trapezoid connectors.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param subtitle Optional subtitle.
#' @param output_file Output file path.
#' @param formats Output formats.
#' @param width,height Figure dimensions in inches.
#' @param verbose Print progress.
#' @return Invisibly, a list with `plot`, `files`, `lead`, `ld_to_lead`,
#'   `suggestive_threshold`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' gwas <- data.frame(SNP = c("snp005", "snp020", "snp045"),
#'                    Chromosome = "1",
#'                    Position = c(1800, 4800, 9800),
#'                    P = c(1e-8, 0.01, 0.5))
#' hbr_plot_superior_haplotype(
#'   gene_id    = "gene1",
#'   output_dir = toy,
#'   gwas_data  = gwas,
#'   trait_col  = "trait_value",
#'   formats    = "pdf",
#'   verbose    = FALSE)
#'
#' @noRd
hbr_plot_superior_haplotype <- function(
    gene_id,
    output_dir,
    gwas_data,
    result        = NULL,
    reference     = NULL,
    gff3_path     = NULL,
    geno_file     = NULL,
    geno_format   = NULL,
    plink_path    = "plink",
    trait_col     = "NT",
    suggestive_threshold = 1e-5,
    window_bp     = 200000,
    zoom_gap_label = NULL,
    palette       = "haplobreed",
    theme         = "haplobreed",
    subtitle      = NULL,
    output_file   = NULL,
    formats       = c("pdf", "png"),
    width         = 13,
    height        = 20,
    verbose       = TRUE) {

  if (!palette %in% hbr_palette_names())
    stop("Unknown palette '", palette, "'.", call. = FALSE)
  if (!theme %in% hbr_theme_names())
    stop("Unknown theme '", theme, "'.", call. = FALSE)

  th     <- hbr_theme(theme, base_size = 10)
  pal    <- .hbr_palette_anchors(palette)
  srf    <- .hbr_theme_surface(theme)
  gwas_l <- hbr_semantic("gwas_lead",   palette, theme)

  reference <- .hbr_resolve_reference(reference = reference, result = result,
                                       gff3_path = gff3_path,
                                       verbose = verbose)

  if (missing(gene_id) || !length(gene_id))
    stop("`gene_id` required.", call. = FALSE)
  gene_id <- as.character(gene_id)[1]

  if (!dir.exists(output_dir))
    stop("`output_dir` does not exist: ", output_dir, call. = FALSE)
  if (missing(gwas_data) || is.null(gwas_data))
    stop("`gwas_data` required.", call. = FALSE)

  gwas <- .hbr_normalize_gwas_data(gwas_data)
  if (is.null(gwas) || !nrow(gwas))
    stop("Could not extract POS/P columns from `gwas_data`.", call. = FALSE)

  chr_col <- intersect(c("Chromosome", "CHR", "Chr", "chr"),
                       names(gwas_data))[1]
  gwas$Chr <- if (!is.na(chr_col)) as.character(gwas_data[[chr_col]])
              else NA_character_

  pieces <- .hbr_load_gene_pieces(
    gene_id = gene_id, per_gene_dir = output_dir,
    summary_csv = file.path(output_dir, "publication_table_functional.csv"),
    geno_file = geno_file, geno_format = geno_format,
    plink_path = plink_path, trait_col = trait_col, verbose = verbose)

  func_snps <- pieces$func_snps
  pos <- as.numeric(func_snps$POS); pos <- pos[!is.na(pos)]
  if (!length(pos))
    stop("No functional SNP positions for ", gene_id, call. = FALSE)
  g_start <- min(pos); g_end <- max(pos)

  snp_chr <- if ("Chr" %in% names(func_snps))
    as.character(func_snps$Chr[1]) else NA_character_
  if (is.na(snp_chr) && !all(is.na(gwas$Chr))) {
    mid <- (g_start + g_end) / 2
    nearest <- gwas[which.min(abs(gwas$POS - mid))]
    snp_chr <- nearest$Chr[1]
  }

  in_signal_window <- function(g) {
    if (is.na(snp_chr)) return(rep(FALSE, nrow(g)))
    as.character(g$Chr) == snp_chr &
      g$POS >= g_start - 50000 & g$POS <= g_end + 50000
  }

  win_gd <- gwas[in_signal_window(gwas) & !is.na(gwas$P) & gwas$P > 0]
  lead_info <- if (nrow(win_gd)) {
    i <- which.min(win_gd$P)
    list(pos     = win_gd$POS[i],
         snp     = win_gd$SNP[i],
         p       = win_gd$P[i],
         neglogP = -log10(win_gd$P[i]),
         label   = sprintf("P = %s", hbr_format_p(win_gd$P[i])))
  } else NULL

  ld_to_lead <- NULL
  if (!is.null(pieces$geno) && !is.null(lead_info) &&
      lead_info$snp %in% colnames(pieces$geno)) {
    ld_to_lead <- .hbr_compute_ld_to_lead(pieces$geno, lead_info$snp)
  }

  gwas_plot <- data.table::copy(gwas)
  if (all(is.na(gwas_plot$Chr))) gwas_plot[, Chr := "1"]
  gwas_plot[, Chr := as.character(Chr)]
  chr_levels <- unique(gwas_plot$Chr)
  chr_num <- suppressWarnings(as.numeric(chr_levels))
  chr_levels <- if (all(!is.na(chr_num))) chr_levels[order(chr_num)]
  else sort(chr_levels)
  gwas_plot[, Chr := factor(Chr, levels = chr_levels)]
  gwas_plot <- .hbr_build_cumulative_positions(gwas_plot)
  chr_mid <- gwas_plot[, .(mid = mean(range(GenomicPos))), by = Chr]

  chrom_alt <- hbr_semantic("chrom_alt", palette, theme)
  chr_fill  <- rep(chrom_alt, length.out = length(chr_levels))
  gwas_plot[, fill_col := chr_fill[as.integer(Chr)]]
  gwas_plot[, is_signal := in_signal_window(gwas_plot)]
  gwas_plot[is_signal == TRUE, fill_col := gwas_l]

  bg_pts <- gwas_plot[is_signal == FALSE]
  fg_pts <- gwas_plot[is_signal == TRUE]

  bonf_p <- 0.05 / nrow(gwas_plot)
  sugg_p <- suggestive_threshold

  p_manhattan <- ggplot2::ggplot() +
    ggplot2::geom_point(data = bg_pts,
                        ggplot2::aes(x = GenomicPos, y = -log10(P)),
                        colour = bg_pts$fill_col,
                        size = 0.9, alpha = 0.75) +
    ggplot2::geom_point(data = fg_pts,
                        ggplot2::aes(x = GenomicPos, y = -log10(P)),
                        colour = gwas_l,
                        size = 3.2, alpha = 0.95) +
    ggplot2::geom_hline(yintercept = -log10(bonf_p), linetype = "dashed",
                        colour = srf$text_muted, linewidth = 0.4) +
    ggplot2::annotate("text", x = Inf, y = -log10(bonf_p),
                      label = sprintf("P = %s", hbr_format_p(bonf_p)),
                      hjust = 1.05, vjust = -0.35,
                      size = 2.6, colour = srf$text_muted)

  if (!is.null(sugg_p) && is.finite(sugg_p)) {
    p_manhattan <- p_manhattan +
      ggplot2::geom_hline(yintercept = -log10(sugg_p),
                          linetype = "dotted",
                          colour = srf$text_muted,
                          linewidth = 0.4) +
      ggplot2::annotate("text", x = Inf, y = -log10(sugg_p),
                        label = sprintf("suggestive P = %s",
                                        hbr_format_p(sugg_p)),
                        hjust = 1.05, vjust = -0.35,
                        size = 2.6, colour = srf$text_muted)
  }

  p_manhattan <- p_manhattan +
    ggplot2::scale_x_continuous(breaks = chr_mid$mid,
                                labels = as.character(chr_mid$Chr),
                                expand = c(0.01, 0.01)) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](P))) +
    th +
    ggplot2::theme(
      plot.margin        = ggplot2::margin(4, 6, 4, 4),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(size = 8))

  p_manhattan <- .hbr_letter_panel(p_manhattan, "A.)",
                                    "Genome-wide association",
                                    colour = srf$text_primary)

  pvals <- gwas$P[!is.na(gwas$P) & gwas$P > 0 & gwas$P <= 1]
  n_q   <- length(pvals)
  if (n_q >= 2) {
    obs_v  <- sort(-log10(pvals))
    exp_v  <- -log10(stats::ppoints(n_q))
    lambda <- stats::median(
      stats::qchisq(pvals, df = 1, lower.tail = FALSE)) /
      stats::qchisq(0.5, df = 1)
    lim_hi <- max(c(obs_v, exp_v), na.rm = TRUE) * 1.05
    if (!is.finite(lim_hi) || lim_hi <= 0) lim_hi <- 1

    sig_mask  <- in_signal_window(gwas) &
      !is.na(gwas$P) & gwas$P > 0 & gwas$P <= 1
    sig_pvals <- gwas$P[sig_mask]

    qq_sig_df <- NULL
    if (length(sig_pvals)) {
      sig_ranks <- vapply(sig_pvals, function(p)
        max(1L, sum(pvals <= p)), integer(1))
      sig_ranks <- pmin(sig_ranks, n_q)
      qq_sig_df <- data.frame(
        Expected = -log10(stats::ppoints(n_q)[sig_ranks]),
        Observed = -log10(sig_pvals))
    }

    p_qq <- ggplot2::ggplot(
      data.frame(Expected = exp_v, Observed = obs_v),
      ggplot2::aes(x = Expected, y = Observed)) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           colour = srf$text_secondary, linewidth = 0.5) +
      ggplot2::geom_line(
        colour = hbr_semantic("qq_line", palette, theme), linewidth = 0.9)

    if (!is.null(qq_sig_df) && nrow(qq_sig_df)) {
      p_qq <- p_qq + ggplot2::geom_point(
        data    = qq_sig_df,
        mapping = ggplot2::aes(x = Expected, y = Observed),
        colour  = gwas_l, size = 2.4, alpha = 0.95,
        inherit.aes = FALSE)
    }

    p_qq <- p_qq +
      ggplot2::annotate(
        "text",
        x = 0.02 * lim_hi, y = 0.98 * lim_hi,
        label = sprintf("lambda == %.3f ~ '|' ~ n == '%s'",
                        lambda, format(n_q, big.mark = ",")),
        parse = TRUE,
        hjust = 0, vjust = 1,
        size = 2.8, colour = srf$text_muted) +
      ggplot2::labs(x = expression(-log[10](P)[expected]),
                    y = expression(-log[10](P)[observed])) +
      ggplot2::coord_cartesian(xlim = c(0, lim_hi), ylim = c(0, lim_hi)) +
      th +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 8),
        axis.text  = ggplot2::element_text(size = 8))
    p_qq <- .hbr_letter_panel(p_qq, "", NULL,
                               colour = srf$text_primary)
  } else {
    p_qq <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0,
                        label = "Q-Q unavailable\n(no P values)",
                        size = 3.5, colour = srf$text_muted)
  }

  panel_a <- if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::patchworkGrob(
      patchwork::wrap_plots(p_manhattan, p_qq,
                            ncol = 2, widths = c(1.6, 1.2)))
  } else {
    gridExtra::arrangeGrob(p_manhattan, p_qq, ncol = 2,
                            widths = c(1.6, 1.2))
  }

  chr_pos <- gwas[!is.na(gwas$Chr) & gwas$Chr == snp_chr, ]
  if (!nrow(chr_pos)) chr_pos <- gwas

  tier1 <- .hbr_build_zoom_tier(
    gwas_norm = chr_pos,
    xlim      = range(chr_pos$POS, na.rm = TRUE),
    chr       = snp_chr,
    title     = "Chromosome view",
    box_range = c(g_start - window_bp, g_end + window_bp),
    lead      = lead_info,
    palette   = palette, theme = theme)

  tier2 <- .hbr_build_zoom_tier(
    gwas_norm       = gwas,
    xlim            = c(max(1, g_start - window_bp), g_end + window_bp),
    chr             = snp_chr,
    title           = "Regional view",
    highlight_range = c(g_start, g_end),
    box_range       = c(g_start - 5000, g_end + 5000),
    ld_to_lead      = ld_to_lead,
    lead            = lead_info,
    palette         = palette, theme = theme)

  stats_box <- NULL
  if (!is.null(lead_info)) {
    stats_box <- sprintf("Lead: %s\nP = %s",
                         lead_info$snp,
                         hbr_format_p(lead_info$p))
  }

  tier3 <- .hbr_build_zoom_tier(
    gwas_norm       = gwas,
    xlim            = c(g_start - 5000, g_end + 5000),
    chr             = snp_chr,
    title           = "Gene view",
    highlight_range = c(g_start, g_end),
    ld_to_lead      = ld_to_lead,
    lead            = lead_info,
    stats_box       = stats_box,
    palette         = palette, theme = theme)

  gap_tier <- .hbr_zoom_gap_grob(label = zoom_gap_label, colour = gwas_l)

  tier1_letter <- .hbr_letter_panel(tier1, "B.)",
                                     "Progressive zoom to the candidate locus",
                                     colour = srf$text_primary)

  p_struct <- NULL
  if (!is.null(reference)) {
    gene_gff_id <- .hbr_resolve_gene_gff_id(gene_id, func_snps, reference)
    if (!is.null(gene_gff_id)) {
      p_struct <- tryCatch(
        .hbr_gs_build_one(
          gene_id = gene_id, gene_gff_id = gene_gff_id,
          reference = reference, variants = func_snps,
          show_variants = TRUE,
          variant_shape_by_region = TRUE,
          aa_full_names = TRUE,
          palette = palette, theme = theme),
        error = function(e) NULL)
      if (!is.null(p_struct)) {
        p_struct <- .hbr_letter_panel(p_struct, "C.)",
                                       "Gene structure",
                                       colour = srf$text_primary) +
          ggplot2::labs(subtitle = NULL) +
          ggplot2::theme(
            plot.subtitle = ggplot2::element_blank(),
            plot.margin   = ggplot2::margin(4, 6, 4, 6))
      }
    }
  }

  p_ld <- NULL
  if (!is.null(pieces$geno) && ncol(pieces$geno) >= 2) {
    ld_mat <- .hbr_compute_ld_matrix(pieces$geno, stat = "r2")
    n_ld   <- nrow(ld_mat)
    cells <- do.call(rbind, lapply(seq_len(n_ld), function(i) {
      js <- seq.int(i, n_ld)
      data.frame(row = i, x = js - (i - 1) / 2,
                 value = ld_mat[cbind(i, js)])
    }))
    apex_x <- (n_ld + 1) / 2
    outline_df <- data.frame(
      x = c(0.5, n_ld + 0.5, apex_x + 0.5, apex_x - 0.5, 0.5),
      y = c(0.5, 0.5, n_ld + 0.5, n_ld + 0.5, 0.5))
    ld_pal <- hbr_semantic("ld_gradient", palette, theme)
    snp_id <- paste0("S", seq_len(n_ld))

    p_heat <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = cells,
                         ggplot2::aes(x = x, y = row, fill = value),
                         width = 1, height = 1,
                         colour = srf$panel_bg, linewidth = 0.35) +
      ggplot2::geom_path(data = outline_df,
                         ggplot2::aes(x = x, y = y),
                         colour = srf$text_primary, linewidth = 0.55,
                         inherit.aes = FALSE) +
      ggplot2::scale_fill_gradientn(
        colours = ld_pal, limits = c(0, 1),
        na.value = srf$grid_major, name = expression(r^2)) +
      ggplot2::scale_x_continuous(breaks = seq_len(n_ld), labels = snp_id,
                                  limits = c(0.5, n_ld + 0.5),
                                  expand = c(0, 0)) +
      ggplot2::scale_y_reverse(breaks = seq_len(n_ld), labels = snp_id,
                               limits = c(n_ld + 0.5, 0.5),
                               expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL) +
      th +
      ggplot2::theme(
        plot.margin       = ggplot2::margin(4, 8, 2, 8),
        panel.background  = ggplot2::element_rect(fill = srf$table_bg,
                                                  colour = NA),
        panel.grid        = ggplot2::element_blank(),
        axis.ticks        = ggplot2::element_blank(),
        axis.text         = ggplot2::element_text(size = 8,
                                                  colour = srf$text_primary),
        legend.key.width  = grid::unit(0.6, "lines"),
        legend.key.height = grid::unit(1.0, "lines"))

    p_heat <- .hbr_letter_panel(p_heat, "D.)", "LD heatmap (r\u00b2)",
                                 colour = srf$text_primary)

    p_cap <- .hbr_snp_caption_strip(pieces$func_snps,
                                     palette = palette, theme = theme)

    p_ld <- gridExtra::arrangeGrob(
      ggplot2::ggplotGrob(p_heat),
      ggplot2::ggplotGrob(p_cap),
      ncol = 1, heights = c(0.82, 0.18))
  }

  panels  <- list(panel_a, tier1_letter, gap_tier, tier2, gap_tier, tier3)
  heights <- c(5.0,     2.1,   0.25,     2.1,   0.25,     2.1)

  if (!is.null(p_struct)) {
    gap_struct <- .hbr_zoom_gap_grob(label = NULL, colour = gwas_l)
    panels  <- c(panels,  list(gap_struct, p_struct))
    heights <- c(heights, 0.20,        1.9)
  }
  if (!is.null(p_ld)) {
    gap_ld  <- .hbr_zoom_gap_grob(label = NULL, colour = gwas_l)
    panels  <- c(panels,  list(gap_ld, p_ld))
    heights <- c(heights, 0.25,     6.0)
  }

  fig <- do.call(gridExtra::arrangeGrob,
                 c(panels, list(ncol = 1, heights = heights)))

  fig <- gtable::gtable_add_rows(fig, grid::unit(2.5, "lines"), 0)
  top_title <- if (is.null(subtitle) || !nzchar(subtitle))
    paste0(gene_id, " \u2014 candidate haplotype locus")
  else paste0(gene_id, " \u2014 candidate haplotype locus\n", subtitle)

  fig <- gtable::gtable_add_grob(
    fig,
    grid::textGrob(top_title,
                   gp = grid::gpar(fontface = "bold", fontsize = 13,
                                   col = srf$text_primary),
                   x = grid::unit(10, "pt"), hjust = 0, vjust = 1),
    1, 1, 1, ncol(fig))

  figures_out <- .hbr_resolve_figures_dir(output_dir, NULL, "figures")
  files <- .hbr_save_fig(fig, gene_id, "candidate_locus",
                         figures_out, width, height, formats)

  if (verbose)
    message("\nCandidate-locus figure written: ",
            paste(unlist(files), collapse = ", "))

  invisible(list(plot = fig, files = files,
                 lead = lead_info, ld_to_lead = ld_to_lead,
                 suggestive_threshold = suggestive_threshold))
}

#' Plot cross ranking by usefulness vs mid-parent index
#'
#' Visualises how the ranking of recommended crosses changes when progeny
#' variance is accounted for. Each point is one cross; the x-axis is the
#' rank by mid-parent index alone, the y-axis is the rank by usefulness
#' (Schnell and Utz, 1975). Points on the diagonal are crosses whose
#' recommendation is identical under the two criteria; points off the
#' diagonal are crosses whose rank depends on progeny variance. The most
#' displaced cross is optionally highlighted and labelled.
#'
#' Use this to document that the genetic map is doing real work in the
#' cross-design step, or to justify why a cross with a lower mid-parent
#' index is nonetheless recommended.
#'
#' @param x Either an `hgsi_advanced` object (from
#'   [hbr_advanced_decision()]) or a data.frame with columns `parent_a`,
#'   `parent_b`, `mid_index`, and `uc`.
#' @param top_n Number of top crosses to display. Crosses are ordered by
#'   `uc`. Default `20L`.
#' @param highlight A cross to annotate. Either a length-2 character
#'   vector `c(parent_a, parent_b)`, or `"auto"` to highlight the cross
#'   with the largest absolute rank displacement (default).
#' @param label_size Text size for the highlight label.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param output_dir Directory for the figure. `NULL` returns the plot
#'   without writing anything.
#' @param file_name Base name for the output file (without extension).
#' @param width,height Figure dimensions in inches.
#' @param formats Character vector of output formats. Any of `"pdf"`,
#'   `"png"`, `"tiff"`, `"svg"`.
#' @param verbose Print progress.
#'
#' @return Invisibly, a list with `plot` (the ggplot object) and `files`
#'   (the paths written, or `character(0)` if `output_dir` is `NULL`).
#'   The `plot` element can be assigned and modified further with
#'   ggplot2.
#'
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' crosses <- data.frame(
#'   parent_a  = sprintf("IRIS_%04d", 1:20),
#'   parent_b  = sprintf("IRIS_%04d", 21:40),
#'   mid_index = runif(20, 5, 12),
#'   uc        = runif(20, 5, 15))
#' hbr_plot_cross_ranking(crosses, top_n = 20)
#'
#' @export
hbr_plot_cross_ranking <- function(
    x,
    top_n         = 20L,
    highlight     = c("auto"),
    label_size    = 3,
    palette       = "haplobreed",
    theme         = "haplobreed",
    output_dir    = NULL,
    file_name     = "cross_ranking",
    width         = 6,
    height        = 6,
    formats       = c("pdf", "png"),
    verbose       = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  # ---- Resolve input -----------------------------------------------------
  if (inherits(x, "hgsi_advanced")) {
    cr <- data.table::as.data.table(x$crosses)
  } else if (is.data.frame(x)) {
    cr <- data.table::as.data.table(x)
  } else {
    stop("`x` must be an hgsi_advanced object or a data.frame of crosses.",
         call. = FALSE)
  }
  required <- c("parent_a", "parent_b", "mid_index", "uc")
  miss <- setdiff(required, names(cr))
  if (length(miss))
    stop("`x` is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (nrow(cr) < 2L)
    stop("`x` must contain at least two crosses.", call. = FALSE)
  
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 2L)
    stop("`top_n` must be an integer >= 2.", call. = FALSE)
  top_n <- min(top_n, nrow(cr))
  
  # ---- Rank -------------------------------------------------------------
  cr <- cr[!is.na(mid_index) & !is.na(uc)]
  if (nrow(cr) < 2L)
    stop("Not enough non-missing crosses to plot.", call. = FALSE)
  
  top <- cr[order(-uc)][seq_len(top_n)]
  top[, rank_uc  := data.table::frank(-uc)]
  top[, rank_mid := data.table::frank(-mid_index)]
  top[, shift    := rank_uc - rank_mid]           # signed
  top[, cross_id := paste(parent_a, parent_b, sep = " \u00d7 ")]
  
  # ---- Highlight --------------------------------------------------------
  if (identical(highlight, "auto")) {
    hi <- top[which.max(abs(shift))]
  } else if (is.character(highlight) && length(highlight) == 2L) {
    hi <- top[parent_a == highlight[1] & parent_b == highlight[2]]
    if (!nrow(hi))
      hi <- top[parent_a == highlight[2] & parent_b == highlight[1]]
    if (!nrow(hi)) {
      warning("Highlighted cross not found in the top ", top_n,
              " crosses; skipping highlight.", call. = FALSE)
      hi <- top[0]
    }
  } else {
    hi <- top[0]
  }
  
  # ---- Summary ----------------------------------------------------------
  n_shifted <- sum(top$shift != 0L)
  max_shift <- max(abs(top$shift))
  if (verbose)
    message(sprintf(
      "  %d of %d top crosses reordered (max |shift| = %d)",
      n_shifted, top_n, max_shift))
  
  # ---- Colours and theme -----------------------------------------------
  pal    <- .hbr_palette_anchors(palette)
  srf    <- .hbr_theme_surface(theme)
  th     <- hbr_theme(theme, base_size = 11)
  
  accent      <- pal$superior
  accent_dark <- pal$superior_dark
  muted       <- srf$text_muted
  txt_primary <- srf$text_primary
  txt_second  <- srf$text_secondary
  
  # ---- Plot -------------------------------------------------------------
  p <- ggplot2::ggplot(top,
                       ggplot2::aes(x = rank_mid, y = rank_uc)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = muted,
                         linewidth = 0.4) +
    ggplot2::geom_point(size = 2.6, colour = accent, alpha = 0.85) +
    ggplot2::geom_point(data = hi,
                        ggplot2::aes(x = rank_mid, y = rank_uc),
                        inherit.aes = FALSE,
                        shape = 21, size = 4.5, stroke = 0.9,
                        fill = NA, colour = accent_dark) +
    ggplot2::scale_x_continuous(
      breaks = pretty(c(1, top_n), n = 5),
      limits = c(0.5, top_n + 0.5),
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = pretty(c(1, top_n), n = 5),
      limits = c(0.5, top_n + 0.5),
      expand = c(0, 0)) +
    ggplot2::labs(
      x = "Rank by mid-parent index",
      y = "Rank by usefulness (UC)",
      title    = "Cross ranking under two criteria",
      subtitle = sprintf(
        "%d of %d top crosses reorder; max shift %d positions",
        n_shifted, top_n, max_shift)) +
    th +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(colour = srf$grid_major,
                                                 linewidth = 0.3),
      plot.title         = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = txt_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(
        size = 9, hjust = 0, colour = txt_second,
        margin = ggplot2::margin(b = 6)),
      axis.title         = ggplot2::element_text(
        face = "bold", size = 10, colour = txt_primary),
      axis.text          = ggplot2::element_text(
        size = 9, colour = txt_primary),
      plot.margin        = ggplot2::margin(8, 12, 8, 8))
  
  if (nrow(hi)) {
    lbl <- sprintf("%s\n(mid rank %d, UC rank %d)",
                   hi$cross_id, hi$rank_mid, hi$rank_uc)
    p <- p + ggplot2::annotate(
      "text",
      x = hi$rank_mid, y = hi$rank_uc,
      label = lbl,
      hjust = -0.05, vjust = -0.6,
      size = label_size,
      colour = accent_dark,
      fontface = "italic",
      lineheight = 0.95)
  }
  
  # ---- Save -------------------------------------------------------------
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files <- .hbr_save_fig(p, file_name, "cross_ranking",
                           output_dir, width, height, formats)
    if (verbose && length(files))
      message("  Figure(s) written: ",
              paste(basename(files), collapse = ", "))
  }
  
  invisible(list(plot = p, files = files))
}

#' Forest plot of cross-validation comparisons
#'
#' Displays the paired CV comparison of M2 vs M0, M3 vs M0 and M3 vs M2
#' per trait, with bootstrap confidence intervals. Points whose interval
#' excludes zero are coloured by the palette's superior colour; the rest
#' are drawn in the muted text colour. Use this as the headline figure
#' for the honest-vs-naive comparison.
#'
#' @param x An `hgsi_advanced` object or the `cv_compare` table.
#' @param metric Metric to display. `"r"` by default.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param output_dir Directory for the figure. `NULL` returns the plot.
#' @param file_name Base file name (without extension).
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param verbose Print progress.
#'
#' @return Invisibly, `list(plot = <ggplot>, files = <character>)`.
#' @examples
#' cmp <- data.frame(
#'   trait     = rep(c("YLD", "DTF"), each = 3),
#'   method    = rep(c("M2", "M3", "M3"), 2),
#'   baseline  = rep(c("M0", "M0", "M2"), 2),
#'   metric    = "r",
#'   mean_diff = c(0.05, 0.02, 0.03, 0.04, 0.01, 0.03),
#'   boot_lo   = c(0.01, -0.01, 0.00, 0.01, -0.02, 0.00),
#'   boot_hi   = c(0.09, 0.05, 0.06, 0.07, 0.04, 0.06),
#'   stringsAsFactors = FALSE)
#' hbr_plot_cv_forest(cmp)
#'
#' @export
hbr_plot_cv_forest <- function(
    x,
    metric      = "r",
    palette     = "haplobreed",
    theme       = "haplobreed",
    output_dir  = NULL,
    file_name   = "cv_forest",
    width       = 8,
    height      = NULL,
    formats     = c("pdf", "png"),
    verbose     = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  if (inherits(x, "hgsi_advanced")) {
    cmp <- data.table::as.data.table(x$cv_compare)
  } else if (is.data.frame(x)) {
    cmp <- data.table::as.data.table(x)
  } else stop("`x` must be an hgsi_advanced object or a data.frame.",
              call. = FALSE)
  
  required <- c("trait", "method", "baseline",
                "mean_diff", "boot_lo", "boot_hi")
  miss <- setdiff(required, names(cmp))
  if (length(miss)) stop("`x` is missing column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  
  if ("metric" %in% names(cmp) && !is.null(metric))
    cmp <- cmp[cmp$metric == metric]
  if (!nrow(cmp)) stop("No rows for metric '", metric, "'.", call. = FALSE)
  
  cmp[, cmp_order := data.table::fcase(
    method == "M2" & baseline == "M0", 1L,
    method == "M3" & baseline == "M0", 2L,
    method == "M3" & baseline == "M2", 3L,
    default = 99L)]
  cmp[, comparison := paste0(method, " \u2212 ", baseline)]
  data.table::setorder(cmp, trait, cmp_order)
  cmp[, row_label := sprintf("%-8s  %s", trait, comparison)]
  cmp[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  cmp[, sig := sign(boot_lo) == sign(boot_hi)]
  
  pal  <- .hbr_palette_anchors(palette)
  srf  <- .hbr_theme_surface(theme)
  th   <- hbr_theme(theme, base_size = 11)
  
  cmp[, col := ifelse(sig, pal$superior, srf$text_muted)]
  
  p <- ggplot2::ggplot(cmp,
                       ggplot2::aes(x = mean_diff, y = row_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = srf$text_muted, linewidth = 0.4) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = boot_lo, xmax = boot_hi, colour = col),
      height = 0.18, linewidth = 0.65) +
    ggplot2::geom_point(ggplot2::aes(colour = col), size = 2.7) +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(
      x = sprintf("\u0394%s  (95%% bootstrap CI)", metric),
      y = NULL,
      title    = "Cross-validation: paired comparisons",
      subtitle = sprintf(
        "%d of %d comparisons have CIs excluding zero",
        sum(cmp$sig), nrow(cmp))) +
    th +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = srf$text_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(
        size = 9.5, hjust = 0, colour = srf$text_secondary,
        margin = ggplot2::margin(b = 6)),
      axis.text.y        = ggplot2::element_text(
        family = "mono", size = 9, colour = srf$text_primary),
      axis.title.x       = ggplot2::element_text(
        face = "bold", size = 10, colour = srf$text_primary),
      plot.margin        = ggplot2::margin(8, 14, 8, 6))
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    h <- height %||% max(3, 0.34 * nrow(cmp) + 1.4)
    files <- .hbr_save_fig(p, file_name, "cv_forest",
                           output_dir, width, h, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = p, files = files))
}

#' Per-trait prediction quality
#'
#' For each trait, draws a scatter of observed vs predicted values. When
#' `haplotype` is supplied for a trait, points are coloured by carrier
#' status of the superior haplotype at the chosen gene (two colours); the
#' panel subtitle reports the carrier mean difference. Otherwise points
#' are coloured by individual reliability and a reliability histogram is
#' drawn underneath.
#'
#' @param x An `hgsi_advanced` object.
#' @param haplotype Optional named list. Each element is a list with
#'   `gene`, optional `trait_dir`, and optional `superior_haplotype`,
#'   used to colour the trait's panel by carrier status. Traits not
#'   listed fall back to the reliability colouring.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param point_size,point_alpha Scatter styling.
#' @param output_dir Directory for the figure.
#' @param file_name Base file name.
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param verbose Print progress.
#'
#' @return Invisibly, `list(plot = <grob>, files = <character>)`.
#'
#' @noRd
hbr_plot_prediction_quality <- function(
    x,
    haplotype    = NULL,
    palette      = "haplobreed",
    theme        = "haplobreed",
    point_size   = 1.8,
    point_alpha  = 0.85,
    output_dir   = NULL,
    file_name    = "prediction_quality",
    width        = NULL,
    height       = 5.5,
    formats      = c("pdf", "png"),
    verbose      = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  if (!inherits(x, "hgsi_advanced"))
    stop("`x` must be an hgsi_advanced object.", call. = FALSE)
  
  traits <- x$traits
  Tn     <- length(traits)
  pal    <- .hbr_palette_anchors(palette)
  srf    <- .hbr_theme_surface(theme)
  th     <- hbr_theme(theme, base_size = 10)
  accent <- pal$superior
  muted  <- srf$text_muted
  
  # Resolve carrier vectors up front so a missing gene fails early.
  carriers <- list()
  if (!is.null(haplotype)) {
    for (tr in names(haplotype)) {
      if (!tr %in% traits) {
        warning("haplotype entry '", tr, "' is not a trait in x; skipping.",
                call. = FALSE)
        next
      }
      spec <- haplotype[[tr]]
      carriers[[tr]] <- hbr_gene_carriers(
        x, trait = tr, gene = spec$gene,
        trait_dir = spec$trait_dir %||% NULL,
        superior_haplotype = spec$superior_haplotype %||% NULL)
    }
  }
  
  panels_sc <- list(); panels_h <- list()
  for (tr in traits) {
    fit <- x$fits[[tr]]
    df  <- data.frame(obs  = x$pheno_mat[, tr],
                      pred = fit$pred,
                      rel  = fit$reliability)
    df  <- df[!is.na(df$obs) & !is.na(df$pred), , drop = FALSE]
    rv  <- if (nrow(df) >= 3) stats::cor(df$obs, df$pred) else NA_real_
    
    use_hap <- tr %in% names(carriers)
    
    if (use_hap) {
      # Align carrier status to df rows
      keep_idx <- !is.na(x$pheno_mat[, tr]) & !is.na(fit$pred)
      df$carrier <- carriers[[tr]]$carrier[keep_idx]
      sup_hap    <- carriers[[tr]]$superior_haplotype
      gene       <- carriers[[tr]]$gene
      n_car      <- sum(df$carrier, na.rm = TRUE)
      mean_car   <- mean(df$obs[df$carrier],     na.rm = TRUE)
      mean_non   <- mean(df$obs[!df$carrier],    na.rm = TRUE)
      diff_mean  <- mean_car - mean_non
      
      # Two colours: muted for non-carrier, superior for carrier
      fills <- c(`FALSE` = muted, `TRUE` = accent)
      
      p_sc <- ggplot2::ggplot(df,
                              ggplot2::aes(x = pred, y = obs,
                                           colour = carrier)) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             linetype = "dashed", colour = muted,
                             linewidth = 0.4) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::scale_colour_manual(
          values = fills, name = NULL,
          labels = c(`FALSE` = "non-carrier", `TRUE` = "carrier")) +
        ggplot2::labs(
          x = "Predicted", y = "Observed",
          title    = tr,
          subtitle = sprintf(
            "%s | %s n=%d | \u0394mean = %+.2f",
            gene, sup_hap, n_car, diff_mean)) +
        th +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(
            face = "bold", size = 10.5, hjust = 0,
            colour = srf$text_primary, margin = ggplot2::margin(b = 1)),
          plot.subtitle = ggplot2::element_text(
            size = 7.8, hjust = 0, colour = srf$text_secondary,
            margin = ggplot2::margin(b = 4)),
          legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = ggplot2::element_rect(
            fill = "white", colour = NA),
          legend.key = ggplot2::element_blank(),
          legend.text  = ggplot2::element_text(size = 8),
          legend.key.size = grid::unit(0.35, "cm"),
          plot.margin = ggplot2::margin(4, 6, 4, 4))
      
    } else {
      # Original reliability-coloured scatter
      p_sc <- ggplot2::ggplot(df,
                              ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             linetype = "dashed", colour = muted,
                             linewidth = 0.4) +
        ggplot2::geom_point(ggplot2::aes(colour = rel),
                            size = point_size, alpha = point_alpha) +
        ggplot2::scale_colour_gradientn(
          colours = c("#F0F0F0", accent),
          limits  = c(0, 1), name = "rel.") +
        ggplot2::labs(
          x = "Predicted", y = "Observed",
          title    = tr,
          subtitle = if (is.na(rv)) "" else sprintf("r = %.3f", rv)) +
        th +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(
            face = "bold", size = 10.5, hjust = 0,
            colour = srf$text_primary, margin = ggplot2::margin(b = 1)),
          plot.subtitle = ggplot2::element_text(
            size = 8.5, hjust = 0, colour = srf$text_secondary,
            margin = ggplot2::margin(b = 4)),
          legend.position = "right",
          legend.key.height = grid::unit(0.8, "cm"),
          legend.title  = ggplot2::element_text(size = 8),
          legend.text   = ggplot2::element_text(size = 7.5),
          plot.margin   = ggplot2::margin(4, 6, 4, 4))
    }
    
    p_h <- ggplot2::ggplot(df, ggplot2::aes(x = rel)) +
      ggplot2::geom_histogram(bins = 20, fill = accent,
                              colour = "white", alpha = 0.85,
                              linewidth = 0.2) +
      ggplot2::labs(x = "Reliability", y = NULL) +
      th +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = 8),
        axis.title = ggplot2::element_text(face = "bold", size = 9,
                                           colour = srf$text_primary),
        plot.margin = ggplot2::margin(2, 6, 4, 4))
    
    panels_sc[[tr]] <- p_sc
    panels_h[[tr]]  <- p_h
  }
  
  fig <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    row_sc <- patchwork::wrap_plots(panels_sc, nrow = 1)
    row_h  <- patchwork::wrap_plots(panels_h,  nrow = 1)
    fig    <- patchwork::wrap_plots(row_sc, row_h,
                                    ncol = 1, heights = c(1, 0.55))
  } else {
    fig <- gridExtra::arrangeGrob(
      do.call(gridExtra::arrangeGrob,
              c(lapply(panels_sc, ggplot2::ggplotGrob), list(nrow = 1))),
      do.call(gridExtra::arrangeGrob,
              c(lapply(panels_h,  ggplot2::ggplotGrob), list(nrow = 1))),
      ncol = 1, heights = c(1, 0.55))
  }
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    w <- width %||% max(4, Tn * 3.2)
    files <- .hbr_save_fig(fig, file_name, "prediction_quality",
                           output_dir, w, height, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = fig, files = files))
}
#' Selected-set coverage against the full panel
#'
#' For each trait and for the multi-trait index, draws box plots of the
#' genetic values (or index) split by membership in the simulated-
#' annealing-selected parent set. Shows whether the selected parents are
#' at the top of each trait and of the index, and by how much the index
#' separates them from the rest of the panel.
#'
#' @param x An `hgsi_advanced` object.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param show_points Overlay individual values as jittered points.
#' @param point_alpha Alpha for the jittered points.
#' @param output_dir Directory for the figure. `NULL` returns the plot.
#' @param file_name Base file name.
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param verbose Print progress.
#'
#' @return Invisibly, `list(plot = <ggplot>, files = <character>)`.
#'
#' @noRd
hbr_plot_selected_set_coverage <- function(
    x,
    palette       = "haplobreed",
    theme         = "haplobreed",
    show_points   = TRUE,
    point_alpha   = 0.35,
    output_dir    = NULL,
    file_name     = "selected_set_coverage",
    width         = 8,
    height        = 4,
    formats       = c("pdf", "png"),
    verbose       = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  if (!inherits(x, "hgsi_advanced"))
    stop("`x` must be an hgsi_advanced object.", call. = FALSE)
  if (is.null(x$selected_set))
    stop("`x` has no selected_set; rerun with run_set_selection = TRUE.",
         call. = FALSE)
  
  sel_ids <- x$selected_set$ids
  if (!length(sel_ids))
    stop("The selected set is empty.", call. = FALSE)
  
  # Build long table: one row per (panel, individual)
  parts <- list()
  for (tr in x$traits) {
    parts[[length(parts) + 1L]] <- data.table::data.table(
      panel    = tr,
      id       = x$ids,
      value    = x$fits[[tr]]$g,
      selected = x$ids %in% sel_ids)
  }
  parts[[length(parts) + 1L]] <- data.table::data.table(
    panel    = "Index",
    id       = x$ids,
    value    = x$index,
    selected = x$ids %in% sel_ids)
  dt <- data.table::rbindlist(parts, use.names = TRUE)
  
  # Panel order: Index first, then traits
  dt[, panel := factor(panel, levels = c("Index", x$traits))]
  
  pal <- .hbr_palette_anchors(palette)
  srf <- .hbr_theme_surface(theme)
  th  <- hbr_theme(theme, base_size = 11)
  fills <- c(`FALSE` = srf$text_muted, `TRUE` = pal$superior)
  
  p <- ggplot2::ggplot(dt,
                       ggplot2::aes(x = selected, y = value,
                                    fill = selected)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA,
                          colour = srf$text_primary,
                          linewidth = 0.5, alpha = 0.75)
  
  if (isTRUE(show_points))
    p <- p + ggplot2::geom_jitter(width = 0.14, size = 0.9,
                                  alpha = point_alpha,
                                  colour = srf$text_primary)
  
  p <- p +
    ggplot2::facet_wrap(~ panel, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_x_discrete(labels = c(`FALSE` = "Panel",
                                         `TRUE`  = "Selected")) +
    ggplot2::labs(
      x = NULL, y = "Genetic value",
      title    = "Selected parents vs the panel",
      subtitle = sprintf(
        "%d selected parents; mean relatedness %.3f; coverage %.2f",
        length(sel_ids), x$selected_set$mean_relatedness,
        x$selected_set$coverage)) +
    th +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = srf$text_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(
        size = 9, hjust = 0, colour = srf$text_secondary,
        margin = ggplot2::margin(b = 6)),
      axis.title.y       = ggplot2::element_text(
        face = "bold", size = 10, colour = srf$text_primary),
      strip.text         = ggplot2::element_text(face = "bold", size = 10),
      plot.margin        = ggplot2::margin(8, 12, 8, 8))
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    w <- width %||% max(6, length(unique(dt$panel)) * 1.8)
    files <- .hbr_save_fig(p, file_name, "selected_set_coverage",
                           output_dir, w, height, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = p, files = files))
}

#' Per-trait GEBV profile of the top ranked parents
#'
#' Draws a heatmap of the top ranked parents across the traits. Each
#' trait's GEBV is sign-flipped so that "higher is always better", then
#' scaled by its own maximum absolute value, so the cells are comparable
#' across traits on the same colour scale. The trait column order is
#' preserved; the parent rows are ordered by the multi-trait index.
#'
#' @param x An `hgsi_advanced` object.
#' @param top_n Number of top parents to display. Default `25L`.
#' @param show_values Print the numeric scaled value inside each cell.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param output_dir Directory for the figure. `NULL` returns the plot.
#' @param file_name Base file name.
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param verbose Print progress.
#'
#' @return Invisibly, `list(plot = <ggplot>, files = <character>)`.
#'
#' @noRd
hbr_plot_parent_profile <- function(
    x,
    top_n       = 25L,
    show_values = FALSE,
    palette     = "haplobreed",
    theme       = "haplobreed",
    output_dir  = NULL,
    file_name   = "parent_profile",
    width       = 7,
    height      = NULL,
    formats     = c("pdf", "png"),
    verbose     = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  if (!inherits(x, "hgsi_advanced"))
    stop("`x` must be an hgsi_advanced object.", call. = FALSE)
  
  traits <- x$traits
  top_n  <- min(as.integer(top_n), nrow(x$parents))
  par    <- x$parents[seq_len(top_n)]
  
  # Matrix: parents x traits, with direction applied
  mat <- as.matrix(par[, paste0(traits, "_gebv"), with = FALSE])
  colnames(mat) <- traits
  for (tr in traits) {
    if (identical(x$directions[[tr]], "lower"))
      mat[, tr] <- -mat[, tr]
    mx <- max(abs(mat[, tr]), na.rm = TRUE)
    if (is.finite(mx) && mx > 0) mat[, tr] <- mat[, tr] / mx
  }
  rownames(mat) <- par$id
  
  long <- data.table::data.table(
    parent = factor(rep(rownames(mat), times = length(traits)),
                    levels = rev(rownames(mat))),
    trait  = factor(rep(colnames(mat), each = nrow(mat)),
                    levels = traits),
    value  = as.vector(mat))
  
  pal <- .hbr_palette_anchors(palette)
  srf <- .hbr_theme_surface(theme)
  th  <- hbr_theme(theme, base_size = 10)
  
  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = trait, y = parent, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low      = srf$text_muted,
      mid      = "white",
      high     = pal$superior,
      midpoint = 0, limits = c(-1, 1),
      name     = "GEBV\n(scaled)") +
    ggplot2::labs(
      x = NULL, y = NULL,
      title    = "Trait GEBV profile of the top parents",
      subtitle = "Direction applied so higher is always better; scaled per trait") +
    th +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 8.5,
                                           colour = srf$text_primary),
      axis.text.x  = ggplot2::element_text(size = 10, face = "bold",
                                           colour = srf$text_primary),
      plot.title   = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = srf$text_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(
        size = 9, hjust = 0, colour = srf$text_secondary,
        margin = ggplot2::margin(b = 6)),
      legend.position = "right",
      legend.key.height = grid::unit(1.3, "cm"),
      legend.key.width  = grid::unit(0.4, "cm"),
      plot.margin = ggplot2::margin(8, 8, 8, 8))
  
  if (isTRUE(show_values))
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", value)),
      size = 2.4, colour = srf$text_primary)
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    h <- height %||% max(4, 0.24 * top_n + 1.6)
    files <- .hbr_save_fig(p, file_name, "parent_profile",
                           output_dir, width, h, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = p, files = files))
}

#' Cross recommendation landscape
#'
#' Draws a symmetric grid of the recommended crosses: the top ranked
#' parents on both axes, with a filled cell where the two parents appear
#' together in the crosses table. Cell colour is usefulness; the two
#' axes are ordered by the multi-trait index. Cells that do not appear in
#' the crosses table are left blank, so the plot shows the sparse
#' structure of the recommended pairs.
#'
#' @param x An `hgsi_advanced` object or a list with elements `crosses`
#'   (data.frame with columns `parent_a`, `parent_b`, `uc`) and
#'   `parents` (data.frame with columns `id`, `rank`).
#' @param top_n Number of parents shown on each axis. Default `12L`.
#' @param palette One of [hbr_palette_names()].
#' @param theme One of [hbr_theme_names()].
#' @param output_dir Directory for the figure. `NULL` returns the plot.
#' @param file_name Base file name.
#' @param width,height Figure dimensions in inches.
#' @param formats Output formats.
#' @param verbose Print progress.
#'
#' @return Invisibly, `list(plot = <ggplot>, files = <character>)`.
#'
#' @noRd
hbr_plot_cross_landscape <- function(
    x,
    top_n      = 12L,
    palette    = "haplobreed",
    theme      = "haplobreed",
    output_dir = NULL,
    file_name  = "cross_landscape",
    width      = 7,
    height     = 7,
    formats    = c("pdf", "png"),
    verbose    = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  if (inherits(x, "hgsi_advanced")) {
    cr <- data.table::as.data.table(x$crosses)
    par <- data.table::as.data.table(x$parents)
  } else if (is.list(x) && all(c("crosses", "parents") %in% names(x))) {
    cr  <- data.table::as.data.table(x$crosses)
    par <- data.table::as.data.table(x$parents)
  } else {
    stop("`x` must be an hgsi_advanced object or a list with ",
         "`crosses` and `parents`.", call. = FALSE)
  }
  
  required_cr <- c("parent_a", "parent_b", "uc")
  miss <- setdiff(required_cr, names(cr))
  if (length(miss)) stop("crosses missing column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  if (!all(c("id", "rank") %in% names(par)))
    stop("parents must have columns `id` and `rank`.", call. = FALSE)
  
  # Restrict to top_n parents by rank
  top_par <- par[order(rank)][seq_len(min(as.integer(top_n), nrow(par)))]
  top_ids <- as.character(top_par$id)
  
  cr <- cr[parent_a %in% top_ids & parent_b %in% top_ids]
  if (!nrow(cr))
    stop("No recommended crosses between the top ", top_n,
         " parents.", call. = FALSE)
  
  # Symmetric long table
  m1 <- data.table::data.table(pa = cr$parent_a, pb = cr$parent_b,
                               uc = cr$uc)
  m2 <- data.table::data.table(pa = cr$parent_b, pb = cr$parent_a,
                               uc = cr$uc)
  mm <- data.table::rbindlist(list(m1, m2))
  
  ranks <- stats::setNames(top_par$rank, top_par$id)
  mm[, x := match(pa, top_ids)]
  mm[, y := match(pb, top_ids)]
  
  # Position labels with rank
  labs <- sprintf("%s\n(#%d)", top_ids, ranks[top_ids])
  names(labs) <- top_ids
  
  pal <- .hbr_palette_anchors(palette)
  srf <- .hbr_theme_surface(theme)
  th  <- hbr_theme(theme, base_size = 10)
  
  # Sequential gradient: muted -> superior
  grad <- grDevices::colorRampPalette(c("#F5F5F5", pal$superior))(100)
  
  p <- ggplot2::ggplot(mm,
                       ggplot2::aes(x = x, y = y, fill = uc)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", uc)),
                       size = 2.6, colour = srf$text_primary) +
    ggplot2::scale_fill_gradientn(
      colours = grad, name = "UC",
      limits  = c(min(mm$uc), max(mm$uc))) +
    ggplot2::scale_x_continuous(breaks = seq_along(top_ids),
                                labels = labs[top_ids],
                                expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_along(top_ids),
                             labels = labs[top_ids],
                             expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = NULL, y = NULL,
      title    = "Cross recommendation landscape",
      subtitle = sprintf(
        "Top %d crosses between the %d highest-ranked parents; cell = UC",
        nrow(cr), length(top_ids))) +
    th +
    ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(size = 7.5, angle = 45,
                                           hjust = 1,
                                           colour = srf$text_primary),
      axis.text.y  = ggplot2::element_text(size = 7.5,
                                           colour = srf$text_primary),
      plot.title   = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = srf$text_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(
        size = 9, hjust = 0, colour = srf$text_secondary,
        margin = ggplot2::margin(b = 6)),
      legend.position = "right",
      legend.key.height = grid::unit(1.4, "cm"),
      plot.margin = ggplot2::margin(8, 12, 8, 8))
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files <- .hbr_save_fig(p, file_name, "cross_landscape",
                           output_dir, width, height, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = p, files = files))
}

#' Optimism absorbed by full-data feature discovery
#'
#' For each trait, draws three bars: the naive improvement (M2 vs M0),
#' the honest improvement (M3 vs M0), and the optimism (M2 vs M3) that
#' the naive number absorbed. The height of the optimism bar is what the
#' M0/M2/M3 design was built to expose.
#'
#' @inheritParams hbr_plot_cv_forest
#' @param bar_alpha Alpha for the bars.
#' @param legend_position Legend position.
#'
#' @return Invisibly, `list(plot = <ggplot>, files = <character>)`.
#'
#' @noRd
hbr_plot_cv_optimism <- function(
    x,
    palette          = "haplobreed",
    theme            = "haplobreed",
    bar_alpha        = 0.9,
    legend_position  = "top",
    output_dir       = NULL,
    file_name        = "cv_optimism",
    width            = 7,
    height           = 4,
    formats          = c("pdf", "png"),
    verbose          = TRUE) {
  
  .hbr_check_palette_theme(palette, theme)
  
  if (inherits(x, "hgsi_advanced")) {
    cmp <- data.table::as.data.table(x$cv_compare)
  } else if (is.data.frame(x)) {
    cmp <- data.table::as.data.table(x)
  } else stop("`x` must be an hgsi_advanced object or a data.frame.",
              call. = FALSE)
  
  required <- c("trait", "method", "baseline", "mean_diff")
  miss <- setdiff(required, names(cmp))
  if (length(miss)) stop("`x` is missing column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  
  naive   <- cmp[method == "M2" & baseline == "M0",
                 .(trait, naive = mean_diff)]
  honest  <- cmp[method == "M3" & baseline == "M0",
                 .(trait, honest = mean_diff)]
  optim   <- cmp[method == "M2" & baseline == "M3",
                 .(trait, optimism = mean_diff)]
  dt      <- Reduce(function(a, b) merge(a, b, by = "trait"),
                    list(naive, honest, optim))
  
  long <- data.table::melt(
    dt, id.vars = "trait",
    measure.vars = c("naive", "honest", "optimism"),
    variable.name = "kind", value.name = "gain")
  long[, kind := factor(kind,
                        levels = c("naive", "honest", "optimism"),
                        labels = c("Naive  (M2 \u2212 M0)",
                                   "Honest (M3 \u2212 M0)",
                                   "Optimism  (M2 \u2212 M3)"))]
  
  pal <- .hbr_palette_anchors(palette)
  srf <- .hbr_theme_surface(theme)
  th  <- hbr_theme(theme, base_size = 11)
  
  bar_cols <- stats::setNames(
    c(srf$text_muted, pal$superior, pal$other_dark),
    levels(long$kind))
  
  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = trait, y = gain, fill = kind)) +
    ggplot2::geom_hline(yintercept = 0, colour = srf$text_muted,
                        linewidth = 0.4) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.82),
                      width = 0.72, alpha = bar_alpha) +
    ggplot2::scale_fill_manual(values = bar_cols, name = NULL) +
    ggplot2::labs(
      x = NULL,
      y = expression(Delta * r),
      title    = "Honest improvement vs absorbed optimism",
      subtitle = "Optimism = what the naive (M2) accuracy would have overstated") +
    th +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(
        face = "bold", size = 12, hjust = 0,
        colour = srf$text_primary, margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(
        size = 9.5, hjust = 0, colour = srf$text_secondary,
        margin = ggplot2::margin(b = 6)),
      legend.position    = legend_position,
      legend.text        = ggplot2::element_text(size = 9),
      legend.key.size    = grid::unit(0.4, "cm"),
      axis.title.y       = ggplot2::element_text(
        face = "bold", size = 10, colour = srf$text_primary),
      plot.margin        = ggplot2::margin(8, 12, 8, 8))
  
  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files <- .hbr_save_fig(p, file_name, "cv_optimism",
                           output_dir, width, height, formats)
    if (verbose && length(files))
      message("  Figure(s): ", paste(basename(files), collapse = ", "))
  }
  invisible(list(plot = p, files = files))
}
