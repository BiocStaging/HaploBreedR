# ============================================================================
# R/plot-functions.r
# All user-facing plotting functions.
#
# Merged from: plot-haplotype.r, plot-gene-list.r, plot-gene-structure.r,
#              plot-haplotype-network.r, plot-haplotype-table.r,
#              plot-ld-heatmap.r, plot-superior.r
# ============================================================================


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
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' hbr_plot_gene_list(
#'   gene_list     = "gene1",
#'   output_dir    = tempdir(),
#'   per_gene_dir  = toy,
#'   trait_col     = "trait_value",
#'   formats       = "pdf",
#'   combined_pdf  = FALSE,
#'   verbose       = FALSE)
#'
#' @export
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
#' @export
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
#' @export
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
