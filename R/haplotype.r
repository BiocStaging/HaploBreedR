# ============================================================================
# R/haplotype.r
# Variant discovery, haplotype matrix construction, and haplotype
# grouping by subpopulation.
#
# Merged from: variant-identification.r, haplotype-build.r, subpopulation-freq.r
# ============================================================================


# ----------------------------------------------------------------------------
# Variant identification helpers
# ----------------------------------------------------------------------------

#' Resolve a RAP-DB ID to a GFF3 gene_id
#'
#' Attempts, in order: a direct match in `cds_dt$gene_id`; the `g` -> `t`
#' transcription-unit substitution; and a coordinate-overlap match
#' against `all_genes` with the majority transcript as the winner.
#'
#' @param rapdb_id RAP-DB gene identifier.
#' @param cds_dt CDS feature table.
#' @param all_genes Gene table.
#'
#' @return A GFF3 gene_id, or `NULL`.
#'
#' @noRd
hbr_resolve_gene_id <- function(rapdb_id, cds_dt, all_genes = NULL) {
  if (rapdb_id %in% cds_dt$gene_id) return(rapdb_id)

  tids <- cds_dt[grepl(sub("g", "t", rapdb_id), transcript_id), transcript_id]
  if (length(tids) > 0) return(cds_dt[transcript_id == tids[1], gene_id][1])

  if (is.null(all_genes) || !nrow(all_genes)) return(NULL)

  rg <- all_genes[gene_id == rapdb_id]
  if (nrow(rg) == 0) return(NULL)
  chr  <- as.character(rg$Chr[1]); gst <- rg$Start[1]; gend <- rg$End[1]

  hits <- cds_dt[Chr == chr & start <= gend & end >= gst]
  if (nrow(hits) == 0) return(NULL)
  best <- hits[, .N, by = gene_id][order(-N)]
  best$gene_id[1]
}

#' GFF3 to FASTA transcript pairs
#'
#' For a gene, returns the table of its transcripts that also exist in
#' the CDS FASTA. Transcripts are ordered by total CDS length.
#'
#' @param gene_gff_id GFF3 gene identifier.
#' @param cds_dt CDS feature table.
#' @param cds_seqs CDS sequence set.
#'
#' @return A `data.table` with `gff_tx` and `fasta_tx` columns.
#'
#' @noRd
hbr_get_transcripts <- function(gene_gff_id, cds_dt, cds_seqs) {
  gene_cds <- cds_dt[gene_id == gene_gff_id]
  if (!nrow(gene_cds))
    return(data.table::data.table(gff_tx = character(0),
                                   fasta_tx = character(0)))

  tx_len <- gene_cds[, .(CDS_bp = sum(width)), by = transcript_id]
  data.table::setorder(tx_len, -CDS_bp)

  out <- data.table::data.table(gff_tx = character(0),
                                fasta_tx = character(0))
  for (tx in tx_len$transcript_id) {
    if (tx %in% names(cds_seqs)) {
      out <- rbind(out, data.table::data.table(gff_tx = tx, fasta_tx = tx))
    } else {
      alt <- c(paste0("transcript:", tx), sub("-\\d+$", ".1", tx))
      hit <- alt[alt %in% names(cds_seqs)]
      if (length(hit))
        out <- rbind(out, data.table::data.table(gff_tx = tx,
                                                 fasta_tx = hit[1]))
    }
  }
  out
}

#' Canonical transcript for a gene
#'
#' Longest-CDS transcript that also exists in the CDS FASTA.
#'
#' @param gene_gff_id GFF3 gene identifier.
#' @param cds_dt CDS feature table.
#' @param cds_seqs CDS sequence set.
#'
#' @return A FASTA transcript name, or `NULL`.
#'
#' @noRd
hbr_get_canonical_transcript <- function(gene_gff_id, cds_dt, cds_seqs) {
  tx_table <- hbr_get_transcripts(gene_gff_id, cds_dt, cds_seqs)
  if (!nrow(tx_table)) return(NULL)
  tx_table$fasta_tx[1]
}

#' Build a CDS coordinate map for a transcript
#'
#' Orders the CDS exons of a transcript by strand and assigns cumulative
#' CDS offsets.
#'
#' @param transcript_id Transcript identifier.
#' @param cds_dt CDS feature table.
#'
#' @return A `data.table` with `Chr`, `Start`, `End`, `Strand`, `Width`,
#'   `CDS_Pos_Start`, `CDS_Pos_End`, `Transcript_ID`, or `NULL`.
#'
#' @noRd
hbr_build_cds_map <- function(transcript_id, cds_dt) {
  target_tx <- transcript_id
  tx_cds <- cds_dt[cds_dt$transcript_id == target_tx,
                   .(Chr = Chr, Start = start, End = end,
                     Strand = strand, Width = width)]
  if (nrow(tx_cds) == 0) return(NULL)

  strand <- unique(tx_cds$Strand)[1]
  if (strand == "+") data.table::setorder(tx_cds, Start)
  else data.table::setorder(tx_cds, -Start)

  tx_cds[, CDS_Pos_Start := cumsum(c(1L, utils::head(Width, -1)))]
  tx_cds[, CDS_Pos_End   := cumsum(Width)]
  tx_cds[, Transcript_ID := target_tx]
  tx_cds
}

#' Build a UTR coordinate map for a transcript
#'
#' @param transcript_id Transcript identifier.
#' @param utr_dt UTR feature table.
#'
#' @return A `data.table` with `Chr`, `Start`, `End`, `Strand`, `Width`,
#'   `UTR_Type`, `Transcript_ID`, or `NULL`.
#'
#' @noRd
hbr_build_utr_map <- function(transcript_id, utr_dt) {
  target_tx <- transcript_id
  tx_utr <- utr_dt[utr_dt$transcript_id == target_tx,
                   .(Chr = Chr, Start = start, End = end,
                     Strand = strand, Width = width, UTR_Type = UTR_Type)]
  if (!nrow(tx_utr)) return(NULL)
  strand <- unique(tx_utr$Strand)[1]
  if (strand == "+") data.table::setorder(tx_utr, Start)
  else data.table::setorder(tx_utr, -Start)
  tx_utr[, Transcript_ID := target_tx]
  tx_utr
}

#' Union of CDS intervals across a gene's transcripts
#'
#' @param cds_dt CDS feature table.
#' @param gene_gff_id GFF3 gene identifier.
#' @param transcripts Character vector of transcript IDs.
#'
#' @return A list with `bp` and `n_exons`.
#'
#' @noRd
hbr_union_cds_intervals <- function(cds_dt, gene_gff_id, transcripts) {
  sub <- cds_dt[gene_id == gene_gff_id & transcript_id %in% transcripts]
  if (!nrow(sub)) return(list(bp = 0L, n_exons = 0L))
  gr   <- GenomicRanges::GRanges(seqnames = sub$Chr,
                                 ranges = IRanges::IRanges(start = sub$start,
                                                           end = sub$end),
                                 strand = sub$strand[1])
  gr_u <- GenomicRanges::reduce(gr)
  list(bp = sum(GenomicRanges::width(gr_u)), n_exons = length(gr_u))
}

#' Union of UTR intervals across a gene's transcripts
#'
#' @param utr_dt UTR feature table.
#' @param gene_gff_id GFF3 gene identifier.
#' @param transcripts Character vector of transcript IDs.
#'
#' @return A list with `bp5` and `bp3`.
#'
#' @noRd
hbr_union_utr_intervals <- function(utr_dt, gene_gff_id, transcripts) {
  sub <- utr_dt[gene_id == gene_gff_id & transcript_id %in% transcripts]
  if (!nrow(sub)) return(list(bp5 = 0L, bp3 = 0L))
  calc_bp <- function(tp) {
    s <- sub[UTR_Type == tp]
    if (!nrow(s)) return(0L)
    gr <- GenomicRanges::GRanges(seqnames = s$Chr,
                                 ranges = IRanges::IRanges(s$start, s$end),
                                 strand = s$strand[1])
    sum(GenomicRanges::width(GenomicRanges::reduce(gr)))
  }
  list(bp5 = calc_bp("5UTR"), bp3 = calc_bp("3UTR"))
}

#' Map SNPs to CDS exons
#'
#' For each CDS exon in `cds_map`, finds all SNPs in `ref_alt` that fall
#' within the exon on the same chromosome and computes their position
#' within the concatenated CDS.
#'
#' @param cds_map A CDS coordinate map.
#' @param ref_alt A SNP table.
#'
#' @return A `data.table` of SNP-to-CDS hits, or `NULL`.
#'
#' @noRd
hbr_map_snps_to_cds <- function(cds_map, ref_alt) {
  tryCatch({
    if (is.null(cds_map) || !nrow(cds_map)) return(NULL)
    if (is.null(ref_alt) || !nrow(ref_alt)) return(NULL)

    find_col <- function(dt, candidates, required = TRUE) {
      hit <- candidates[candidates %in% names(dt)][1]
      if (is.na(hit) && required)
        stop("Required column not found. Expected one of: ",
             paste(candidates, collapse = ", "), call. = FALSE)
      hit
    }

    chr_col <- find_col(cds_map,
                        c("Chr", "CHR", "Chromosome", "chromosome", "seqid"))
    start_col <- find_col(cds_map, c("Start", "start", "START"))
    end_col <- find_col(cds_map, c("End", "end", "END"))
    strand_col <- find_col(cds_map, c("Strand", "strand"))
    cds_start_col <- find_col(cds_map,
                              c("CDS_Pos_Start", "CDS_Start", "cds_start",
                                "CDS_Position_Start"))
    cds_end_col <- find_col(cds_map,
                            c("CDS_Pos_End", "CDS_End", "cds_end",
                              "CDS_Position_End"))

    snp_chr_col <- find_col(ref_alt,
                            c("CHR", "Chr", "chr", "Chromosome",
                              "chromosome", "seqid"))
    pos_col <- find_col(ref_alt,
                        c("POS", "Pos", "pos", "Position", "position",
                          "BP", "bp"))
    snp_id_col <- find_col(ref_alt,
                           c("SNP", "Snp", "snp", "Marker", "marker",
                             "ID", "id"))
    ref_col <- find_col(ref_alt, c("REF", "Ref", "ref", "Reference",
                                    "reference"))
    alt_col <- find_col(ref_alt, c("ALT", "Alt", "alt", "Alternative",
                                    "alternative"))

    normalize_chr <- function(x) {
      x <- trimws(as.character(x))
      x <- sub("^(chr|Chr|CHR)", "", x)
      x <- sub("^chromosome[_-]?", "", x, ignore.case = TRUE)
      x <- sub("^0+([0-9]+)$", "\\1", x)
      toupper(x)
    }

    cds_chr <- normalize_chr(cds_map[[chr_col]])
    snp_chr <- normalize_chr(ref_alt[[snp_chr_col]])

    snps <- data.table::data.table(
      SNP = as.character(ref_alt[[snp_id_col]]),
      POS = as.numeric(ref_alt[[pos_col]]),
      REF = toupper(as.character(ref_alt[[ref_col]])),
      ALT = toupper(as.character(ref_alt[[alt_col]])),
      Chr = snp_chr)

    target_chr <- unique(cds_chr[!is.na(cds_chr)])
    if (length(target_chr) == 1L) {
      snps_chr <- snps[Chr == target_chr]
    } else {
      snps_chr <- snps[Chr %in% target_chr]
    }
    if (!nrow(snps_chr)) return(NULL)

    data.table::setorder(snps_chr, POS)

    hits_list <- vector("list", nrow(cds_map))
    h <- 0L

    for (i in seq_len(nrow(cds_map))) {
      chr_i      <- cds_chr[i]
      start_i    <- as.numeric(cds_map[[start_col]][i])
      end_i      <- as.numeric(cds_map[[end_col]][i])
      strand_i   <- as.character(cds_map[[strand_col]][i])
      cds_start_i <- as.numeric(cds_map[[cds_start_col]][i])
      cds_end_i  <- as.numeric(cds_map[[cds_end_col]][i])

      if (is.na(start_i) || is.na(end_i) || is.na(cds_start_i)) next

      if (!is.na(chr_i) && chr_i != target_chr[1]) next

      lo <- findInterval(start_i, snps_chr$POS) + 1L
      hi <- findInterval(end_i,   snps_chr$POS)
      if (hi < lo) next

      s <- data.table::copy(snps_chr[lo:hi])

      h <- h + 1L
      s[, `:=`(Strand = strand_i, Exon_Start = start_i, Exon_End = end_i,
               Exon_CDS_Start = cds_start_i, Exon_CDS_End = cds_end_i)]
      hits_list[[h]] <- s
    }

    if (!h) return(NULL)

    hits <- data.table::rbindlist(hits_list[seq_len(h)], fill = TRUE)
    hits <- unique(hits, by = c("SNP", "POS", "Chr",
                                "Exon_Start", "Exon_End"))

    hits[, Pos_in_CDS := ifelse(Strand == "+",
                                Exon_CDS_Start + (POS - Exon_Start),
                                Exon_CDS_Start + (Exon_End - POS))]
    hits[, Codon_Num := ceiling(Pos_in_CDS / 3)]
    hits[, Codon_Pos := ((Pos_in_CDS - 1L) %% 3L) + 1L]
    data.table::setorder(hits, Chr, Pos_in_CDS, SNP)

    hits
  }, error = function(e) NULL)
}

#' Map SNPs to UTR intervals
#'
#' @param utr_map A UTR coordinate map.
#' @param ref_alt A SNP table.
#'
#' @return A `data.table` of SNP-to-UTR hits, or `NULL`.
#'
#' @noRd
hbr_map_snps_to_utr <- function(utr_map, ref_alt) {
  tryCatch({
    if (is.null(utr_map) || !nrow(utr_map)) return(NULL)
    if (is.null(ref_alt) || !nrow(ref_alt)) return(NULL)

    find_col <- function(dt, candidates, required = TRUE) {
      hit <- candidates[candidates %in% names(dt)][1]
      if (is.na(hit) && required)
        stop("Required column not found. Expected one of: ",
             paste(candidates, collapse = ", "), call. = FALSE)
      hit
    }

    chr_col <- find_col(utr_map,
                        c("Chr", "CHR", "Chromosome", "chromosome", "seqid"))
    start_col <- find_col(utr_map, c("Start", "start", "START"))
    end_col <- find_col(utr_map, c("End", "end", "END"))
    strand_col <- find_col(utr_map, c("Strand", "strand"))
    region_col <- find_col(utr_map,
                           c("UTR_Type", "Region", "region", "Type", "type"))

    snp_chr_col <- find_col(ref_alt,
                            c("CHR", "Chr", "chr", "Chromosome",
                              "chromosome", "seqid"))
    pos_col <- find_col(ref_alt,
                        c("POS", "Pos", "pos", "Position", "position",
                          "BP", "bp"))
    snp_id_col <- find_col(ref_alt,
                           c("SNP", "Snp", "snp", "Marker", "marker",
                             "ID", "id"))
    ref_col <- find_col(ref_alt, c("REF", "Ref", "ref"))
    alt_col <- find_col(ref_alt, c("ALT", "Alt", "alt"))

    normalize_chr <- function(x) {
      x <- trimws(as.character(x))
      x <- sub("^(chr|Chr|CHR)", "", x)
      x <- sub("^chromosome[_-]?", "", x, ignore.case = TRUE)
      x <- sub("^0+([0-9]+)$", "\\1", x)
      toupper(x)
    }

    utr_chr <- normalize_chr(utr_map[[chr_col]])
    snp_chr <- normalize_chr(ref_alt[[snp_chr_col]])

    snps <- data.table::data.table(
      SNP = as.character(ref_alt[[snp_id_col]]),
      POS = as.numeric(ref_alt[[pos_col]]),
      REF = toupper(as.character(ref_alt[[ref_col]])),
      ALT = toupper(as.character(ref_alt[[alt_col]])),
      Chr = snp_chr)

    hits_list <- vector("list", nrow(utr_map))
    h <- 0L

    for (i in seq_len(nrow(utr_map))) {
      chr_i <- utr_chr[i]
      start_i <- as.numeric(utr_map[[start_col]][i])
      end_i <- as.numeric(utr_map[[end_col]][i])
      strand_i <- as.character(utr_map[[strand_col]][i])
      region_i <- as.character(utr_map[[region_col]][i])

      if (is.na(start_i) || is.na(end_i)) next

      s <- snps[Chr == chr_i & POS >= start_i & POS <= end_i]
      if (!nrow(s)) next

      h <- h + 1L
      s[, `:=`(Strand = strand_i, Region = region_i)]
      hits_list[[h]] <- s
    }

    if (!h) return(NULL)
    unique(data.table::rbindlist(hits_list[seq_len(h)], fill = TRUE))
  }, error = function(e) NULL)
}

#' Validate a transcript's frame and translation
#'
#' @param transcript_id Transcript identifier.
#' @param fasta_tx FASTA transcript name.
#' @param cds_seqs CDS sequence set.
#' @param pep_seqs Peptide sequence set.
#'
#' @return A list with `ok` and `reason`.
#'
#' @noRd
hbr_validate_transcript_frame <- function(transcript_id, fasta_tx,
                                           cds_seqs, pep_seqs) {
  tx_seq <- cds_seqs[[fasta_tx]]
  if (is.null(tx_seq))
    return(list(ok = FALSE, reason = "CDS sequence not found"))
  cds_len <- length(tx_seq)
  if (cds_len %% 3 != 0)
    return(list(ok = FALSE,
                reason = sprintf("CDS length (%d bp) not a multiple of 3",
                                 cds_len)))

  pep_seq <- pep_seqs[[fasta_tx]]
  if (is.null(pep_seq)) return(list(ok = TRUE, reason = NA_character_))

  pep_len <- length(pep_seq)
  if (pep_len > 0L) {
    pep_last <- as.character(Biostrings::subseq(pep_seq, pep_len, pep_len))
    if (pep_last == "*") pep_len <- pep_len - 1L
  }

  translated <- tryCatch(
    Biostrings::translate(tx_seq, if.fuzzy.codon = "solve"),
    error = function(e) NULL)
  if (is.null(translated))
    return(list(ok = FALSE, reason = "translation failed"))

  translated_str <- as.character(translated)
  ends_in_stop   <- grepl("\\*$", translated_str)
  core_len       <- nchar(translated_str) - as.integer(ends_in_stop)

  if (core_len != pep_len) {
    return(list(ok = FALSE, reason = sprintf(
      "translated length (%d) does not match peptide length (%d) for %s",
      core_len, pep_len, transcript_id)))
  }
  list(ok = TRUE, reason = NA_character_)
}

#' Translate variant codons to amino-acid changes
#'
#' @param cds_snps A SNP-to-CDS hit table.
#' @param fasta_tx FASTA transcript name.
#' @param cds_seqs CDS sequence set.
#' @param pep_seqs Optional peptide sequence set.
#'
#' @return The input table with added translation columns.
#'
#' @noRd
hbr_translate_variant <- function(cds_snps, fasta_tx, cds_seqs,
                                   pep_seqs = NULL) {

  if (is.null(cds_snps) || !nrow(cds_snps)) return(NULL)

  tx_seq <- cds_seqs[[fasta_tx]]
  if (is.null(tx_seq)) {
    message("    CDS FASTA entry not found for: ", fasta_tx)
    return(NULL)
  }
  tx_seq  <- Biostrings::DNAString(as.character(tx_seq))
  cds_len <- nchar(as.character(tx_seq))

  pep_seq <- NULL
  if (!is.null(pep_seqs) && !is.null(pep_seqs[[fasta_tx]])) {
    pep_seq <- as.character(pep_seqs[[fasta_tx]])
    pep_seq <- sub("\\*$", "", pep_seq)
    pep_seq <- Biostrings::AAString(pep_seq)
  }

  cds_snps <- data.table::copy(cds_snps)
  cds_snps[, `:=`(
    Codon_REF = NA_character_, Codon_ALT = NA_character_,
    AA_REF = NA_character_, AA_ALT = NA_character_,
    AA_Change = NA_character_, Is_Nonsynonymous = NA,
    AA_Validated = NA, REF_Matches_FASTA = NA)]

  for (k in seq_len(nrow(cds_snps))) {
    codon_n <- as.integer(cds_snps$Codon_Num[k])
    codon_p <- as.integer(cds_snps$Codon_Pos[k])
    strand  <- as.character(cds_snps$Strand[k])
    ref_nt  <- toupper(as.character(cds_snps$REF[k]))
    alt_nt  <- toupper(as.character(cds_snps$ALT[k]))

    if (is.na(codon_n) || is.na(codon_p) ||
        is.na(ref_nt) || is.na(alt_nt)) next
    if (nchar(ref_nt) != 1 || nchar(alt_nt) != 1) next
    if (!grepl("^[ACGT]$", ref_nt) || !grepl("^[ACGT]$", alt_nt)) next

    codon_start <- (codon_n - 1L) * 3L + 1L
    codon_end   <- codon_start + 2L
    if (codon_start < 1 || codon_end > cds_len) next

    ref_codon <- as.character(Biostrings::subseq(tx_seq,
                                                  start = codon_start,
                                                  end = codon_end))

    if (strand == "-") {
      ref_cds <- as.character(Biostrings::complement(
        Biostrings::DNAString(ref_nt)))
      alt_cds <- as.character(Biostrings::complement(
        Biostrings::DNAString(alt_nt)))
    } else {
      ref_cds <- ref_nt
      alt_cds <- alt_nt
    }

    fasta_base <- substr(ref_codon, codon_p, codon_p)
    ref_ok     <- identical(ref_cds, fasta_base)
    cds_snps[k, REF_Matches_FASTA := ref_ok]
    if (!ref_ok) next

    alt_chars <- strsplit(ref_codon, "", fixed = TRUE)[[1]]
    alt_chars[codon_p] <- alt_cds
    alt_codon <- paste0(alt_chars, collapse = "")

    aa_ref <- tryCatch(
      as.character(Biostrings::translate(Biostrings::DNAString(ref_codon),
                                          if.fuzzy.codon = "X")),
      error = function(e) NA_character_)
    aa_alt <- tryCatch(
      as.character(Biostrings::translate(Biostrings::DNAString(alt_codon),
                                          if.fuzzy.codon = "X")),
      error = function(e) NA_character_)
    if (is.na(aa_ref) || is.na(aa_alt)) next

    cds_snps[k, `:=`(Codon_REF = ref_codon, Codon_ALT = alt_codon,
                     AA_REF = aa_ref, AA_ALT = aa_alt)]
    is_ns <- !identical(aa_ref, aa_alt)
    cds_snps[k, Is_Nonsynonymous := is_ns]

    aa_valid <- NA
    if (!is.null(pep_seq)) {
      pep_len <- nchar(as.character(pep_seq))
      if (codon_n >= 1 && codon_n <= pep_len) {
        pep_aa <- as.character(Biostrings::subseq(pep_seq,
                                                   start = codon_n,
                                                   end = codon_n))
        aa_valid <- identical(aa_ref, pep_aa)
      }
    }
    cds_snps[k, AA_Validated := aa_valid]
    cds_snps[k, AA_Change := hbr_aa_change_label(aa_ref, aa_alt, codon_n)]
  }

  cds_snps
}

#' Identify functional variants for one gene
#'
#' @param rapdb_id RAP-DB or matching gene identifier.
#' @param cds_dt CDS feature table.
#' @param utr_dt UTR feature table.
#' @param cds_seqs CDS sequence set.
#' @param pep_seqs Peptide sequence set.
#' @param all_genes Gene table.
#' @param ref_alt SNP table.
#' @param use_all_transcripts Use all FASTA-matched transcripts.
#' @param regions_to_include Regions retained in `functional_snps`.
#'
#' @return A list describing the gene's functional variants, or `NULL`.
#'
#' @noRd
hbr_get_functional_variants <- function(rapdb_id, cds_dt, utr_dt,
                                         cds_seqs, pep_seqs, all_genes,
                                         ref_alt,
                                         use_all_transcripts = TRUE,
                                         regions_to_include =
                                           c("CDS_nonsynonymous",
                                             "5UTR", "3UTR")) {

  message("\n  --- Functional variant identification ---")

  gene_gff_id <- hbr_resolve_gene_id(rapdb_id, cds_dt, all_genes)
  if (is.null(gene_gff_id)) {
    message("  Gene not found in GFF3: ", rapdb_id); return(NULL)
  }
  message("  GFF3 gene_id: ", gene_gff_id)

  tx_table <- hbr_get_transcripts(gene_gff_id, cds_dt, cds_seqs)
  if (!nrow(tx_table)) {
    cat("  No FASTA-matched transcript for:", gene_gff_id, "\n"); return(NULL)
  }
  if (!use_all_transcripts) tx_table <- tx_table[1]
  message("  Transcripts considered: ", nrow(tx_table),
          " (", paste(tx_table$gff_tx, collapse = ", "), ")")

  strand_g <- unique(cds_dt[gene_id == gene_gff_id, strand])[1]

  frame_warnings <- character(0)
  cds_per_tx <- list()
  for (i in seq_len(nrow(tx_table))) {
    gff_tx <- tx_table$gff_tx[i]; fasta_tx <- tx_table$fasta_tx[i]

    fcheck <- hbr_validate_transcript_frame(gff_tx, fasta_tx,
                                             cds_seqs, pep_seqs)
    if (!fcheck$ok)
      frame_warnings <- c(frame_warnings,
                          sprintf("%s: %s", gff_tx, fcheck$reason))

    cds_map_tx <- hbr_build_cds_map(gff_tx, cds_dt)
    if (is.null(cds_map_tx)) next
    snps_tx <- hbr_map_snps_to_cds(cds_map_tx, ref_alt)
    if (is.null(snps_tx) || !nrow(snps_tx)) next
    snps_tx <- hbr_translate_variant(snps_tx, fasta_tx, cds_seqs, pep_seqs)
    if (is.null(snps_tx)) next

    snps_tx[, Transcript_ID := gff_tx]
    snps_tx[, Region := ifelse(Is_Nonsynonymous,
                                "CDS_nonsynonymous", "CDS_synonymous")]
    snps_tx[, Variant_Label := ifelse(Is_Nonsynonymous, AA_Change,
                                       paste0("syn:", REF, ">", ALT))]
    cds_per_tx[[gff_tx]] <- snps_tx
  }
  if (length(frame_warnings))
    message("  Frame-validation warning(s): ",
            paste(frame_warnings, collapse = "; "))

  if (length(cds_per_tx)) {
    cds_all <- data.table::rbindlist(cds_per_tx, fill = TRUE)
  } else {
    cds_all <- data.table::data.table()
  }
  if (!"Region" %in% names(cds_all)) cds_all[, Region := character(.N)]

  utr_per_tx <- list()
  if (!is.null(utr_dt) && nrow(utr_dt)) {
    for (i in seq_len(nrow(tx_table))) {
      gff_tx <- tx_table$gff_tx[i]
      utr_map_tx <- hbr_build_utr_map(gff_tx, utr_dt)
      if (is.null(utr_map_tx)) next
      snps_tx <- hbr_map_snps_to_utr(utr_map_tx, ref_alt)
      if (is.null(snps_tx) || !nrow(snps_tx)) next
      snps_tx[, Transcript_ID := gff_tx]
      snps_tx[, AA_Change := NA_character_]
      snps_tx[, Variant_Label := paste0(Region, ":", REF, ">", ALT)]
      snps_tx[, Is_Nonsynonymous := NA]
      snps_tx[, AA_Validated := NA]
      snps_tx[, REF_Matches_FASTA := NA]
      utr_per_tx[[gff_tx]] <- snps_tx
    }
  }
  if (length(utr_per_tx)) {
    utr_all <- data.table::rbindlist(utr_per_tx, fill = TRUE)
  } else {
    utr_all <- data.table::data.table()
  }
  if (!"Region" %in% names(utr_all)) utr_all[, Region := character(.N)]

  if (!nrow(cds_all) && !nrow(utr_all)) {
    cat("  No SNPs overlap CDS or UTR (any transcript)\n")
    return(NULL)
  }

  combined <- data.table::rbindlist(list(cds_all, utr_all), fill = TRUE)
  if (!"Region" %in% names(combined)) combined[, Region := character(.N)]
  combined[, Region := as.character(Region)]

  region_priority <- c(CDS_nonsynonymous = 1L,
                       `5UTR`            = 2L,
                       `3UTR`            = 2L,
                       CDS_synonymous    = 3L)
  combined$.rank <- unname(region_priority[combined$Region])
  data.table::setorderv(combined, c("SNP", ".rank"))
  combined_dedup <- combined[!duplicated(SNP)]
  combined_dedup[, .rank := NULL]

  ns_tx_map <- combined[Is_Nonsynonymous == TRUE,
                        .(NS_Transcripts = paste(unique(Transcript_ID),
                                                  collapse = ";")),
                        by = SNP]
  combined_dedup <- merge(combined_dedup, ns_tx_map, by = "SNP", all.x = TRUE)

  combined_dedup[, Genomic_Order := if (strand_g == "+") POS else -POS]
  data.table::setorder(combined_dedup, Genomic_Order)

  cds_only <- combined_dedup[Region %in% c("CDS_nonsynonymous",
                                            "CDS_synonymous")]
  cat("  CDS SNPs total       :", nrow(cds_only), "\n")
  cat("    Synonymous         :", sum(cds_only$Region == "CDS_synonymous"),
      "\n")
  cat("    Non-synonymous     :",
      sum(cds_only$Region == "CDS_nonsynonymous"), "\n")
  cat("  5'UTR SNPs           :", sum(combined_dedup$Region == "5UTR"),
      "\n")
  cat("  3'UTR SNPs           :", sum(combined_dedup$Region == "3UTR"),
      "\n")

  functional_snps <- combined_dedup[Region %in% regions_to_include]
  if (!nrow(functional_snps)) {
    cat("  No SNPs in requested regions\n"); return(NULL)
  }
  cat("  Functional SNPs kept (",
      paste(regions_to_include, collapse = ","), "): ",
      nrow(functional_snps), "\n", sep = "")
  cat("  Changes:",
      paste(functional_snps$Variant_Label, collapse = ", "), "\n")

  n_ref_mismatch <- sum(!functional_snps$REF_Matches_FASTA, na.rm = TRUE)
  if (n_ref_mismatch > 0)
    cat("  WARNING:", n_ref_mismatch,
        "CDS SNP(s) have a REF allele that does not match the reference ",
        "CDS FASTA base.\n")

  cds_union <- hbr_union_cds_intervals(cds_dt, gene_gff_id, tx_table$gff_tx)
  utr_union <- if (!is.null(utr_dt) && nrow(utr_dt))
    hbr_union_utr_intervals(utr_dt, gene_gff_id, tx_table$gff_tx)
  else list(bp5 = 0L, bp3 = 0L)

  list(gene_gff_id             = gene_gff_id,
       transcript_id           = tx_table$gff_tx[1],
       transcripts_considered  = paste(tx_table$gff_tx, collapse = ";"),
       n_transcripts_considered = nrow(tx_table),
       total_cds_bp            = cds_union$bp,
       n_cds_exons             = cds_union$n_exons,
       utr5_bp                 = utr_union$bp5,
       utr3_bp                 = utr_union$bp3,
       all_cds_snps            = cds_only,
       all_utr_snps            = combined_dedup[Region %in% c("5UTR",
                                                               "3UTR")],
       functional_snps         = functional_snps,
       n_ref_mismatch          = n_ref_mismatch,
       frame_warnings          = frame_warnings)
}


# ----------------------------------------------------------------------------
# Haplotype matrix construction
# ----------------------------------------------------------------------------

#' Build a haplotype character matrix
#'
#' Converts the numeric dosage matrix of a gene's functional SNPs into a
#' matrix of nucleotide characters. Missing calls are marked `"N"`.
#' Columns and rows with too much missingness are dropped, and a
#' concatenated NT string and a per-individual AA string (nonsynonymous
#' CDS SNPs only) are produced.
#'
#' @param geno_data Output of [hbr_extract_genotypes()].
#'
#' @return A list with `char_matrix`, `hap_nt`, `hap_aa`,
#'   `individuals`, `n_individuals`, `n_snps`, `snp_names`,
#'   `functional_snps_used`. Returns `NULL` if fewer than 10 individuals
#'   remain after filtering.
#'
#' @noRd
hbr_build_haplotype_matrix <- function(geno_data) {
  G <- geno_data$geno
  snps <- geno_data$functional_snps
  n_ind <- nrow(G); n_snps <- ncol(G)

  char_mat <- matrix("N", n_ind, n_snps)
  colnames(char_mat) <- colnames(G); rownames(char_mat) <- rownames(G)

  for (j in seq_len(n_snps)) {
    snp_row <- snps[j]
    ref_nt  <- snp_row$REF
    alt_nt  <- gsub(",.*", "", snp_row$ALT)
    g <- G[, j]
    char_mat[, j] <- ifelse(is.na(g), "N",
                             ifelse(g == 0, ref_nt,
                                    ifelse(g == 2, alt_nt, "N")))
  }

  miss_snp <- apply(char_mat == "N", 2, mean)
  char_mat <- char_mat[, miss_snp <= 0.05, drop = FALSE]
  if (!ncol(char_mat)) return(NULL)

  miss_ind <- apply(char_mat == "N", 1, mean)
  char_mat <- char_mat[miss_ind <= 0.10, , drop = FALSE]
  if (nrow(char_mat) < 10) return(NULL)

  all_n <- apply(char_mat, 1, function(x) all(x == "N"))
  char_mat <- char_mat[!all_n, , drop = FALSE]
  if (nrow(char_mat) < 10) return(NULL)

  hap_nt <- apply(char_mat, 1, paste0, collapse = "-")

  snps_kept <- snps[SNP %in% colnames(char_mat)]
  snps_kept <- snps_kept[order(Genomic_Order)]

  hap_aa <- apply(char_mat, 1, function(row_alleles) {
    vals <- vapply(seq_along(row_alleles), function(j) {
      nt <- row_alleles[j]
      if (nt == "N") return("?")
      is_cds <- !is.na(snps_kept$Region[j]) &&
        snps_kept$Region[j] == "CDS_nonsynonymous"
      if (!is_cds) return(nt)
      if (nt == snps_kept$REF[j]) snps_kept$AA_REF[j]
      else snps_kept$AA_ALT[j]
    }, character(1))
    paste(vals, collapse = "-")
  })

  list(
    char_matrix = char_mat,
    hap_nt = hap_nt,
    hap_aa = hap_aa,
    individuals = rownames(char_mat),
    n_individuals = nrow(char_mat),
    n_snps = ncol(char_mat),
    snp_names = colnames(char_mat),
    functional_snps_used = snps_kept)
}

#' Group individuals into haplotypes
#'
#' Assigns each individual to a haplotype group based on its
#' concatenated NT string, applies minimum-size and minimum-frequency
#' filters, and computes diversity statistics (Nei's Hd, effective
#' number of haplotypes, PIC). For 3+ haplotypes a pegas network is
#' built when `pegas` and `ape` are installed.
#'
#' @param hap_mat_data Output of [hbr_build_haplotype_matrix()].
#' @param min_individuals Minimum group size.
#' @param min_frequency Minimum haplotype frequency.
#'
#' @return A list with `haplotypes`, `haplotype_map`, `n_haplotypes`,
#'   `n_individuals`, `network_data`, `diversity`, `is_single_snp`.
#'   Returns `NULL` if fewer than 2 haplotypes pass the filters.
#'
#' @noRd
hbr_define_haplotypes <- function(hap_mat_data, min_individuals = 3,
                                   min_frequency = 0) {
  seqs  <- hap_mat_data$hap_nt
  aas   <- hap_mat_data$hap_aa
  inds  <- hap_mat_data$individuals
  n_snp <- hap_mat_data$n_snps

  counts <- table(seqs)
  total  <- length(seqs)
  keep   <- names(counts)[!grepl("^(N-?)+$", names(counts))]
  keep   <- keep[counts[keep] >= min_individuals &
                   (counts[keep] / total) >= min_frequency]
  if (length(keep) < 2) return(NULL)

  aa_rep <- vapply(keep, function(nt) aas[which(seqs == nt)[1]],
                   character(1))

  hap_map <- data.table::data.table(
    NT_Sequence = keep,
    AA_Sequence = aa_rep,
    Count       = as.integer(counts[keep]),
    Frequency   = as.numeric(counts[keep]) / total)
  data.table::setorder(hap_map, -Count)
  hap_map$Haplotype <- paste0("Hap", seq_len(nrow(hap_map)))

  hap_assign <- hap_map$Haplotype[match(seqs, hap_map$NT_Sequence)]
  haplo_df <- data.table::data.table(ID = inds, NT_Seq = seqs,
                                     AA_Seq = aas, Haplotype = hap_assign)
  haplo_df <- haplo_df[!is.na(haplo_df$Haplotype), ]

  p  <- hap_map$Frequency
  n  <- sum(hap_map$Count)
  Hd <- if (n > 1) (n / (n - 1)) * (1 - sum(p^2)) else NA
  Neff <- 1 / sum(p^2)
  PIC  <- 1 - sum(p^2) -
    sum(outer(p, p, function(a, b) 2 * a^2 * b^2)[
      upper.tri(matrix(0, length(p), length(p)))])

  net_data <- NULL
  if (length(keep) >= 3 &&
      requireNamespace("ape", quietly = TRUE) &&
      requireNamespace("pegas", quietly = TRUE)) {
    char_mat_haps <- do.call(rbind, lapply(keep, function(s)
      strsplit(s, "-")[[1]]))
    rownames(char_mat_haps) <- keep

    for (j in seq_len(ncol(char_mat_haps))) {
      col_j <- char_mat_haps[, j]
      is_n  <- toupper(col_j) == "N"
      if (any(is_n)) {
        non_n <- col_j[!is_n]
        if (length(non_n) > 0) {
          modal_base <- names(sort(table(non_n), decreasing = TRUE))[1]
          if (!is.na(modal_base))
            char_mat_haps[is_n, j] <- modal_base
        }
      }
    }

    net_data <- tryCatch({
      dna <- ape::as.DNAbin(char_mat_haps)
      haps_obj <- suppressWarnings(pegas::haplotype(dna))
      net <- pegas::haploNet(haps_obj)
      list(network = net,
           frequencies = stats::setNames(hap_map$Frequency,
                                          hap_map$NT_Sequence))
    }, error = function(e) NULL)
  }

  list(
    haplotypes    = haplo_df,
    haplotype_map = hap_map,
    n_haplotypes  = nrow(hap_map),
    n_individuals = nrow(haplo_df),
    is_single_snp = (n_snp == 1),
    network_data  = net_data,
    diversity     = list(Hd = Hd, Neff = Neff, PIC = PIC))
}


# ----------------------------------------------------------------------------
# Sub-population frequency
# ----------------------------------------------------------------------------

#' Compute haplotype frequencies by subpopulation
#'
#' @param hap_data Output of [hbr_define_haplotypes()].
#' @param pheno_data Phenotype table with an ID and subpopulation column.
#' @param subpop_col Subpopulation column name.
#'
#' @return A list with `long` and `wide`, or `NULL`.
#'
#' @noRd
hbr_subpop_freq <- function(hap_data, pheno_data, subpop_col) {

  if (is.null(subpop_col) || !subpop_col %in% names(pheno_data))
    return(NULL)

  mg <- merge(hap_data$haplotypes[, .(ID, Haplotype)],
              pheno_data[, c("ID", subpop_col), with = FALSE],
              by = "ID")
  mg <- mg[!is.na(get(subpop_col))]
  if (!nrow(mg)) return(NULL)

  ft <- mg[, .(N = .N), by = c("Haplotype", subpop_col)]
  ft[, Total := sum(N), by = subpop_col]
  ft[, Freq_pct := round(N / Total * 100, 1)]

  wide <- data.table::dcast(ft, Haplotype ~ get(subpop_col),
                             value.var = "Freq_pct", fill = 0)

  list(long = ft, wide = wide)
}
