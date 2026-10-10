# ============================================================================
# R/reference.r
# Reference object construction, remote fetch (MSU7, RAP-DB, Ensembl Plants),
# sequence QC, functional annotation, and gene characterization.
#
# ----------------------------------------------------------------------------
# Crop lookup
# ----------------------------------------------------------------------------

#' Default crop lookup table
#'
#' Built-in mapping from short crop names to Ensembl Plants species and
#' assembly identifiers, used when no cached lookup is available.
#'
#' @return A `data.table` with columns `crop`, `species`, `assembly`.
#'
#' @noRd
.hbr_crop_lookup_default <- function() {
  data.table::data.table(
    crop = c(
      "rice_japonica", "rice_indica", "maize", "wheat", "barley",
      "sorghum", "foxtail_millet", "brachypodium", "soybean",
      "common_bean", "peanut", "potato", "cassava", "sugar_beet",
      "brassica_rapa", "grape", "peach", "sunflower", "poplar",
      "eucalyptus", "arabidopsis"),
    species = c(
      "oryza_sativa", "oryza_indica", "zea_mays",
      "triticum_aestivum", "hordeum_vulgare",
      "sorghum_bicolor", "setaria_italica", "brachypodium_distachyon",
      "glycine_max", "phaseolus_vulgaris", "arachis_hypogaea",
      "solanum_tuberosum", "manihot_esculenta", "beta_vulgaris",
      "brassica_rapa", "vitis_vinifera", "prunus_persica",
      "helianthus_annuus", "populus_trichocarpa",
      "eucalyptus_grandis", "arabidopsis_thaliana"),
    assembly = c(
      "IRGSP-1.0", "ASM465v1", "Zm-B73-REFERENCE-NAM-5.0", "IWGSC",
      "MorexV3_pseudomolecules_assembly",
      "Sorghum_bicolor_NCBIv3", "Setaria_italica_v2.0",
      "Brachypodium_distachyon_v3.0", "Glycine_max_v2.1",
      "PhaVulg1_0", "arahy.Tifrunner.gnm2.J5K5", "SolTub_3.0",
      "M.esculenta_v8", "RefBeet-1.2.2", "Brapa_1.0",
      "ASM3070453v1", "Prunus_persica_NCBIv2",
      "HanXRQr2.0-SUNRISE", "Pop_tri_v4", "Egrandis1_0", "TAIR10"))
}

#' Common crops mapped to their Ensembl Plants species and assembly
#'
#' Returns the crop-to-species-to-assembly lookup table. If a cached
#' table exists under `tools::R_user_dir("haplobreed", "cache")`, it is
#' returned; otherwise the built-in default from
#' [`.hbr_crop_lookup_default()`] is used.
#' @return A `data.table` with columns `crop`, `species`, `assembly`.
#' @examples
#' crop <- hbr_crop_lookup()
#' head(crop)
#' crop$crop
#'
#' @export
hbr_crop_lookup <- function() {
  cached <- file.path(tools::R_user_dir("haplobreed", "cache"),
                      "crop_lookup.rds")
  if (file.exists(cached)) {
    out <- tryCatch(readRDS(cached), error = function(e) NULL)
    if (!is.null(out) && is.data.frame(out) && nrow(out) > 0L)
      return(data.table::as.data.table(out))
  }
  .hbr_crop_lookup_default()
}

#' Refresh the crop lookup table by querying Ensembl Plants
#'
#' Probes the Ensembl Plants FTP directory listing for each crop in the
#' default lookup table and replaces any stale assembly name with the
#' one actually present. Optionally writes the corrected table to the
#' package cache.
#'
#' @param release Ensembl Plants release number.
#' @param write_cache Write the corrected table to the package cache.
#' @param verbose Print per-crop progress.
#'
#' @return The corrected crop table, invisibly.
#'
#' @noRd
hbr_refresh_crop_lookup <- function(release = 63L, write_cache = TRUE,
                                     verbose = TRUE) {
  crops   <- .hbr_crop_lookup_default()
  n_fixed <- 0L

  for (i in seq_len(nrow(crops))) {
    files <- .hbr_ensembl_gff3_listing(crops$species[i], release)
    if (!length(files)) {
      if (verbose)
        cat(sprintf("%-20s no listing (species renamed?)\n",
                    crops$crop[i]))
      next
    }
    asm <- .hbr_assembly_from_gff3(files[1], crops$species[i], release)
    if (!identical(asm, crops$assembly[i])) {
      if (verbose)
        cat(sprintf("%-20s %s  ->  %s\n",
                    crops$crop[i], crops$assembly[i], asm))
      crops$assembly[i] <- asm
      n_fixed <- n_fixed + 1L
    } else if (verbose) {
      cat(sprintf("%-20s unchanged (%s)\n", crops$crop[i], asm))
    }
  }

  if (verbose) cat("\nFixed:", n_fixed, "of", nrow(crops), "rows\n")

  if (isTRUE(write_cache)) {
    cache <- file.path(tools::R_user_dir("haplobreed", "cache"),
                       "crop_lookup.rds")
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    saveRDS(crops, cache)
    if (verbose) cat("Saved corrected table:", cache, "\n")
  }

  invisible(crops)
}


# ----------------------------------------------------------------------------
# Reference cache directory
# ----------------------------------------------------------------------------

#' Resolve the reference cache directory
#'
#' Returns a valid writable cache directory, creating it if necessary.
#' Uses the supplied `cache` argument when given, otherwise
#' `tools::R_user_dir("haplobreed", "cache")`.
#'
#' @param cache Optional explicit cache directory.
#'
#' @return A single directory path.
#'
#' @noRd
.hbr_reference_cache_dir <- function(cache = NULL) {
  if (!is.null(cache) && nzchar(cache)) {
    if (!dir.exists(cache))
      dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    return(cache)
  }
  d <- tryCatch(
    tools::R_user_dir("haplobreed", "cache"),
    error = function(e) file.path(tempdir(), "haplobreed_cache"))
  if (!dir.exists(d))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


# ----------------------------------------------------------------------------
# Remote fetch
# ----------------------------------------------------------------------------

#' Current Ensembl Plants release number
#'
#' Scrapes the Ensembl Plants FTP index page for the current release
#' number. Falls back to 63L if the page cannot be read.
#'
#' @return Integer scalar.
#'
#' @noRd
.hbr_ensembl_current_release <- function() {
  url <- "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/plants/"
  out <- tryCatch({
    con <- url(url, "r")
    on.exit(close(con), add = TRUE)
    html <- paste(readLines(con, n = 30L, warn = FALSE), collapse = " ")
    m <- regmatches(html, regexpr("release-[0-9]+", html))
    if (length(m) == 1L && nzchar(m)) as.integer(sub("release-", "", m))
    else 63L
  }, error = function(e) 63L)
  out
}

#' Ensembl Plants URLs for one species and assembly
#'
#' @param species Ensembl species identifier (e.g. `"oryza_sativa"`).
#' @param assembly Assembly identifier (e.g. `"IRGSP-1.0"`).
#' @param release Ensembl Plants release number.
#'
#' @return A named list with `gff3`, `cds`, `pep`, `genome_fasta`.
#'
#' @noRd
.hbr_ensembl_urls <- function(species, assembly, release = 63L) {
  base <- sprintf(
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-%d",
    as.integer(release))
  sp_cap <- paste0(toupper(substr(species, 1, 1)),
                   substr(species, 2, nchar(species)))
  gff3_name    <- sprintf("%s.%s.%d.gff3.gz", sp_cap, assembly, release)
  fasta_prefix <- sprintf("%s.%s", sp_cap, assembly)
  list(
    gff3         = sprintf("%s/gff3/%s/%s", base, species, gff3_name),
    cds          = sprintf("%s/fasta/%s/cds/%s.cds.all.fa.gz",
                           base, species, fasta_prefix),
    pep          = sprintf("%s/fasta/%s/pep/%s.pep.all.fa.gz",
                           base, species, fasta_prefix),
    genome_fasta = sprintf("%s/fasta/%s/dna/%s.dna.toplevel.fa.gz",
                           base, species, fasta_prefix))
}

#' MSU7 download URLs
#'
#' @return A named list with `gff3`, `cds`, `pep`, `genome_fasta`.
#'
#' @noRd
.hbr_msu7_urls <- function() {
  base <- "https://rice.uga.edu/osa1r7_download"
  list(
    gff3         = paste0(base, "/osa1_r7.all_models.gff3.gz"),
    cds          = paste0(base, "/osa1_r7.all_models.cds.fa.gz"),
    pep          = paste0(base, "/osa1_r7.all_models.pep.fa.gz"),
    genome_fasta = paste0(base, "/osa1_r7.asm.fa.gz"))
}

#' RAP-DB download URLs
#'
#' @param release Ensembl Plants release number.
#'
#' @return A named list with `gff3`, `cds`, `pep`, `genome_fasta`.
#'
#' @noRd
.hbr_rapdb_urls <- function(release = 63L) {
  .hbr_ensembl_urls("oryza_sativa", "IRGSP-1.0", release)
}

#' Download a URL to the cache, with retries
#'
#' If the file already exists in `cache_dir`, returns it immediately.
#' Otherwise downloads it (three attempts, one-hour timeout, method
#' `libcurl` when available) and returns the cached path. Fails softly:
#' returns `NULL` when offline or after all retries fail.
#'
#' @param url URL to fetch.
#' @param cache_dir Cache directory.
#' @param verbose Print progress.
#'
#' @return Path to the cached file, or `NULL`.
#'
#' @noRd
.hbr_fetch_cached <- function(url, cache_dir, verbose = TRUE) {

  if (is.null(url) || !nzchar(url)) return(NULL)
  if (file.exists(url)) return(url)

  fname <- basename(sub("\\?.*$", "", url))
  dest  <- file.path(cache_dir, fname)

  if (file.exists(dest) && file.info(dest)$size > 0) {
    if (verbose)
      cat("    cache hit:", fname,
          sprintf("(%.1f MB)\n", file.info(dest)$size / 1e6))
    return(dest)
  }

  if (verbose) cat("    downloading:", fname, "...\n")

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(3600L, old_timeout))

  method <- if (isTRUE(capabilities("libcurl"))) "libcurl" else "auto"

  ok <- FALSE
  for (attempt in seq_len(3L)) {
    ok <- tryCatch({
      utils::download.file(url, destfile = dest, mode = "wb",
                           quiet = !verbose, method = method)
      file.exists(dest) && file.info(dest)$size > 0
    }, error = function(e) {
      if (verbose)
        cat("      attempt ", attempt, " failed: ",
            conditionMessage(e), "\n", sep = "")
      if (file.exists(dest)) unlink(dest)
      FALSE
    })
    if (ok) break
    if (attempt < 3L) Sys.sleep(5L)
  }

  if (!ok) {
    message("    could not download ", fname,
            " after 3 attempts; continuing without it.")
    return(NULL)
  }

  if (verbose)
    cat(sprintf("    saved: %s (%.1f MB)\n",
                fname, file.info(dest)$size / 1e6))
  dest
}


# ----------------------------------------------------------------------------
# Reference path resolution
# ----------------------------------------------------------------------------

#' Resolve reference file paths
#'
#' @param gff3,cds,pep,genome_fasta Optional local file paths.
#' @param source One of `"local"`, `"msu7"`, `"rapdb"`, `"ensembl"`.
#' @param crop_species,crop_assembly Ensembl Plants species and assembly.
#' @param ensembl_release Ensembl Plants release.
#' @param cache_dir Optional cache directory.
#' @param verbose Print progress.
#'
#' @return A list of resolved paths.
#'
#' @noRd
.hbr_resolve_reference_paths <- function(gff3, cds, pep, genome_fasta,
                                         source = c("local", "msu7",
                                                    "rapdb", "ensembl"),
                                         crop_species    = NULL,
                                         crop_assembly   = NULL,
                                         ensembl_release = NULL,
                                         cache_dir       = NULL,
                                         verbose         = TRUE) {

  source <- match.arg(source)

  if (source == "local") {
    return(list(gff3 = gff3, cds = cds, pep = pep,
                genome_fasta = genome_fasta,
                cache_dir = NULL, source = "local", downloaded = FALSE,
                crop_species = NULL, crop_assembly = NULL,
                ensembl_release = NULL))
  }

  cache_dir <- .hbr_reference_cache_dir(cache_dir)

  if (source == "ensembl") {
    if (is.null(crop_species) || is.null(crop_assembly)) {
      lookup <- data.table::as.data.table(hbr_crop_lookup())
      if (!is.null(crop_species)) {
        hit <- lookup[species == crop_species]
        if (nrow(hit) && is.null(crop_assembly))
          crop_assembly <- as.character(hit$assembly[1])
      }
      if (is.null(crop_species) || is.null(crop_assembly))
        stop("For reference_source = 'ensembl', supply `crop_species` ",
             "and `crop_assembly`.", call. = FALSE)
    }
    if (is.null(ensembl_release))
      ensembl_release <- .hbr_ensembl_current_release()
    urls <- .hbr_ensembl_urls(crop_species, crop_assembly,
                              release = ensembl_release)
  } else if (source == "msu7") {
    urls <- .hbr_msu7_urls()
  } else if (source == "rapdb") {
    if (is.null(ensembl_release))
      ensembl_release <- .hbr_ensembl_current_release()
    urls <- .hbr_rapdb_urls(ensembl_release)
  }

  if (verbose)
    cat("    reference_source =", source,
        if (!is.null(crop_species))
          paste0("| ", crop_species, " / ", crop_assembly) else "",
        if (!is.null(ensembl_release))
          paste0("| release ", ensembl_release) else "",
        "\n    cache:", cache_dir, "\n")

  downloaded_any <- FALSE
  fetch_one <- function(user_path, url) {
    if (!is.null(user_path) && nzchar(user_path) && file.exists(user_path))
      return(user_path)
    downloaded_any <<- TRUE
    .hbr_fetch_cached(url, cache_dir, verbose = verbose)
  }

  list(
    gff3            = fetch_one(gff3,         urls$gff3),
    cds             = fetch_one(cds,          urls$cds),
    pep             = fetch_one(pep,          urls$pep),
    genome_fasta    = fetch_one(genome_fasta, urls$genome_fasta),
    cache_dir       = cache_dir,
    source          = source,
    downloaded      = downloaded_any,
    crop_species    = crop_species,
    crop_assembly   = crop_assembly,
    ensembl_release = ensembl_release)
}

#' Extract a `.tar.gz` GFF3 archive
#'
#' @param path Path to a file.
#' @param cache_dir Cache directory for the extracted file.
#'
#' @return Path to the extracted `.gff3` file, or `path`.
#'
#' @noRd
.hbr_maybe_untar_gff3 <- function(path, cache_dir) {
  if (is.null(path) || !file.exists(path)) return(path)
  if (!grepl("\\.tar\\.gz$", path, ignore.case = TRUE)) return(path)

  base   <- sub("\\.tar\\.gz$", "", basename(path))
  target <- file.path(cache_dir, paste0(base, ".gff3"))
  if (file.exists(target) && file.info(target)$size > 0) return(target)

  tmp <- file.path(cache_dir, paste0(base, "_untar"))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  utils::untar(path, exdir = tmp)

  hit <- list.files(tmp, pattern = "\\.gff3(\\.gz)?$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(hit))
    stop("No .gff3 found inside ", basename(path), call. = FALSE)

  file.copy(hit[1], target, overwrite = TRUE)
  unlink(tmp, recursive = TRUE)
  target
}

#' List GFF3 files in an Ensembl Plants species directory
#'
#' @param species Ensembl species identifier.
#' @param release Ensembl Plants release number.
#'
#' @return Character vector of GFF3 file names (may be empty).
#'
#' @noRd
.hbr_ensembl_gff3_listing <- function(species, release = 63L) {

  url <- sprintf(
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-%d/gff3/%s/",
    as.integer(release), species)

  html <- tryCatch({
    con <- url(url, "r")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    txt <- readLines(con, warn = FALSE)
    if (!length(txt)) "" else paste(txt, collapse = "\n")
  }, error = function(e) "")

  if (!is.character(html) || !nzchar(html)) return(character(0))

  m <- gregexpr('href="[^"]+\\.gff3\\.gz"', html)
  hits <- unlist(regmatches(html, m), use.names = FALSE)
  if (!length(hits)) return(character(0))

  files <- sub('href="', "", hits, fixed = TRUE)
  files <- sub('"$', "", files)
  files <- sub("^.*/", "", files)
  files <- unique(files)

  files <- files[!grepl("\\.chr\\.gff3\\.gz$", files)]
  files <- files[!grepl("\\.chromosome\\.[^\\.]+\\.gff3\\.gz$", files)]

  prefer <- sprintf("^[A-Z][a-z]+_[a-z_]+\\.%d\\.gff3\\.gz$", release)
  toplevel <- files[grepl(prefer, files)]
  if (length(toplevel)) return(toplevel)

  files
}

#' Extract an assembly name from a GFF3 filename
#'
#' @param filename GFF3 file name.
#' @param species Ensembl species identifier.
#' @param release Ensembl Plants release number.
#'
#' @return Character scalar.
#'
#' @noRd
.hbr_assembly_from_gff3 <- function(filename, species, release = 63L) {
  base <- sub("\\.gff3\\.gz$", "", filename)
  sp_cap <- paste0(toupper(substr(species, 1, 1)),
                   substr(species, 2, nchar(species)))
  base <- sub(paste0("^", sp_cap, "\\."), "", base)
  base <- sub(paste0("\\.", release, "$"), "", base)
  base
}


# ----------------------------------------------------------------------------
# Reference loading
# ----------------------------------------------------------------------------

#' Load a standardized reference object for HaploBreedR
#'
#' Parses a GFF3 file into genes, transcripts, CDS, and UTR features and
#' optionally attaches CDS, peptide, and genome FASTA sequences.
#' Chromosome names on the reference are standardized via
#' [hbr_standardize_chr()] so they match the genotype convention applied
#' by [hbr_normalize_genotype()].
#'
#' All console output is gated by `verbose`. `TRUE` (default) reports the
#' paths being read, gene/transcript counts, and the sequence-QC status;
#' `FALSE` is silent (errors and warnings still fire).
#'
#' @param genome Genome build label.
#' @param annotation One of `"MSU7"`, `"RAPDB"`, `"CUSTOM"`.
#' @param gff3 Path to a GFF3 file.
#' @param cds Optional CDS FASTA path.
#' @param pep Optional peptide FASTA path.
#' @param genome_fasta Optional whole-genome FASTA path.
#' @param verbose Print progress. `TRUE` (default) reports loading steps;
#'   `FALSE` silences everything except errors.
#'
#' @return A list of class `hbr_reference`.
#'
#' @noRd
hbr_load_reference <- function(genome,
                               annotation = c("MSU7", "RAPDB", "CUSTOM"),
                               gff3,
                               cds = NULL,
                               pep = NULL,
                               genome_fasta = NULL,
                               verbose = TRUE) {

  annotation <- match.arg(annotation)
  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  if (missing(gff3) || is.null(gff3) || !file.exists(gff3))
    stop("hbr_load_reference(): `gff3` is required and must exist. Got: ",
         if (missing(gff3)) "<missing>" else gff3, call. = FALSE)
  if (!is.null(cds) && !file.exists(cds))
    stop("hbr_load_reference(): `cds` file not found: ", cds, call. = FALSE)
  if (!is.null(pep) && !file.exists(pep))
    stop("hbr_load_reference(): `pep` file not found: ", pep, call. = FALSE)
  if (!is.null(genome_fasta) && !file.exists(genome_fasta))
    stop("hbr_load_reference(): `genome_fasta` file not found: ",
         genome_fasta, call. = FALSE)

  vcat("Loading reference:", genome, "/", annotation, "\n")
  vcat("  GFF3:", gff3, "\n")

  gff_raw <- rtracklayer::import(gff3)
  gff_df  <- data.table::as.data.table(as.data.frame(gff_raw))
  data.table::setnames(gff_df, "seqnames", "Chr")
  gff_df[, Chr    := as.character(Chr)]
  gff_df[, type   := as.character(type)]
  gff_df[, strand := as.character(strand)]

  first_chr    <- hbr_first_chr
  strip_prefix <- hbr_strip_prefix

  # ---- genes ----
  genes_raw <- gff_df[type == "gene"]
  if (!nrow(genes_raw))
    stop("hbr_load_reference(): no 'gene' type features found in GFF3.",
         call. = FALSE)
  if ("gene_id" %in% names(genes_raw)) genes_raw[, gene_id := NULL]
  genes_raw[, gene_id := strip_prefix(first_chr(ID))]
  desc_col <- intersect(c("description", "Description", "Note", "note"),
                        names(genes_raw))
  genes_raw[, gene_function := if (length(desc_col) > 0L)
    first_chr(get(desc_col[1])) else NA_character_]
  genes <- unique(genes_raw[!is.na(gene_id) & nzchar(gene_id),
                            .(gene_id, Chr, Start = start, End = end,
                              Strand = strand, gene_function)])

  # ---- transcripts ----
  tx_raw <- gff_df[type %in% c("mRNA", "transcript")]
  if (!nrow(tx_raw))
    stop("hbr_load_reference(): no mRNA/transcript features found in GFF3.",
         call. = FALSE)
  if ("gene_id" %in% names(tx_raw)) tx_raw[, gene_id := NULL]
  tx_raw[, transcript_id := strip_prefix(first_chr(ID))]
  tx_raw[, gene_id       := strip_prefix(first_chr(Parent))]
  transcripts <- unique(tx_raw[!is.na(transcript_id) & !is.na(gene_id),
                               .(transcript_id, gene_id, Chr, Strand = strand)])

  # ---- CDS ----
  cds_dt <- gff_df[type == "CDS"]
  if (nrow(cds_dt)) {
    if ("gene_id" %in% names(cds_dt)) cds_dt[, gene_id := NULL]
    cds_dt[, transcript_id := strip_prefix(first_chr(Parent))]
    cds_dt <- merge(cds_dt,
                    transcripts[, .(transcript_id, gene_id)],
                    by = "transcript_id", all.x = TRUE)
    cds_dt <- cds_dt[!is.na(gene_id)]
  } else {
    cds_dt <- NULL
  }

  # ---- UTR ----
  type_upper  <- toupper(gff_df$type)
  utr5_labels <- toupper(c("five_prime_UTR", "5UTR", "5'UTR"))
  utr3_labels <- toupper(c("three_prime_UTR", "3UTR", "3'UTR"))
  utr_dt <- gff_df[type_upper %in% c(utr5_labels, utr3_labels)]
  if (nrow(utr_dt)) {
    if ("gene_id" %in% names(utr_dt)) utr_dt[, gene_id := NULL]
    utr_dt[, UTR_Type := ifelse(toupper(type) %in% utr5_labels,
                                "5UTR", "3UTR")]
    utr_dt[, transcript_id := strip_prefix(first_chr(Parent))]
    utr_dt <- merge(utr_dt,
                    transcripts[, .(transcript_id, gene_id)],
                    by = "transcript_id", all.x = TRUE)
    utr_dt <- utr_dt[!is.na(gene_id)]
  } else {
    utr_dt <- NULL
  }

  # ---- Optional sequence resources ----
  cds_seqs <- pep_seqs <- genome_seqs <- NULL
  if (!is.null(cds)) {
    vcat("  CDS FASTA:", cds, "\n")
    cds_seqs <- Biostrings::readDNAStringSet(cds)
    names(cds_seqs) <- gsub(
      "^transcript:", "",
      vapply(strsplit(names(cds_seqs), " "),
             function(x) x[1],
             character(1)))
  }
  if (!is.null(pep)) {
    vcat("  Peptide FASTA:", pep, "\n")
    pep_seqs <- Biostrings::readAAStringSet(pep)
    names(pep_seqs) <- gsub(
      "^peptide:|^protein:", "",
      vapply(strsplit(names(pep_seqs), " "),
             function(x) x[1],
             character(1)))
  }
  if (!is.null(genome_fasta)) {
    vcat("  Genome FASTA:", genome_fasta, "\n")
    genome_seqs <- Biostrings::readDNAStringSet(genome_fasta)
    names(genome_seqs) <- vapply(strsplit(names(genome_seqs), " "),
                                 function(x) x[1],
                                 character(1))
  }

  # ---- Chromosome lengths ----
  if (!is.null(genome_seqs)) {
    chromosomes <- data.table::data.table(
      Chr = names(genome_seqs),
      Length_Approx = Biostrings::width(genome_seqs),
      Length_Source = "genome_fasta")
  } else {
    chromosomes <- gff_df[, .(Length_Approx = max(end, na.rm = TRUE)), by = Chr]
    chromosomes[, Length_Source := "gff3_max_feature"]
    vcat("  Note: no genome_fasta supplied -- chromosome lengths are ",
         "approximated from the GFF3's\n",
         "  own maximum feature coordinate per chromosome.\n", sep = "")
  }

  # ---- Standardize chromosome names on every slot ----
  chromosomes$Chr  <- hbr_standardize_chr(chromosomes$Chr)
  genes$Chr        <- hbr_standardize_chr(genes$Chr)
  transcripts$Chr  <- hbr_standardize_chr(transcripts$Chr)
  if (!is.null(cds_dt)) cds_dt$Chr <- hbr_standardize_chr(cds_dt$Chr)
  if (!is.null(utr_dt)) utr_dt$Chr <- hbr_standardize_chr(utr_dt$Chr)
  if (!is.null(genome_seqs))
    names(genome_seqs) <- hbr_standardize_chr(names(genome_seqs))

  ref <- list(
    genome = genome, annotation = annotation,
    gff3_path = gff3, cds_path = cds, pep_path = pep,
    genome_fasta_path = genome_fasta,
    chromosomes = chromosomes[], genes = genes, transcripts = transcripts,
    cds = cds_dt, utr = utr_dt,
    cds_seqs = cds_seqs, pep_seqs = pep_seqs, genome_seqs = genome_seqs,
    has_cds = !is.null(cds_seqs),
    has_pep = !is.null(pep_seqs),
    has_genome_fasta = !is.null(genome_seqs))
  class(ref) <- c("hbr_reference", "list")

  vcat("  Genes:", nrow(genes), "| Transcripts:", nrow(transcripts),
       "| CDS features:", if (!is.null(cds_dt)) nrow(cds_dt) else 0L,
       "| UTR features:", if (!is.null(utr_dt)) nrow(utr_dt) else 0L, "\n")

  ref$validation <- hbr_validate_sequence_resources(ref, verbose = verbose)
  vcat("  Sequence-resource QC:", ref$validation$status, "\n")

  ref
}

#' Print an hbr_reference object
#'
#' @param x An object of class `hbr_reference`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' ref <- structure(
#'   list(genome = "TOY", annotation = "CUSTOM",
#'        genes       = data.frame(gene_id = "g1"),
#'        transcripts = data.frame(gene_id = "g1",
#'                                 transcript_id = "g1.1"),
#'        has_cds = FALSE, has_pep = FALSE, has_genome_fasta = FALSE,
#'        validation = list(status = "PASS")),
#'   class = "hbr_reference")
#' print(ref)
#'
#' @export
print.hbr_reference <- function(x, ...) {
  cat("<hbr_reference>", x$genome, "/", x$annotation, "\n")
  cat("  Genes:", nrow(x$genes), " Transcripts:", nrow(x$transcripts), "\n")
  cat("  CDS available:", x$has_cds, " PEP available:", x$has_pep,
      " Genome FASTA available:", x$has_genome_fasta, "\n")
  cat("  Sequence-resource QC:", x$validation$status, "\n")
  invisible(x)
}


# ----------------------------------------------------------------------------
# Sequence-resource QC
# ----------------------------------------------------------------------------

#' Validate consistency between GFF3 transcripts, CDS FASTA, and peptide FASTA
#'
#' Cross-checks the reference's transcript table against its CDS and
#' peptide sequence sets. By default the check is length-only: each CDS
#' length is required to be a multiple of three, and the CDS length
#' divided by three (with and without a stop codon) is compared to the
#' peptide length. When `check_translation = TRUE`, the CDS sequences
#' are bulk-translated and the actual translation length is compared to
#' the peptide length.
#'
#' @param ref A `hbr_reference` object from [hbr_load_reference()].
#' @param check_translation If `FALSE` (default), QC is length-only.
#' @param verbose Print progress.
#'
#' @return A list of class `hbr_sequence_qc`.
#'
#' @noRd
hbr_validate_sequence_resources <- function(ref,
                                            check_translation = FALSE,
                                            verbose = TRUE) {

  if (!inherits(ref, "hbr_reference"))
    stop("hbr_validate_sequence_resources(): `ref` must be a ",
         "hbr_reference object.", call. = FALSE)

  t_start <- Sys.time()
  tx      <- data.table::copy(ref$transcripts)
  n_tx    <- nrow(tx)

  if (!ref$has_cds && !ref$has_pep) {
    out <- list(
      status  = "NO SEQUENCE DATA",
      summary = data.table::data.table(
        GFF3_transcripts     = n_tx,
        CDS_sequences        = 0L,
        PEP_sequences        = 0L,
        Complete_CDS_PEP     = 0L,
        CDS_only             = 0L,
        PEP_only             = 0L,
        Translation_mismatch = 0L,
        Invalid_CDS_length   = 0L),
      per_transcript = tx[, .(transcript_id, gene_id)])
    class(out) <- c("hbr_sequence_qc", "list")
    return(out)
  }

  idx_of <- function(xset) {
    if (is.null(xset)) return(integer(0))
    stats::setNames(seq_along(xset), names(xset))
  }
  cds_idx <- idx_of(ref$cds_seqs)
  pep_idx <- idx_of(ref$pep_seqs)

  tx[, has_cds_fasta := FALSE]
  tx[, has_pep_fasta := FALSE]
  if (ref$has_cds) {
    ci <- unname(cds_idx[tx$transcript_id])
    tx[, has_cds_fasta := !is.na(ci)]
  }
  if (ref$has_pep) {
    pi <- unname(pep_idx[tx$transcript_id])
    tx[, has_pep_fasta := !is.na(pi)]
  }

  tx[, cds_length          := NA_integer_]
  tx[, cds_length_ok       := NA]
  tx[, translated_length   := NA_integer_]
  tx[, pep_length          := NA_integer_]
  tx[, translation_matches := NA]

  if (ref$has_cds) {
    idx_cds <- which(tx$has_cds_fasta)
    if (verbose)
      cat(sprintf("  [seq-QC] checking %d CDS sequences (%s)...\n",
                  length(idx_cds),
                  if (check_translation) "length + translation"
                  else                   "length only"))

    tx_ids  <- tx$transcript_id[idx_cds]
    cds_int <- unname(cds_idx[tx_ids])

    cds_subset <- ref$cds_seqs[cds_int]
    widths_v   <- Biostrings::width(cds_subset)

    tx$cds_length[idx_cds]    <- widths_v
    tx$cds_length_ok[idx_cds] <- (widths_v %% 3L == 0L)

    pep_widths_v <- rep(NA_integer_, length(idx_cds))
    if (ref$has_pep) {
      pep_int <- unname(pep_idx[tx_ids])
      has_pep <- !is.na(pep_int)
      if (any(has_pep)) {
        pep_widths_v[has_pep] <-
          Biostrings::width(ref$pep_seqs[pep_int[has_pep]])
      }
      tx$pep_length[idx_cds] <- pep_widths_v
    }

    core_len_v <- rep(NA_integer_, length(idx_cds))
    ok_pos     <- which(widths_v %% 3L == 0L)

    if (check_translation && length(ok_pos)) {
      t0_bulk <- Sys.time()
      if (verbose)
        cat(sprintf("  [seq-QC] bulk-translating %d in-frame CDS...\n",
                    length(ok_pos)))

      translated_str <- tryCatch({
        ts <- Biostrings::translate(cds_subset[ok_pos])
        v  <- as.character(ts)
        v[!nzchar(v) | v == "<NA>"] <- NA_character_
        v
      }, error = function(e) NULL)

      if (is.null(translated_str)) {
        translated_str <- character(length(ok_pos))
        for (m in seq_along(ok_pos)) {
          seq_m <- cds_subset[[ok_pos[m]]]
          v <- tryCatch(
            as.character(Biostrings::translate(seq_m)),
            error = function(e) tryCatch(
              as.character(Biostrings::translate(
                seq_m, if.fuzzy.codon = "solve")),
              error = function(e) NA_character_))
          translated_str[m] <- v
        }
      }

      if (verbose)
        cat(sprintf("  [seq-QC] translation done in %.1fs\n",
                    as.numeric(difftime(Sys.time(), t0_bulk,
                                        units = "secs"))))

      ends_stop                 <- grepl("\\*$", translated_str)
      cl                        <- nchar(translated_str) -
        as.integer(ends_stop)
      cl[is.na(translated_str)] <- NA_integer_
      core_len_v[ok_pos]        <- cl
      tx$translated_length[idx_cds[ok_pos]] <- cl

      if (ref$has_pep) {
        has_pep_ok <- !is.na(pep_widths_v[ok_pos])
        tx$translation_matches[idx_cds[ok_pos][has_pep_ok]] <-
          (cl[has_pep_ok] == pep_widths_v[ok_pos][has_pep_ok])
      }
    }

    if (!check_translation && ref$has_pep) {
      with_stop    <- widths_v %/% 3L - 1L
      without_stop <- widths_v %/% 3L
      has_pep_ok   <- !is.na(pep_widths_v)
      actual_pep   <- pep_widths_v[has_pep_ok]
      tx$translation_matches[idx_cds[has_pep_ok]] <-
        (actual_pep == with_stop[has_pep_ok]) |
        (actual_pep == without_stop[has_pep_ok])
    }
  }

  if (ref$has_pep) {
    idx_pep_only <- which(tx$has_pep_fasta & !tx$has_cds_fasta)
    if (length(idx_pep_only)) {
      tx_ids  <- tx$transcript_id[idx_pep_only]
      pep_int <- unname(pep_idx[tx_ids])
      tx$pep_length[idx_pep_only] <-
        Biostrings::width(ref$pep_seqs[pep_int])
    }
  }

  n_complete <- sum(tx$has_cds_fasta & tx$has_pep_fasta, na.rm = TRUE)
  n_cds_only <- sum(tx$has_cds_fasta & !tx$has_pep_fasta, na.rm = TRUE)
  n_pep_only <- sum(tx$has_pep_fasta & !tx$has_cds_fasta, na.rm = TRUE)
  n_mismatch <- sum(!is.na(tx$translation_matches) &
                      !tx$translation_matches)
  n_bad_len  <- sum(!is.na(tx$cds_length_ok) & !tx$cds_length_ok)

  summary_dt <- data.table::data.table(
    GFF3_transcripts     = n_tx,
    CDS_sequences        = if (ref$has_cds) length(ref$cds_seqs) else 0L,
    PEP_sequences        = if (ref$has_pep) length(ref$pep_seqs) else 0L,
    Complete_CDS_PEP     = n_complete,
    CDS_only             = n_cds_only,
    PEP_only             = n_pep_only,
    Translation_mismatch = n_mismatch,
    Invalid_CDS_length   = n_bad_len)

  status <- if (n_mismatch == 0L && n_bad_len == 0L) "PASS"
  else "PASS WITH WARNINGS"

  if (verbose)
    cat(sprintf("  [seq-QC] finished in %.1fs  (status: %s)\n",
                as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                status))

  out <- list(status         = status,
              summary        = summary_dt,
              per_transcript = tx[])
  class(out) <- c("hbr_sequence_qc", "list")
  out
}

#' Print a sequence-resource QC report
#'
#' @param x An object of class `hbr_sequence_qc`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' x <- structure(
#'   list(status = "PASS",
#'        summary = data.table::data.table(
#'          GFF3_transcripts = 0L, CDS_sequences = 0L,
#'          PEP_sequences = 0L, Complete_CDS_PEP = 0L,
#'          CDS_only = 0L, PEP_only = 0L,
#'          Translation_mismatch = 0L, Invalid_CDS_length = 0L)),
#'   class = "hbr_sequence_qc")
#' print(x)
#'
#' @export
print.hbr_sequence_qc <- function(x, ...) {
  cat("SEQUENCE RESOURCE QC\n====================\n\n")
  s <- x$summary
  cat(sprintf("GFF3 transcripts       : %d\n", s$GFF3_transcripts))
  cat(sprintf("CDS sequences          : %d\n", s$CDS_sequences))
  cat(sprintf("PEP sequences          : %d\n\n", s$PEP_sequences))
  cat(sprintf("Complete CDS+PEP       : %d\n", s$Complete_CDS_PEP))
  cat(sprintf("CDS only               : %d\n", s$CDS_only))
  cat(sprintf("PEP only               : %d\n", s$PEP_only))
  cat(sprintf("Translation mismatch   : %d\n", s$Translation_mismatch))
  cat(sprintf("Invalid CDS length     : %d\n\n", s$Invalid_CDS_length))
  cat("STATUS:", x$status, "\n")
  invisible(x)
}


# ----------------------------------------------------------------------------
# Annotation paths and loading
# ----------------------------------------------------------------------------

#' Annotation file URLs
#'
#' @param source One of `"msu7"`, `"rapdb"`.
#'
#' @return A named list of URLs.
#'
#' @noRd
.hbr_annotation_urls <- function(source = c("msu7", "rapdb")) {
  source <- match.arg(source)

  if (source == "msu7") {
    base <- "https://rice.uga.edu/osa1r7_download"
    return(list(
      functional   = paste0(base, "/osa1_r7.all_models.functional_annotation.txt.gz"),
      goslim       = paste0(base, "/osa1_r7.all_models.GOSlim.txt.gz"),
      interpro     = paste0(base, "/osa1_r7.all_models.interpro.txt.gz"),
      pfam         = paste0(base, "/osa1_r7.all_models.pfam.txt.gz"),
      locus_brief  = paste0(base, "/osa1_r7.locus_brief_info.txt.gz"),
      expression   = paste0(base, "/osa1_r7.all_models.gene_exp_matrix.txt.gz"),
      coexpression = paste0(base, "/osa1_r7.all_models.coexpression_modules.txt.gz"),
      module_list  = paste0(base, "/osa1_r7.coexpression_module_list.txt.gz")))
  }

  base <- "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1"
  list(
    functional   = paste0(base, "/IRGSP-1.0_representative_function.txt.gz"),
    goslim       = paste0(base, "/IRGSP-1.0_representative_GO.txt.gz"),
    interpro     = paste0(base, "/IRGSP-1.0_representative_interpro.txt.gz"),
    pfam         = paste0(base, "/IRGSP-1.0_representative_pfam.txt.gz"),
    locus_brief  = NULL,
    expression   = NULL,
    coexpression = NULL,
    module_list  = NULL)
}

#' Resolve annotation file paths, downloading when necessary
#'
#' @param functional,goslim,interpro,pfam,locus_brief,expression,coexpression,module_list
#'   User-supplied paths.
#' @param source One of `"msu7"`, `"rapdb"`.
#' @param cache_dir Cache directory.
#' @param verbose Print progress.
#'
#' @return A named list of local file paths.
#'
#' @noRd
hbr_resolve_annotation_paths <- function(
    functional   = NULL,
    goslim       = NULL,
    interpro     = NULL,
    pfam         = NULL,
    locus_brief  = NULL,
    expression   = NULL,
    coexpression = NULL,
    module_list  = NULL,
    source       = c("msu7", "rapdb"),
    cache_dir    = NULL,
    verbose      = TRUE) {

  source    <- match.arg(source)
  cache_dir <- .hbr_reference_cache_dir(cache_dir)
  urls      <- .hbr_annotation_urls(source)

  fetch_one <- function(user_path, url) {
    if (!is.null(user_path) && nzchar(user_path) && file.exists(user_path))
      return(user_path)
    if (is.null(url)) return(NULL)
    if (verbose) cat("    fetching:", basename(url), "...\n")
    tryCatch(.hbr_fetch_cached(url, cache_dir, verbose = verbose),
             error = function(e) {
               if (verbose)
                 cat("      failed:", conditionMessage(e), "\n")
               NULL
             })
  }

  list(
    functional   = fetch_one(functional,   urls$functional),
    goslim       = fetch_one(goslim,       urls$goslim),
    interpro     = fetch_one(interpro,     urls$interpro),
    pfam         = fetch_one(pfam,         urls$pfam),
    locus_brief  = fetch_one(locus_brief,  urls$locus_brief),
    expression   = fetch_one(expression,   urls$expression),
    coexpression = fetch_one(coexpression, urls$coexpression),
    module_list  = fetch_one(module_list,  urls$module_list),
    cache_dir    = cache_dir,
    source       = source)
}


# ----------------------------------------------------------------------------
# Functional annotation loading
# ----------------------------------------------------------------------------

#' Read Oryza sativa functional annotation tables
#'
#' @param functional_annotation_file,goslim_file,interpro_file,pfam_file
#'   Local paths.
#' @param locus_brief_file,expression_file,coexpression_file,module_list_file
#'   Optional files.
#' @param source One of `"msu7"`, `"rapdb"`.
#' @param cache_dir Cache directory for downloads.
#' @param verbose Print progress.
#'
#' @return A list with the loaded tables.
#'
#' @noRd
hbr_read_osa1_r7_annotations <- function(
    functional_annotation_file = NULL,
    goslim_file                = NULL,
    interpro_file              = NULL,
    pfam_file                  = NULL,
    locus_brief_file           = NULL,
    expression_file            = NULL,
    coexpression_file          = NULL,
    module_list_file           = NULL,
    source                     = c("msu7", "rapdb"),
    cache_dir                  = NULL,
    verbose                    = TRUE) {

  source <- match.arg(source)

  paths <- hbr_resolve_annotation_paths(
    functional   = functional_annotation_file,
    goslim       = goslim_file,
    interpro     = interpro_file,
    pfam         = pfam_file,
    locus_brief  = locus_brief_file,
    expression   = expression_file,
    coexpression = coexpression_file,
    module_list  = module_list_file,
    source       = source,
    cache_dir    = cache_dir,
    verbose      = verbose)

  read_tsv <- function(path, header = TRUE) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    dt <- data.table::fread(path, sep = "\t", header = FALSE,
                            data.table = TRUE, showProgress = FALSE,
                            na.strings = c("", "NA", "N/A"), quote = "")
    if (!nrow(dt) || !isTRUE(header)) return(dt)
    first <- as.character(dt[[1]][1])
    if (!isTRUE(grepl("^(LOC_)?Os[0-9]{2}[gt][0-9]+", first))) {
      nm <- trimws(as.character(unlist(dt[1], use.names = FALSE)))
      bad <- is.na(nm) | nm == ""
      nm[bad] <- paste0("V", which(bad))
      nm <- make.unique(nm)
      dt <- dt[-1]
      data.table::setnames(dt, nm)
      for (j in names(dt))
        data.table::set(dt, j = j,
                        value = utils::type.convert(dt[[j]], as.is = TRUE))
    }
    dt
  }

  read_maybe_gz <- function(path, header = TRUE) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
      if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
      data.table::as.data.table(readxl::read_excel(path))
    } else {
      read_tsv(path, header = header)
    }
  }

  if (verbose) message("[1/8] Functional annotation...")
  functional <- read_tsv(paths$functional)
  if (!is.null(functional) && nrow(functional)) {
    data.table::setnames(functional, old = names(functional)[1],
                         new = "model", skip_absent = TRUE)
    if (!"annotation" %in% names(functional) && ncol(functional) >= 2)
      data.table::setnames(functional, old = names(functional)[2],
                           new = "annotation", skip_absent = TRUE)
    functional[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
  }

  if (verbose) message("[2/8] GOSlim annotation...")
  goslim <- read_tsv(paths$goslim, header = FALSE)
  if (!is.null(goslim) && nrow(goslim)) {
    if (ncol(goslim) >= 6) {
      data.table::setnames(goslim, old = names(goslim)[seq_len(6L)],
                           new = c("model", "GO_ID", "GO_Aspect",
                                   "GO_Term", "Evidence", "Reference"))
      goslim[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
    } else goslim <- NULL
  }

  if (verbose) message("[3/8] InterPro annotation...")
  interpro <- read_tsv(paths$interpro)
  if (!is.null(interpro) && nrow(interpro)) {
    data.table::setnames(interpro, old = names(interpro)[1],
                         new = "model", skip_absent = TRUE)
    interpro[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
  }

  if (verbose) message("[4/8] PFAM annotation...")
  pfam <- read_tsv(paths$pfam)
  if (!is.null(pfam) && nrow(pfam)) {
    data.table::setnames(pfam, old = names(pfam)[1],
                         new = "model", skip_absent = TRUE)
    pfam[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
  }

  if (verbose) message("[5/8] Locus brief info...")
  locus_brief <- read_tsv(paths$locus_brief)
  if (!is.null(locus_brief) && nrow(locus_brief)) {
    data.table::setnames(locus_brief, old = names(locus_brief)[1],
                         new = "model", skip_absent = TRUE)
    locus_brief[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
  }

  if (verbose) message("[6/8] Expression matrix...")
  expression <- read_maybe_gz(paths$expression)
  if (!is.null(expression) && nrow(expression)) {
    id_col <- intersect(c("transcript_id", "model", "gene_id",
                          "Gene_ID", "ID"), names(expression))[1]
    if (!is.na(id_col)) {
      if (id_col != "model")
        data.table::setnames(expression, id_col, "model",
                             skip_absent = TRUE)
      expression[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
    }
  }

  if (verbose) message("[7/8] Coexpression modules...")
  coexpression <- read_maybe_gz(paths$coexpression)
  if (!is.null(coexpression) && nrow(coexpression)) {
    id_col <- intersect(c("gene_model", "model", "gene_id",
                          "Gene_ID", "transcript_id"),
                        names(coexpression))[1]
    if (!is.na(id_col)) {
      if (id_col != "model")
        data.table::setnames(coexpression, id_col, "model",
                             skip_absent = TRUE)
      coexpression[, Gene_ID := hbr_normalize_osa1_gene_id(model)]
    }
  }

  if (verbose) message("[8/8] Module catalog...")
  module_list <- read_maybe_gz(paths$module_list)

  if (verbose) {
    message("")
    message("=== OSA1 R7 ANNOTATION LOADED ===")
    for (nm in c("functional", "goslim", "interpro", "pfam",
                 "locus_brief", "expression", "coexpression",
                 "module_list")) {
      obj <- get(nm)
      if (!is.null(obj))
        message(sprintf("  %-14s %d rows x %d cols",
                        nm, nrow(obj), ncol(obj)))
    }
  }

  list(functional   = functional,
       goslim       = goslim,
       interpro     = interpro,
       pfam         = pfam,
       locus_brief  = locus_brief,
       expression   = expression,
       coexpression = coexpression,
       module_list  = module_list,
       paths        = paths)
}

#' Detect a column by content: the character column best matching `regex`
#'
#' @param dt A data.table.
#' @param regex Regular expression the values should match.
#' @param min_frac Minimum fraction of matching values.
#' @param exclude Columns to skip.
#'
#' @return A column name or `NA_character_`.
#'
#' @noRd
.hbr_detect_acc_col <- function(dt, regex, min_frac = 0.05,
                                exclude = c("model", "Gene_ID")) {
  cols <- setdiff(names(dt), exclude)
  best <- NA_character_; best_f <- 0
  for (cn in cols) {
    v <- dt[[cn]]
    if (!is.character(v)) next
    v <- v[!is.na(v)]
    if (!length(v)) next
    if (length(v) > 20000L) v <- v[seq_len(20000L)]
    f <- mean(grepl(regex, v))
    if (f >= min_frac && f > best_f) { best <- cn; best_f <- f }
  }
  best
}

#' Column after `acc_col` if it is a character column (description)
#'
#' @param dt A data.table.
#' @param acc_col Column name returned by [`.hbr_detect_acc_col()`].
#'
#' @return A column name or `NA_character_`.
#'
#' @noRd
.hbr_next_desc_col <- function(dt, acc_col) {
  if (is.na(acc_col)) return(NA_character_)
  k <- match(acc_col, names(dt)) + 1L
  if (is.na(k) || k > ncol(dt)) return(NA_character_)
  if (is.character(dt[[k]]) && !names(dt)[k] %in% c("Gene_ID", "model"))
    names(dt)[k] else NA_character_
}


# ----------------------------------------------------------------------------
# Functional characterization of superior candidate genes
# ----------------------------------------------------------------------------

#' Functional characterization of superior candidate genes
#'
#' @param candidate_genes Candidate gene table or path.
#' @param functional_annotation_file,goslim_file,interpro_file,pfam_file
#'   Local annotation paths.
#' @param expression_file,coexpression_file Optional XLSX files.
#' @param coexpression_module_list_file,locus_brief_info_file Optional
#'   supplementary files.
#' @param source One of `"msu7"`, `"rapdb"`.
#' @param cache_dir Cache directory for downloads.
#' @param output_dir Output directory.
#' @param output_prefix Output file prefix.
#' @param verbose Print progress.
#' @return Invisibly, a list with `table` and file paths.
#' @examples
#' \donttest{
#' # Downloads the MSU7 reference from the network; run interactively only.
#' if (interactive()) {
#'   cand <- data.frame(
#'     Gene_ID = c("LOC_Os01g0100100", "LOC_Os01g0100200"))
#'   hbr_superior_candidate_genes_functional_characterization(
#'     candidate_genes = cand,
#'     source          = "msu7",
#'     output_dir      = tempdir())
#' }
#' }
#' @noRd
hbr_superior_candidate_genes_functional_characterization <- function(
    candidate_genes,
    functional_annotation_file = NULL,
    goslim_file                = NULL,
    interpro_file              = NULL,
    pfam_file                  = NULL,
    expression_file            = NULL,
    coexpression_file          = NULL,
    coexpression_module_list_file = NULL,
    locus_brief_info_file      = NULL,
    source = c("msu7", "rapdb"),
    cache_dir                  = NULL,
    output_dir                 = "HaploBreedR_results/downstream",
    output_prefix              = "superior_haplotype_publication_table",
    verbose                    = TRUE) {

  source <- match.arg(source)

  read_txt <- function(file, header = TRUE) {
    if (is.null(file)) return(NULL)
    if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
    data.table::fread(file, sep = "\t", header = header,
                      data.table = TRUE, showProgress = FALSE)
  }

  read_candidate_input <- function(x) {
    if (is.data.frame(x)) return(data.table::as.data.table(x))
    if (is.character(x) && length(x) == 1) {
      if (!file.exists(x))
        stop("Candidate genes file not found: ", x, call. = FALSE)
      if (grepl("\\.rds$", x, ignore.case = TRUE))
        return(data.table::as.data.table(readRDS(x)))
      return(data.table::fread(x, data.table = TRUE))
    }
    stop("'candidate_genes' must be a data.frame, path, or RDS.",
         call. = FALSE)
  }

  find_column <- function(dt, candidates) {
    hit <- candidates[candidates %in% names(dt)]
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  safe_paste <- function(x, sep = "; ") {
    x <- as.character(x); x <- trimws(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(NA_character_)
    paste(unique(x), collapse = sep)
  }

  if (verbose) message("[1/5] Loading candidate genes...")
  master <- read_candidate_input(candidate_genes)

  gene_col <- find_column(master,
                           c("Gene_ID", "Gene", "gene_id", "gene",
                             "Locus", "MSU7_Gene_ID", "MSU_Gene_ID",
                             "RAPDB_Gene_ID"))
  if (is.null(gene_col))
    stop("Candidate genes must contain one of: Gene_ID, Gene, gene_id, ",
         "gene, Locus.", call. = FALSE)
  if (gene_col != "Gene_ID")
    data.table::setnames(master, gene_col, "Gene_ID")

  master[, .Gene_ID_Source := as.character(Gene_ID)]
  master[, Gene_ID := hbr_normalize_osa1_gene_id(Gene_ID)]
  master <- master[!is.na(Gene_ID) & Gene_ID != ""]
  master <- unique(master, by = "Gene_ID")

  if (verbose) message("[2/5] Loading OSA1 R7 annotations...")
  annotations <- hbr_read_osa1_r7_annotations(
    functional_annotation_file = functional_annotation_file,
    goslim_file                = goslim_file,
    interpro_file              = interpro_file,
    pfam_file                  = pfam_file,
    expression_file            = expression_file,
    coexpression_file          = coexpression_file,
    source                     = source,
    cache_dir                  = cache_dir,
    verbose                    = verbose)

  if (!is.null(annotations$functional)) {
    fn <- annotations$functional[, .(
      Gene_ID = Gene_ID,
      Functional_Description = annotation)]
    fn <- fn[!duplicated(Gene_ID)]
    master <- merge(master, fn, by = "Gene_ID", all.x = TRUE)
  }

  if (!is.null(annotations$locus_brief)) {
    lb <- annotations$locus_brief
    desc_col <- intersect(c("annotation", "description",
                            "Description", "locus_description"),
                          names(lb))[1]
    if (!is.na(desc_col)) {
      lb <- lb[, .(Gene_ID = Gene_ID,
                   Locus_Description = get(desc_col))]
      lb <- lb[!duplicated(Gene_ID)]
      master <- merge(master, lb, by = "Gene_ID", all.x = TRUE)
    }
  }

  if (!is.null(annotations$goslim) && nrow(annotations$goslim)) {
    go <- annotations$goslim
    go_bp <- go[toupper(GO_Aspect) == "P" | toupper(GO_Aspect) == "BP",
                .(GO_BP = paste(unique(GO_Term), collapse = "; ")),
                by = Gene_ID]
    go_mf <- go[toupper(GO_Aspect) == "F" | toupper(GO_Aspect) == "MF",
                .(GO_MF = paste(unique(GO_Term), collapse = "; ")),
                by = Gene_ID]
    go_cc <- go[toupper(GO_Aspect) == "C" | toupper(GO_Aspect) == "CC",
                .(GO_CC = paste(unique(GO_Term), collapse = "; ")),
                by = Gene_ID]
    go_all <- go[, .(GO_IDs = paste(unique(GO_ID), collapse = "; "),
                     GO_Term_Count = data.table::uniqueN(GO_ID)),
                 by = Gene_ID]
    for (g in list(go_bp, go_mf, go_cc, go_all))
      master <- merge(master, g, by = "Gene_ID", all.x = TRUE)
  }

  if (!is.null(annotations$interpro) && nrow(annotations$interpro)) {
    ip <- annotations$interpro
    acc_col  <- intersect(c("interpro_id", "interpro_accession",
                            "accession"), names(ip))[1]
    desc_col <- intersect(c("description", "name"), names(ip))[1]
    if (is.na(acc_col)) {
      acc_col <- .hbr_detect_acc_col(ip, "^IPR[0-9]{6}$", min_frac = 0.01)
      if (is.na(acc_col))
        acc_col <- .hbr_detect_acc_col(ip, "^[A-Za-z]+[0-9]{4,}$")
      desc_col <- .hbr_next_desc_col(ip, acc_col)
    }
    if (verbose)
      message(sprintf("    InterPro columns used: accession = %s | description = %s",
                      acc_col, desc_col))
    ip_agg <- ip[, .(
      InterPro_Accessions = if (!is.na(acc_col))
        paste(unique(get(acc_col)), collapse = "; ") else NA_character_,
      InterPro_Descriptions = if (!is.na(desc_col))
        paste(unique(get(desc_col)), collapse = "; ") else NA_character_,
      InterPro_Domain_Count = if (!is.na(acc_col))
        data.table::uniqueN(get(acc_col)) else NA_integer_),
      by = Gene_ID]
    master <- merge(master, ip_agg, by = "Gene_ID", all.x = TRUE)
  }

  if (!is.null(annotations$pfam) && nrow(annotations$pfam)) {
    pf <- annotations$pfam
    acc_col  <- intersect(c("pfam_id", "pfam_accession", "accession"),
                          names(pf))[1]
    desc_col <- intersect(c("description", "name"), names(pf))[1]
    if (is.na(acc_col)) {
      acc_col  <- .hbr_detect_acc_col(pf, "^PF[0-9]{5}")
      desc_col <- .hbr_next_desc_col(pf, acc_col)
    }
    if (verbose)
      message(sprintf("    PFAM columns used: accession = %s | description = %s",
                      acc_col, desc_col))
    pf_agg <- pf[, .(
      PFAM_Accessions = if (!is.na(acc_col))
        paste(unique(get(acc_col)), collapse = "; ") else NA_character_,
      PFAM_Descriptions = if (!is.na(desc_col))
        paste(unique(get(desc_col)), collapse = "; ") else NA_character_,
      PFAM_Domain_Count = if (!is.na(acc_col))
        data.table::uniqueN(get(acc_col)) else NA_integer_),
      by = Gene_ID]
    master <- merge(master, pf_agg, by = "Gene_ID", all.x = TRUE)
  }

  if (!is.null(annotations$expression) && nrow(annotations$expression)) {
    ex <- annotations$expression
    sample_cols <- setdiff(names(ex),
                           c("model", "Gene_ID", "transcript_id", "id"))
    sample_cols <- sample_cols[vapply(ex[, ..sample_cols],
                                      is.numeric, logical(1))]
    if (length(sample_cols)) {
      ex_agg <- ex[, {
        vals <- as.matrix(.SD)
        storage.mode(vals) <- "numeric"
        row_mean <- rowMeans(vals, na.rm = TRUE)
        row_max  <- apply(vals, 1, max, na.rm = TRUE)
        breadth  <- rowSums(vals > 0.5, na.rm = TRUE)
        .(Mean_TPM = row_mean, Max_TPM = row_max,
          N_Tissues_Expressed = breadth)
      }, by = Gene_ID, .SDcols = sample_cols]
      ex_agg <- suppressWarnings(ex_agg[, .(
        Mean_TPM = mean(Mean_TPM, na.rm = TRUE),
        Max_TPM = max(Max_TPM, na.rm = TRUE),
        N_Tissues_Expressed = max(N_Tissues_Expressed, na.rm = TRUE)),
        by = Gene_ID])
      ex_agg[, Expression_Breadth :=
               cut(N_Tissues_Expressed,
                   breaks = c(-Inf, 0, 3, 6, Inf),
                   labels = c("None", "Narrow", "Broad", "Ubiquitous"))]
      master <- merge(master, ex_agg, by = "Gene_ID", all.x = TRUE)
    }
  }

  if (!is.null(annotations$coexpression) && nrow(annotations$coexpression)) {
    co <- annotations$coexpression
    mod_col <- intersect(c("module", "Module", "module_id", "Module_ID",
                           "module_name", "Module_Name"),
                         names(co))[1]
    if (!is.na(mod_col)) {
      co_agg <- co[, .(Coexpression_Module =
                         paste(unique(get(mod_col)), collapse = "; ")),
                   by = Gene_ID]
      master <- merge(master, co_agg, by = "Gene_ID", all.x = TRUE)
    }
  }

  if (verbose) message("[3/5] Loading locus brief info...")
  locus <- read_txt(locus_brief_info_file, header = TRUE)
  if (!is.null(locus)) {
    locus_col <- find_column(locus, c("model", "gene_model",
                                       "transcript_id", "Gene_ID",
                                       "gene_id"))
    if (!is.null(locus_col)) {
      if (locus_col != "Gene_ID")
        data.table::setnames(locus, locus_col, "Gene_ID")
      locus[, Gene_ID := hbr_normalize_osa1_gene_id(Gene_ID)]
      locus <- locus[, lapply(.SD, safe_paste), by = Gene_ID]
      master <- merge(master, locus, by = "Gene_ID", all.x = TRUE)
    }
  }

  if (verbose) message("[4/5] Loading coexpression module list...")
  module_list <- read_txt(coexpression_module_list_file, header = TRUE)
  if (!is.null(module_list)) {
    mgc <- find_column(module_list, c("Gene_ID", "gene_id", "gene_model",
                                       "model", "Gene", "gene"))
    if (!is.null(mgc)) {
      if (mgc != "Gene_ID")
        data.table::setnames(module_list, mgc, "Gene_ID")
      module_list[, Gene_ID := hbr_normalize_osa1_gene_id(Gene_ID)]
      mcol <- find_column(module_list,
                           c("module", "Module", "module_id",
                             "Module_ID", "module_name", "Module_Name"))
      if (!is.null(mcol) && mcol != "Coexpression_Module")
        data.table::setnames(module_list, mcol, "Coexpression_Module")
      if ("Coexpression_Module" %in% names(module_list)) {
        module_list <- module_list[
          , .(Coexpression_Module = safe_paste(Coexpression_Module)),
          by = Gene_ID]
        master <- merge(master, module_list, by = "Gene_ID", all.x = TRUE)
      }
    }
  }

  if (verbose) message("[5/5] Preparing publication table...")

  for (col in c("Mean_Expression_TPM", "Maximum_Expression_TPM",
                "Expressed_Sample_Count", "Expression_Breadth",
                "Expression_Category", "Coexpression_Function",
                "Module_Peak_Expression", "Evidence_Score",
                "Evidence_Class"))
    if (!col %in% names(master)) master[, (col) := NA]

  evidence_cols <- c("Functional_Evidence", "GOSlim_Evidence",
                      "PFAM_Evidence", "InterPro_Evidence",
                      "Expression_Evidence", "Coexpression_Evidence")
  for (col in evidence_cols)
    if (!col %in% names(master)) master[, (col) := NA]

  publication_columns <- c(
    "Gene_ID",
    "Chromosome", "Start_Position", "End_Position", "Strand",
    "Functional_Description", "Locus_Description",
    "GO_BP", "GO_MF", "GO_CC",
    "GO_IDs", "GO_Term_Count",
    "InterPro_Domain_Count", "InterPro_Descriptions",
    "InterPro_Accessions",
    "PFAM_Domain_Count", "PFAM_Descriptions", "PFAM_Accessions",
    "Mean_TPM", "Max_TPM",
    "N_Tissues_Expressed", "Expression_Breadth",
    "Coexpression_Module",
    "Superior_Haplotype", "Superior_NT_Seq", "Superior_AA_Seq",
    "Superior_Mean", "Superior_SE", "Cohens_d",
    "P_value", "Permutation_P", "Adjusted_P",
    "Is_Superior", "Improvement_Pct",
    "LD_Region_ID")
  keep_cols <- c(intersect(publication_columns, names(master)),
                 setdiff(names(master), publication_columns))
  master <- master[, ..keep_cols]
  master <- unique(master, by = "Gene_ID")
  if (".Gene_ID_Source" %in% names(master)) {
    master[, Gene_ID := .Gene_ID_Source][, .Gene_ID_Source := NULL]
  }

  master <- .hbr_drop_blank_cols(master, keep_always = c("Gene_ID"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tsv_file <- file.path(output_dir, paste0(output_prefix, ".tsv"))
  csv_file <- file.path(output_dir, paste0(output_prefix, ".csv"))
  rds_file <- file.path(output_dir, paste0(output_prefix, ".rds"))

  data.table::fwrite(master, tsv_file, sep = "\t", na = "")
  data.table::fwrite(master, csv_file, sep = ",",  na = "")
  saveRDS(master, rds_file)

  if (verbose) {
    message("")
    message("==============================================")
    message("FUNCTIONAL ANNOTATION TABLE CREATED")
    message("==============================================")
    message("Genes   : ", nrow(master))
    message("Columns : ", ncol(master))
    message("TSV     : ", tsv_file)
    message("==============================================")
  }

  invisible(list(table = master, tsv = tsv_file,
                 csv = csv_file, rds = rds_file))
}