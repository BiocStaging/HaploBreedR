# HaploBreedR 0.99.1

## Initial submission

First release. The package accompanies the manuscript describing
haplotype-based breeding decisions in rice.

## Superior-haplotype discovery

* `haplobreed_analysis()`: single-trait pipeline from GWAS hits or a
  gene list to a superior haplotype per candidate gene. Builds LD decay
  and LD regions from the genotype, discovers functional variants
  (nonsynonymous CDS and UTR) in every tested gene, tests
  haplotype-trait associations (Student's t-test / one-way ANOVA, optional
  PC correction, Tukey and Duncan post-hoc), and applies effect-size,
  permutation and false-discovery gates. Supports local, MSU7, RAP-DB
  and Ensembl Plants reference sources.

* `hbr_generate_input_qc_report()`: unified input-QC report covering
  sequence-resource quality, genotype-reference compatibility and
  individual-ID overlap, written as a text report and three CSV tables.

* `hbr_check_inputs()`: pre-flight validation of PLINK files, phenotype
  records, GFF3 features and sample-ID overlap before launching the
  full pipeline.

* `hbr_go_enrichment()`: GO over-representation per aspect (BP, MF, CC)
  against the tested-gene universe, with BH-FDR correction.

* `hbr_read_results()`: reads the output of `haplobreed_analysis()`
  (in-memory result, per-trait output directory, or a named list of
  either) into a classed `hbr_run` object. The returned object carries
  the trait direction, effective gate settings, the PCs used, tested
  and superior gene lists, functional SNPs per gene, the phenotype
  table, the publication-style results table, and the pipeline-summary
  row.

* `hbr_regate()`: re-applies the superior-haplotype rule with new
  thresholds. Accepts either the legacy summary-table input (a
  data.frame, returned with the significance and superior flags
  recomputed) or an `hbr_run` object, in which case all four gates -
  including the permutation gate and, when stored, the trait
  threshold - are recomputed and the genes gained or lost versus the
  stored call are reported.

* `hbr_gene_carriers()`: extracts per-individual carrier status of the
  pipeline's superior haplotype at a given gene, used to colour
  prediction-quality panels.

## Genomic prediction and breeding decisions

* `hbr_advanced_decision()`: multi-trait parent ranking and cross design
  from one or more `haplobreed_analysis()` runs. Two-kernel genomic
  BLUP (VanRaden 2008; Kang et al. 2008) with the pipeline's superior
  haplotypes as the feature kernel; honest nested cross-validation with
  in-fold discovery replay (M0 = GBLUP, M2 = features chosen on all
  data, M3 = features re-discovered inside each training fold) and a
  paired comparison of the three (Nadeau and Bengio 2003 corrected
  t-test plus a paired bootstrap over individuals).

* `hbr_fit_models()`: the phenotype-dependent half of the pipeline
  (loading, GRM, haplotype codes, gate replay, cross-validation, final
  models, multi-trait index). Returns an `hgsi_fit` object that can be
  saved once and reused.

* `hbr_cross_design()`: the second half of the pipeline. Consumes an
  `hgsi_fit` object, groups superior haplotypes into haplotype blocks,
  computes cross usefulness (Schnell and Utz 1975; Zhong and Jannink
  2007), pyramiding probabilities and required population sizes
  (Sedcole 1977), and selects a diverse parent set by block coverage
  with a relatedness ceiling.

* `hbr_cross_plan()`: plans any cross between named parents in a
  decision result - not only the top-ranked ones. Reports which
  haplotype blocks are fixed, segregating, split between the parents
  or missing, the exact probability that a doubled-haploid (or RIL)
  progeny carries every segregating target haplotype, the F2
  probabilities, and the number of plants needed to recover at least
  one carrier at a chosen confidence.

* Exported modular building blocks: `hbr_build_grm()`,
  `hbr_reml_two_kernel()`, `hbr_cross_usefulness()`,
  `hbr_pyramid_probability()`, `hbr_population_size()`,
  `hbr_trait_conflicts()`.

## Haplotype blocks

* `hbr_haplotype_blocks()`: groups target loci into haplotype blocks.
  Accepts the original low-level inputs (`units`, `cm_snp`, `chr_snp`)
  and returns integer block IDs, or a pipeline run (`hbr_run`), in
  which case the superior haplotypes are rebuilt, the discovery gate
  replayed, and a classed `hbr_blocks` object returned with the block
  table, per-line carrier matrix and SNP windows. Supports the
  map-based rules (`"recomb"`, `"ld"`, `"either"`), the SNP-level
  methods (`"fourgamete"`, `"gabriel"`) and the two combinations
  (`"recomb+gabriel"`, `"recomb+fourgamete"`). Gabriel is implemented
  with the Gabriel et al. (2002) D' confidence interval and vectorised
  pairwise gamete counting.

* `hbr_haplotype_block_table()`: the block summary table and the
  per-locus assignment in one call, for users who want the block table
  from a low-level locus set without a pipeline run.

## Genetic map

* `hbr_rice_genetic_map()`: the Nagamura et al. (1997) rice RFLP
  linkage-map anchors (n approx. 1,650 markers, total 1,510 cM) as a
  data frame with `chr`, `bp`, `cm`. Interpolated at run time onto the
  physical positions of any PLINK `.bim` with a monotone cubic Hermite
  spline.

## Visualisation

* Publication-quality figures for haplotype structure, gene models and
  cross design: `hbr_plot_haplotype()`, `hbr_plot_gene_structure()`,
  `hbr_plot_haplotype_network()`, `hbr_plot_cross_ranking()` and
  `hbr_plot_cv_forest()`.

* Ten built-in colour palettes (`hbr_palette()`,
  `hbr_palette_names()`) and ten ggplot2 themes (`hbr_theme()`,
  `hbr_theme_names()`); per-gene marker colours derived from the
  IUPAC allele set.

## Manuscript tables

* `hbr_build_manuscript_tables()`: a multi-sheet Excel workbook from a
  pipeline run, with study summary, gene results, donor tables and
  breeding decisions.

## Reference layer

* `hbr_crop_lookup()`: built-in crop-to-species-to-assembly lookup for
  Ensembl Plants; refreshable against the live FTP listing.

* `hbr_annotate_regions_with_genes()` and `hbr_read_gff3_genes()`:
  parse a GFF3 annotation and overlap the resulting gene table with a
  set of genomic regions, labelling each region-gene pair with its
  relationship (contained, overlaps left/right, contains region).

* Local, MSU7, RAP-DB and Ensembl Plants reference sources; chromosome
  names standardised on every slot so they match the genotype
  convention applied by `hbr_normalize_genotype()`.

## S3 methods

* `print()`, `summary()` and `plot()` methods for the classed results
  of the pipeline: `hbr_run` (single-trait discovery), `hgsi_fit`
  (fitted genomic prediction models), `hgsi_advanced` (cross-design
  result), `hbr_regate` (re-gated call) and `hbr_blocks` (haplotype
  block assignment).

* `print()` methods for the reference and QC objects:
  `hbr_reference`, `hbr_sequence_qc`, `hbr_compat_report`,
  `hbr_id_report` and `hbr_genotype`.

## Reproducibility and performance

* Every run writes `pipeline_summary.csv`, `phenotype_data_used.csv`,
  `functional_snps_all_tested_genes.csv`, `gene_attrition_log.csv`
  and `pipeline_result_for_breeding.rds`; the RDS carries the gate
  parameters and the exact phenotype table used, so downstream runs
  can be replayed without any external state.

* Genotype reads are streamed through `BEDMatrix`, so peak memory is
  one column block rather than the full genotype matrix; the RDS
  genotype packing is streamed from the same handle.

* `hbr_advanced_decision()` records `provenance`, `warnings` and the
  evaluation label per trait (`honest_discovery_in_fold`,
  `conditional_on_gwas_candidates`, `conditional_on_discovered_set`,
  `approximate_replay`).

* `verbose` in every function selects a single-line-per-phase console
  summary; the detailed narration is available by passing
  `verbose = TRUE` where noted.