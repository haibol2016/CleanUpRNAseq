#' A package for detecting and correcting for DNA contamination in RNA-seq data
#'
#' The CleanUpRNAseq package provides two categories of functions for detecting
#' DNA contamination and correcting for contamination from RNA-seq count data
#'
#' @section functions for detecting DNA contamination: check_expr_distribution,
#' ccheck_expressed_gene_percentage, check_read_distribution,
#' check_sample_correlation, exploratory_analysis, summarize_reads,
#' get_feature_saf, granges_to_saf, make_ensdb
#' @section functions for correcting for DNA contamination: global_correction,
#' filter_unexpressed_genes, get_unexpressed_spliced_genes, salmon_res

#'
#' @docType package
#' @name CleanUpRNAseq
globalVariables(c("sample_name", "region_type",
                  "percent",
                  "Sample", "Statistics", "Status",
                  "log2CPM", "log2Count", "CPB",
                  "percent_gene_tpm_gt1", "gc_content",
                  "PC1", "PC2", "assigned_percent",
                  "percent_gene_cpm_gt1", "group"))
