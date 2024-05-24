#' Correction for gDNA contamination in RNA-seq data
#'
#' Correction for gDNA contamination in stranded RNA-seq data based on the
#' strandedness information or in the unstranded RNA-seq data using the "GC%",
#' "IR%" or "Global" method.
#'
#' @param is_stranded_library A logical(1), specifying whether the RNA-seq
#'   library is stranded or not.
#' @param salmon_summary A list of matrices containing gene-level abundances,
#'   counts, lengths, such as the output of the [salmon_res()] function. For
#'   more details, See [tximport::tximport()]. It can be an element named
#'   *salmon_summary* of the list returned by the [create_diagnostic_plot()]
#'   function.
#' @param correction_method A character(1), specifying a method for correcting
#'   for gDNA contamination in unstranded RNA-seq data. The options are "GC%",
#'   "Global", "IR%". This parameter is ignored for stranded RNA-seq data.
#' @param featurecounts_summary It can be an element named
#'   *featurecounts_summary* of the list returned by the
#'   [create_diagnostic_plot()] function. It is only
#'   needed for correcting for gDNA contamination in unstranded RNA-seq data
#'   using the "Global" or the "GC%" method.
#' @param saf_list A list of data frames containing annotation in the SAF
#'   format. It can be an element named saf_list* of the list returned by the
#'   [create_diagnostic_plot()] function. It is only needed for correcting for
#'   gDNA contamination in unstranded RNA-seq data using the "GC%" method.
#' @param unstranded_metadata  A data frame, with column names: sample_name,
#'   group, IR_rate (the percentages of reads mapping to the intergenic region
#'   for all samples), batch if any, and other possible covariates for linear
#'   models used to correct for gDNA contamination. The *sample_name* and
#'   *group* columns contain the unique sample labels and the experimental
#'   conditions for each sample, respectively. The order of the columns
#'   doesn't matter. It can be an element named *metadata* of the list returned
#'   by the [create_diagnostic_plot()] function. It is only needed for
#'   unstranded RNA-seq data.
#' @param covariates A character(n), specifying the names of covariates for
#'   linear models used to correcting for gDNA contamination and other
#'   nuisance paramenters. The names of covariates must be contained in the
#'   column names of the `metadata` data frame.
#' @param covariate_types A character(n), specifying the data types of
#'   covariates. "numeric" for a quantitative covariate; "categorical" for a
#'   non-quantitative covariate, such as batch information. The length of
#'   `covariate_types` must be the same as that of `covariates`.
#' @param Ensembl_BSgenome An object of the [BSgenome::BSgenome-class] in the
#'   UCSCstyle, which is not avaible from Bioconductor. You can build such a
#'   package using a multi-fasta file for the reference genome of interest
#'   using the [make_BSgenome()] function.
#' @param UCSC_BSgenome An object of the [BSgenome::BSgenome-class] in the UCSC
#'   style, such as BSgenome.Hsapiens.UCSC.hg38.It is only needed for correcting
#'   for gDNA contamination in unstranded RNA-seq data using the "GC%" method
#'   and a UCSC-style BSgenome is provided as an alternative to building a
#'   BSgenome package from scratch.
#' @param genome_version A character(1), specifying the genome version, such as
#'   "GRCh38". Caution: make sure the genome_version here is the same as that
#'   for the [make_ensdb()]. It is only needed for correcting for
#'   gDNA contamination in unstranded RNA-seq data using the "GC%" method.
#' @param seqname_alias A data frame or a tab-delimited file to a data frame,
#'   with two columns: ucsc and ensembl for UCSC-style seqnames
#'   (chromosome/scaffold names) and Ensembl-style seqnames, respectively. It is
#'   only needed for correcting for gDNA contamination in unstranded RNA-seq
#'   data using the "GC%" method and a UCSC-style BSgenome is provided. All the
#'   chromosomes/scaffolds in the primary reference genome assembly should be
#'   provided in the data frame.
#' @param stranded_metadata A data frame or a path to a tab-delimited
#'   text file with column names: `salmon_quant_file_strand`,
#'   `salmon_quant_file_opposite_strand`, `sample_name`, and `group`. The
#'   `salmon_quant_file_strand` and `salmon_quant_file_opposite_strand` column
#'   contains the full paths to `quant.sf` files for quantitation of gene
#'   expression by Salmon pseudo-alignment with the library strandedness type
#'   set to the true and opposite orientations. See Salmon's description of
#'   library type (<https://salmon.readthedocs.io/en/latest/library_type.html>)
#'   for details. The `sample_name` columns contain the unique sample labels
#'   for each sample. The order of the three columns doesn't matter. It is only
#'   needed for stranded RNA-seq data.
#' @param ensdb_sqlite An EnsDb object or a character(1) vector, specifying a
#'   path to an SQLite database for an object of the [ensembldb::EnsDb-class].
#'   It can be an element named `ensdb_sqlite_file` of the list returned by the
#'   [create_diagnostic_plot()] function. It is needed for the "GC%" method-
#'   based correction for unstranded RNA-seq data and correction for stranded
#'   RNA-seq data.
#'
#' @return For unstranded RNA-seq datat, it returns a gDNA contamination
#'   corrected count matrix if the "GC%" or "Global" correction method is
#'   chosen, or a log2CPM matrix if the "IR%" correction method is used.
#'   For stranded RNA-seq data, it returns a gDNA contamination corrected
#'   count matrix.
#'
#' @importFrom limma voom lmFit
#' @importFrom stats as.formula model.matrix
#' @export
#'
#' @examples
#' \dontrun{
#' metadata <- read.delim(system.file("extdata",
#'     "CD1A.RNAseq.metadata.txt",
#'     package = "CleanUpRNAseq"
#' ))
#' options(timeout = max(3000, getOption("timeout")))
#' gtf_url <- paste0(
#'     "https://ftp.ensembl.org/pub/release-110/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
#' )
#' gtf <- basename(gtf_url)
#' tmp_dir <- tempdir()
#' download.file(
#'     url = gtf_url,
#'     destfile = file.path(tmp_dir, gtf),
#'     mode = "wb"
#' )
#' gtf <- file.path(tmp_dir, gtf)
#'
#' diagnosis_res <-
#'     create_diagnostic_plot(
#'         gtf = gtf,
#'         metadata = metadata,
#'         normalization = "DESeq2",
#'         ensdb_sqlite =
#'             file.path(tmp_dir, "GRCh38.GTF.EnsDb.sqlite"),
#'         out_dir = tmp_dir,
#'         threads = 1
#'     )
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' UCSC_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#' sequences_to_keep <- seqnames(UCSC_BSgenome)[!grepl(
#'     "_(alt|fix|hap\\d+)$",
#'     seqnames(UCSC_BSgenome)
#' )]
#' ensembl_chr <- gsub(
#'     "v", ".",
#'     gsub(
#'         "_random$", "",
#'         gsub("^chr.+?_", "", sequences_to_keep)
#'     )
#' )
#' ensembl_chr <- gsub("^chr", "", gsub("chrM", "chrMT", ensembl_chr))
#'
#'
#' corrected_expression <-
#'     correct_for_contamination(
#'         is_stranded_library = FALSE,
#'         salmon_summary = diagnosis_res$salmon_summary,
#'         correction_method = "GC%",
#'         featurecounts_summary =
#'             diagnosis_res$featurecounts_summary,
#'         saf_list = diagnosis_res$saf_list,
#'         unstranded_metadata = diagnosis_res$metadata,
#'         covariates = c("batch", "IR_rate"),
#'         covariate_types = c("categorical", "numeric"),
#'         UCSC_BSgenome = UCSC_BSgenome,
#'         seqname_alias =
#'             data.frame(
#'                 ucsc = sequences_to_keep,
#'                 ensembl = ensembl_chr
#'             ),
#'         ensdb_sqlite = diagnosis_res$ensdb_sqlite_file
#'     )
#' }
#'
correct_for_contamination <-
    function(is_stranded_library = FALSE,
             salmon_summary = NULL,
             correction_method = c("GC%", "Global", "IR%"),
             featurecounts_summary = NULL,
             saf_list = NULL,
             unstranded_metadata =
                 data.frame(
                     sample_name = vector(mode = "character"),
                     group = vector(mode = "character"),
                     IR_rate = vector(mode = "numeric"),
                     batch = vector(mode = "character")
                 ),
             covariates = c("batch", "IR_rate"),
             covariate_types = c("categorical", "numeric"),
             Ensembl_BSgenome = NULL,
             UCSC_BSgenome = NULL,
             genome_version = "GRCh38",
             seqname_alias =
                 data.frame(
                     ucsc = paste0("chr", c(seq_len(22), "X", "Y", "M")),
                     ensembl = c(seq_len(22), "X", "Y", "MT")
                 ),
             stranded_metadata =
                 data.frame(
                     salmon_quant_file_strand = vector(mode = "character"),
                     salmon_quant_file_opposite_strand =
                         vector(mode = "character"),
                     sample_name = vector(mode = "character")
                 ),
             ensdb_sqlite = NULL) {
        if (!is_stranded_library) {
            correction_method <- match.arg(correction_method)
            if (correction_method == "Global") {
                corrected_counts <-
                    global_correction(
                        intergenic_featureCounts_res =
                            featurecounts_summary$intergenic_region,
                        salmon_res = salmon_summary,
                        lambda = 1
                    )
                corrected_counts
            } else if (correction_method == "GC%") {
                if (!is(ensdb_sqlite, "EnsDb") && !is.character(ensdb_sqlite)) {
                    stop("Please provide a single EnsDb object or a EnsDb ",
                         "SQLite database file!")
                } else if (is.character(ensdb_sqlite)) {
                    if (!file.exists(ensdb_sqlite)) {
                        stop("ensdb_sqlite is a non-existing SQLite ",
                             "database file")
                    }
                }

                if (is.null(Ensembl_BSgenome) && is.null(UCSC_BSgenome)) {
                    stop("You must provide a Ensembl-style BSgenome, ",
                         " or a UCSC-style BSgenome!")
                } else if (!is.null(UCSC_BSgenome) &&
                    is(UCSC_BSgenome, "BSgenome")) {
                    bsgenome <- UCSC2Ensembl(
                        UCSC_BSgenome = UCSC_BSgenome,
                        genome_version = genome_version,
                        seqname_alias = seqname_alias
                    )
                } else if (!is.null(Ensembl_BSgenome) &&
                    is(Ensembl_BSgenome, "BSgenome")) {
                    bsgenome <- Ensembl_BSgenome
                }

                intergenic_gc <-
                    calculate_region_gc(
                        region = saf_list$intergenic_region,
                        BSgenome = bsgenome,
                        batch_size = 2000
                    )

                gene_gc <- calculate_gene_gc(
                    ensdb_sqlite = ensdb_sqlite,
                    BSgenome = bsgenome,
                    batch_size = 2000
                )

                corrected_counts <-
                    gc_bias_correction(
                        salmon_res = salmon_summary,
                        gene_gc = gene_gc,
                        intergenic_counts =
                            featurecounts_summary$intergenic_region$counts,
                        intergenic_gc = intergenic_gc,
                        plot = FALSE
                    )

                corrected_counts
            } else {
                if (!is.data.frame(unstranded_metadata) ||
                    any(
                        !c(
                            "sample_name", "group",
                            "IR_rate", "batch"
                        ) %in% colnames(unstranded_metadata)
                    )) {
                    stop(
                        "'metadata' must be a data frame containing ",
                        "all columns: ",
                        "sample_name, group, IR_rate, and batch."
                    )
                }
                metadata <- unstranded_metadata
                if (any(!covariates %in% colnames(metadata))) {
                    stop("Some covariates are not provided in ",
                         "the unstranded_metadata!")
                }

                metadata$group <- factor(metadata$group)
                for (i in seq_along(covariates))
                {
                    if (covariate_types[i] == "categorical") {
                        metadata[, covariates[i]] <-
                          factor(metadata[, covariates[i]])
                    }
                }
                model <- model.matrix(
                    as.formula(paste(
                        "~ 0 + group",
                        paste(covariates, collapse = " + "),
                        sep = " + "
                    )),
                    data = metadata
                )

                if (!all(sort(as.character(metadata$sample_name)) ==
                        sort(colnames(salmon_summary$counts)))) {
                    stop(
                        "Sample names in metadata are not the same ",
                        "as those in the Salmon count table!"
                    )
                }
                counts <- salmon_summary$counts[, metadata$sample_name]

                adj.exp <- IR_percent_correction(
                    expr.data = counts,
                    design = model,
                    meta.data = metadata
                )
                adj.exp
            }
        } else {
            if (!is(ensdb_sqlite, "EnsDb") && !is.character(ensdb_sqlite)) {
                stop("Please provide a single EnsDb object or a EnsDb ",
                     "SQLite database file!")
            } else if (is.character(ensdb_sqlite)) {
                if (!file.exists(ensdb_sqlite)) {
                    stop("ensdb_sqlite is a non-existing SQLite ",
                         "database file")
                }
            }
            corrected_counts <-
                correct_stranded_lib(
                    metadata = stranded_metadata,
                    ensdb_sqlite = ensdb_sqlite
                )
            corrected_counts
        }
    }

## helper function-using voom and limma linear model to adjust gene expression
IR_percent_correction <- function(expr.data = NULL,
                                  design = NULL,
                                  meta.data = NULL) {
    ## adjusted expression using voom()
    v <- voom(expr.data, design = design)
    fit <- lmFit(v, design = design)
    group_num <- nlevels(as.factor(meta.data$group))
    if (ncol(design) > group_num) {
        col.covariate <- (group_num + 1):ncol(design)
        adj.exp <- v$E - fit$coefficients[, col.covariate] %*%
            t(design[, col.covariate])
    } else {
        adj.exp <- v$E
    }
    adj.exp
}
