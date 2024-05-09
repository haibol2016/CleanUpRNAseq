#' Import Salmon quantification output into R
#'
#' Import Salmon quantification output into R using tximport
#'
#' @param metadata A data frame with column names: salmon_quant_file,
#'   sample_name, and group. The *salmon_quant_file* column contains the full
#'   paths to *quant.sf* files for quantitation of gene expression by Salmon
#'   pseudo-alignment. See <https://salmon.readthedocs.io/en/latest/>. The
#'   *sample_name* and *group* columns contain the unique sample labels, and
#'   the experimental conditions for each sample. The order of the three
#'   columns doesn't matter.
#' @param ensdb_sqlite An object of the [ensembldb::EnsDb-class] or a
#'   character(1) vector, specifying a path to an SQLite database for an
#'   object of the [ensembldb::EnsDb-class].
#' @param filtered_gene_biotypes A character(n) vector,specifying the biotypes
#'   of genes which will not be considered for downstream gene expression
#'   analysis. By default, genes of the following biotypes are excluded:
#'   "artifact","TEC","miRNA","tRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "rRNA",
#'   "rRNA_pseudogene", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA",
#'   "vault_RNA", "TR_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene",
#'   "IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
#'   because their transcripts are too short (< 200 nt) or not likely expressed
#'   at all.
#'
#' @return A list of of matrices of gene-level abundance, counts, and length.
#'   See [tximport::tximport()].
#'   \describe{
#'   \item{abundance}{A numeric matrix containing abundance (TPM) for each gene
#'                    of each sample}
#'   \item{counts}{A numeric matrix containing read count (fraction) for each
#'                 gene of each sample}
#'   \item{length}{A numeric matrix containing length (bp) for each gene of
#'                 each sample}
#' }
#' @export
#' @importFrom tximport tximport
#' @importFrom ensembldb transcripts
#' @importFrom AnnotationFilter AnnotationFilter
#' @importFrom utils write.table
#'
#' @examples
#' options(timeout = max(3000, getOption("timeout")))
#' ## using trucated quant.sf files for demonstration purpose
#' quant_files <- system.file("extdata",
#'     "K084CD7PCD1N.quant.sf",
#'     package = "CleanUpRNAseq"
#' )
#'
#' metadata <- data.frame(
#'     sample_name = "CD1AN_m2_1",
#'     salmon_quant_file = quant_files,
#'     group = "CD1AN"
#' )
#'
#' # download the EnsDb SQLite database: GRCh38.V110.ensdb.sqlite.zip
#' tmp_dir <- tempdir()
#' ensdb_url <- paste0(
#'     "https://www.dropbox.com/scl/fi/5i8mwuxe7mcbxz4qvf43g/",
#'     "GRCh38.V110.ensdb.sqlite.zip?rlkey=udcpv8yafbdujvom628u6aucb&dl=1"
#' )
#' download.file(
#'     url = ensdb_url,
#'     destfile = file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite.zip"),
#'     mode = "wb"
#' )
#' unzip(file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite.zip"),
#'     exdir = tmp_dir
#' )
#' hs_ensdb_sqlite <- file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite")
#'
#' salmon_counts <- salmon_res(
#'     metadata = metadata,
#'     ensdb_sqlite = hs_ensdb_sqlite
#' )
#'
salmon_res <- function(metadata =
                           data.frame(
                               salmon_quant_file = vector(mode = "character"),
                               sample_name = vector(mode = "character")
                           ),
                       ensdb_sqlite = NULL,
                       filtered_gene_biotypes = c(
                           "artifact",
                           "TEC",
                           "miRNA",
                           "tRNA",
                           "misc_RNA",
                           "Mt_rRNA",
                           "Mt_tRNA",
                           "rRNA",
                           "rRNA_pseudogene",
                           "scaRNA",
                           "scRNA",
                           "snoRNA",
                           "snRNA",
                           "sRNA",
                           "vault_RNA",
                           "TR_V_pseudogene",
                           "IG_C_pseudogene",
                           "IG_J_pseudogene",
                           "IG_V_pseudogene",
                           "IG_pseudogene",
                           "TR_J_pseudogene",
                           "TR_V_pseudogene"
                       )) {
    if (!is(ensdb_sqlite, "EnsDb") && !is.character(ensdb_sqlite)) {
        stop("Please provide a single EnsDb object or a EnsDb SQLite",
             " database file!")
    } else if (is.character(ensdb_sqlite)) {
        if (file.exists(ensdb_sqlite)) {
            ensdb <- EnsDb(ensdb_sqlite)
        } else {
            stop("ensdb_sqlite is a non-existing SQLite database file")
        }
    } else {
        ensdb <- ensdb_sqlite
    }

    if (!is.data.frame(metadata)) {
        stop("metadata should be a data frame.")
    } else {
        if (any(is.na(metadata)) || any(metadata == "") ||
            any(metadata == " ")) {
            stop("Missing value(s) found in the metadata.")
        }
        quant_files <- metadata$salmon_quant_file
        sample_name <- metadata$sample_name
        if (any(duplicated(quant_files)) ||
            any(duplicated(sample_name))) {
            stop("Some Salmon quant.sf files or sample names are not unique.")
        }
    }

    if (any(!(
        c("sample_name", "salmon_quant_file") %in%
            colnames(metadata)
    ))) {
        stop(
            "sample_name or salmon_quant_file is not included in the metadata ",
            "columns."
        )
    }

    if (length(quant_files) < 1 || any(!file.exists(quant_files))) {
        stop("At least one existing Salmon quant.sf file should be provided!")
    }


    transcripts <- transcripts(
        ensdb,
        columns = c("gene_id", "tx_id"),
        AnnotationFilter(~ gene_biotype !=
            filtered_gene_biotypes)
    )
    ## gene-level count
    tx2gene <- as.data.frame(transcripts)[, c("tx_id", "gene_id")]

    names(quant_files) <- sample_name

    null <- lapply(quant_files, function(.x) {
        quant_sf <- read.delim(.x, header = TRUE, as.is = TRUE)
        if (any(grepl("[._].+$", quant_sf$Name))) {
            quant_sf$Name <- gsub("[._].+$", "", quant_sf$Name)
        }
        write.table(
            quant_sf,
            file = .x,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    })

    txi <- tximport(
        quant_files,
        type = "salmon",
        txOut = FALSE,
        tx2gene = tx2gene
    )
    txi
}

#' Global correction for DNA contamination
#'
#' Correct for DNA contamination in RNA-seq data using the 2.5 times of median
#' counts per base of intergenic regions with at least one count.
#'
#' @param intergenic_featureCounts_res A sublist of the function
#'   summarize_reads, containing read summary for intergenic regions.
#' @param salmon_res A list of matrices containing gene-level abundances,
#'   counts, lengths. See [tximport::tximport()].
#' @param lambda A positive number specifying how many times of the median read
#'   coverage of non-zero count intergenic regions used as an estimate of DNA
#'   contamination. The default of *lambda* is 1,but it could be adjusted based
#'   on the gene-level count distributions of the resulting corrected count
#'   table output by the [check_read_distribution()] function. A value between
#'   1 and 3 can be tried. Ideally the distributions of all samples from a
#'   given condition should be very similar.
#'
#' @return A matrix containing an RNA-seq count table corrected for DNA
#'   contamination, with rows for genes and columns for samples.
#' @export
#' @examples
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://www.dropbox.com/scl/fi/lyvh6bsljnqxtq85nnugq/",
#'     "read_count_summary.RData?rlkey=e0tmpehpxtnr1fdx4fz0h8sa0&dl=1"
#' )
#' download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#'
#' # download the Salmon quantification results
#' salmon_url <- paste0(
#'     "https://drive.google.com/",
#'     "uc?export=download&id=13vsobXnENoiYOBo-Vf_e_00xbR8ekc4Z"
#' )
#' destfile <- file.path(tmp_dir, "salmon_quant_summary.RData")
#' download.file(
#'     url = salmon_url,
#'     destfile = destfile,
#'     mode = "wb"
#' )
#'
#' ## load the salmon_quant object
#' load(destfile)
#'
#' corrected_counts <- global_correction(
#'     intergenic_featureCounts_res =
#'         counts_summary$intergenic_region,
#'     salmon_res = salmon_quant
#' )

global_correction <- function(intergenic_featureCounts_res = NULL,
                              salmon_res = NULL,
                              lambda = 1) {
    if (!is.list(salmon_res) ||
        any(!c("abundance", "counts", "length") %in%
            names(salmon_res))) {
        stop(
            deparse(substitute(salmon_res)), " is not a valid ",
            "tximport output!"
        )
    }

    if (!is.list(intergenic_featureCounts_res) ||
        any(!c("counts", "annotation") %in%
            names(intergenic_featureCounts_res))) {
        stop(
            deparse(substitute(intergenic_featureCounts_res)),
            " is not a valid ",
            "featureCounts output!"
        )
    }

    if (any(
        !colnames(intergenic_featureCounts_res$counts) %in%
            colnames(salmon_res$counts)
    ) ||
        any(
            !colnames(salmon_res$counts) %in%
                colnames(intergenic_featureCounts_res$counts)
        )) {
        stop(
            "Column names of the 2 count tables generated by featureCounts",
            " and Salmon are different!"
        )
    }

    if (lambda < 1) {
        stop("lambda should not be less than 1.")
    }

    ## intergenic read count per base
    intergenic_counts <-
        as.data.frame(intergenic_featureCounts_res$counts)
    ## make sure colnames match bwteen two count tables
    intergenic_counts <-
        intergenic_counts[, colnames(salmon_res$counts)]

    intergenic_length <-
        intergenic_featureCounts_res$annotation[, 6, drop = FALSE]
    intergenic_fpb <- mapply(function(.x, .y) {
        .x / .y
    }, intergenic_counts, intergenic_length, SIMPLIFY = FALSE)

    intergenic_fpb <- as.data.frame(do.call("cbind", intergenic_fpb))
    colnames(intergenic_fpb) <- colnames(intergenic_counts)

    ## remove rows from intergenic_fpb, where all samples have 0
    intergenic_fpb <- intergenic_fpb[rowSums(intergenic_fpb) != 0, ]
    median_fpb <- vapply(intergenic_fpb, function(x) {
        10^(median(log10(x + 1))) - 1
    }, numeric(1))

    gene_len <- as.data.frame(salmon_res$length)
    dna_contamination_count <-
        do.call(
            cbind,
            mapply(function(len, cov) {
                len * cov
            }, gene_len, lambda * median_fpb, SIMPLIFY = FALSE)
        )

    colnames(dna_contamination_count) <- colnames(intergenic_counts)
    dna_contamination_count <-
        dna_contamination_count[, colnames(salmon_res$counts)]
    counts <- as.matrix(salmon_res$counts - dna_contamination_count)

    # set negative values to 0
    counts <- ifelse(counts < 0, 0, counts)
    counts <- round(counts)
    counts <- counts[rowSums(counts) != 0, ]
    counts
}

#' Correct for gDNA contamination in stranded libraries
#'
#' Correct for gDNA contamination in stranded libraries based on Salmon
#' quantitation using the real and opposite library strandedness information
#'
#' @param metadata A data frame with column names: BAM_file, sample_name, and
#'   group. The *salmon_quant_file_strand* and
#'   *salmon_quant_file_reverse_strand* column contains the full paths to
#'   *quant.sf* files for quantitation of gene expression by Salmon
#'   pseudo-alignment with the library strandedness type set to the true and
#'   reverse orientations. See Salmon's description of library type
#'   (<https://salmon.readthedocs.io/en/latest/library_type.html>)
#'   for details. The *sample_name* columns contain the unique sample labels
#'   for each sample. The order of the three columns doesn't matter.
#' @param ensdb_sqlite A character(1), a path to an SQLite database for
#'   an object of the [ensembldb::EnsDb-class].
#'
#' @return A list of of matrices of gene-level abundance, counts, and length.
#'   See [tximport::tximport()]. The count matrix is corrected for DNA
#'   contamination and rounded into integers.
#'   \describe{
#'   \item{abundance}{A numeric matrix containing corrected abundance (TPM) for
#'                    each gene of each sample}
#'   \item{counts}{An *integer* matrix containing *corrected* read count for
#'                 each gene of each sample}
#'   \item{length}{A numeric matrix containing length (bp) for each gene of
#'                 each sample}
#' }
#' @export
#' @importFrom tximport tximport
#' @importFrom ensembldb transcripts
#'
#' @examples
#' ## using made-up salmon quant.sf files, NOT RUN. For real situation,
#' ## replace the quant.sf files and sample names with real ones.
#' quant_files <-
#'     system.file("extdata",
#'                 c("K084CD7PCD1N.quant.sf", "K084CD7PCD1P.quant.sf"),
#'                 package = "CleanUpRNAseq")
#' metadata <- data.frame(
#'     salmon_quant_file_strand = quant_files,
#'     salmon_quant_file_reverse_strand = quant_files,
#'     sample_name = c("CD1AN_m2_1", "CD1AP_m2_1")
#' )
#'
#' # download the EnsDb SQLite database: GRCh38.V110.ensdb.sqlite.zip
#' tmp_dir <- tempdir()
#' ensdb_url <- paste0(
#'     "https://www.dropbox.com/scl/fi/5i8mwuxe7mcbxz4qvf43g/",
#'     "GRCh38.V110.ensdb.sqlite.zip?rlkey=udcpv8yafbdujvom628u6aucb&dl=1"
#' )
#' download.file(
#'     url = ensdb_url,
#'     destfile = file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite.zip"),
#'     mode = "wb"
#' )
#' unzip(file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite.zip"),
#'     exdir = tmp_dir
#' )
#' hs_ensdb_sqlite <- file.path(tmp_dir, "GRCh38.V110.ensdb.sqlite")
#'
#' corrected_counts <- correct_stranded_lib(
#'     metadata = metadata,
#'     ensdb_sqlite = hs_ensdb_sqlite
#' )
#'
correct_stranded_lib <-
    function(metadata =
                 data.frame(
                     salmon_quant_file_strand = vector(mode = "character"),
                     salmon_quant_file_reverse_strand =
                         vector(mode = "character"),
                     sample_name = vector(mode = "character")
                 ),
             ensdb_sqlite = NULL) {
        if (!is(ensdb_sqlite, "EnsDb") && !is.character(ensdb_sqlite)) {
            stop("Please provide a single EnsDb object or a EnsDb ",
                 "SQLite database file!")
        } else if (is.character(ensdb_sqlite)) {
            if (!file.exists(ensdb_sqlite)) {
                stop("ensdb_sqlite is a non-existing SQLite ",
                     "database file")
            }
        }

        if (!is.data.frame(metadata)) {
            stop("metadata must be a data frame.")
        } else if (any(is.na(metadata)) || any(metadata == "") ||
                   any(metadata == " ")) {
            stop("Missing value(s) found in the metadata.")
        }

        if (any(
            !c(
                "salmon_quant_file_strand",
                "salmon_quant_file_reverse_strand",
                "sample_name"
            ) %in% colnames(metadata)
        )) {
            stop(
                "metadata must contain the three columns: ",
                "salmon_quant_file_strand",
                "salmon_quant_file_reverse_strand",
                " and sample_name."
            )
        }
        quant_files_strand <- metadata$salmon_quant_file_strand
        quant_files_reverse_strand <-
            metadata$salmon_quant_file_reverse_strand
        sample_name <- metadata$sample_name
        if (any(duplicated(quant_files_strand)) ||
            any(duplicated(quant_files_reverse_strand)) ||
            any(duplicated(sample_name)))
        {
            stop("Some Salmon quant.sf files or sample names are ",
                 "not unique.")
        }

        if (length(quant_files_strand) < 1 ||
            length(quant_files_reverse_strand) < 1 ||
            any(!file.exists(quant_files_reverse_strand)) ||
            any(!file.exists(quant_files_strand))) {
            stop("At least one existing Salmon quant.sf file should be",
                 " provided!")
        }

        ## real strandedness
        metadata_strand <-
            metadata[, c("salmon_quant_file_strand", "sample_name")]
        colnames(metadata_strand) <-
            c("salmon_quant_file", "sample_name")
        salmon_strand <- salmon_res(
            metadata = metadata_strand,
            ensdb_sqlite = ensdb_sqlite
        )

        ## opposite strandedness
        metadata_reverse_strand <-
            metadata[, c("salmon_quant_file_reverse_strand", "sample_name")]
        colnames(metadata_reverse_strand) <-
            c("salmon_quant_file", "sample_name")
        salmon_reverse_strand <- salmon_res(
            metadata = metadata_reverse_strand ,
            ensdb_sqlite = ensdb_sqlite
        )

        txi <- mapply(
            function(.strand, .reverse_strand) {
                val <- .strand - .reverse_strand
                val[val < 0] <- 0
                val
            },
            salmon_strand[c("abundance", "counts")],
            salmon_reverse_strand[c("abundance", "counts")],
            SIMPLIFY = FALSE
        )

        txi$length <- salmon_strand$length
        txi$counts <- round(txi$counts)
        txi
    }
