#' Calculate GC content of genomic regions
#'
#' @param region A data frame containing the columns: "Chr", "Start", "End",
#'   and "Strand", such as a data frame contained in a sublist named
#'   *intergenic_region* from the output of the [get_feature_saf()] function,
#'   or a object of [GenomicRanges::GRangesList-class] or
#'   [GenomicRanges::GRanges-class].
#' @param BSgenome An object of [BSgenome::BSgenome-class]. Make sure the
#'   chromosome names (aka seqnames) in the BSgenome object match those in the
#'   *region_saf* data frame and the *region_gr* object.
#' @param batch_size An integer(1) vector, specifying how many regions are
#'   processed each batch.
#' @return A data.frame contains two columns: gc_content and width.
#'   \describe{
#'   \item{gc_content}{GC contents (proportion) of genomic regions}
#'   \item{width}{widths of genomic regions}
#'   }
#' @importFrom Biostrings letterFrequency
#' @importFrom BSgenome getSeq
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics width width<- strand<-
#' @export
#'
#' @examples
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' ucsc_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#' ucsc_seqnames <-
#'     seqnames(ucsc_BSgenome)[!grepl(
#'         "_alt|fix|hap\\d+",
#'         seqnames(ucsc_BSgenome)
#'     )]
#'
#' ## Wierd: BSgenome.Hsapiens.UCSC.hg19 has both chrM and chrMT for
#' ## mitochondrial genome. renomve chrMT.
#'
#' if (all(c("chrM", "chrMT") %in% ucsc_seqnames)) {
#'     ucsc_seqnames <- ucsc_seqnames[!ucsc_seqnames %in% "chrMT"]
#' }
#'
#' ensembl_seqnames <- gsub(
#'     "^chr", "",
#'     gsub(
#'         "chrM$", "MT",
#'         gsub(
#'             "v", ".",
#'             gsub("_random|chr[^_]+_", "", ucsc_seqnames)
#'         )
#'     )
#' )
#' ## for BSgenome.Hsapiens.UCSC.hg19, scaffold seqnames start with lower case,
#' ## should be changed to upper case. For example, change "gl000231" to
#' ## "GL000231". Add a suffix ".1" to these scaffold seqnames.
#'
#' seqname_alias <- data.frame(ucsc = ucsc_seqnames,
#'                             ensembl = ensembl_seqnames)
#' bsgenome <- UCSC2Ensembl(
#'     UCSC_BSgenome = BSgenome.Hsapiens.UCSC.hg38,
#'     genome_version = "GRCh38",
#'     seqname_alias = seqname_alias
#' )
#' region <- data.frame(
#'     GeneID = as.character(seq_len(10)),
#'     Chr = rep("1", 10),
#'     Start = 100000 * (seq_len(10)),
#'     End = 100000 * (seq_len(10)) + 1000,
#'     Strand = rep("+", 10)
#' )
#'
#' gc_contents <- calculate_region_gc(
#'     region = region,
#'     BSgenome = bsgenome,
#'     batch_size = 2000
#' )
#'
calculate_region_gc <- function(region = NULL,
                                BSgenome = NULL,
                                batch_size = 2000) {
    if (!is(BSgenome, "BSgenome")) {
        stop("Please provide a valid BSgenome object.")
    }
    if (is.null(region)) {
        stop("region must be specified.")
    }

    if (batch_size <= 0 || batch_size != as.integer(batch_size)) {
        stop("batch_size must be a positive integer.")
    }

    if (!is.null(region)) {
        if (is.data.frame(region) &&
            any(c("GeneID", "Chr", "Start", "End", "Strand") !=
                colnames(region))) {
            stop(
                "region is a data farme, but not in a valid format. ",
                'It must contain the following columns: "GeneID"',
                '"Chr", "Start", "End", and "Strand".'
            )
        } else if (!is.data.frame(region) && !is(region, "GRanges") &&
            !is(region, "GRangesList")) {
            stop(
                "region is not a SAF, thus it must be providedas a object of ",
                "GRanges or GRangesList."
            )
        }

        if (is(region, "GRanges") && is.null(names(region))) {
            stop("region is a GRanges, but it doesn't have names.")
        }
    }

    if (is.data.frame(region)) {
        saf_seqnames <- unique(region$Chr)
        if (any(!saf_seqnames %in% seqnames(BSgenome))) {
            stop(
                "Some chromosome names in the intergenic SAF data frame ",
                "is not in the BSgenome object."
            )
        }
        grouping <- rep(seq_len(ceiling(nrow(region) / batch_size)),
            each = batch_size
        )[seq_len(nrow(region))]
        region_list <- split(region, f = grouping)

        gc_contents <- do.call(rbind, lapply(region_list, function(.x) {
            getSeq(
                BSgenome,
                names = .x$Chr,
                start = .x$Start,
                end = .x$End,
                strand = "+",
                as.character = FALSE
            ) %>%
                letterFrequency("GC", as.prob = TRUE) %>%
                as.data.frame()
        }))
        rownames(gc_contents) <- region$GeneID
        colnames(gc_contents) <- "gc_content"
        gc_contents$width <- region$End - region$Start + 1
    } else {
        grouping_f <- rep(seq_len(ceiling(length(region) / batch_size)),
            each = batch_size
        )[seq_len(length(region))]
        grouping <- split(seq_len(length(region)), grouping_f)
        if (is(region, "GRangesList")) {
            gc_contents <- do.call(rbind, lapply(grouping, function(.x) {
                sequences <- getSeq(BSgenome,
                    region[.x],
                    as.character = FALSE
                )
                gc_contents <-
                    do.call(rbind, lapply(sequences, function(.y) {
                        gene_size <- sum(width(.y))
                        gene_gc <-
                            sum(letterFrequency(.y, "GC",
                                                as.prob = FALSE)[, "G|C"]) /
                            gene_size
                        data.frame(gc_content = gene_gc, width = gene_size)
                    }))
                gc_contents
            }))
        } else if (is(region, "GRanges")) {
            strand(region) <- "+"
            gc_contents <- do.call(rbind, lapply(grouping, function(.g) {
                gc_contents <-
                    getSeq(BSgenome, region[.g], as.character = FALSE) %>%
                    letterFrequency("GC", as.prob = TRUE) %>%
                    as.data.frame()
                colnames(gc_contents) <- "gc_content"
                gc_contents$width <- width(region[.g])
                gc_contents
            }))
        }
        rownames(gc_contents) <- names(region)
    }
    gc_contents
}


#' Calculate GC content of genes
#'
#' Calculate GC content of all genes based on sequences of collapsed exons of
#' each gene.
#'
#' @param ensdb_sqlite A character(1) specifying a path to an SQLite file to
#'   store the an object of the [ensembldb::EnsDb-class] or an object of the
#'   [ensembldb::EnsDb-class].
#' @param BSgenome An object of [BSgenome::BSgenome-class]. Make sure the
#'   chromosome names (aka seqnames) in the BSgenome object match those in the
#'   [ensembldb::EnsDb-class] specified by *ensdb_sqlite*.
#' @param batch_size An integer(1) vector, specifying how many regions are
#'   processed each batch.
#' @return A data frame contains two columns: gc_content and width.
#'   \describe{
#'   \item{gc_content}{GC contents (proportion) of genes}
#'   \item{width}{widths of genes}
#'   }
#' @importFrom ensembldb exonsBy seqlevels
#' @importFrom GenomicRanges reduce
#' @importFrom AnnotationFilter SeqNameFilter
#' @export
#'
#' @examples
#' if (interactive()) {
#'     library("BSgenome.Hsapiens.UCSC.hg38")
#'     library("EnsDb.Hsapiens.v86")
#'     edb <- EnsDb.Hsapiens.v86
#'
#'     ucsc_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#'     ucsc_seqnames <-
#'         seqnames(ucsc_BSgenome)[!grepl("_alt|fix|hap\\d+",
#'         seqnames(ucsc_BSgenome))]
#'
#'     ## Wierd: BSgenome.Hsapiens.UCSC.hg19 has both chrM and chrMT for
#'     ## mitochondrial
#'     ## genome. renomve chrMT.
#'
#'     if (all(c("chrM", "chrMT") %in% ucsc_seqnames)) {
#'         ucsc_seqnames <- ucsc_seqnames[!ucsc_seqnames %in% "chrMT"]
#'     }
#'
#'     ensembl_seqnames <- gsub(
#'         "^chr", "",
#'         gsub(
#'             "chrM$", "MT",
#'             gsub(
#'                 "v", ".",
#'                 gsub("_random|chr[^_]+_", "", ucsc_seqnames)
#'             )
#'         )
#'     )
#'     ## for BSgenome.Hsapiens.UCSC.hg19, scaffold seqnames start with
#'     ## lower case, should be changed to upper case. For example,
#'     ## change "gl000231" to "GL000231". Add a suffix ".1" to these
#'     ## scaffold seqnames.
#'
#'     seqname_alias <- data.frame(ucsc = ucsc_seqnames,
#'                                 ensembl = ensembl_seqnames)
#'
#'     bsgenome <- UCSC2Ensembl(
#'         UCSC_BSgenome = ucsc_BSgenome,
#'         genome_version = "GRCh38",
#'         seqname_alias = seqname_alias
#'     )
#'     gene_contents <- calculate_gene_gc(
#'         ensdb_sqlite = edb,
#'         BSgenome = bsgenome,
#'         batch_size = 2000
#'     )
#' }
calculate_gene_gc <- function(ensdb_sqlite = NULL,
                              BSgenome = NULL,
                              batch_size = 2000) {
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

    if (!is(BSgenome, "BSgenome")) {
        stop("Please provide a valid BSgenome object.")
    }
    if (batch_size <= 0 || batch_size != as.integer(batch_size)) {
        stop("batch_size must be a positive integer.")
    }
    common_seqnames <- intersect(seqnames(BSgenome), seqlevels(ensdb))
    exons_grlist <- exonsBy(
        ensdb,
        by = "gene",
        filter = SeqNameFilter(value = common_seqnames),
        use.names = FALSE
    )
    exons_grlist <- reduce(exons_grlist)
    gene_gc <- calculate_region_gc(
        region = exons_grlist,
        BSgenome = BSgenome,
        batch_size = batch_size
    )
    gene_gc
}


#' Correct for DNA contamination in consideratio of GC-bias
#'
#'
#' @param intergenic_counts A data frame or matrix containing counts assigned to
#'   each intergenic region of each sample, such as the sublist, *counts*  from
#'   sublist named *intergenic_region* of an output from the [summarize_reads()]
#'   function.
#' @param intergenic_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and lengths of each intergenic region, such as an output from
#'   the [calculate_region_gc()] function.
#' @param plot A logical(1) vector, specifying whether to output a panel of
#'   scatter plot showing a bi-variate distribution of GC content and count per
#'   base of intergenic regions in each sample. Default is TRUE. A loess
#'   regression line is displayed in each subplot.
#'
#' @return A matrix containing estimated DNA contamination in the count per base
#'   of intergenic region binned by GC content, with a bin size of 5% of each
#'   sample, 20 bins. Column names are sample names, rownames are GC content
#'   bins in the form of '(0,0.05]'.
#'
#' @importFrom stats quantile loess predict
#' @import ggplot2
#'
gc_bin_contamination <- function(intergenic_counts,
                                 intergenic_gc,
                                 plot = TRUE) {
    estimate_contamination <- function(.count, .sample_name, .gc) {
        ## remove zero count intergenic entries
        .rm <- .count == 0
        .count <- .count[!.rm]
        .gc <- .gc[!.rm, ]

        count_per_base <- .count / .gc$width

        # remove first and last quantile
        keep <-
            # count_per_base >= quantile(count_per_base, probs = 0.01) &
            count_per_base <= quantile(count_per_base, probs = 0.99)

        .gc <- .gc[keep, ]
        count_per_base <- count_per_base[keep]

        ## loess fitting
        l <- loess(count_per_base ~ .gc$gc_content,
            family = "symmetric",
            span = 0.25
        )

        data <- cbind(
            .gc, count_per_base,
            predict(l, newdata = .gc$gc_content)
        )
        colnames(data) <- c("gc_content", "width", "CPB", "predict")
        data$sample_name <- .sample_name

        temp_data <- data[order(data$gc_content), ]
        temp_data <- split(temp_data, f = cut(temp_data$gc_content,
            breaks = seq(0, 1, 0.05)
        ))
        perc5_bin_cpb <- vapply(
            temp_data,
            function(.x) {
                quantile(.x$CPB, prob = 0.5)
            },
            numeric(1)
        )
        perc5_bin_cpb[10:20] <- c(perc5_bin_cpb[9:1], 0, 0)
        list(contamination = perc5_bin_cpb, plot_data = data)
    }

    contamination <- mapply(
        function(.intergenic_count,
                 .sample_name,
                 .intergenic_gc) {
            contamination_level <- estimate_contamination(.intergenic_count,
                .sample_name,
                .gc = .intergenic_gc
            )
        },
        intergenic_counts,
        colnames(intergenic_counts),
        MoreArgs = list(.intergenic_gc = intergenic_gc),
        SIMPLIFY = FALSE
    )
    plot_data <- do.call("rbind", lapply(contamination, "[[", 2))
    plot_data$sample_name <- factor(plot_data$sample_name,
        levels = unique(plot_data$sample_name)
    )
    gc_bin_contamination <-
        do.call(cbind, lapply(contamination, "[[", 1))
    rownames(gc_bin_contamination) <- gsub(
        "\\.\\d+%", "",
        rownames(gc_bin_contamination)
    )

    p <- ggplot(plot_data, aes(x = gc_content, y = CPB)) +
        geom_point(
            color = "gray",
            alpha = 0.3,
            size = 0.3
        ) +
        geom_line(aes(x = gc_content, y = predict),
            color = "blue"
        ) +
        xlab("GC%") +
        ylab("Count per base") +
        facet_wrap(. ~ sample_name, ncol = 3)
    if (plot && interactive()) {
        p
    }
    gc_bin_contamination
}


#' Correct DNA contamination considering GC-bias effect
#'
#' Correct DNA contamination considering GC-bias effect on fragment
#' amplification. Intergenic regions are binned based on their GC content
#' ranging from 0 to 100%, with a bin size of 5%. Per gene DNA contamination
#' is estimated as the product of count per base in a GC content matching bin
#' of intergenic regions and the total collapsed exons of a gene and is
#' subtracted away from the gene count matrix.
#'
#'
#' @param salmon_res A list of matrices containing gene-level abundances,
#'   counts, lengths, such as the output of the [salmon_res()] function.
#'   For more details, See [tximport::tximport()].
#' @param gene_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and total exon lengths of each gene. An output of the
#'   [calculate_gene_gc()] function.
#' @param intergenic_counts A data frame or matrix containing counts assigned to
#'   each intergenic region of each sample, such as the sublist, *counts*  from
#'   sublist named *intergenic_region* of an output from the [summarize_reads()]
#'   function.
#' @param intergenic_gc A data frame or matrix with two columns: gc_content in
#'   proportion between 0 and 1, and width in basepair, containing
#'   GC-content and lengths of each intergenic region, such as an output from
#'   the [calculate_region_gc()] function.
#' @param plot A logical(1) vector, specifying whether to output a panel of
#'   scatter plot showing a bi-variate distribution of GC content and count per
#'   base of intergenic regions in each sample. Default is TRUE.
#' @return A data frame containing corrected count for each gene (row) of each
#'   sample (column).
#' @export
#'
#' @examples
#' if (interactive()) {
#'    options(timeout = max(3000, getOption("timeout")))
#'    tmp_dir <- tempdir()
#'    ## download feaureCounts results
#'    count_url <- paste0(
#'        "https://www.dropbox.com/scl/fi/lyvh6bsljnqxtq85nnugq/",
#'        "read_count_summary.RData?rlkey=e0tmpehpxtnr1fdx4fz0h8sa0&dl=1"
#'    )
#'    download.file(
#'        url = count_url,
#'        destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'        mode = "wb"
#'    )
#'    load(file.path(tmp_dir, "read_count_summary.RData"))
#'
#'    # download the Salmon quantification results
#'    salmon_url <- paste0(
#'        "https://drive.google.com/",
#'        "uc?export=download&id=13vsobXnENoiYOBo-Vf_e_00xbR8ekc4Z"
#'    )
#'    salmon_destfile <- file.path(tmp_dir, "salmon_quant_summary.RData")
#'    download.file(url = salmon_url,
#'                  destfile = salmon_destfile,
#'                  mode = "wb")
#'
#'    ## load the salmon_quant object
#'    load(salmon_destfile)
#'
#'    gene_gc_content_url <-
#'        paste0(
#'            "https://www.dropbox.com/scl/fi/taoj18w5t4up22d04xihk/GRCh38.gene.",
#'            "exon.collapsed.GC.content.RDS?rlkey=zj31ul2slckrn5zm8vq6tgq31&dl=1"
#'        )
#'    gene_gc_destfile <- file.path(tmp_dir,
#'                                  "GRCh38.gene.exon.collapsed.GC.content.RDS")
#'    download.file(url = gene_gc_content_url,
#'                  destfile = gene_gc_destfile,
#'                  mode = "wb")
#'
#'    ## load the gene_ object
#'    gene_gc <- readRDS(gene_gc_destfile)
#'
#'    intergenic_gc_content_url <-
#'        paste0(
#'            "https://www.dropbox.com/scl/fi/ik5z5etbr63nnv43jna8c/GRCh38.",
#'            "intergenic.GC.content.RDS?rlkey=kxrn6ixdnvhpba9h7w2r265xp&dl=1"
#'        )
#'    intergenic_gc_destfile <- file.path(
#'        tmp_dir,
#'        "GRCh38.intergenic.GC.content.RDS"
#'    )
#'    download.file(
#'        url = intergenic_gc_content_url,
#'        destfile = intergenic_gc_destfile,
#'        mode = "wb"
#'    )
#'
#'    ## load the salmon_quant object
#'    intergenic_gc <- readRDS(intergenic_gc_destfile)
#'    gc_bias_corrected_count <-
#'        gc_bias_correction(
#'            salmon_res = salmon_quant,
#'            gene_gc = gene_gc,
#'            intergenic_counts =
#'                counts_summary$intergenic_region$counts,
#'            intergenic_gc = intergenic_gc,
#'            plot = FALSE
#'        )
#' }


gc_bias_correction <- function(salmon_res = NULL,
                               gene_gc = NULL,
                               intergenic_counts = NULL,
                               intergenic_gc = NULL,
                               plot = FALSE) {
    if (!is.list(salmon_res) ||
        any(!c("abundance", "counts", "length") %in%
            names(salmon_res))) {
        stop(
            deparse(substitute(salmon_res)), " is not a valid ",
            "tximport output!"
        )
    }

    if (!is.data.frame(intergenic_counts) &&
        !is.matrix(intergenic_counts)) {
        stop("intergenic_counts must be a data.frame or matrix.")
    }

    if (!is.data.frame(gene_gc) && !is.matrix(gene_gc)) {
        stop("intergenic_counts must be a data.frame or matrix.")
    }

    if (!is.data.frame(intergenic_gc) && !is.matrix(intergenic_gc)) {
        stop("intergenic_counts must be a data.frame or matrix.")
    }

    if (!all(colnames(salmon_res$counts) == colnames(intergenic_counts))) {
        stop(
            "colnames of the Salmon count table is not exactly the same ",
            "as those of the intergenic_counts table."
        )
    }

    if (!all(rownames(salmon_res$counts) %in% rownames(gene_gc)) ||
        !all(rownames(gene_gc) %in% rownames(salmon_res$counts))) {
        warning(
            "Gene names in the Salmon count table is not exactly the same as\n",
            "those of the gene_gc data frame. The intersection of gene ",
            "names from both\n",
            "sources will be used."
        )

        common_genes <-
            intersect(rownames(salmon_res$counts), rownames(gene_gc))
        if (length(common_genes) < 1000) {
            stop(
                "There are only",
                length(common_genes),
                "genes shared between the ",
                "Salmon count table and the gene_gc data frame. Please double ",
                "check the rownames."
            )
        }
        salmon_res <- lapply(salmon_res[seq_len(3)], function(.x) {
            .x <- .x[common_genes, ]
        })

        gene_gc <- gene_gc[common_genes, ]
    }

    ## intergenic gc and counts
    if (any(!rownames(intergenic_counts) %in% rownames(intergenic_gc)) ||
        any(!rownames(intergenic_gc) %in% rownames(intergenic_counts))) {
        stop("rownames of intergenic_gc DO NOT match those of",
             " intergenic_counts!")
    }

    intergenic_counts <-
        as.data.frame(intergenic_counts[rownames(intergenic_gc), ])

    gc_bin_cpb <- gc_bin_contamination(intergenic_counts,
        intergenic_gc,
        plot = plot
    )

    gc_bin_cpb <- as.data.frame(gc_bin_cpb)
    gc_bin_cpb$bin <- rownames(gc_bin_cpb)
    rownames(gc_bin_cpb) <- NULL

    gene_gc <- gene_gc[order(gene_gc$gc_content), ]
    gene_gc$bin <- cut(gene_gc$gc_content, breaks = seq(0, 1, 0.05))
    gene_gc$width <- NULL
    gene_gc$GeneID <- rownames(gene_gc)

    gene_gc_sample <- merge(gene_gc, gc_bin_cpb,
        by = "bin",
        all.x = TRUE,
        sort = FALSE
    )
    rownames(gene_gc_sample) <- gene_gc_sample$GeneID
    gene_gc_sample$GeneID <- NULL
    gene_gc_sample <- gene_gc_sample[rownames(salmon_res$length), ]
    gene_contamination <-
        gene_gc_sample[, -c(1, 2)] * salmon_res$length
    corrected_counts <- salmon_res$counts - gene_contamination
    corrected_counts[corrected_counts < 0] <- 0
    corrected_counts <- round(corrected_counts)
}
