#' Generate SAF files for genomic features
#'
#' Generate SAF (simplified annotation format) files for genomic features: genic
#' regions,intergenic regions, exonic regions, intronic regions, rRNA genes,
#' mitochondrial genome, chloroplast genome (only for plants).
#'
#' @param ensdb_sqlite A character(1) specifying a path to an SQLite file to
#'   store the an object of the [ensembldb::EnsDb-class] or an object of the
#'   [ensembldb::EnsDb-class].
#' @param bamfile A character(1), a path to a BAM file for the experiment of
#'   interest. The BAM file is used to extract genome information:
#'   chromosome/scaffold names and lengths.
#' @param mitochondrial_genome A character(1), mitochondrial genome name in the
#'   EnsDb database (ie, the mitochondrial name in column 1 of the GTF file
#'   used to generate the EnsDb SQLite database).
#' @param chloroplast_genome A character(1), chloroplast genome name in the
#'   EnsDb database (ie, the chloroplast name in column 1 of the GTF file
#'   used to generate the EnsDb SQLite database). This is only relevant
#'   for plants.
#'
#' @importFrom ensembldb genes exons EnsDb
#' @importFrom AnnotationFilter GeneBiotypeFilter SeqNameFilter
#' @importFrom plyranges reduce_ranges setdiff_ranges as_granges
#' @importFrom Rsamtools seqinfo BamFile
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqnames
#'
#' @return A list of data frames containing SAF for genomic features:
#'   genic regions, intergenic regions, exonic regions, intronic regions,
#'   rRNA genes, mitochondrial genome, chloroplast genome (only for plants),
#'   respectively.
#'   \describe{
#'   \item{gene}{a data frame containing a SAF for genes}
#'   \item{exon}{a data frame containing  a SAF for exons}
#'   \item{intergenic_region}{a data frame containing  a SAF for intergenic
#'                            regions}
#'   \item{intronic_region}{a data frame containing a SAF for intronic region}
#'   \item{rRNA}{a data frame containing a SAF for rRNA exons}
#'   \item{mitochondrion}{a data frame containing a SAF for the mitochodrion}
#'   \item{chloroplast (optional)}{a data frame containing a SAF for
#'                                 chloroplast, plnat only}
#'  }
#'
#' @export
#' @examples
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "CleanUpRNAseq")
#' tmp_dir <- tempdir()
#' hs_ensdb_sqlite <- make_ensdb(
#'     gtf = gtf,
#'     outfile = file.path(tmp_dir, "EnsDb.hs.v110.sqlite"),
#'     organism_latin_name = "Homo_Sapiens",
#'     genome_version = "GRCh38",
#'     Ensembl_release_version = 110
#' )
#' bam_file <- system.file("extdata", "K084CD7PCD1N.srt.bam",
#'     package = "CleanUpRNAseq"
#' )
#' saf_list <- get_feature_saf(
#'     ensdb_sqlite = hs_ensdb_sqlite,
#'     bamfile = bam_file,
#'     mitochondrial_genome = "MT"
#' )

get_feature_saf <- function(ensdb_sqlite = NULL,
                            bamfile = NULL,
                            mitochondrial_genome = c("MT", "chrM"),
                            chloroplast_genome = c("chrPltd", "Pltd")) {
    if (!file.exists(bamfile)) {
        stop("A single BAM file should be specified via bamfile.")
    }

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

    genome_info <- as.data.frame(seqinfo(BamFile(bamfile)))

    gene_GR <-
        genes(ensdb, filter = SeqNameFilter(rownames(genome_info))) %>%
        reduce_ranges()
    gene_saf <- granges_to_saf(gene_GR)

    exon_GR <-
        exons(ensdb, filter = SeqNameFilter(rownames(genome_info))) %>%
        reduce_ranges()
    exon_saf <- granges_to_saf(exon_GR)

    ## chromosome/scaffolds GRanges

    genome_GR <- data.frame(
        seqnames = rownames(genome_info),
        start = 1L,
        end = genome_info$seqlengths,
        strand = "*"
    ) %>%
        as_granges()

    # intergenic regions
    intergenic_GR <- setdiff_ranges(genome_GR, gene_GR)
    intergenic_saf <- granges_to_saf(intergenic_GR)

    ## intronic regions
    intronic_GR <- setdiff_ranges(gene_GR, exon_GR)
    intronic_saf <- granges_to_saf(intronic_GR)

    ## rRNA transcripts
    rrna_GR <- exons(ensdb, filter = GeneBiotypeFilter("rRNA"))
    rrna_saf <- granges_to_saf(rrna_GR)

    saf_list <- list(
        gene = gene_saf,
        exon = exon_saf,
        intergenic_region = intergenic_saf,
        intronic_region = intronic_saf,
        rRNA = rrna_saf
    )
    if (any(mitochondrial_genome %in% rownames(genome_info))) {
        mitochondrion <- genome_GR[as.character(seqnames(genome_GR)) %in%
            mitochondrial_genome]
        saf_list$mitochondrion <- granges_to_saf(mitochondrion)
    }

    if (any(chloroplast_genome %in% rownames(genome_info))) {
        chloroplast <- genome_GR[as.character(seqnames(genome_GR)) %in%
            chloroplast_genome]
        saf_list$chloroplast <- granges_to_saf(chloroplast)
    }
    saf_list
}


#' Summarize reads for different genomic features
#'
#' Summarize reads in alignment files, SAM or BAM, to different genomic regions,
#' such as genic regions, intergenic regions, exonic regions, intronic regions,
#' rRNA genes, mitochrondrial genome, chloroplast genome (only for plants), and
#' gene-level exonic regions using [Rsubread::featureCounts()].
#'
#' @param metadata A data frame with column names: BAM_file, sample_name, and
#'   group. The *BAM_file* column contains the full paths to BAM files;
#'   The *sample_name* and *group* columns contain the unique sample labels,
#'   and the experimental conditions for each sample. The order of the three
#'   columns doesn't matter.
#' @param isPairedEnd A logical(1), indicating whether the sequencing data is
#'   paired-end or not.
#' @param strandSpecific An integer(1), specifying the strandedness of the
#'   RNA-seq libraries. The possible options are 0 (unstranded), 1 (stranded,
#'   read 1 or single-end read comes from the forward/sense strand) and 2
#'   (reversely stranded, read 1 or single-end read comes from the
#'   reverse/antisense strand).
#'   See <https://salmon.readthedocs.io/en/latest/library_type.html> and
#'   <https://chipster.csc.fi/manual/library-type-summary.html> for more details
#'   about types of stranded libraries.
#' @param allowMultiOverlap A logical(1), indicating if a read is allowed to be
#'   assigned to more than one feature (or meta-feature) if it is found to
#'   overlap with more than one feature (or meta-feature). FALSE by default. A
#'   read (or read pair) will be assigned to the feature (or meta-feature) that
#'   has the largest number of overlapping bases, if the read (or read pair)
#'   overlaps with multiple features (or meta-features).
#' @param countMultiMappingReads A logical(1), indicating if multi-mapping
#'   reads/fragments should be counted, TRUE by default. ‘NH’ tag is used to
#'   located multi-mapping reads in the input BAM/SAM files.
#' @param fraction A logical(1) indicating if fractional counts are produced for
#'   multi-mapping reads and/or multi-overlapping reads. FALSE by default.
#' @param minMQS An integer(1), giving the minimum mapping quality score a read
#'   must satisfy in order to be counted. For paired-end reads, at least one end
#'   should satisfy this criteria. 0 by default.
#' @param saf_list A list of data frames containing annotation in the SAF
#'   format, such as the output of the [get_feature_saf()] function.
#' @param gtf A character(1), specifying a path to a GTF file.
#' @param threads An integer(1), number of threads for
#'   [Rsubread::featureCounts()] calling.
#' @param verbose A logical(1) vector, indicating if verbose information for
#' debugging will be generated. This may include information such as unmatched
#' chromosomes/contigs between reads and annotation.
#'
#' @return A list of lists, each of the sublist contains an output of the
#'   [Rsubread::featureCounts()] function for each type of genomic features.
#'   For more details of each sublist, see the *value* section of the
#'   documentation for the [Rsubread::featureCounts()] function.
#'   \describe{
#'   \item{gene}{featureCounts output with a SAF for genes as annotation}
#'   \item{exon}{featureCounts output with a SAF for exons as annotation}
#'   \item{intergenic_region}{featureCounts output with a SAF for intergenic
#'                            regions as annotation}
#'   \item{intronic_region}{featureCounts output with a SAF for introns
#'                          as annotation}
#'   \item{rRNA}{featureCounts output with a SAF for rRNA exons as annotation}
#'   \item{mitochondrion}{featureCounts output with a SAF for the mitochondrion
#'                        as annotation}
#'   \item{chloroplast (optional)}{featureCounts output with a SAF for the
#'                                 chloroplast as annotation, plant only}
#'   \item{gtf}{featureCounts output with a GTF as annotation}
#' }
#' @importFrom Rsubread featureCounts
#' @export
#' @examples
#' if (interactive()){
#' library(R.utils)
#' options(timeout = max(3000, getOption("timeout")))
#' bams <- system.file("extdata",
#'     "K084CD7PCD1N.srt.bam",
#'     package = "CleanUpRNAseq"
#' )
#'
#' metadata <- data.frame(
#'     sample_name = "CD1AN_m2_1",
#'     BAM_file = bams,
#'     group = "CD1AN"
#' )
#'
#' # download the GTF file
#' gtf_url <- paste0(
#'     "https://ftp.ensembl.org/pub/release-110/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
#' )
#' gtf <- basename(gtf_url)
#' tmp_dir <- tempdir()
#' retry_download({download.file(
#'     url = gtf_url,
#'     destfile = file.path(tmp_dir, gtf),
#'     mode = "wb"
#' )})
#' gunzip(file.path(tmp_dir, gtf), remove = TRUE)
#'
#' # download the SAF file
#' saf_url <- paste0("https://zenodo.org/records/11458839/files/",
#'                   "hs_saf_list.RData?download=1")
#' retry_download({download.file(
#'     url = saf_url,
#'     destfile = file.path(tmp_dir, "hs_saf_list.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "hs_saf_list.RData"))
#'
#' counts_list <- summarize_reads(
#'     metadata = metadata,
#'     isPaired = TRUE,
#'     strandSpecific = 0,
#'     saf_list = saf_list,
#'     gtf = gsub(".gz", "", file.path(tmp_dir, gtf)),
#'     threads = 1
#' )
#' }
#'
summarize_reads <- function(metadata =
                                data.frame(
                                    sample_name = vector(mode = "character"),
                                    BAM_file = vector(mode = "character"),
                                    group = vector(mode = "character")
                                ),
                            isPairedEnd = TRUE,
                            strandSpecific = 0,
                            allowMultiOverlap = FALSE,
                            countMultiMappingReads = TRUE,
                            fraction = TRUE,
                            minMQS = 0,
                            saf_list = NULL,
                            gtf = NULL,
                            threads = 1,
                            verbose = FALSE) {
    if (!is.data.frame(metadata)) {
        stop("metadata must be a data frame.")
    }
    if (any(is.na(metadata)) ||
        any(metadata == "") || any(metadata == " ")) {
        stop("Missing value(s) found in the metadata.")
    }
    if (any(!c("sample_name", "BAM_file", "group") %in% colnames(metadata))) {
        stop(
            "Not all column names: sample_name, BAM_file, group are included ",
            "in the metadata."
        )
    }
    bamfiles <- metadata$BAM_file
    sample_name <- metadata$sample_name
    if (any(duplicated(bamfiles)) || any(duplicated(sample_name))) {
        stop("Some BAM files or sample names are not unique.")
    }

    if (length(bamfiles) < 1 || any(!file.exists(bamfiles))) {
        stop("At least one existing BAM file should be provided!")
    }

    if (length(saf_list) < 5 ||
        any(
            !c(
                "gene",
                "exon",
                "intergenic_region",
                "intronic_region",
                "rRNA"
            ) %in%
                names(saf_list)
        )) {
        stop(
            "A valid SAF list for gene, exon, intergenic region,",
            "intronic region, rRNA genes is needed!"
        )
    }

    if (!file.exists(gtf)) {
        stop("A single GTF file is needed!")
    }
    if (!strandSpecific %in% c(0:2)) {
        stop("strandSpecific must be 0, 1, or 2.")
    }

    saf_list$gtf <- gtf
    count_res <- mapply(
        function(annotation, isGTF, junction) {
            res <- featureCounts(
                files = bamfiles,
                # annotation
                annot.ext = annotation,
                isGTFAnnotationFile = isGTF,
                GTF.featureType = "exon",
                GTF.attrType = "gene_id",

                # overlap between reads and features
                allowMultiOverlap = allowMultiOverlap,
                minOverlap = 1,
                fracOverlap = 0,
                fracOverlapFeature = 0,
                largestOverlap = TRUE,
                nonOverlap = NULL,
                nonOverlapFeature = NULL,
                # multi-mapping reads
                countMultiMappingReads = countMultiMappingReads,
                # fractional counting
                fraction = fraction,
                # read filtering
                minMQS = minMQS,
                splitOnly = FALSE,
                nonSplitOnly = FALSE,
                primaryOnly = FALSE,
                ignoreDup = FALSE,
                # strandness
                strandSpecific = strandSpecific,
                # exon-exon junctions
                juncCounts = junction,
                genome = NULL,
                # parameters specific to paired end reads
                isPairedEnd = isPairedEnd,
                countReadPairs = TRUE,
                requireBothEndsMapped = FALSE,
                checkFragLength = FALSE,
                countChimericFragments = TRUE,
                autosort = TRUE,
                # number of CPU threads
                nthreads = threads,
                # read group
                byReadGroup = FALSE,
                # report assignment result for each read
                reportReads = NULL,
                reportReadsPath = NULL,
                # miscellaneous
                tmpDir = ".",
                verbose = verbose
            )
            colnames(res$counts) <- sample_name
            colnames(res$stat) <- c("Status", sample_name)
            if ("counts_junction" %in% names(res)) {
                colnames(res$counts_junction)[9:ncol(res$counts_junction)] <-
                    sample_name
            }

            res
        },
        saf_list,
        c(rep(FALSE, length(saf_list) - 1), TRUE),
        c(rep(FALSE, length(saf_list) - 1), TRUE),
        SIMPLIFY = FALSE
    )
    count_res
}

#' Visualize assignment statistics of reads/fragments by featureCounts
#'
#' @param assignment_stat A data frame containing the assignment statistics with
#'   GTF as annotation by the [Rsubread::featureCounts()] function, such as the
#'   `stat` element in the sublist `gtf` in the output of the
#'   [summarize_reads()] function. You can access it via
#'   `featureCounts()$gtf$stat`.
#'
#' @return A ggplot object, showing percentages and number of fragments in each
#'   assignment category as determined by [Rsubread::featureCounts()].
#'
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2
#' @examples
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#'
#' p <- check_read_assignment_stat(
#'     assignment_stat =
#'         counts_summary$gtf$stat
#' )
#' p
#'
check_read_assignment_stat <- function(assignment_stat = NULL) {
    if (is.null(assignment_stat)) {
        stop("Read assignment statistics is no provided.")
    }
    if (!is.data.frame(assignment_stat) ||
        (is.data.frame(assignment_stat) &&
            colnames(assignment_stat)[1] != "Status")) {
        stop(
            "assignment_stat must be a data frame with column",
            " 1 named as 'Status'"
        )
    }

    assignment_stat <-
        assignment_stat[rowSums(assignment_stat[, -1]) != 0, ]
    assignment_stat_pct <- assignment_stat
    assignment_stat_pct[, 2:ncol(assignment_stat_pct)] <-
        mapply(
            function(.x, .y) {
                .x / .y * 100
            },
            assignment_stat[, -1],
            colSums(assignment_stat[, -1]),
            SIMPLIFY = FALSE
        )

    assignment_stat[, 2:length(assignment_stat)] <-
        assignment_stat[, 2:length(assignment_stat)] / 10^6

    assignment_count_long <- melt(
        assignment_stat,
        id.vars = "Status",
        variable.name = "Sample",
        value.name = "Statistics"
    )
    assignment_count_long$Stat <- "Count"
    assignment_pct_long <- melt(
        assignment_stat_pct,
        id.vars = "Status",
        variable.name = "Sample",
        value.name = "Statistics"
    )
    assignment_pct_long$Stat <- "Percentage"
    assignment_long <-
        rbind(assignment_count_long, assignment_pct_long)

    assignment_long$Sample <- factor(assignment_long$Sample,
        levels = unique(assignment_long$Sample)
    )

    p <- ggplot(assignment_long, aes(
        x = Sample, y = Statistics,
        fill = Status
    )) +
        geom_bar(stat = "identity", color = "white") +
        ylab("Reads") +
        facet_wrap(~Stat, scales = "free_y") +
        theme(
            text = element_text(size = 8),
            axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 0,
                size = 8
            ),
            legend.key.size = unit(0.5, "cm"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6)
        )
    p
}

#' Check reads distribution over different types of genomic features
#'
#' Generate plots showing percentages of fragments assigned to different
#' genomic features.
#'
#' @param featurecounts_list A list of lists, containing output of the
#'   [summarize_reads()] function.
#' @param metadata  A data frame with column names: sample_name and group for
#'   sample names and experimental groups for each sample, respectively.
#'
#' @return A list of two elements.
#'   \describe{
#'   \item{p}{A ggplot object, showing percentages of fragments assigned to
#'   different genomic features, such as genic regions, intergenic regions,
#'   exonic regions, intronic regions, rRNA genes, mitochondrial genome,
#'   chloroplast genome (only for plants)}
#'   \item{IR_rate}{a data frame containing the percentages of reads mapping to
#'   intergenic regions, which can be used for the "IR%" method for correcting
#'   for gDNA contamination}
#'   }
#'
#' @export
#' @importFrom graphics pairs par smoothScatter strwidth text
#' @importFrom grDevices dev.off pdf
#' @importFrom stats cor dist median prcomp
#' @examples
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#' metadata <- counts_summary$metadata
#' counts_summary$metadata <- NULL
#' p <- check_read_distribution(
#'     featurecounts_list = counts_summary,
#'     metadata = metadata
#' )
#' p$p

check_read_distribution <-
    function(featurecounts_list = NULL,
             metadata = data.frame(
                 sample_name = vector(mode = "character"),
                 group = vector(mode = "character")
             )) {
        if (!is.list(featurecounts_list) ||
            any(!c("counts", "annotation", "stat") %in%
                names(featurecounts_list[[1]]))) {
            stop(
                "featurecounts_list must be a list of lists, each sublist",
                " containing featureCounts output for different genomic ",
                "features, such as output ",
                "of the summarize_reads function."
            )
        }
        if (!is.data.frame(metadata)) {
            stop("metadata must be a data frame.")
        } else if (any(!c("sample_name", "group") %in% colnames(metadata))) {
            stop(
                "metadata is a data frame, but it does not contain columns:",
                " sample_name or group."
            )
        }


        assigned_per_region <- mapply(
            function(.x, .y) {
                assigned <- .x$stat[.x$stat$Status == "Assigned", ,
                                    drop = FALSE]
                assigned <- assigned[, -1, drop = FALSE]
                assigned <- as.data.frame(t(assigned))
                colnames(assigned) <- "assigned_count"
                assigned$sample_name <- rownames(assigned)
                rownames(assigned) <- NULL
                assigned$region_type <- .y
                assigned
            },
            featurecounts_list[-length(featurecounts_list)],
            ## no GTF-base summary
            names(featurecounts_list)[-length(featurecounts_list)],
            SIMPLIFY = FALSE
        )

        assigned_per_region <- do.call("rbind", assigned_per_region)
        total_mapped_frags <-
            colSums(featurecounts_list$gtf$stat[, -1, drop = FALSE])
        total_mapped_frags <-
            rep(total_mapped_frags, length(featurecounts_list) - 1)
        assigned_per_region$assigned_percent <-
            assigned_per_region$assigned_count / total_mapped_frags * 100

        levels <- c(
            "gene",
            "exon",
            "intronic_region",
            "intergenic_region",
            "rRNA",
            "mitochondrion"
        )
        labels <- c(
            "Gene",
            "Exon",
            "Intron",
            "Intergenic region",
            "rRNA",
            "Mitochondrion"
        )
        if (length(unique(assigned_per_region$region_type)) == 7) {
            levels <- c(levels, "chloroplast")
            labels <- c(labels, "Chloroplast")
        }
        assigned_per_region$region_type <-
            factor(assigned_per_region$region_type,
                levels = levels,
                labels = labels
            )
        if (!all(sort(as.character(metadata$sample_name)) ==
            sort(unique(assigned_per_region$sample_name)))) {
            stop(
                "sample names in metadata are not consistent with those in ",
                "featurecounts_list."
            )
        }

        assigned_per_region <- merge(assigned_per_region,
            metadata[, c("sample_name", "group")],
            by = "sample_name",
            all.x = TRUE,
            sort = FALSE
        )

        assigned_per_region <-
            assigned_per_region[order(
                assigned_per_region$group,
                assigned_per_region$sample_name
            ), ]
        assigned_per_region$group <- factor(assigned_per_region$group,
            levels = unique(assigned_per_region$group)
        )
        assigned_per_region$sample_name <-
            factor(assigned_per_region$sample_name,
                levels = unique(assigned_per_region$sample_name)
            )


        ## read distribution among different genomic regions
        p <- ggplot(
            assigned_per_region,
            aes(
                x = sample_name,
                y = assigned_percent,
                color = group
            )
        ) +
            geom_point() +
            xlab("Sample") +
            ylab("Assigned reads (%)") +
            facet_wrap(~region_type) +
            guides(color = guide_legend(title = "Group")) +
            theme(
                text = element_text(size = 8),
                axis.text.x = element_text(
                    angle = 90,
                    vjust = 0.5,
                    hjust = 0,
                    size = 8
                )
            )
        IR_rate <-
            assigned_per_region[assigned_per_region$region_type ==
                "Intergenic region", ]
        rownames(IR_rate) <- IR_rate$sample_name
        IR_rate <- IR_rate[metadata$sample_name, ]
        IR_rate <- IR_rate[, "assigned_percent", drop = FALSE]
        colnames(IR_rate) <- "IR_rate"
        list(p = p, IR_rate = IR_rate)
    }

#' Check sample correlation
#'
#' Generate a panel of plots based on a count table, with the diagonal
#' showing the sample names, the lower triangle showing smoothed scatterplots
#' for gene expression of pairwise samples, and the upper triangle showing
#' Pearson correlation coefficients of gene expression of pairwise samples.
#'
#' @param counts A numeric matrix or data frame containing gene expression count
#'   data from an RNA-seq experiment, with row names for gene ID and column
#'   names for sample names. If values are not integers, it will be rounded.
#' @return NULL. A plot with pairwise scatter plots and Pearson correlation
#'   coefficients is generated. When the sample size is big, save it as tiff
#'   file.
#'
#' @export
#' @importFrom methods is
#' @importFrom graphics smoothScatter
#' @import KernSmooth
#' @examples
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#'
#' check_sample_correlation(counts = counts_summary$gtf$counts)
#'

check_sample_correlation <- function(counts = NULL) {
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop("A count table must be provied as a matrix or data frame")
    }
    if (is.data.frame(counts)) {
        counts <- as.matrix(counts)
    }
    if (!is.numeric(counts)) {
        stop("The count table must be numeric")
    }

    ## scatter plot and correlation plot
    counts <- counts[rowSums(counts) > 0, ]
    counts <- log2(counts + 1)
    counts <- counts[, order(colnames(counts))]

    panel.cor <-
        function(x,
                 y,
                 digits = 2,
                 prefix = "",
                 cex.cor,
                 ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            Cor <- cor(x, y) # Remove abs function if desired
            txt <-
                paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
            if (missing(cex.cor)) {
                cex.cor <- 0.4 / strwidth(txt)
            }
            text(0.5, 0.5, txt,
                cex = 1 + cex.cor * Cor
            )
        }

    panel.smth <- function(x, y) {
        smoothScatter(x, y, add = TRUE)
    }

    pairs(
        counts,
        lower.panel = panel.smth,
        upper.panel = panel.cor,
        font.labels = 1,
        row1attop = TRUE,
        gap = 0.2,
        log = ""
    )
}


#' Expression distribution
#'
#' Compare expression distribution by boxplot, density plot and empirical
#' cumulative distribution plot
#'
#' @param counts A numeric matrix or data frame containing gene expression count
#'   data for an RNA-seq experiment, with row names for gene ID and column names
#'   for sample names. If values are not integers, it will be rounded. You can
#'   use the count table outputted by the [summarize_reads()] function with
#'   a GTF file as annotation.
#' @param normalization A character(1) vector, specifying a between-sample
#'   normalization methods: DESeq2's  median of ratios method, smooth quantile
#'   normalization method [qsmooth::qsmooth()], or `none`.
#' @param metadata A data frame including column names: sample_name and group
#'   for sample names and experimental groups of samples, respectively. The
#'   order of the sample name in the rows must match those in the count table
#'   specified by `counts`.
#' @return A list of 3 ggplot objects.
#' \describe{
#'   \item{box_plot}{boxplots showing DESeq2-normalized gene-level count
#'                  distribution, on a log scale}
#'   \item{density_plot}{density plots showing DESeq2-normalized gene-level
#'                       count distribution, on a log scale}
#'   \item{ecd_plot}{plots showing empricial cumulative distribution of
#'                  fraction  of genes with CPM greater than or equal to a
#'                  given CPM on a log scale}
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom edgeR filterByExpr cpm
#' @importFrom qsmooth qsmooth qsmoothData
#' @importFrom reshape2 melt
#' @examples
#' library(patchwork)
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#' metadata <- counts_summary$metadata
#' metadata$group <- gsub("CD1A\\(-\\)", "CD1A_N", metadata$group)
#' metadata$group <- gsub("CD1A\\(\\+\\)", "CD1A_P", metadata$group)
#'
#' p <- check_expr_distribution(
#'     counts = counts_summary$gtf$counts,
#'     normalization = "DESeq2",
#'     metadata = metadata
#' )
#' wrap_plots(p, ncol = 1, nrow = 3, widths = 6, heights = 8)

check_expr_distribution <-
    function(counts = NULL,
             normalization = c("DESeq2", "qsmooth", "none"),
             metadata =
                 data.frame(
                     sample_name = vector(mode = "character"),
                     group = vector(mode = "character")
                 )) {
        if (!is.matrix(counts) && !is.data.frame(counts)) {
            stop("A count table must be provied as a matrix or data frame")
        }
        if (is.data.frame(counts)) {
            counts <- as.matrix(counts)
        }

        normalization <- match.arg(normalization)

        # verify the metadata file
        if (!is.data.frame(metadata)) {
            stop("metadata must be a data frame.")
        } else if (any(!c("sample_name", "group") %in%
            colnames(metadata))) {
            stop("metadata must contain columns: sample_name, and group!")
        }

        if (!is.numeric(counts)) {
            stop("The count table must be numeric.")
        }
        if (!all(sort(colnames(counts)) ==
                 sort(as.character(metadata$sample_name)))) {
            stop(
                "Column names of the raw count matrix DO NOT match ",
                "sample names in the metadata!"
            )
        } else {
            counts <- counts[, metadata$sample_name]
        }
        metadata <-
            metadata[order(metadata$group, metadata$sample_name), ]
        metadata$group <- factor(metadata$group,
            levels = unique(metadata$group)
        )
        metadata$sample_name <- factor(metadata$sample_name,
            levels = metadata$sample_name
        )

        keep <- filterByExpr(counts, group = metadata$group)
        counts <- counts[keep, ]

        if (normalization == "qsmooth") {
            # library size normalization
            lib_sizes <- colSums(counts)
            geometric_mean <- exp(mean(log(lib_sizes)))
            size_factors <- lib_sizes / geometric_mean
            counts <- sweep(counts, 2, size_factors, FUN = "/")
            counts_qs <- qsmooth(
                object = counts,
                group_factor = metadata$group
            )
            counts_1 <- qsmoothData(counts_qs)
        } else if (normalization == "DESeq2") {
            ## use DESeq2 normalization method instead
            dds <-
                DESeqDataSetFromMatrix(
                    countData = as.matrix(round(counts)),
                    colData = metadata,
                    design = ~ 0 + group
                )
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds,
                fitType = "local",
                maxit = 1000
            )

            counts_1 <- counts(dds, normalized = TRUE)
        } else {
            counts_1 <- counts
        }

        counts_1_log <- as.data.frame(log2(counts_1 + 1))
        counts_1_log$GeneID <- rownames(counts_1_log)
        counts_1_log_long <- melt(
            counts_1_log,
            id.vars = "GeneID",
            variable.name = "Sample",
            value.name = "log2Count"
        )
        counts_1_log_long <- merge(
            counts_1_log_long,
            metadata[, c("sample_name", "group")],
            by.x = "Sample",
            by.y = "sample_name",
            all.x = TRUE,
            sort = FALSE
        )
        counts_1_log_long <-
            counts_1_log_long[order(
                counts_1_log_long$group,
                counts_1_log_long$Sample
            ), ]
        counts_1_log_long$group <- factor(counts_1_log_long$group,
            levels = unique(counts_1_log_long$group)
        )
        counts_1_log_long$Sample <- factor(counts_1_log_long$Sample,
            levels = unique(counts_1_log_long$Sample)
        )
        p_box <- ggplot(
            counts_1_log_long,
            aes(
                x = Sample,
                y = log2Count,
                fill = group
            )
        ) +
            geom_boxplot() +
            ylab(expression(log[2](counts + 1))) +
            xlab("Sample") +
            guides(fill = guide_legend(title = "Group")) +
            ggtitle("Boxplot-normalized gene expression") +
            theme(
                text = element_text(size = 8),
                axis.text.x = element_text(
                    angle = 90,
                    vjust = 0.5,
                    hjust = 0,
                    size = 8
                ),
                plot.title = element_text(hjust = 0.5),
                legend.key.size = unit(0.5, "cm"),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )
        if (length(unique(metadata$sample_name)) > 12) {
            p_density <- ggplot(
                counts_1_log_long,
                aes(
                    x = log2Count,
                    group = Sample,
                    color = group
                )
            ) +
                guides(color = guide_legend(title = "Group"))
        } else {
            p_density <- ggplot(
                counts_1_log_long,
                aes(
                    x = log2Count,
                    group = Sample,
                    color = Sample
                )
            ) +
                guides(color = guide_legend(title = "Sample"))
        }

        p_density <- p_density +
            geom_density(linewidth = 0.8) +
            ylab("Density") +
            xlab(expression(log[2](counts + 1))) +
            ggtitle("Density plot-normalized counts") +
            theme(
                text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )

        cpms_0 <- as.data.frame(cpm(counts))
        cpms_log <- log2(cpms_0 + 1 / 10 * min(cpms_0[cpms_0 != 0]))
        cpms_log$GeneID <- rownames(cpms_log)
        cpms_log_long <- melt(
            cpms_log,
            id.vars = "GeneID",
            variable.name = "Sample",
            value.name = "log2CPM"
        )
        cpms_log_long <- merge(
            cpms_log_long,
            metadata[, c("sample_name", "group")],
            by.x = "Sample",
            by.y = "sample_name",
            all.x = TRUE,
            sort = FALSE
        )
        cpms_log_long <- cpms_log_long[order(
            cpms_log_long$group,
            cpms_log_long$Sample
        ), ]
        cpms_log_long$group <- factor(cpms_log_long$group,
            levels = unique(cpms_log_long$group)
        )
        cpms_log_long$Sample <- factor(cpms_log_long$Sample,
            levels = unique(cpms_log_long$Sample)
        )

        if (length(unique(metadata$sample_name)) > 12) {
            p_ecdf <- ggplot(
                cpms_log_long,
                aes(
                    x = log2CPM,
                    group = Sample,
                    colour = group
                )
            ) +
                guides(color = guide_legend(title = "Group"))
        } else {
            p_ecdf <- ggplot(
                cpms_log_long,
                aes(
                    x = log2CPM,
                    group = Sample,
                    colour = Sample
                )
            ) +
                guides(color = guide_legend(title = "Sample"))
        }
        p_ecdf <- p_ecdf +
            stat_ecdf(geom = "step", linewidth = 0.8) +
            xlab(bquote(log[2] * CPM)) +
            ylab("Proportion") +
            ggtitle("Empirical cumulative distributions") +
            theme(
                text = element_text(size = 8),
                legend.key.size = unit(0.5, "cm"),
                plot.title = element_text(hjust = 0.5),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6)
            )

        list(
            box_pot = p_box,
            density_plot = p_density,
            ecd_plot = p_ecdf
        )
    }


#' Check the percentage of genes with counts greater than minimal CPM
#' @param metadata A data frame with column names: sample_name, and group.
#'   The *sample_name* and *group* columns contain the unique sample labels,
#'   and the experimental conditions for each sample. The order of the two
#'   columns doesn't matter.
#' @param counts A numeric matrix or data frame containing gene expression count
#'   data for an RNA-seq experiment, with row names for gene ID and column names
#'   for sample names. If values are not integers, it will be rounded.
#' @param abundance A numeric matrix or data frame containing gene expression
#'   data in TPM. for an RNA-seq experiment, with row names for gene ID and
#'   column names for sample names. It can be the sublist of the [salmon_res()]
#'   output.
#' @param min_cpm A numeric(1), minimal CPM threshold.
#' @param min_tpm A numeric(1), minimal TPM threshold.
#' @return A ggplot object if `counts` is not `NULL`, showing percentages of
#'   genes with counts above the user-specified minimal CPM (count per million)
#'   in each sample. Or a ggplot object of two panels if both `counts` and
#'   `abundance` are not `NULL`, showing percentages of genes  with counts
#'   above the user-specified minimal CPM and minimal TPM (transcript per
#'   million) in each sample.
#' @details
#'  The axis title contains unicode, so please output the plot in the svg format
#'  using the svglite package, instead of [grDevices::pdf()], for high-
#'  resolution plot.
#' @return A ggplot object.
#'
#' @export
#' @examples
#' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})
#' load(file.path(tmp_dir, "read_count_summary.RData"))
#'
#' p <- check_expressed_gene_percentage(
#'     metadata = counts_summary$metadata,
#'     counts = counts_summary$gtf$counts,
#'     min_cpm = 1
#' )
#' p

check_expressed_gene_percentage <-
    function(metadata =
                 data.frame(
                     sample_name = vector(mode = "character"),
                     group = vector(mode = "character")
                 ),
             counts = NULL,
             min_cpm = 1,
             abundance = NULL,
             min_tpm = 1) {
        if (!is.data.frame(metadata)) {
            stop("metadata must be a data frame.")
        } else if (any(!c("sample_name", "group") %in% colnames(metadata))) {
            stop(
                "Not all column names: sample_name and group are included ",
                "in the metadata."
            )
        }

        if (!is.null(counts) &&
            !is.matrix(counts) && !is.data.frame(counts)) {
            stop(
                "A count table must be provied as a matrix or",
                "data frame if not NULL."
            )
        }

        if (!is.null(abundance) && !is.matrix(abundance) &&
            !is.data.frame(abundance)) {
            stop(
                "An abundance table must be provied as a",
                "matrix or data frame if it not NULL."
            )
        }
        if (is.null(counts) && is.null(abundance)) {
            stop("Either counts or abundance must be provided.")
        }
        if (is.data.frame(counts)) {
            counts <- as.matrix(counts)
        }

        if (!is.numeric(min_cpm) || min_cpm <= 0) {
            stop("min_cpm must be a single positive number")
        }
        if (!is.numeric(min_tpm) || min_tpm <= 0) {
            stop("min_tpm must be a single positive number")
        }

        ## filter all zero genes
        # counts <- counts[rowSums(counts)> 0, ]
        if ((!is.null(counts) &&
            any(!colnames(counts) %in% metadata$sample_name)) ||
            (!is.null(abundance) &&
                any(!colnames(abundance) %in% metadata$sample_name))) {
            stop("colnames of expression data DON'T match those in metadata!")
        }

        if (!is.null(counts)) {
            cpms <- as.data.frame(cpm(counts))
            pct_cpm_gt1 <- vapply(
                cpms,
                function(.x) {
                    sum(.x >= min_cpm)
                },
                numeric(1)
            ) / nrow(counts) * 100

            pct_cpm_gt1_df <- data.frame(
                sample_name = names(pct_cpm_gt1),
                percent = pct_cpm_gt1
            )
            pct_cpm_gt1_df <- merge(pct_cpm_gt1_df,
                metadata,
                by = "sample_name",
                all.x = TRUE,
                sort = FALSE
            )
            pct_cpm_gt1_df <- pct_cpm_gt1_df[order(
                pct_cpm_gt1_df$group,
                pct_cpm_gt1_df$sample
            ), ]
            pct_cpm_gt1_df$group <- factor(pct_cpm_gt1_df$group,
                levels = unique(pct_cpm_gt1_df$group)
            )
            pct_cpm_gt1_df$sample_name <-
                factor(pct_cpm_gt1_df$sample_name,
                    levels = unique(pct_cpm_gt1_df$sample_name)
                )
        }

        if (!is.null(abundance)) {
            if (is.matrix(abundance)) {
                abundance <- as.data.frame(abundance)
            }
            ## filter out all zero abundance
            # abundance <- abundance[rowSums(abundance)> 0, ]
            pct_tpm_gt1 <- vapply(
                abundance,
                function(.x) {
                    sum(.x >=
                        min_tpm)
                },
                numeric(1)
            ) / nrow(abundance) * 100

            pct_tpm_gt1_df <- data.frame(
                sample_name = names(pct_tpm_gt1),
                percent = pct_tpm_gt1
            )

            pct_tpm_gt1_df <- merge(pct_tpm_gt1_df,
                metadata,
                by = "sample_name",
                all.x = TRUE,
                sort = FALSE
            )
            pct_tpm_gt1_df <- pct_tpm_gt1_df[order(
                pct_tpm_gt1_df$group,
                pct_tpm_gt1_df$sample
            ), ]
            pct_tpm_gt1_df$group <- factor(pct_tpm_gt1_df$group,
                levels = unique(pct_tpm_gt1_df$group)
            )
            pct_tpm_gt1_df$sample_name <-
                factor(pct_tpm_gt1_df$sample_name,
                    levels = unique(pct_tpm_gt1_df$sample_name)
                )
        }

        if (!is.null(counts) && !is.null(abundance)) {
            pct_tpm_gt1_df$type <- "TPM"
            pct_cpm_gt1_df$type <- "CPM"
            cpm_tpm <- rbind(pct_cpm_gt1_df, pct_tpm_gt1_df)

            cpm_tpm$sample_name <-
                factor(cpm_tpm$sample_name,
                    levels = unique(cpm_tpm$sample_name)
                )


            p <- ggplot(cpm_tpm, aes(
                x = sample_name,
                y = percent,
                color = group
            )) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 threshold)")) +
                facet_wrap(~type) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        } else if (!is.null(counts)) {
            p <- ggplot(
                pct_cpm_gt1_df,
                aes(
                    x = sample_name,
                    y = percent,
                    color = group
                )
            ) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 ", min_cpm, " CPM)")) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        } else {
            p <- ggplot(
                pct_tpm_gt1_df,
                aes(
                    x = sample_name,
                    y = percent,
                    color = group
                )
            ) +
                geom_point() +
                xlab("Sample") +
                ylab(paste0("%Gene (\U2265 ", min_cpm, " TPM)")) +
                theme(
                    text = element_text(size = 8),
                    axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 0,
                        size = 6
                    )
                )
        }
        return(p)
    }


#' Convert genomic features from GRanges to Simplified Annotation Format (SAF)
#'
#' @param granges An object of [GenomicRanges::GRanges-class]
#' @export
#' @return A data frame with columns: "GeneID", "Chr", "Start", "End","Strand".
#' @examples
#' library("GenomicRanges")
#' gr0 <- GRanges(Rle(
#'     c("chr2", "chr2", "chr1", "chr3"),
#'     c(1, 3, 2, 4)
#' ), IRanges(seq_len(10), width = 10:1))
#' saf <- granges_to_saf(gr0)
#'
granges_to_saf <- function(granges) {
    if (!is(granges, "GRanges")) {
        stop(deparse(substitute(granges)), " is not a GRanges object.")
    }
    saf <- as.data.frame(granges)[, c(seq_len(3), 5)]
    saf$GeneID <- rownames(saf)
    rownames(saf) <- NULL
    colnames(saf) <- c("Chr", "Start", "End", "Strand", "GeneID")
    saf_colnames <- c("GeneID", "Chr", "Start", "End", "Strand")
    saf <- saf[, saf_colnames]
    saf
}


#' Download files via HTTP with retry if failure occurs
#'
#' @details This function is adopted from David Weber's post on the
#' StackOverflow with modification:
#' https://stackoverflow.com/questions/63340463/download-files-until-it-works.
#'
#' @param expr A R expression using [utils::download.file()] to download a file.
#' @param isError A function to determine if a *try-error* is returned
#' @param maxErrors An integer(1), specifying the maximal number of tries before
#'   stopping retry.
#' @param sleep A numeric(1), specifying the time period in second to wait before
#'   retry.
#' @return 0 if success, otherwise a *try-error*.
#' @export
#'
#' @examples
#' #' tmp_dir <- tempdir()
#' options(timeout = max(3000, getOption("timeout")))
#' ## download feaureCounts results
#' count_url <- paste0(
#'     "https://zenodo.org/records/11458839/files/",
#'     "read_count_summary.RData?download=1"
#' )
#' retry_download({download.file(
#'     url = count_url,
#'     destfile = file.path(tmp_dir, "read_count_summary.RData"),
#'     mode = "wb"
#' )})

retry_download <- function(expr,
                           isError = function(x) {"try-error" %in% class(x)},
                           maxErrors = 10, sleep = 5) {
    attempts <- 0
    retval <- try(eval(expr))
    while (isError(retval)) {
        attempts <- attempts + 1
        if (attempts >= maxErrors) {
            stop("Already tried too many times. ",
                 "The network is too busy to download the file!")
        }
        if (sleep > 0) Sys.sleep(sleep)
        retval <- try(eval(expr))
    }
    invisible(retval)
}
