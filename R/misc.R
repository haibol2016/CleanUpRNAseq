#' @include R6_class.R
NULL

#' Convert a BSgenome from the UCSC to the Ensembl style
#'
#' Convert a BSgenome from the UCSC style to the Ensembl style, changing
#' "chrM" to "MT", removing the "chr" prefix from chromosome/scaffold seqnames,
#' only keeping seqnames in the primary genome assembly, and changing the genome
#' version name form the UCSC name to the Ensembl name, such as from "hg38" to
#' "GRCh38".
#'
#' @param UCSC_BSgenome An object of the [BSgenome::BSgenome-class] in the UCSC
#'   style, such as BSgenome.Hsapiens.UCSC.hg38.
#' @param genome_version A character(1), specifying the genome version, such as
#'   "GRCh38".
#' @param seqname_alias A data frame or a tab-delimited file to a data frame,
#'   with two columns: `ucsc` and `ensembl` for UCSC-style seqnames and
#'   Ensembl-style seqnames, respectively.
#'
#' @return An object of the [BSgenome::BSgenome-class] in the Ensembl style,
#'   such as BSgenome.Dvirilis.Ensembl.dvircaf1.
#' @importFrom  stats complete.cases setNames
#' @importFrom GenomeInfoDb seqnames seqnames<-
#' @export
#'
#' @examplesIf require("BSgenome.Hsapiens.UCSC.hg38")
#' ucsc_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#' ucsc_seqnames <-
#'     seqnames(ucsc_BSgenome)[!grepl("_alt|fix|hap\\d+",
#'                             seqnames(ucsc_BSgenome))]
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
#'
#' bsgenome <- style_BSgenome(
#'     UCSC_BSgenome = ucsc_BSgenome,
#'     genome_version = "GRCh38",
#'     seqname_alias = seqname_alias
#' )

style_BSgenome <-
    function(
        UCSC_BSgenome = NULL,
        genome_version = "GRCh38",
        seqname_alias =
            data.frame(
                ucsc = paste0("chr", c(seq_len(22), "X", "Y", "M")),
                ensembl = c(seq_len(22), "X", "Y", "MT")
            )) {
        if (!is(UCSC_BSgenome, "BSgenome")) {
            stop("UCSC_BSgenome is not a BSgenome!")
        }

        stopifnot(length(genome_version) == 1, is.character(genome_version))

        if (!is.data.frame(seqname_alias)) {
            if (!file.exists(seqname_alias)) {
                stop("seqname_alias is provided by a non-existing file!")
            }
            seqname_alias <-
                read.delim(seqname_alias, header = TRUE, as.is = TRUE)
        }

        if (!all(c("ucsc", "ensembl") %in% colnames(seqname_alias))) {
            stop('column names must be "ucsc" and "ensembl"!')
        }
        if (any(!complete.cases(seqname_alias)) ||
            any(duplicated(seqname_alias$ucsc)) ||
            any(duplicated(seqname_alias$ensembl))) {
            stop("missing or duplicated values are not allowed in",
                 " seqname_alias!")
        }
        if (any(!seqname_alias$ucsc %in% seqnames(UCSC_BSgenome))) {
            stop("some seqname in seqname_alias is not in UCSC_BSgenome!")
        }

        keep_BSgenome_Sequences <- function(genome, seqnames) {
            stopifnot(all(seqnames %in% seqnames(genome)))
            genome@user_seqnames <- setNames(seqnames, seqnames)
            genome@seqinfo <- genome@seqinfo[seqnames]
            genome
        }

        ## remove seqnames not in primary assembly
        sequences_to_keep <-
            seqnames(UCSC_BSgenome)[!grepl(
                "_(alt|fix|hap\\d+)$",
                seqnames(UCSC_BSgenome)
            )]
        seqname_alias <-
            seqname_alias[seqname_alias$ucsc %in% sequences_to_keep, ]
        UCSC_BSgenome <-
            keep_BSgenome_Sequences(UCSC_BSgenome, seqname_alias$ucsc)
        seqnames(UCSC_BSgenome) <- seqname_alias$ensembl
        UCSC_BSgenome@seqinfo@seqnames <- seqname_alias$ensembl

        # Don't change the names of user_seqnames, otherwise cause bug when
        # call BSgenome::getSeq.
        # names(UCSC_BSgenome@user_seqnames) <- seqname_alias$ensembl
        UCSC_BSgenome@seqinfo@genome <-
            rep(genome_version, length(UCSC_BSgenome@seqinfo@genome))
        UCSC_BSgenome@metadata$genome <- genome_version

        UCSC_BSgenome
    }

#' Convert genomic features from GRanges to Simplified Annotation Format (SAF)
#'
#' @param granges An object of [GenomicRanges::GRanges-class]
#' @return A data frame with columns: "GeneID", "Chr", "Start", "End","Strand".
#' @noRd
#' @examples
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


#' Generate a SAF file for genomic features
#'
#' Generate a SAF (simplified annotation format) file for genomic features:
#' genicregions, intergenic regions, exonic regions, intronic regions, rRNA
#' genes, mitochondrial genome, and chloroplast genome (only for plants).
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
#' @importFrom GenomicRanges reduce setdiff
#' @importFrom methods as
#' @importFrom Rsamtools seqinfo BamFile
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
#' require("ensembldb")
#' tmp_dir <- tempdir()
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "CleanUpRNAseq")
#' hs_ensdb_sqlite <-
#'     ensembldb::ensDbFromGtf(
#'         gtf = gtf,
#'         outfile = file.path(tmp_dir, "EnsDb.hs.v110.sqlite"),
#'         organism = "Homo_Sapiens",
#'         genomeVersion = "GRCh38",
#'         version = 110
#'     )
#' bam_file <- system.file("extdata", "K084CD7PCD1N.srt.bam",
#'     package = "CleanUpRNAseq"
#' )
#' saf_list <- get_saf(
#'     ensdb_sqlite = hs_ensdb_sqlite,
#'     bamfile = bam_file,
#'     mitochondrial_genome = "MT"
#' )

get_saf <- function(ensdb_sqlite = NULL,
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
        reduce(genes(ensdb, filter = SeqNameFilter(rownames(genome_info))))
    gene_saf <- granges_to_saf(gene_GR)

    exon_GR <-
        reduce(exons(ensdb, filter = SeqNameFilter(rownames(genome_info))))
    exon_saf <- granges_to_saf(exon_GR)

    ## chromosome/scaffolds GRanges

    genome_GR <- as(data.frame(
        seqnames = rownames(genome_info),
        start = 1L,
        end = genome_info$seqlengths,
        strand = "*"
    ), "GRanges")

    ## intronic regions
    intronic_GR <- setdiff(gene_GR, exon_GR)
    intronic_saf <- granges_to_saf(intronic_GR)
    
    ## intergenic regions
    intergenic_GR <- setdiff(genome_GR, gene_GR, ignore.strand = TRUE)
    intergenic_saf <- granges_to_saf(intergenic_GR)

    ## rRNA transcripts
    rrna_GR <- exons(ensdb, filter = GeneBiotypeFilter("rRNA"))

    ## in case no rRNA annotation in the GTF
    if (length(rrna_GR) > 0) {
        rrna_saf <- granges_to_saf(rrna_GR)
        saf_list <- list(
            gene = gene_saf,
            exon = exon_saf,
            intergenic_region = intergenic_saf,
            intronic_region = intronic_saf,
            rRNA = rrna_saf
        )
    } else {
        saf_list <- list(
            gene = gene_saf,
            exon = exon_saf,
            intergenic_region = intergenic_saf
        )
    }

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
#' @param SummarizedCounts An object of [SummarizedCounts].
#' @param isPairedEnd A logical(1), whether the RNA-seq data is paired-end.
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
#'   format, such as the output of the [get_saf()] function.
#' @param gtf A character(1), specifying a path to a GTF file.
#' @param threads An integer(1), number of threads for
#'   [Rsubread::featureCounts()] calling.
#' @param verbose A logical(1) vector, indicating if verbose information for
#' debugging will be generated. This may include information such as unmatched
#' chromosomes/contigs between reads and annotation.
#'
#' @return An object of [SummarizedCounts] with the feature_counts field
#'   populated by modified output from the [Rsubread::featureCounts()]
#'   function for each type of genomic features.
#'
#' @importFrom Rsubread featureCounts
#' @export
#' @examplesIf require("R.utils") && require("ensembldb")
#' tmp_dir <- tempdir()
#' in_dir <- system.file("extdata", package = "CleanUpRNAseq")
#' gtf.gz <- dir(in_dir, ".gtf.gz$", full.name = TRUE)
#' gtf <- file.path(tmp_dir, gsub("\\.gz", "", basename(gtf.gz)))
#' R.utils::gunzip(gtf.gz, destname= gtf,
#'                 overwrite = TRUE, remove = FALSE)
#'
#' in_dir <- system.file("extdata", package = "CleanUpRNAseq")
#' BAM_file <- dir(in_dir, ".bam$", full.name = TRUE)
#' salmon_quant_file <- dir(in_dir, ".sf$", full.name = TRUE)
#' sample_name = gsub(".+/(.+?).srt.bam", "\\1", BAM_file)
#' salmon_quant_file_opposite_strand <- salmon_quant_file
#' col_data <- data.frame(sample_name = sample_name,
#'                        BAM_file = BAM_file,
#'                        salmon_quant_file = salmon_quant_file,
#'                        salmon_quant_file_opposite_strand =
#'                            salmon_quant_file_opposite_strand,
#'                        group = c("CD1N", "CD1P"))
#'
#' sc <- create_summarizedcounts(lib_strand = 0, colData = col_data)
#'
#' bam_file <- system.file("extdata", "K084CD7PCD1N.srt.bam",
#'     package = "CleanUpRNAseq"
#' )
#'
#' tmp_dir <- tempdir()
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "CleanUpRNAseq")
#' hs_ensdb_sqlite <-
#'     ensembldb::ensDbFromGtf(
#'         gtf = gtf,
#'         outfile = file.path(tmp_dir, "EnsDb.hs.v110.sqlite"),
#'         organism = "Homo_Sapiens",
#'         genomeVersion = "GRCh38",
#'         version = 110
#'     )
#' saf_list <- get_saf(
#'     ensdb_sqlite = hs_ensdb_sqlite,
#'     bamfile = bam_file,
#'     mitochondrial_genome = "MT"
#' )
#'
#' counts_list <- summarize_reads(
#'     SummarizedCounts = sc,
#'     saf_list = saf_list,
#'     gtf = gtf
#' )
#'
#'
summarize_reads <- function(SummarizedCounts = NULL,
                            isPairedEnd = TRUE,
                            allowMultiOverlap = FALSE,
                            countMultiMappingReads = TRUE,
                            fraction = TRUE,
                            minMQS = 0,
                            saf_list = NULL,
                            gtf = NULL,
                            threads = 1,
                            verbose = FALSE) {
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))
    if (length(saf_list) < 5 ||
        any(
            !c("gene",
                "exon",
                "intergenic_region",
                "intronic_region"
            ) %in%
            names(saf_list)
        )) {
        stop(
            "A valid SAF list for gene, exon, intergenic region,",
            "and intronic region is needed"
        )
    }

    if (!file.exists(gtf)) {
        stop("A single GTF file is needed")
    }

    bamfiles <- SummarizedCounts$col_data$BAM_file
    sample_name <- SummarizedCounts$col_data$sample_name
    strandSpecific <- SummarizedCounts$lib_strand
    saf_list$gtf <- gtf

    count_res <- mapply(
        function(annotation, isGTF, saflist_name, sample_name) {
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
                juncCounts = FALSE,
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
                tmpDir = tempdir(),
                verbose = verbose
            )
            colnames(res$counts) <- sample_name
            colnames(res$stat) <- c("Status", sample_name)
            if ("counts_junction" %in% names(res)) {
                colnames(res$counts_junction)[9:ncol(res$counts_junction)] <-
                    sample_name
            }
            res$targets <- NULL
            if (!saflist_name %in% c("intergenic_region", "gtf")) {
                res$counts <- NULL
                res$annotation <- NULL
                return(res$stat)
            }
            if (saflist_name == "gtf") {
                res$annotation <- NULL
            }
            return(res)
        },
        saf_list,
        c(rep(FALSE, length(saf_list) - 1), TRUE),
        names(saf_list),
        MoreArgs = list(sample_name = sample_name),
        SIMPLIFY = FALSE
    )
    names(count_res) <- names(saf_list)
    SummarizedCounts$set_feature_counts(count_res)
    count_res
}


#' Import Salmon quantification output into R
#'
#' Import Salmon quantification output into R using tximport
#'
#' @inheritParams summarize_reads
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
#' @return An object of [SummarizedCounts] with the feature_counts field
#'   populated by a list of of matrices of gene-level abundance, counts, and
#'   length outputted by [tximport::tximport()], as described in the folowing:
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
#' @importFrom utils write.table read.delim
#'
#' @examples
#' require("ensembldb")
#' in_dir <- system.file("extdata", package = "CleanUpRNAseq")
#' BAM_file <- dir(in_dir, ".bam$", full.name = TRUE)
#' salmon_quant_file <- dir(in_dir, ".sf$", full.name = TRUE)
#' sample_name = gsub(".+/(.+?).srt.bam", "\\1", BAM_file)
#' salmon_quant_file_opposite_strand <- salmon_quant_file
#' col_data <- data.frame(sample_name = sample_name,
#'                        BAM_file = BAM_file,
#'                        salmon_quant_file = salmon_quant_file,
#'                        salmon_quant_file_opposite_strand =
#'                            salmon_quant_file_opposite_strand,
#'                        group = c("CD1N", "CD1P"))
#'
#' sc <- create_summarizedcounts(lib_strand = 0, colData = col_data)
#'
#' tmp_dir <- tempdir()
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "CleanUpRNAseq")
#' hs_ensdb_sqlite <-
#'     ensembldb::ensDbFromGtf(
#'         gtf = gtf,
#'         outfile = file.path(tmp_dir, "EnsDb.hs.v110.sqlite"),
#'         organism = "Homo_Sapiens",
#'         genomeVersion = "GRCh38",
#'         version = 110
#'     )
#' salmon_counts <- salmon_tximport(
#'     SummarizedCounts = sc,
#'     ensdb_sqlite = hs_ensdb_sqlite
#' )
#'
salmon_tximport <- function(SummarizedCounts = NULL,
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
    stopifnot(is(SummarizedCounts, "SummarizedCounts"))

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

    transcripts <- transcripts(
        ensdb,
        columns = c("gene_id", "tx_id"),
        AnnotationFilter(~ gene_biotype !=
                             filtered_gene_biotypes)
    )
    ## gene-level count
    tx2gene <- as.data.frame(transcripts)[, c("tx_id", "gene_id")]
    col_data <- SummarizedCounts$col_data
    quant_files <- col_data$salmon_quant_file
    sample_name <- col_data$sample_name
    names(quant_files) <- sample_name
    txi <- .import_salmon(quant_files, tx2gene)
    SummarizedCounts$set_salmon_quant(txi)

    ## for stranded RNA-seq data
    if ("salmon_quant_file_opposite_strand" %in%
        colnames(col_data)) {
        opposite_strand_quant_files <-
            col_data$salmon_quant_file_opposite_strand
        names(quant_files) <- sample_name
        txi_opposite <- .import_salmon(opposite_strand_quant_files, tx2gene)
        SummarizedCounts$set_salmon_quant_opposite(txi_opposite)
    }
}

#' Import salmon output to R by collapsing quantification of tx to gene's
#'
#' @param quant_files Salmon `quant.sf` files
#' @param tx2gene transcript-to-gene mapping
#'
#' @return A list of three matrices: `counts`, `abundance`, and `length`
#' @noRd
#'
.import_salmon <- function(quant_files, tx2gene) {
    null <- lapply(quant_files, function(.x) {
        quant_sf <- read.delim(.x, header = TRUE, as.is = TRUE)
        if (any(grepl("[._].+$", quant_sf$Name))) {
            quant_sf$Name <- gsub("[._].+$", "", quant_sf$Name)
            write.table(
                quant_sf,
                file = .x,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE
            )
        }
    })

    txi <- tximport(
        quant_files,
        type = "salmon",
        txOut = FALSE,
        tx2gene = tx2gene
    )
    txi$countsFromAbundance <- NULL
    txi
}

#' Calculate GC content of genomic regions
#'
#' @param region A data frame containing the columns: "Chr", "Start", "End",
#'   and "Strand", such as a data frame contained in a sublist named
#'   *intergenic_region* from the output of the [get_saf()] function,
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
#' @examplesIf require("BSgenome.Hsapiens.UCSC.hg38")
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
#' bsgenome <- style_BSgenome(
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
#' gc_contents <- calc_region_gc(
#'     region = region,
#'     BSgenome = bsgenome,
#'     batch_size = 2000
#' )
#'
calc_region_gc <- function(region = NULL,
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
            ) |>
                letterFrequency("GC", as.prob = TRUE) |>
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
                    getSeq(BSgenome, region[.g], as.character = FALSE) |>
                    letterFrequency("GC", as.prob = TRUE) |>
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
#' @examplesIf interactive() && require("BSgenome.Hsapiens.UCSC.hg38")
#' require("ensembldb")
#' ucsc_BSgenome <- BSgenome.Hsapiens.UCSC.hg38
#' ucsc_seqnames <-
#'     seqnames(ucsc_BSgenome)[!grepl("_alt|fix|hap\\d+",
#'                             seqnames(ucsc_BSgenome))]
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
#'
#' bsgenome <- style_BSgenome(
#'     UCSC_BSgenome = ucsc_BSgenome,
#'     genome_version = "GRCh38",
#'     seqname_alias = seqname_alias
#' )
#'
#'
#' tmp_dir <- tempdir()
#' gtf <- system.file("extdata", "example.gtf.gz",
#'                    package = "CleanUpRNAseq")
#' hs_ensdb_sqlite <-
#'     ensembldb::ensDbFromGtf(
#'         gtf = gtf,
#'         outfile = file.path(tmp_dir, "EnsDb.hs.v110.sqlite"),
#'         organism = "Homo_Sapiens",
#'         genomeVersion = "GRCh38",
#'         version = 110
#'     )
#' gene_gc <- calc_gene_gc(
#'     ensdb_sqlite = hs_ensdb_sqlite,
#'     BSgenome = bsgenome,
#'     batch_size = 2000
#' )
#'
calc_gene_gc <- function(ensdb_sqlite = NULL,
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
    gene_gc <- calc_region_gc(
        region = exons_grlist,
        BSgenome = BSgenome,
        batch_size = batch_size
    )
    gene_gc
}

