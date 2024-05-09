#' Generate diagnostic plots for detecting gDNA contamination in RNA-seq data
#'
#' @param gtf A character(1), a path to an Ensembl GTF file.
#' @param metadata A data frame  with column
#'   names: sample_name, group, BAM_file, and salmon_quant_file. The
#'   *sample_name* and *group* columns contain the unique sample labels and
#'   the experimental conditions for each sample, respectively. The *BAM_file*
#'   column contains the full paths to BAM files. The *salmon_quant_file*
#'   column contains the full paths to *quant.sf* files for quantitation of
#'   gene expression by Salmon pseudo-alignment.
#'   See <https://salmon.readthedocs.io/en/latest/>.The order of the columns
#'   doesn't matter.
#' @param normalization A character(1) vector, specifying a between-sample
#'   normalization methods: DESeq2's  median of ratios method, smooth quantile
#'   normalization method [qsmooth::qsmooth()], or `none`.
#' @param ensdb_sqlite An EnsDb object or a character(1) vector, specifying a
#'   path to a sqlite file to store an object of the [ensembldb::EnsDb-class].
#' @param organism_latin_name A character(1), Latin name for the organism, with
#'   spaces replaced by underscores "_". It is only require if an EnsDb SQLite
#'   database has to be created from a Ensembl GTF file.
#' @param genome_version A character(1), the assembly version of the reference
#'   genome, such as "GRCh38" for the human reference genome version 38. It is
#'   only require if an EnsDb SQLite database has to be created from a Ensembl
#'   GTF file.
#' @param Ensembl_release_version An integer(1), the Ensembl release version of
#'   the reference genome with patches, such as 110 for Ensembl Release 110
#'   (July 2023). It is only require if an EnsDb SQLite database has to be
#'   created from a Ensembl GTF file.
#' @param mitochondrial_genome A character(1), mitochondrial genome name in the
#'   EnsDb database (ie, the mitochondrial name in column 1 of the GTF file
#'   used to generate the EnsDb SQLite database).
#' @param chloroplast_genome A character(1), chloroplast genome name in the
#'   EnsDb database (ie, the chloroplast name in column 1 of the GTF file
#'   used to generate the EnsDb SQLite database). This is only relevant
#'   for plants.
#' @param out_dir A character(1), specifying an output directory for diagnostic
#'   plots and the related R objects. If it does not exist, it will be created.
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
#' @param fraction A logical(1) indicating if fractional counts are produced
#'   for multi-mapping reads and/or multi-overlapping reads. TRUE by default.
#' @param minMQS An integer(1), giving the minimum mapping quality score a read
#'   must satisfy in order to be counted. For paired-end reads, at least one end
#'   should satisfy this criteria. 0 by default.
#' @param threads An integer(1), number of threads for
#'   [Rsubread::featureCounts()] calling.
#' @param verbose A logical(1) vector, indicating if verbose information for
#'   debugging will be generated. This may include information such as unmatched
#'   chromosomes/contigs between reads and annotation.
#' @param filtered_gene_biotypes A character(n) vector,specifying the biotypes
#'   of genes which will not be considered for downstream gene expression
#'   analysis. By default, genes of the following biotypes are excluded:
#'   "artifact","TEC","miRNA","tRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "rRNA",
#'   "rRNA_pseudogene", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA",
#'   "vault_RNA","TR_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene",
#'   "IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
#'   because their transcripts are too short (< 200 nt) or not likely expressed
#'   at all.
#'
#' @details A set of diagnostic plots will be generated and saved to the
#'   specified directory. ggplot objects for regenerating these plots as
#'   necessary are saved to the same directory as "Diagnostic.plots.objects.RDS"
#'   along with other R objects for correcting for gDNA contamination as a list.
#'   In addition, the metadata augmented with IR rate (IR%, percentage of reads
#'   mapping to intergenic regions) for each sample is saved to the same
#'   directory as "metadata.with.IR.rates.RDS", which can be used for
#'   generalized linear model-based differential expression analysis.
#'
#' @return A list containing the ggplot objects for generating the diagnostic
#'   plots and other objects for gDNA contamination correction.
#' \describe{
#'   \item{p1_featurecount_read_assignment}{A ggplot object of bar plots showing
#'         read assignment statistics by featureCounts}
#'   \item{p2_mapping_rate_by_features}{A ggplot object of dot plots showing the
#'         percentages of reads mapping to different genomic features: Genes,
#'         Exons, Introns, Intergenic regions, rRNA exons, organelle genome(s)}
#'   \item{p4_boxplot}{A ggplot object of boxplot showing normalized gene
#'         expression distributions}
#'   \item{p4_density_plot}{A ggplot object of density plots showing normalized
#'         gene expression distributions}
#'   \item{p4_ecd_plot}{A ggplot object of empirical cumulative distribution
#'         plots showing normalized gene expression distributions}
#'   \item{p5_percent_expressed_genes}{A ggplot object of dot plots showing the
#'         percentages of genes with observed expression level above a given
#'         threshold}
#'   \item{p6_pca}{A ggplot object of PCA score plots showing sample varibility}
#'   \item{p6_heatmap}{a gtable object containing the heatmap, can be used for
#'         combining the heatmap with other plots}
#'   \item{metadata}{the metadata provided by the user augmented with the
#'         percentages of reads mapping to intergenic regions, IR_rate}
#'   \item{ensdb_sqlite_file}{A path to an EnsDb SQLite database or an EnsDb
#'         object}
#'   \item{saf_list}{an output of the [get_feature_saf()] function}
#'   \item{featurecounts_summary}{an output of the [summarize_reads()] function}
#'   \item{salmon_summary}{an output of the [salmon_res()] function}
#' }
#' @importFrom grDevices tiff pdf
#' @importFrom patchwork wrap_plots
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
#' }
create_diagnostic_plot <-
    function(gtf = NULL,
             metadata =
                 data.frame(
                     sample_name = vector(mode = "character"),
                     group = vector(mode = "character"),
                     batch = vector(mode = "character"),
                     BAM_file = vector(mode = "character"),
                     salmon_quant_file = vector(mode = "character")
                 ),
             normalization = c("DESeq2", "qsmooth", "none"),
             ensdb_sqlite = NULL,
             organism_latin_name = "Homo_sapiens",
             genome_version = "GRCh38",
             Ensembl_release_version = 110,
             mitochondrial_genome = "MT",
             chloroplast_genome = "Pltd",
             out_dir = NULL,
             threads = 1,
             isPairedEnd = TRUE,
             strandSpecific = 0,
             allowMultiOverlap = FALSE,
             countMultiMappingReads = TRUE,
             fraction = TRUE,
             minMQS = 0,
             verbose = FALSE,
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
        # verify the metadata file
        if (!is.data.frame(metadata))
        {
            stop("metadata must be a data frame.")
        }
        if (any(
            !c(
                "sample_name",
                "group",
                "batch",
                "BAM_file",
                "salmon_quant_file"
            ) %in% colnames(metadata)
        )) {
            stop(
                "metadata must contain columns: sample_name, group, ",
                "batch, BAM_file, and salmon_quant_file!"
            )
        }

        if (any(is.na(metadata)) ||
            any(metadata == "") || any(metadata == " ")) {
            stop("Missing value(s) found in the metadata.")
        }

        if (any(duplicated(metadata$salmon_quant_file)) ||
            any(duplicated(metadata$BAM_file)) ||
            any(duplicated(metadata$sample_name))) {
            stop(
                "Some BAM files, Salmon quant.sf files or sample names are ",
                "not unique."
            )
        }
        if (any(!file.exists(metadata$BAM_file))) {
            stop(
                "Some BAM files do not exist:\n",
                paste(metadata$BAM_file[!file.exists(metadata$BAM_file)],
                    collapse = "\n"
                )
            )
        }

        if (any(!file.exists(metadata$salmon_quant_file))) {
            stop(
                "Some Salmon quantification files do not exist:\n",
                paste(metadata$salmon_quant_file[
                  !file.exists(metadata$salmon_quant_file)],
                    collapse = "\n"
                )
            )
        }

        if (any(grepl("[^A-Za-z0-9_.]", metadata$group)) ||
            any(grepl("^_", metadata$group))) {
            stop(
                "column 'group' in metadata contains special symbols, which ",
                "is not letter, number, dot, or underscore!"
            )
        }

        if (!file.exists(gtf)) {
            stop("A single GTF file should be specified via gtf.")
        }

        if (!dir.exists(out_dir)) {
            dir.create(out_dir, recursive = TRUE)
        }

        # check EnsDb availability
        # if not available, create one.
        if (!is(ensdb_sqlite, "EnsDb") && !is.character(ensdb_sqlite)) {
            stop("Please provide a single EnsDb object or a EnsDb SQLite",
                 " database file!")
        } else if (is.character(ensdb_sqlite)) {
            if (file.exists(ensdb_sqlite)) {
                ensdb <- EnsDb(ensdb_sqlite)
            } else {
                ensdb <-
                    make_ensdb(
                        gtf = gtf,
                        outfile = ensdb_sqlite,
                        organism_latin_name = organism_latin_name,
                        genome_version = genome_version,
                        Ensembl_release_version = Ensembl_release_version
                    )
            }
        } else {
            ensdb <- ensdb_sqlite
        }

        # get SAF for genomic regions: genes, intron, exon, intergenic regions,
        # rRNA, organelle genome
        bamfile <- metadata$BAM_file[1]

        saf_list <- get_feature_saf(
            ensdb_sqlite = ensdb,
            bamfile = bamfile,
            mitochondrial_genome = mitochondrial_genome,
            chloroplast_genome = chloroplast_genome
        )

        # summarize reads by featureCount
        count_res <- summarize_reads(
            metadata = metadata,
            isPairedEnd = isPairedEnd,
            strandSpecific = strandSpecific,
            allowMultiOverlap = allowMultiOverlap,
            countMultiMappingReads = countMultiMappingReads,
            fraction = fraction,
            minMQS = minMQS,
            saf_list = saf_list,
            gtf = gtf,
            threads = threads,
            verbose = verbose
        )
        # combine Salmon quantification results
        salmon_summary <- salmon_res(
            metadata = metadata,
            ensdb_sqlite = ensdb,
            filtered_gene_biotypes = filtered_gene_biotypes
        )


        # check read assignment statistics
        p1_featurecount_read_assignment <-
            check_read_assignment_stat(assignment_stat = count_res$gtf$stat)
        pdf(
            file.path(
                out_dir,
                "Fig1.featureCount.read.assignment.statistics.pdf"
            ),
            height = 3,
            width = 5
        )
        p1_featurecount_read_assignment
        dev.off()

        # check read mapping rate across different genomic features:
        # genes, introns, exons, intergenic regions, rRNA,
        # mitochondrial/plastid genome
        p2_mapping_rate_by_features_IR_rate <-
            check_read_distribution(
                featurecounts_list = count_res,
                metadata = metadata
            )
        pdf(
            file.path(
                out_dir,
                "Fig2.read.mapping.rates.across.different.genomic.features.pdf"
            ),
            height = 6,
            width = 8
        )
        p2_mapping_rate_by_features_IR_rate$p
        dev.off()

        ## add IR_rate to metadata
        IR_rate <- p2_mapping_rate_by_features_IR_rate$IR_rate
        ir <- IR_rate$IR_rate
        names(ir) <- rownames(IR_rate)
        metadata$IR_rate <- ir[metadata$sample_name]

        # metadata <- merge(
        #     metadata,
        #     IR_rate,
        #     by.x = "sample_name",
        #     by.y = "row.names",
        #     all.x = TRUE,
        #     sort = FALSE
        # )

        # check sample correlation
        tiff(
            filename = file.path(out_dir, "Fig3.sample.correlation.tiff"),
            units = "in",
            width = 8,
            height = 8,
            res = 600
        )
        check_sample_correlation(counts = salmon_summary$counts)
        dev.off()

        # check expression distribution
        p4_expression_distribution <-
            check_expr_distribution(
                counts = salmon_summary$counts,
                normalization = normalization,
                metadata = metadata
            )

        pdf(
            file.path(out_dir, "Fig4.expression.distribution.pdf"),
            height = 9,
            width = 6
        )
        wrap_plots(p4_expression_distribution,
            nrow = 3,
            ncol = 1
        )
        dev.off()


        # check perentages of genes with expression above a threshold
        p5_percent_expressed_genes <-
            check_expressed_gene_percentage(
                metadata = metadata,
                counts = salmon_summary$counts,
                min_cpm = 1,
                abundance = salmon_summary$abundance,
                min_tpm = 1
            )
        pdf(
            file.path(out_dir, "Fig5.percentages.of.genes.expressed.pdf"),
            height = 6,
            width = 8
        )
        p5_percent_expressed_genes
        dev.off()

        # check sample similarity and variation by hierarchical clustering
        # and PCA
        p6_pca_heatmap <-
            exploratory_analysis(
                counts = salmon_summary$counts,
                metadata = metadata
            )

        pdf(
            file.path(out_dir, "Fig6.heatmap.showing.sample.similarity.pdf"),
            height = 6,
            width = 8
        )
        p6_pca_heatmap$heatmap
        dev.off()
        pdf(
            file.path(out_dir, "Fig7.PCA.showing.sample.variability.pdf"),
            height = 6,
            width = 8
        )
        p6_pca_heatmap$pca
        dev.off()

        diagnosis_res <-
            list(
                p1_featurecount_read_assignment =
                  p1_featurecount_read_assignment,
                p2_mapping_rate_by_features =
                  p2_mapping_rate_by_features_IR_rate$p,
                p4_boxplot = p4_expression_distribution$box_pot,
                p4_density_plot = p4_expression_distribution$density_plot,
                p4_ecd_plot = p4_expression_distribution$ecd_plot,
                p5_percent_expressed_genes = p5_percent_expressed_genes,
                p6_pca = p6_pca_heatmap$pca,
                p6_heatmap = p6_pca_heatmap$heatmap,
                metadata = metadata,
                ensdb_sqlite_file = ensdb,
                saf_list = saf_list,
                featurecounts_summary = count_res,
                salmon_summary = salmon_summary
            )
        saveRDS(diagnosis_res,
            file = file.path(out_dir, "Diagnostic.plots.objects.RDS")
        )
        saveRDS(metadata, file = file.path(out_dir,
                                           "metadata.with.IR.rates.RDS"))
        diagnosis_res
    }
