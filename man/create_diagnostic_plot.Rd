% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_diagnostic_plots.R
\name{create_diagnostic_plot}
\alias{create_diagnostic_plot}
\title{Generate diagnostic plots for detecting gDNA contamination in RNA-seq data}
\usage{
create_diagnostic_plot(
  gtf = NULL,
  metadata = data.frame(sample_name = vector(mode = "character"), group = vector(mode =
    "character"), batch = vector(mode = "character"), BAM_file = vector(mode =
    "character"), salmon_quant_file = vector(mode = "character")),
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
  filtered_gene_biotypes = c("artifact", "TEC", "miRNA", "tRNA", "misc_RNA", "Mt_rRNA",
    "Mt_tRNA", "rRNA", "rRNA_pseudogene", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA",
    "vault_RNA", "TR_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene",
    "IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene")
)
}
\arguments{
\item{gtf}{A character(1), a path to an Ensembl GTF file.}

\item{metadata}{A data frame  with column
names: sample_name, group, BAM_file, and salmon_quant_file. The
\emph{sample_name} and \emph{group} columns contain the unique sample labels and
the experimental conditions for each sample, respectively. The \emph{BAM_file}
column contains the full paths to BAM files. The \emph{salmon_quant_file}
column contains the full paths to \emph{quant.sf} files for quantitation of
gene expression by Salmon pseudo-alignment.
See \url{https://salmon.readthedocs.io/en/latest/}.The order of the columns
doesn't matter.}

\item{normalization}{A character(1) vector, specifying a between-sample
normalization methods: DESeq2's  median of ratios method, smooth quantile
normalization method \code{\link[qsmooth:qsmooth]{qsmooth::qsmooth()}}, or \code{none}.}

\item{ensdb_sqlite}{An EnsDb object or a character(1) vector, specifying a
path to a sqlite file to store an object of the \link[ensembldb:EnsDb-class]{ensembldb::EnsDb}.}

\item{organism_latin_name}{A character(1), Latin name for the organism, with
spaces replaced by underscores "_". It is only require if an EnsDb SQLite
database has to be created from a Ensembl GTF file.}

\item{genome_version}{A character(1), the assembly version of the reference
genome, such as "GRCh38" for the human reference genome version 38. It is
only require if an EnsDb SQLite database has to be created from a Ensembl
GTF file.}

\item{Ensembl_release_version}{An integer(1), the Ensembl release version of
the reference genome with patches, such as 110 for Ensembl Release 110
(July 2023). It is only require if an EnsDb SQLite database has to be
created from a Ensembl GTF file.}

\item{mitochondrial_genome}{A character(1), mitochondrial genome name in the
EnsDb database (ie, the mitochondrial name in column 1 of the GTF file
used to generate the EnsDb SQLite database).}

\item{chloroplast_genome}{A character(1), chloroplast genome name in the
EnsDb database (ie, the chloroplast name in column 1 of the GTF file
used to generate the EnsDb SQLite database). This is only relevant
for plants.}

\item{out_dir}{A character(1), specifying an output directory for diagnostic
plots and the related R objects. If it does not exist, it will be created.}

\item{threads}{An integer(1), number of threads for
\code{\link[Rsubread:featureCounts]{Rsubread::featureCounts()}} calling.}

\item{isPairedEnd}{A logical(1), indicating whether the sequencing data is
paired-end or not.}

\item{strandSpecific}{An integer(1), specifying the strandedness of the
RNA-seq libraries. The possible options are 0 (unstranded), 1 (stranded,
read 1 or single-end read comes from the forward/sense strand) and 2
(reversely stranded, read 1 or single-end read comes from the
reverse/antisense strand).
See \url{https://salmon.readthedocs.io/en/latest/library_type.html} and
\url{https://chipster.csc.fi/manual/library-type-summary.html} for more details
about types of stranded libraries.}

\item{allowMultiOverlap}{A logical(1), indicating if a read is allowed to be
assigned to more than one feature (or meta-feature) if it is found to
overlap with more than one feature (or meta-feature). FALSE by default. A
read (or read pair) will be assigned to the feature (or meta-feature) that
has the largest number of overlapping bases, if the read (or read pair)
overlaps with multiple features (or meta-features).}

\item{countMultiMappingReads}{A logical(1), indicating if multi-mapping
reads/fragments should be counted, TRUE by default. ‘NH’ tag is used to
located multi-mapping reads in the input BAM/SAM files.}

\item{fraction}{A logical(1) indicating if fractional counts are produced
for multi-mapping reads and/or multi-overlapping reads. TRUE by default.}

\item{minMQS}{An integer(1), giving the minimum mapping quality score a read
must satisfy in order to be counted. For paired-end reads, at least one end
should satisfy this criteria. 0 by default.}

\item{verbose}{A logical(1) vector, indicating if verbose information for
debugging will be generated. This may include information such as unmatched
chromosomes/contigs between reads and annotation.}

\item{filtered_gene_biotypes}{A character(n) vector,specifying the biotypes
of genes which will not be considered for downstream gene expression
analysis. By default, genes of the following biotypes are excluded:
"artifact","TEC","miRNA","tRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "rRNA",
"rRNA_pseudogene", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA",
"vault_RNA","TR_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene",
"IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
because their transcripts are too short (< 200 nt) or not likely expressed
at all.}
}
\value{
A list containing the ggplot objects for generating the diagnostic
plots and other objects for gDNA contamination correction.
\describe{
\item{p1_featurecount_read_assignment}{A ggplot object of bar plots showing
read assignment statistics by featureCounts}
\item{p2_mapping_rate_by_features}{A ggplot object of dot plots showing the
percentages of reads mapping to different genomic features: Genes,
Exons, Introns, Intergenic regions, rRNA exons, organelle genome(s)}
\item{p4_boxplot}{A ggplot object of boxplot showing normalized gene
expression distributions}
\item{p4_density_plot}{A ggplot object of density plots showing normalized
gene expression distributions}
\item{p4_ecd_plot}{A ggplot object of empirical cumulative distribution
plots showing normalized gene expression distributions}
\item{p5_percent_expressed_genes}{A ggplot object of dot plots showing the
percentages of genes with observed expression level above a given
threshold}
\item{p6_pca}{A ggplot object of PCA score plots showing sample varibility}
\item{p6_heatmap}{a gtable object containing the heatmap, can be used for
combining the heatmap with other plots}
\item{metadata}{the metadata provided by the user augmented with the
percentages of reads mapping to intergenic regions, IR_rate}
\item{ensdb_sqlite_file}{A path to an EnsDb SQLite database or an EnsDb
object}
\item{saf_list}{an output of the \code{\link[=get_feature_saf]{get_feature_saf()}} function}
\item{featurecounts_summary}{an output of the \code{\link[=summarize_reads]{summarize_reads()}} function}
\item{salmon_summary}{an output of the \code{\link[=salmon_res]{salmon_res()}} function}
}
}
\description{
Generate diagnostic plots for detecting gDNA contamination in RNA-seq data
}
\details{
A set of diagnostic plots will be generated and saved to the
specified directory. ggplot objects for regenerating these plots as
necessary are saved to the same directory as "Diagnostic.plots.objects.RDS"
along with other R objects for correcting for gDNA contamination as a list.
In addition, the metadata augmented with IR rate (IR\%, percentage of reads
mapping to intergenic regions) for each sample is saved to the same
directory as "metadata.with.IR.rates.RDS", which can be used for
generalized linear model-based differential expression analysis.
}
\examples{
\dontrun{
metadata <- read.delim(system.file("extdata",
    "CD1A.RNAseq.metadata.txt",
    package = "CleanUpRNAseq"
))
options(timeout = max(3000, getOption("timeout")))
gtf_url <- paste0(
    "https://ftp.ensembl.org/pub/release-110/gtf/",
    "homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
)
gtf <- basename(gtf_url)
tmp_dir <- tempdir()
retry_download({download.file(
    url = gtf_url,
    destfile = file.path(tmp_dir, gtf),
    mode = "wb"
)})
gtf <- file.path(tmp_dir, gtf)

diagnosis_res <-
    create_diagnostic_plot(
        gtf = gtf,
        metadata = metadata,
        normalization = "DESeq2",
        ensdb_sqlite =
            file.path(tmp_dir, "GRCh38.GTF.EnsDb.sqlite"),
        out_dir = tmp_dir,
        threads = 1
    )
}
}
