#' Exploratory analysis of RNA-seq data
#'
#' Perform sample-level exploratory analysis of RNA-seq data, generating heatmap
#' showing sample distances and PCA plot showing sample variations. Internally,
#' DESeq2 is used for vst transformation of count data.
#'
#' @param counts A numeric matrix of data frame containing gene expression count
#'   data for an RNA-seq experiment, with row names for gene ID and column names
#'   for sample names. If values are not integers, it will be rounded.
#' @param metadata A data frame with column names: sample_name and group for
#'   sample names and experimental groups for each sample, respectively. The
#'   order of the sample name in the rows must match those in the count table
#'   specified by `counts`.
#' @param silent A logical(1), specify whether to draw the plot. It is useful
#'   to set it to FALSE useful when using the gtable output.
#' @return A list of a ggplot object and a [gtable::gtable()] object.
#' \describe{
#'   \item{pca}{A *ggplot* object containing the PCA score plot showing sample
#'              similarity}
#'   \item{heatmap}{A *gtable* object containing the heatmap showing
#'                  pairwise sample distances}
#' }
#' @export
#' @importFrom edgeR cpm
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment assay
#' @import DESeq2 ggplot2
#' @importFrom grDevices rainbow
#'
#' @examples
#' library(patchwork)
#' options(timeout = max(3000, getOption("timeout")))
#' tmp_dir <- tempdir()
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
#' metadata <- counts_summary$metadata
#' metadata$group <- gsub("CD1A\\(-\\)", "CD1AN", metadata$group)
#' metadata$group <- gsub("CD1A\\(\\+\\)", "CD1AP", metadata$group)
#' p <- exploratory_analysis(
#'     counts = counts_summary$gtf$counts,
#'     metadata = metadata
#' )
#'
exploratory_analysis <-
    function(counts = NULL,
             metadata = data.frame(
                 sample_name = vector(mode = "character"),
                 group = vector(mode = "character")
             ),
             silent = FALSE) {
        if (!is.matrix(counts) && !is.data.frame(counts)) {
            stop("A count table must be provied as a matrix or data frame")
        }
        if (is.data.frame(counts)) {
            counts <- as.matrix(counts)
        }
        if (!is.data.frame(metadata)) {
            stop("metadata must be a data frame.")
        } else if (any(!c("sample_name", "group") %in% colnames(metadata))) {
            stop(
                "Not all column names: sample_name and group are included ",
                "in the metadata."
            )
        }

        if (!is.numeric(counts)) {
            stop("The count table must be numeric.")
        }
        if (!all(sort(colnames(counts)) == sort(metadata$sample_name))) {
            stop(
                "colnames of the count table don't match the sample names",
                "in the metadata."
            )
        } else {
            counts <- counts[, metadata$sample_name]
        }

        keep <- filterByExpr(counts, group = metadata$group)
        counts <- counts[keep, ]

        metadata$group <- factor(metadata$group,
            levels = unique(metadata$group)
        )

        dds <-
            DESeqDataSetFromMatrix(
                countData = as.matrix(round(counts)),
                colData = metadata,
                design = ~ 0 + group
            )
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds,
            fitType = "parametric",
            maxit = 1000
        )
        ## exploratory analysis
        vsd <- vst(dds, blind = TRUE)
        sampleDists <- dist(t(assay(vsd)))

        ## Heatmap showing sample distances
        distancePlot <- function(sampleDists, sampleNames, metadata) {
            sampleDistMatrix <- as.matrix(sampleDists)
            rownames(sampleDistMatrix) <- sampleNames
            colnames(sampleDistMatrix) <- sampleNames
            colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

            # add group annotation
            anno <- metadata[, c("sample_name", "group")]
            rownames(anno) <- anno$sample_name
            anno <- anno[, "group", drop = FALSE]
            colnames(anno) <- "Group"
            n_groups <- nlevels(anno$Group)
            if (n_groups <= 20) {
                distinct_cols <- c(
                    "#e6194b",
                    "#3cb44b",
                    "#ffe119",
                    "#4363d8",
                    "#f58231",
                    "#911eb4",
                    "#46f0f0",
                    "#f032e6",
                    "#bcf60c",
                    "#fabebe",
                    "#008080",
                    "#e6beff",
                    "#9a6324",
                    "#fffac8",
                    "#800000",
                    "#aaffc3",
                    "#808000",
                    "#ffd8b1",
                    "#000075",
                    "#808080",
                    "#000000",
                    "#ffffff"
                )[seq_len(n_groups)]
            } else {
                distinct_cols <- rainbow(n_groups)
            }
            names(distinct_cols) <- levels(anno$Group)
            # define the colours
            anno_col <- list(Group = distinct_cols)

            p <- pheatmap(
                sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                annotation_row = anno,
                annotation_col = anno,
                annotation_colors = anno_col,
                fontsize = 6,
                silent = silent,
                main = "Heatmap showing sample distances",
                color = colors
            )
        }

        p_pheatmap <- distancePlot(
            sampleDists = sampleDists,
            sampleNames = vsd$sample_name,
            metadata = metadata
        )

        ## PCA plot showing PC1 and PC2 only
        pca <-
            prcomp(t(assay(vsd)),
                scale = TRUE,
                center = TRUE,
                retx = TRUE
            )
        pc12 <- as.data.frame(pca$x[, seq_len(2)])
        colnames(pc12) <- c("PC1", "PC2")
        pc12 <- cbind(pc12, metadata)
        pc12_var <-
            round(pca$sdev[seq_len(2)]^2 / (sum(pca$sdev^2)) * 100, digits = 2)
        pc12$group <- factor(pc12$group, levels = unique(pc12$group))
        pc12$sample_name <- factor(pc12$sample_name,
            levels = unique(pc12$sample_name)
        )
        p_pca <- ggplot(pc12, aes(
            x = PC1,
            y = PC2,
            color = group,
            label = sample_name
        )) +
            geom_text_repel(size = 2.5, show.legend = FALSE) +
            geom_point() +
            xlab(paste0("PC1 (", pc12_var[1], "%)")) +
            guides(color = guide_legend(title = "Group")) +
            ylab(paste0("PC2 (", pc12_var[2], "%)")) +
            ggtitle("PCA score plot") +
            theme(
                plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = 8),
                axis.title = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8)
            )

        list(pca = p_pca, heatmap = p_pheatmap)
    }
