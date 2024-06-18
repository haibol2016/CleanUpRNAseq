#' Convert a BSgenome from the UCSC to Ensembl style
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
#'   "GRCh38". Caution: make sure the genome_version here is the same as that
#'   for the [make_ensdb()].
#' @param seqname_alias A data frame or a tab-delimited file to a data frame,
#'   with two columns: ucsc and ensembl for UCSC-style seqnames and Ensembl-
#'   style seqnames, respectively.
#'
#' @return An object of the [BSgenome::BSgenome-class] in the Ensembl style,
#'   such as BSgenome.Dvirilis.Ensembl.dvircaf1.
#' @importFrom  stats complete.cases setNames
#' @importFrom GenomeInfoDb seqnames seqnames<-
#' @export
#'
#' @examples
#' library("BSgenome.Hsapiens.UCSC.hg38")
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
#' bsgenome <- UCSC2Ensembl(
#'     UCSC_BSgenome = ucsc_BSgenome,
#'     genome_version = "GRCh38",
#'     seqname_alias = seqname_alias
#' )
#'
UCSC2Ensembl <-
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

        if (!is.data.frame(seqname_alias)) {
            if (!file.exists(seqname_alias)) {
                stop("seqname_alias is provided by a non-existing file!")
            }
            seqname_alias <-
                read.delim(seqname_alias, header = TRUE, as.is = TRUE)
        }

        if (!all(c("ucsc", "ensembl") %in% colnames(seqname_alias))) {
            stop('column names must be "ucsc", "ensembl"!')
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

#' Split a mutli-fasta file by chromosome/scaffold into individual fasta file
#'
#' Given a multi-fasta file for a reference genome, create a gzip-compressed
#' fasta file for each chromosome/scaffold, such as chr1.fa.gz, chr2.fa.gz,...
#' A prefix "chr" may be added if necessary.
#'
#' @param genome_fasta A character(1) specifying a path for a multi-fasta file
#'   for a reference genome. It can be a gzip-compressed or uncompressed fasta
#'   file.
#' @param out_dir A character(1) specifying a path for output split, compressed
#'   fasta files
#' @param prefix A character(1), specifying chromosome prefix. It can be "none"
#'   or "chr". Default is "none".
#'
#' @return A character string, which is th path to the directory where split,
#'   compressed fasta files are located.
#' @export
#'
#' @examples
#' out_dir <- tempdir()
#' genome_fasta <- file.path(out_dir, "toy.example.fa")
#' in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#' close(in_fasta)
#'
#' fa_dir <- generate_multifasta(
#'     genome_fasta = genome_fasta,
#'     out_dir = out_dir,
#'     prefix = "none"
#' )
#'
generate_multifasta <- function(genome_fasta = NULL,
                                out_dir = tempdir(),
                                prefix = c("none", "chr")) {
    if (!file.exists(genome_fasta)) {
        stop(
            "A path to a reference genome fasta file, genome_fasta,",
            " are required, but it doesn't exist!"
        )
    }

    if (is.null(out_dir)) {
        stop("out_dir is required!")
    }

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    prefix <-
        match.arg(prefix,
            choices = c("none", "chr"),
            several.ok = FALSE
        )
    if (prefix == "none") {
        prefix <- ""
    }

    if (grepl(".(fa|fasta).gz$", genome_fasta)) {
        in_fasta <- gzfile(genome_fasta, open = "rt")
    } else if (grepl(".(fa|fasta)$", genome_fasta)) {
        in_fasta <- file(genome_fasta, open = "r")
    } else {
        stop(
            "It seems the genome sequence file is not a valid fasta file ",
            "which should with an extension .fa, .fasta, .fa.gz or .fasta.gz!"
        )
    }

    f <- ""
    line_num <- 1
    while (length({
        line <- readLines(in_fasta, n = 1, warn = FALSE)
    }) > 0) {
        line <- trimws(line) # remove leading and trailing white spaces
        if (length(line) == 0) {
            next
        }
        if (grepl("^>", line)) {
            if (line_num > 1) {
                close(f)
            }
            f <- gzfile(
                file.path(
                    out_dir,
                    gsub(
                        "^>(chr)?([^\\s]+).*",
                        paste0(prefix, "\\2.fa.gz"),
                        line,
                        perl = TRUE
                    )
                ),
                "w"
            )
            writeLines(gsub(
                "^>(chr)?([^\\s]+).*",
                paste0(">", prefix, "\\2"),
                line,
                perl = TRUE
            ), f)
        } else {
            writeLines(line, f)
        }
        line_num <- line_num + 1
    }
    close(f)
    close(in_fasta)
    out_dir
}

#' Prepare a seed file for building a BSgenome package
#'
#'
#' @param multifasta_path A character(1), specifying the path to a directory
#'   containing fasta.gz files for individual chromosome/scaffolds.
#' @param latin_name A character(1), specifying the Latin name of the species of
#'   the reference genome, such as "Homo sapiens" for the human.
#' @param common_name A character(1), specifying the common name of the species
#'   of the reference genome, such as "Human".
#' @param genome_version A character(1), specifying the genome build, such as
#'   "GRCh38".
#' @param seed_file_name A character(1), specifying the path to a seed file to
#'   be created for building a BSgenome package.
#' @param fasta_url A character(1), specifying the URL from where the reference
#'   genome fasta file is downloaded.
#' @param release_date  A character(1), specifying the genome assembly release
#'   date, such as "August 2020".
#' @param source A character(1), specifying the source of the genome fasta
#'   file, such as "Ensembl".
#' @param version A character(1), specifying the version of a BSgenome package
#'   to be build, such as "1.0.0".
#'
#' @return A character(2), the package name and the path to the seed file.
#' @export
#'
#' @examples
#' out_dir <- tempdir()
#' genome_fasta <- file.path(out_dir, "toy.example.fa")
#' in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'              ">chr2\nCCCCCCCCCCCGGAAA",
#'              ">chrM\nATACTACTGGA"), in_fasta)
#' close(in_fasta)
#'
#' fa_dir <- generate_multifasta(
#'     genome_fasta = genome_fasta,
#'     out_dir = out_dir,
#'     prefix = "none"
#' )
#'
#' seed_file <- generate_seed_file(
#'     multifasta_path = fa_dir,
#'     latin_name = "Homo sapiens",
#'     common_name = "Human",
#'     genome_version = "GRCh38",
#'     seed_file_name = file.path(
#'         out_dir,
#'         "human.genome.seed.txt"
#'     ),
#'     fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'     release_date = "August 2007",
#'     source = "Ensembl",
#'     version = "1.0.0"
#' )
generate_seed_file <-
    function(
        multifasta_path = NULL,
        latin_name = NULL,
        common_name = NULL,
        genome_version = NULL,
        seed_file_name = NULL,
        fasta_url = NULL,
        release_date = "August 2007",
        source = "Ensembl",
        version = "1.0.0") {
    if (is.null(multifasta_path) || is.null(latin_name) ||
        is.null(common_name) || is.null(genome_version) ||
        is.null(seed_file_name) || is.null(fasta_url)) {
        stop("All arguments except source and version are required!")
    }
    if (!dir.exists(multifasta_path)) {
        stop("Path to multifasta ", multifasta_path, " doesn't exist!")
    }
    chr_fa_files <- dir(multifasta_path, ".fa.gz$")
    if (length(chr_fa_files) < 1) {
        stop(
            "There are not multiple fasta files in the directory ",
            multifasta_path,
            "!"
        )
    }
    seed_dir <- dirname(seed_file_name)
    if (!dir.exists(seed_dir)) {
        dir.create(seed_dir, recursive = TRUE)
    }
    seed_fh <- file(seed_file_name, open = "w")
    BSgenomeObjname <- gsub("^(.).*\\s+(.+)", "\\1\\2", latin_name)
    package_name <- paste("BSgenome", BSgenomeObjname, source,
        genome_version,
        sep = "."
    )
    writeLines(paste0("Package: ", package_name), con = seed_fh)
    writeLines(paste(
        "Title: Full genome sequences for",
        latin_name,
        paste0("(", source),
        " version ",
        paste0(genome_version, ")")
    ), con = seed_fh)
    writeLines(
        paste(
            "Description: Full genome sequences for",
            latin_name,
            paste0("(", common_name, ")"),
            "as provided by",
            source,
            paste0("(", genome_version, ")"),
            "and stored in Biostrings objects."
        ), con = seed_fh
    )
    writeLines(paste0("Version: ", version), con = seed_fh)
    writeLines(paste0("organism: ", latin_name), con = seed_fh)
    writeLines(paste0("common_name: ", common_name), con = seed_fh)
    writeLines(paste0("provider: ", source), con = seed_fh)
    writeLines(paste0("release_date: ", release_date), con = seed_fh)
    writeLines(paste0("genome: ", genome_version), con = seed_fh)
    writeLines(paste0("source_url: ", fasta_url), con = seed_fh)
    writeLines(paste0("BSgenomeObjname: ", BSgenomeObjname), con = seed_fh)
    writeLines(paste0(
        "organism_biocview: ",
        gsub("\\s+", "_", latin_name, perl = TRUE)), con = seed_fh)
    chromosome_names <-
        gsub(".fa.gz$", "", dir(multifasta_path, "fa.gz$"))
    writeLines(paste0(
        'seqnames: c("',
        paste(chromosome_names, collapse = '","'),
        '")'
    ), con = seed_fh)

    circ_seqs <- c("chrM", "MT", "Pltd", "chrPltd")
    circ_seqs <- circ_seqs[circ_seqs %in% chromosome_names]
    if (length(circ_seqs) >=1)
    {
        writeLines(paste0('circ_seqs: c("',
                            paste(circ_seqs, collapse = '","'), '")'),
                    con = seed_fh)
    }

    writeLines(paste0("seqs_srcdir: ", multifasta_path),
                con = seed_fh)
    writeLines(paste0("seqfiles_suffix: .fa.gz"), con = seed_fh)
    close(seed_fh)
    c(package_name, seed_file_name)
}

#' Create and install a BSgenome package
#'
#' Based on a seed file, a BSgenome package is created and installed
#'
#' @param seed_file A character(1), specifying a path to a seed file
#' @param dest_dir A character(1), specifying a directory where a tar-ball
#'   of a BSgenome package is created.
#'
#' @return A BSgenome package name
#' @importFrom BSgenomeForge forgeBSgenomeDataPkg
#' @export
#'
#' @examples
#' if (TRUE) {
#'     out_dir <- tempdir()
#'     genome_fasta <- file.path(out_dir, "toy.example.fa")
#'     in_fasta <- file(genome_fasta, open = "w")
#' writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#'     close(in_fasta)
#'
#'     fa_dir <- generate_multifasta(
#'         genome_fasta = genome_fasta,
#'         out_dir = out_dir,
#'         prefix = "none"
#'     )
#'
#'     seed_file <- generate_seed_file(
#'         multifasta_path = fa_dir,
#'         latin_name = "Homo sapiens",
#'         common_name = "Human",
#'         genome_version = "GRCh38",
#'         seed_file_name = file.path(
#'             out_dir,
#'             "human.genome.seed.txt"
#'         ),
#'         fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'         release_date = "August 2007",
#'         source = "Ensembl",
#'         version = "1.0.0"
#'     )
#'
#'     forge_BSgenome(seed_file[2], dest_dir = out_dir)
#' }
forge_BSgenome <-
    function(
        seed_file = NULL,
        dest_dir = tempdir()) {
    if (!file.exists(seed_file)) {
        stop("seed_file are required!")
    }
    if (!dir.exists(dest_dir)) {
        dir.create(dest_dir, recursive = TRUE)
    }

    pkgname <- gsub(
        "Package:\\s*([^\\s]+)",
        "\\1",
        trimws(readLines(seed_file, n = 1))
    )
    if (dir.exists(file.path(dest_dir, pkgname))) {
        unlink(file.path(dest_dir, pkgname), recursive = TRUE)
    }

    ## crate a BSgenome package from a seed file
    BSgenomeForge::forgeBSgenomeDataPkg(seed_file, destdir = dest_dir)

    ## check, build and install package
    ## OR install using command line
    # R CMD build BSgenome.Hsapiens.Ensembl.GRCh38
    # R CMD check BSgenome.Hsapiens.Ensembl.GRCh38_1.0.0.tar.gz
    # devtools::check(file.path(dest_dir, pkgname))
    # devtools::build(file.path(dest_dir, pkgname), path = dest_dir)
    # devtools::install(file.path(dest_dir, pkgname))
    pkgname
}



#' Create and install a BSgenome package
#'
#' Starting with a multi-fasta file for a reference genome, a BSgenome package
#' is created and installed. This function is a wrapper function of three
#' functions: [generate_multifasta()], [generate_seed_file()], and
#' [forge_BSgenome()].
#'
#' @param genome_fasta A character(1) specifying a path for a multi-fasta file
#'   for a reference genome. It can be a gzip-compressed or uncompressed fasta
#'   file.
#' @param out_dir A character(1) specifying a path for output split, compressed
#'   fasta files
#' @param prefix A character(1), specifying chromosome prefix. It can be an
#'   empty string or "chr". Default is "".
#' @param latin_name A character(1), specifying the Latin name of the species of
#'   the reference genome, such as "Homo sapiens" for the human.
#' @param common_name A character(1), specifying the common name of the species
#'   of the reference genome, such as "Human".
#' @param genome_version A character(1), specifying the genome build, such as
#'   "GRCh38".
#' @param seed_file_name A character(1), specifying the path to a seed file to
#'   be created for building a BSgenome package.
#' @param fasta_url A character(1), specifying the URL from where the reference
#'   genome fasta file is downloaded.
#' @param release_date  A character(1), specifying the genome assembly release
#'   date, such as "August 2020".
#' @param source A character(1), specifying the source of the genome fasta
#'   file, such as "Ensembl".
#' @param version A character(1), specifying the version of a BSgenome package
#'   to be build, such as "1.0.0".
#'
#' @return A BSgenome package name. Users have to install the package
#'   to use it.
#' @export
#'
#' @examples
#' if (TRUE) {
#'     out_dir <- tempdir()
#'     genome_fasta <- file.path(out_dir, "toy.example.fa")
#'     in_fasta <- file(genome_fasta, open = "w")
#'     writeLines(c(">chr1\nATCGCTGCGGATGCGG",
#'               ">chr2\nCCCCCCCCCCCGGAAA",
#'               ">chrM\nATACTACTGGA"), in_fasta)
#'     close(in_fasta)
#'
#'     make_BSgenome(
#'         genome_fasta = genome_fasta,
#'         out_dir = file.path(out_dir, "multifasta"),
#'         prefix = "none",
#'         latin_name = "Homo sapiens",
#'         common_name = "Human",
#'         genome_version = "GRCh38",
#'         seed_file_name = file.path(
#'             out_dir,
#'             "human.genome.seed.txt"
#'         ),
#'         fasta_url = "http://ftp.ensembl.org/xxx.fa",
#'         release_date = "August 2007",
#'         source = "Ensembl",
#'         version = "1.0.0"
#'     )
#' }
make_BSgenome <-
    function(
        genome_fasta = NULL,
        out_dir = tempdir(),
        prefix = c("", "chr"),
        latin_name = NULL,
        common_name = NULL,
        genome_version = NULL,
        seed_file_name = NULL,
        fasta_url = NULL,
        release_date = "August 2007",
        source = "Ensembl",
        version = "1.0.0") {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    multifasta_path <-
        generate_multifasta(
            genome_fasta = genome_fasta,
            out_dir = out_dir,
            prefix = prefix
        )

    seed_file <- generate_seed_file(
        multifasta_path = multifasta_path,
        latin_name = latin_name,
        common_name = common_name,
        genome_version = genome_version,
        seed_file_name = seed_file_name,
        fasta_url = fasta_url,
        release_date = release_date,
        source = source,
        version = version
    )

    pkgname <- forge_BSgenome(
        seed_file = seed_file[2],
        dest_dir = out_dir
    )
    pkgname
}
