% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_BSgenome.R
\name{forge_BSgenome}
\alias{forge_BSgenome}
\title{Create and install a BSgenome package}
\usage{
forge_BSgenome(seed_file = NULL, dest_dir = tempdir())
}
\arguments{
\item{seed_file}{A character(1), specifying a path to a seed file}

\item{dest_dir}{A character(1), specifying a directory where a tar-ball
of a BSgenome package is created.}
}
\value{
A BSgenome package name
}
\description{
Based on a seed file, a BSgenome package is created and installed
}
\examples{
if (TRUE) {
    out_dir <- tempdir()
    genome_fasta <- file.path(out_dir, "toy.example.fa")
    in_fasta <- file(genome_fasta, open = "w")
writeLines(c(">chr1\nATCGCTGCGGATGCGG",
              ">chr2\nCCCCCCCCCCCGGAAA",
              ">chrM\nATACTACTGGA"), in_fasta)
    close(in_fasta)

    fa_dir <- generate_multifasta(
        genome_fasta = genome_fasta,
        out_dir = out_dir,
        prefix = "none"
    )

    seed_file <- generate_seed_file(
        multifasta_path = fa_dir,
        latin_name = "Homo sapiens",
        common_name = "Human",
        genome_version = "GRCh38",
        seed_file_name = file.path(
            out_dir,
            "human.genome.seed.txt"
        ),
        fasta_url = "http://ftp.ensembl.org/xxx.fa",
        release_date = "August 2007",
        source = "Ensembl",
        version = "1.0.0"
    )

    forge_BSgenome(seed_file[2], dest_dir = out_dir)
}
}
