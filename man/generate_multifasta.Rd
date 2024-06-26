% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_BSgenome.R
\name{generate_multifasta}
\alias{generate_multifasta}
\title{Split a mutli-fasta file by chromosome/scaffold into individual fasta file}
\usage{
generate_multifasta(
  genome_fasta = NULL,
  out_dir = tempdir(),
  prefix = c("none", "chr")
)
}
\arguments{
\item{genome_fasta}{A character(1) specifying a path for a multi-fasta file
for a reference genome. It can be a gzip-compressed or uncompressed fasta
file.}

\item{out_dir}{A character(1) specifying a path for output split, compressed
fasta files}

\item{prefix}{A character(1), specifying chromosome prefix. It can be "none"
or "chr". Default is "none".}
}
\value{
A character string, which is th path to the directory where split,
compressed fasta files are located.
}
\description{
Given a multi-fasta file for a reference genome, create a gzip-compressed
fasta file for each chromosome/scaffold, such as chr1.fa.gz, chr2.fa.gz,...
A prefix "chr" may be added if necessary.
}
\examples{
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

}
