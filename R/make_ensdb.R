#' Make an EnsDb from a GTF file
#'
#' A wrapper function from the ensembldb package to make an EnsDb SQLite
#' database from an Ensembl GTF file.
#'
#' @param gtf A character(1), a path to an Ensembl GTF file.
#' @param outfile A character(1), a path to a SQLite file to store an object of
#'   the [ensembldb::EnsDb-class].
#' @param organism_latin_name A character(1), Latin name for the organism, with
#'   spaces replaced by underscores "_".
#' @param genome_version A character(1), the assembly version of the reference
#'   genome, such as "GRCh38" for the human reference genome version 38.
#' @param Ensembl_release_version An integer(1), the Ensembl release version of
#'   the reference genome with patches, such as 110 for Ensembl Release 110
#'   (July 2023).
#' @importFrom ensembldb ensDbFromGtf
#'
#' @return A name of an SQLite file for the built EnsDb. Use the
#'   [ensembldb::EnsDb()] function from the package ensembldb package to
#'   load the SQLite file to the R environment.
#' @export
#'
#' @examples
#' \dontrun{
#' options(timeout = max(3000, getOption("timeout")))
#' gtf_url <- paste0("https://ftp.ensembl.org/pub/release-110/gtf/",
#'                   "homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz")
#' gtf <- basename(gtf_url)
#' tmp_dir <- tempdir()
#' download.file(url = gtf_url,
#'               destfile = file.path(tmp_dir, gtf),
#'               mode = "wb")
#'
#' hs_ensdb_sqlite <- make_ensdb(gtf = file.path(tmp_dir, gtf),
#'                               outfile = "EnsDb.hs.v110.sqlite",
#'                               organism_latin_name = "Homo_Sapiens",
#'                               genome_version = "GRCh38",
#'                               Ensemble_release_version = 110)
#' }


make_ensdb <- function(gtf = NULL,
                       outfile = NULL,
                       organism_latin_name = "Homo_sapiens",
                       genome_version = "GRCh38",
                       Ensembl_release_version = 110)
{
  if (length(gtf)!= 1 || !file.exists(gtf)) {
    stop("A GTF file is required!")
  }
  ensdb_sqlite_file <- ensDbFromGtf(gtf = gtf,
                                    outfile = outfile,
                                    organism = organism_latin_name,
                                    genomeVersion = genome_version,
                                    version = Ensembl_release_version)
  ensdb_sqlite_file
}
