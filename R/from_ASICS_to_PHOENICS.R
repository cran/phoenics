#' @title Prepare quantification data from ASICS outputs
#' 
#' @description Prepare quantification data from \code{ASICS} outputs for
#' \code{test_pathway}. In short, it replaces metabolite names by metabolites
#' KEGG codes and transposes the matrix to have samples in rows and metabolites
#' in columns.
#' 
#' @param quantif output matrix of \code{ASICS} quantification
#' 
#' @return A matrix of quantification with samples in rows and metabolites in
#' columns, properly formatted for \code{test_pathway}
#' 
#' @importFrom utils read.table
#'
#' @examples 
#' data("MTBLS422")
#' quantif <- from_ASICS_to_PHOENICS(quantif)
#' 
#' @author Camille Guilmineau <camille.guilmineau@inrae.fr>\cr
#' Remi Servien <remi.servien@inrae.fr>\cr
#' Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>
#' 
#' @references
#' Lefort G., Liaubet L., Canlet C., Tardivel P., P\`ere M.C., Quesnel H., 
#' Paris A., Iannuccelli N., Vialaneix N. Servien R. (2019). ASICS: an R
#' package for a whole analysis workflow of 1D 1H NMR spectra.
#' \emph{Bioinformatics}, \strong{35}(21): 4356--4363.
#' \doi{10.1093/bioinformatics/btz248}
#' 
#' Tardivel P., Canlet C., Lefort G., Tremblay-Franco M., Debrauwer L.,
#' Concordet D., Servien R. (2017). ASICS: an automatic method for
#' identification and quantification of metabolites in complex 1D 1H NMR
#' spectra. \emph{Metabolomics}, \strong{13}(10): 109.
#' \doi{10.1007/s11306-017-1244-5}
#'
#' @export

from_ASICS_to_PHOENICS <- function(quantif) {
  
  asics_lib <- read.table(system.file("adddata", "ASICS_library_codes.txt",
                                      package = "phoenics"), header = TRUE)
  
  if ("1.3-Diaminopropane" %in% rownames(quantif)) {
    rownames(quantif) <- gsub("1.3-Diaminopropane", "1,3-Diaminopropane",
                              rownames(quantif))
  }
  
  quantif <- t(quantif)
  
  colnames(quantif) <- sapply(colnames(quantif), function(col) {
    asics_lib[which(asics_lib$Query == col), "KEGG"]
  })
  
  n <- ncol(quantif[, is.na(colnames(quantif))])
  if (n != 0) {
    warning("No correspondance found between name and KEGG code for ", n, 
            " metabolites. They have been removed. ", ncol(quantif) - n, 
            " metabolites left.")
    quantif <- quantif[, !is.na(colnames(quantif))]
  }
  
  return(quantif)
}

