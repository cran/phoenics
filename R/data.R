## Documentation for datasets

#' @title Dataset "MTBLS422"
#' @name MTBLS422
#' @aliases quantif
#' @aliases pathways
#' @aliases design
#' 
#' @description  Metabolites quantifications, associated metabolic pathways and
#' experimental design to characterize the effects of two clinically important
#' antibiotic treatments, ciprofloxacin and vancomycin-imipenem on mice.
#' 
#' @docType data
#' 
#' @format 3 datasets are provided: \itemize{
#' \item \code{quantif}: a data frame with 10 rows (metabolites name) and 11
#' columns (samples id), which corresponds to the metabolites quantifications
#' in the samples.
#' \item \code{pathways}: a data frame with 11 rows and 4 columns, which
#' contains informations about pathways and their metabolites.
#' \item \code{design}: a data frame with 11 rows (samples id) and 4 columns
#' (id and effects to be added in the model).
#' }
#' 
#' @details The raw dataset have been made available on MetaboLights (with the
#' id MTBLS422 \url{https://www.ebi.ac.uk/metabolights/editor/MTBLS422}) by
#' [Choo \emph{et al.}, 2017]. Metabolite quantifications were obtained
#' based on the raw signal using \code{ASICS} package. Pathways were obtained
#' using \code{KEGGREST} package. The datasets provided for the example
#' are a subset of the original dataset.
#'  
#' @references Choo J. M., Kanno T., Zain N. M. M., Leong L. E. X.,
#' Abell G. C. J., Keeble J. E., Bruce K. D., Mason A. J., Rogers G. B. (2017).
#' Divergent relationships between fecal microbiota and metabolome following
#' distinct antibiotic-induced disruptions. \emph{mSphere}, \strong{2}(1)
#' \doi{10.1128/msphere.00005-17}
#'  
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
#' Tenenbaum D., Maintainer B. (2022). KEGGREST: Client-side REST access to the
#' Kyoto Encyclopedia of Genes and Genomes (KEGG). R package version 1.38.0.
#' 
#' @examples
#' data(MTBLS422)
#' \donttest{
#' design[1:5, ]
#' pathways[1:5, ]
#' quantif[1:5, 1:5]
#' }
#'
NULL


