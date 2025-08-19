#' @title Calculate overlap coefficient between pathways
#'
#' @param pathwayA a character string of pathway name or pathway code or a
#' character vector of metabolite names or metabolite codes
#' @param pathwayB a character string of pathway name or pathway code or a
#' character vector of metabolite names or metabolite codes
#' @param pathways data.frame or matrix with metabolites in rows and the
#' following information in columns: \itemize{
#'  \item \code{metabolite_code} metabolite code
#'  \item \code{metabolite_name} metabolite name
#'  \item \code{pathway_code} pathway code (identifier)
#'  \item \code{pathway_name} name of the pathway
#' }
#' Used if \code{pathwayA} and \code{pathwayB} are pathway names or pathway
#' codes.
#' @param organism organism code in KEGG database. Required if
#' \code{pathways = NULL} and \code{pathwayA} and \code{pathwayB} are pathway 
#' names and ignored otherwise.
#' 
#' @return A value between 0 and 1 calculated with the formula:
#' \deqn{\text{OC}(A,B) = \frac{\vert A \cap B \vert}{\min(\vert A \vert,
#' \vert B \vert)}}
#' An overlap coefficient of 1 means that one pathway is included in the other.
#' An overlap coefficient of 0 means that there is no overlap between the
#' pathways.
#' 
#' @examples 
#' data("MTBLS422")
#' pathwayA <- "Galactose metabolism"
#' pathwayB <- "Vitamin digestion and absorption"
#' overlap_coefficient(pathwayA, pathwayB, pathways)
#' 
#' if (requireNamespace("KEGGREST", quietly = TRUE)) {
#'   pathwayA <- "Galactose metabolism"
#'   pathwayB <- "Vitamin digestion and absorption"
#'   \donttest{
#'   overlap_coefficient(pathwayA, pathwayB, organism = "mmu")
#'   }
#' }
#' 
#' pathwayA <- "mmu00052"
#' pathwayB <- "mmu00562"
#' overlap_coefficient(pathwayA, pathwayB, pathways)
#' 
#' pathwayA <- c("C00029", "C00116", "C00137", "C00794", "C00984", "C01697")
#' pathwayB <- c("C00191", "C00092", "C00137")
#' overlap_coefficient(pathwayA, pathwayB)
#' 
#' @references
#' Wieder C., Lai R.P.J., Ebbels T.M.D. (2022). Single sample pathway analysis
#' in metabolomics: performance evaluation and application.
#' \emph{BMC Bioinformatics}, \strong{23}, 481.
#' \doi{10.1186/s12859-022-05005-1}
#' 
#' @export

overlap_coefficient <- function(pathwayA, pathwayB, pathways = NULL, organism = NULL) {
  
  if (!is.character(pathwayA) & !is.character(pathwayB)) {
    stop("`pathwayA` and `pathwayB` must be character or vector of characters.")
  }
  
  if (length(pathwayA) == 1 & length(pathwayB) == 1) {
    if (!is.null(pathways)) {
      if (!is.matrix(pathways) & !is.data.frame(pathways)) {
        stop("`pathways` must be NULL or a data.frame or a matrix.")
      }
      
      if (!all(c("metabolite_code", "metabolite_name", "pathway_code",
                 "pathway_name") %in% colnames(pathways))) {
        stop("The column names of `pathways` must be \"metabolite_code\", ",
             "\"metabolite_name\", \"pathway_code\", \"pathway_name\".")
      }
      
      if (!(pathwayA %in% pathways$pathway_code |
            pathwayA %in% pathways$pathway_name) |
          !(pathwayB %in% pathways$pathway_code |
            pathwayB %in% pathways$pathway_name)) {
        stop("`pathwayA` and `pathwayB` must be pathway names or pathway codes",
             " from `pathways`.")
      }
      
      A <- pathways[pathways$pathway_code == pathwayA | pathways$pathway_name ==
                      pathwayA, "metabolite_code"]
      B <- pathways[pathways$pathway_code == pathwayB | pathways$pathway_name ==
                      pathwayB, "metabolite_code"]
    }
    
    if (is.null(pathways)) {
      message("`pathways` is NULL. Pathway composition is searched in KEGG",
              " database.")
      
      if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("Package 'KEGGREST' not available. Automatic pathway searching.", 
             "cannot be performed.")
      }
      
      code <- substr(pathwayA, start = 1, stop = nchar(pathwayA) - 5)
      kegg_org <- use_KEGGREST(KEGGREST::keggList("organism"))
      
      if (!code %in% kegg_org) {
        pathlist <- use_KEGGREST(KEGGREST::keggList("pathway", organism))
        pathlist <- as.data.frame(pathlist)
        pathlist$pathway_code <- rownames(pathlist)
        
        kegg_org <- as.data.frame(kegg_org)
        org_name <- kegg_org[kegg_org$organism == organism, "species"]
        org_name <- paste0(" - ", org_name)
        pathlist$pathlist <- gsub(org_name, "", pathlist$pathlist,
                                      fixed = TRUE)
        pathwayA <- pathlist[pathlist$pathlist == pathwayA, "pathway_code"]
        pathwayB <- pathlist[pathlist$pathlist == pathwayB, "pathway_code"]
      }
      
      pathinfo <- use_KEGGREST(KEGGREST::keggGet(c(pathwayA, pathwayB)))
      A <- pathinfo[[1]]$COMPOUND
      B <- pathinfo[[2]]$COMPOUND
    }
  }
  
  if (length(pathwayA) > 1 & length(pathwayB) > 1) {
    A <- pathwayA
    B <- pathwayB
  }
  
  intersection <- length(intersect(A, B))
  smaller_set <- min(length(A), length(B))
  OC <- intersection / smaller_set
  
  return(OC)
}
