#' @title Query KEGG pathways for a given set of metabolites
#' 
#' @param metab vector of metabolite KEGG codes
#' @param organism organism code in KEGG database
#' @param min_size minimal number of metabolites required for a pathway to be
#' returned
#' 
#' @return a data.frame with metabolites in rows and the following information
#' in columns: \itemize{
#'  \item \code{metabolite_code} metabolite code
#'  \item \code{metabolite_name} metabolite name
#'  \item \code{pathway_code} pathway code (identifier)
#'  \item \code{pathway_name} name of the pathway
#' }
#' 
#' @importFrom tibble rownames_to_column remove_rownames
#' @importFrom tidyr separate
#' 
#' @examples
#' if (requireNamespace("KEGGREST", quietly = TRUE)) {
#'   data("MTBLS422")
#'   quantif <- from_ASICS_to_PHOENICS(quantif)
#'   \donttest{
#'   pathways <- pathway_search(metab = colnames(quantif), organism = "mmu")
#'   }
#' }
#' 
#' @author Camille Guilmineau <camille.guilmineau@inrae.fr>\cr
#' Remi Servien <remi.servien@inrae.fr>\cr
#' Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>
#' 
#' @references
#' Kanehisa M., Goto S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes,
#' \emph{Nucleic Acids Research}, Volume 28, Issue 1, Pages 27--30,
#' \doi{10.1093/nar/28.1.27}
#' 
#' Tenenbaum D., Maintainer B. (2022). KEGGREST: Client-side REST access to the
#' Kyoto Encyclopedia of Genes and Genomes (KEGG). R package version 1.38.0.
#' 
#' @export

pathway_search <- function(metab, organism, min_size = 2) {
  
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("Package 'KEGGREST' not available. This function cannot be used.")
  }
  if (round(min_size) != min_size) {
    stop("`min_size` must be an integer.")
  }
  if (!is.character(metab)) {
    stop("`metab` must be a character vector.")
  }
  if (!is.character(organism) & !length(organism) == 1) {
    stop("`organism` must be a character.")
  }
  kegg_org <- use_KEGGREST(KEGGREST::keggList("organism"))
  kegg_org <- as.data.frame(kegg_org)
  
  if (!(organism %in% kegg_org$organism)) {
    stop("`organism` must be a KEGG organism code.")
  }
  
  message("Starts searching for pathways...")
  
  # search KEGG pathways
  pathlist <- use_KEGGREST(KEGGREST::keggList("pathway", organism))
  pathlist <- as.data.frame(pathlist)
  colnames(pathlist) <- "pathway_name"
  pathlist <- rownames_to_column(pathlist, "pathway_code")
  
  org_name <- kegg_org[kegg_org$organism == organism, "species"]
  org_name <- paste0(" - ", org_name)
  pathlist$pathway_name <- gsub(org_name, "", pathlist$pathway_name,
                                fixed = TRUE)
  
  comp <- sapply(1:nrow(pathlist), function(i) {
    path <- use_KEGGREST(KEGGREST::keggGet(paste0("path:",
                                                  pathlist$pathway_code[i])))
    path <- as.data.frame(path[[1]]$COMPOUND)
  })
  names(comp) <- pathlist$pathway_code
  comp <- comp[sapply(comp, nrow) > 0]
  
  comp_df <- do.call("rbind", comp)
  colnames(comp_df) <- "metabolite_name"
  comp_df <- rownames_to_column(comp_df, "pathway.metabolite")
  comp_df <- separate(comp_df, col = "pathway.metabolite",
                      into = c("pathway_code", "metabolite_code"),
                      remove = FALSE, fill = "right")
  comp_df$pathway.metabolite <- NULL
  
  pathlist <- merge(pathlist, comp_df, by = "pathway_code")
  
  # subset the pathways with the requested metabolites
  pathlist_sub <- pathlist[pathlist$metabolite_code %in% metab, ]
  tab <- as.data.frame(table(pathlist_sub$pathway_code))
  tab <- tab[tab$Freq >= min_size, ]
  pathlist_sub <- pathlist_sub[pathlist_sub$pathway_code %in% tab$Var1, ]
  
  pathlist_sub <- pathlist_sub[, c("metabolite_code", "metabolite_name", 
                                   "pathway_code", "pathway_name")]
  
  message("Done.")
  
  return(pathlist_sub)
}
