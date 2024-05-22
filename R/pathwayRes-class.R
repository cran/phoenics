#' Class pathwayRes
#' 
#' S3 class for pathway differential analysis results
#' 
#' Methods for the class pathwayRes
#' 
#' @name pathwayRes
#' @param object,x object of class \code{pathwayRes}
#' @param ... not used
NULL

#' @rdname pathwayRes
#' @aliases summary.pathwayRes
#' @examples 
#' data("MTBLS422")
#' quantif <- from_ASICS_to_PHOENICS(quantif)
#' out_test <- test_pathway(quantif, design, pathways, 
#'                          fixed = c("Age", "Treatment"), random = "Mouse", 
#'                          npc = 2, model = "blmer")
#' summary(out_test)
#' @exportS3Method
summary.pathwayRes <- function(object, ...) {
  cat("Number of pathways: ", length(object), "\n")
  cat("Tested effects: ",  object[[1]]$test_pathway$Fixed_effect)
}

#' @rdname pathwayRes
#' @aliases print.pathwayRes
#' @examples 
#' print(out_test)
#' @exportS3Method 
print.pathwayRes <- function(x, ...) {
  cat("A list of pathways which contains for each pathway: \n")
  path_res <- x[[1]]
  cat("$", names(path_res)[1], ": ", "Pathway name \n", sep = "")
  cat("$", names(path_res)[2], ": ", "Pathway code \n", sep = "")
  cat("$", names(path_res)[3], ": ", "Metabolites in the pathway \n", sep = "")
  cat("$", names(path_res)[4], ": ", "PCA results \n", sep = "")
  cat("$", names(path_res)[5], ": ", "Model results \n", sep = "")
  cat("$", names(path_res)[6], ": ", "Pathway test results \n", sep = "")
}

#' @rdname pathwayRes
#' @aliases plot.pathwayRes
#' @param pathway_id a character string or vector of pathway codes or names.
#' Default to \code{NULL} in which case the first pathway is displayed
#' @param plot a character string indicating the type of plot to return. Default
#' to \code{"eig"} (the screegraph of the PCA is displayed)
#' @param habillage a character string indicating the column of the design used 
#' to color the individuals. Only used when \code{plot = "ind"}. Default to 
#' \code{"none"} (no color)
#' @importFrom factoextra fviz_pca_ind fviz_eig fviz_pca_var
#' @examples 
#' plot(out_test)
#' plot(out_test, "mmu00052", plot = "var")
#' plot(out_test, "mmu00052", plot = "ind", habillage = "Age")
#' @exportS3Method
plot.pathwayRes <- function(x, ..., pathway_id = NULL, 
                            plot = c("eig", "var", "ind"), habillage = "none") {
  
  plot <- match.arg(plot)
  
  if (is.null(pathway_id)) {
    message("`pathway_id` not specified... defaulting to the first pathway")
    pathway_id <- names(x)[1]
  }
  
  object_sub <- extract.pathwayRes(x, pathway_id)
  lapply(object_sub, function (path_res) {
    if (plot == "var") {
      p <- factoextra::fviz_pca_var(path_res$PCA, title = path_res$pathway_name,
                                    repel = TRUE)
    }
    if (plot == "ind") {
      p <- factoextra::fviz_pca_ind(path_res$PCA, habillage = habillage, 
                                    title = path_res$pathway_name, 
                                    geom = "point", pointsize = 2, 
                                    mean.point = FALSE)
    }
    if (plot == "eig") {
      p <- factoextra::fviz_eig(path_res$PCA, addlabels = TRUE,
                                title = path_res$pathway_name)
    }
    return(p)
  })
}

#' @rdname pathwayRes
#' @aliases head
#' @aliases head.pathwayRes
#' @examples 
#' head(out_test)
#' @export
head <- function(object) {
  UseMethod("head")
}

#' @exportS3Method
head.pathwayRes <- function(object){
  cat("Pathways: \n")
  names(object)
}

#' @rdname pathwayRes
#' @aliases extract
#' @aliases extract.pathwayRes
#' @param pathway_id a character string or vector of pathway codes or names
#' @return The function \code{extract} returns an object of class
#' \code{pathwayRes} which is a list of pathway results, containing only the
#' pathways in \code{pathway_id}.
#' @examples 
#' extract(out_test, "mmu00562")
#' @export
extract <- function(object, pathway_id) {
  UseMethod("extract")
}

#' @exportS3Method
extract.pathwayRes <- function(object, pathway_id) {
  cond <- lapply(object, function(path) path$pathway_code %in% pathway_id |
                   path$pathway_name %in% pathway_id)
  object_sub <- object[unlist(cond)]
  class(object_sub) <- "pathwayRes"
  return(object_sub)
}

#' @rdname pathwayRes
#' @aliases adjust_pval
#' @aliases adjust_pval.pathwayRes
#' @param method a character string indicating the correction method to be used
#' for multiple testing correction (authorized values are those of 
#' \code{\link{p.adjust.methods}})
#' @param n number of comparisons for multiple testing correction
#' @return The function \code{adjust_pval} returns a data.frame with pathways
#' in rows and the following information in columns:
#' \item{pathway_name}{name of the pathway}
#' \item{pathway_code}{pathway code (identifier)}
#' \item{Fixed_effect}{tested effect}
#' \item{pval}{raw p-value of the pathway}
#' \item{adjusted_pval}{adjusted p-value of the pathway}
#' @importFrom stats p.adjust p.adjust.methods
#' @examples 
#' adj_pval <- adjust_pval(out_test)
#' @export
adjust_pval <- function(object, method = p.adjust.methods, n = length(object)) {
  UseMethod("adjust_pval")
}

#' @exportS3Method 
adjust_pval.pathwayRes <- function(object, method = p.adjust.methods,
                                   n = length(object)) {
  pval_df <- lapply(object, function(path) {
    pathway_name <- path$pathway_name
    pathway_code <- path$pathway_code
    pathway_pval <- path$test_pathway
    pathway_pval <- cbind(pathway_name, pathway_code, pathway_pval)
  })
  pval_df <- do.call("rbind", pval_df)
  rownames(pval_df) <- NULL
  
  res <- lapply(unique(pval_df$Fixed_effect), function(effect) {
    sub <- pval_df[pval_df$Fixed_effect == effect, ]
    sub$adjusted_pval <- p.adjust(sub$pval, method, n)
    return(sub)
  })
  res <- do.call("rbind", res)
  return(res)
}
