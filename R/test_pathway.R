#' @title Pathway differential analysis based on longitudinal metabolomics data
#' 
#' @description Perform a differential analysis at pathway level based on 
#' metabolite quantifications and information on pathway metabolite composition.
#' The method relies on a PCA step.
#' 
#' @param quantif data.frame or matrix of the metabolite quantification with 
#' samples in rows (sample identifiers must be row names) and metabolites in 
#' columns (metabolite codes must be column names)
#' @param design data.frame or matrix with samples in rows (sample identifiers
#' must be row names) and the different effects to be included in the model
#' in columns. Column names must be provided and are used as arguments for
#' \code{fixed} and \code{random}
#' @param pathways data.frame or matrix with metabolites in rows (all 
#' metabolites in columns of \code{quantif} must have a row in this input) and the
#' following information in columns: \itemize{
#'  \item \code{metabolite_code} metabolite code
#'  \item \code{metabolite_name} metabolite name
#'  \item \code{pathway_code} pathway code (identifier)
#'  \item \code{pathway_name} name of the pathway
#' }
#' @param npc number of PCs for the PCA step
#' @param model a character string indicating if the model has to be fitted
#' using \link[lme4]{lmer} or \link[blme]{blmer}. Default to \code{"lmer"}
#' @param fixed character vector of fixed effects to be included in the model.
#' They must correspond to column names of \code{design}
#' @param random character vector of random effects to be included in the model.
#' They must correspond to column names of \code{design}
#' @param organism organism code in KEGG database. Required if
#' \code{pathways = "auto"} and ignored otherwise
#' @param min_size minimal number of metabolites in a pathway. Required if
#' \code{pathways = "auto"} and ignored otherwise. Default to \code{2}
#' 
#' @return an object of class \code{pathwayRes} which is a list of pathway 
#' results. Each element of the list contain the following entries:
#' \item{pathway_name}{a character corresponding to the pathway name}
#' \item{pathway_code}{a character corresponding to the pathway code}
#' \item{metabolites}{a data.frame with the names and codes of the quantified
#' metabolites in the pathway}
#' \item{PCA}{the result of the pathway PCA (a \code{PCA} object as obtained 
#' from \link[FactoMineR]{PCA})}
#' \item{model}{the output of the mixed model fit}
#' \item{test_pathway}{a data.frame with the p-values for each tested fixed
#' effect}
#' 
#' @importFrom FactoMineR PCA
#' @importFrom lme4 lmer
#' @importFrom blme blmer
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom stats as.formula anova
#' @importFrom tidyr separate
#' 
#' @examples 
#' data("MTBLS422")
#' quantif <- from_ASICS_to_PHOENICS(quantif)
#' out_test <- test_pathway(quantif, design, pathways, 
#'                          fixed = c("Age", "Treatment"), random = "Mouse", 
#'                          npc = 2, model = "blmer")
#' out_test
#' 
#' if (requireNamespace("KEGGREST", quietly = TRUE)) { 
#' \donttest{
#'   out_test2 <- test_pathway(quantif, design, pathways = "auto", 
#'                             fixed = c("Age", "Treatment"), random = "Mouse", 
#'                             npc = 2, model = "blmer", organism = "mmu")
#'   out_test2
#'  }
#' }
#' 
#' @seealso
#' \link{pathway_search} for the automatic annotation of metabolites in KEGG 
#' pathways.
#' 
#' @details
#' If \code{pathways = "auto"}, information on pathways in which metabolites are
#' present is automatically obtained by the function \link{pathway_search} using
#' \code{KEGGREST} that queries KEGG database. Results are specific to a given
#' organism (passed in \code{organism}). Pathways containing less than
#' \code{min_size} metabolites are filtered out.
#' 
#' @author Camille Guilmineau <camille.guilmineau@inrae.fr>\cr
#' Remi Servien <remi.servien@inrae.fr>\cr
#' Nathalie Vialaneix <nathalie.vialaneix@inrae.fr>
#' 
#' @export

test_pathway <- function(quantif, design, pathways = "auto", fixed, random, 
                         npc = 1L, model = c("lmer", "blmer"),
                         organism = NULL, min_size = 2) {
  
  # Inputs tests
  if (!is.character(fixed)) {
    stop("`fixed` must be a character vector.")
  }
  
  if (!is.character(random) ) {
    stop("`random` must be a character vector.")
  }
  
  if (!all(fixed %in% colnames(design))) {
    stop("Fixed effects must be in `design` column names.")
  }
  
  if (!all(random %in% colnames(design))) {
    stop("Random effect must be in `design` column names.")
  }
  
  
  if (!identical(sort(rownames(quantif)), sort(rownames(design)))) {
    stop("`quantif` and `design` row names must be identical.")
  } 
  
  if (!identical(rownames(quantif), rownames(design))) {
    order_names <- match(rownames(quantif), rownames(design))
    quantif <- quantif[order_names, ]
    message("`quantif` and `design` did not have the same order for row names:",
            " changing row order for `quantif` to match `design`.")
  }
  
  model <- match.arg(model)
  
  if (round(npc) != npc) {
    stop("`npc` must be an integer.")
  }
  
  if (is.character(pathways)) {
    if (pathways != "auto") {
      stop("`pathways` must be \"auto\" or a data.frame or a matrix.")
    }
    if (is.null(organism)) {
      stop("`organism` is required when `pathways = \"auto\"`")
    }
    if (round(min_size) != min_size) {
      stop("`min_size` must be an integer.")
    }
    
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
      stop("Package 'KEGGREST' not available. Automatic pathway searching", 
           " cannot be performed.")
    }
    kegg_org <- KEGGREST::keggList("organism")
    kegg_org <- as.data.frame(kegg_org)
    kegg_org <- unique(kegg_org$organism)
    
    if (!(organism %in% kegg_org)) {
      stop("`organism` must be a KEGG organism code.")
    }
    
    metab <- as.vector(colnames(quantif))
    pathways <- pathway_search(metab, organism, min_size)
    
  } else {
    if (!is.matrix(pathways) & !is.data.frame(pathways)) {
      stop("`pathways` must be \"auto\" or a data.frame or a matrix.")
    }
    if (!all(c("metabolite_code", "metabolite_name", "pathway_code",
               "pathway_name") %in% colnames(pathways))) {
      stop("The column names of `pathways` must be \"metabolite_code\", ",
           "\"metabolite_name\", \"pathway_code\", \"pathway_name\".")
    }
    if (!all(pathways$metabolite_code %in% colnames(quantif))) {
      stop("Metabolites codes from `pathways` must be in the column names of", 
           " `quantif`.")
    }
  }
  
  if (npc > min(table(pathways$pathway_code))) {
    stop(paste0("`npc` is larger than the number of metabolites in a pathway.",
                "\nChoose a smaller value."))
  }
  
  if (!all(colnames(quantif) %in% pathways$metabolite_code)) {
    quantif <- quantif[, pathways$metabolite_code]
    warning(paste("Metabolites in quantifications have been filtered based on",
                  "metabolites in pathways.\n",
                  "Quantification matrix now contains:", ncol(quantif),
                  "metabolites."))
  }
  
  # Main part: loop over pathways to perform test
  res_pathway <- lapply(unique(pathways$pathway_code), function(path) {
    
    # Calculate pathway score matrix
    pathways_sub <- pathways[pathways$pathway_code == path, ]
    metab_path <- pathways_sub$metabolite_code
    
    quantif_path <- quantif[, metab_path]
    quantif_path <- merge(quantif_path, design, by = 0)
    quantif_path <- column_to_rownames(quantif_path, "Row.names")
    
    res_pca <- PCA(quantif_path, graph = FALSE, ncp = npc,
                   quali.sup = colnames(design))
    
    score <- data.frame(res_pca$ind$coord[, 1:npc, drop = FALSE])
    colnames(score) <- paste0(path, "_PC", 1:ncol(score))
    
    # Mixed model
    score_path <- merge(score, design, by = 0)
    score_path <- column_to_rownames(score_path, "Row.names")
    
    form <- paste0("~ ", paste(fixed, collapse = " + "), " + ", 
                   paste(paste0("(1|", random, ")"), collapse = " + "))
    fmla <- sapply(paste0(colnames(score), form), as.formula)
    
    model_function <- eval(as.name(model))
    res_lmm <- lapply(fmla, function(path_form) {
      mod <- model_function(formula = path_form, data = score_path)
    })
    names(res_lmm) <- colnames(score)
    
    # Test fixed effects
    null_formulas <- lapply(fixed, function(eff) {
      in_formula <- fixed[-which(fixed %in% eff)]
      form <- paste0("~ ", paste(in_formula, collapse = " + "), " + ", 
                     paste(paste0("(1|", random, ")"), collapse = " + "))
    })
    names(null_formulas) <- fixed
    
    res_test <- lapply(null_formulas, function(null_form) {
      fmla_null <- sapply(paste0(colnames(score), null_form), as.formula)
      names(fmla_null) <- colnames(score)
      
      res_lmm_null <- lapply(fmla_null, function(path_form) {
        path_name <- as.character(path_form[2])
        mod <- model_function(formula = path_form, data = score_path)
      })
      
      res_aov <- mapply(function(model1, model2) anova(model1, model2,
                                                       test = "LRT"), 
                        res_lmm, res_lmm_null, SIMPLIFY = FALSE)
      names(res_aov) <- colnames(score)
      
      res_aov_df <- do.call("rbind", res_aov)
      res_aov_df <- res_aov_df[!is.na(res_aov_df$`Pr(>Chisq)`), ]
      rownames(res_aov_df) <- names(res_aov)
      
      pval <- res_aov_df$`Pr(>Chisq)`
      
      signif_pathway <- data.frame(rownames(res_aov_df), pval)
      colnames(signif_pathway) <- c("pathway_PC", "pval")
      signif_pathway <- separate(signif_pathway, col = "pathway_PC",
                                        into = c("pathway", "PC"), "_")
      signif_pathway$pathway_PC <- paste0(signif_pathway$pathway, "_",
                                          signif_pathway$PC)
      
      # Simes procedure
      simes <- sapply(unique(signif_pathway$pathway), function(path) {
        p_vals <- signif_pathway[signif_pathway$pathway == path, "pval"]
        pvals_sort <- p_vals[order(p_vals)]
        out <- min(1, min(length(pvals_sort) * pvals_sort / rank(pvals_sort)))
        return(out)
      })
      simes <- as.data.frame(simes)
      colnames(simes)[1] <- "pval"
      return(simes)
    })
    
    res_test <- do.call("rbind", res_test)
    res_test <- rownames_to_column(res_test, "Fixed_effect")
    
    path_name <- unique(pathways_sub[, "pathway_name"])
    path_code <- unique(pathways_sub[, "pathway_code"])
    out_met <- unique(pathways_sub[, c("metabolite_code", "metabolite_name")])
    
    res <- list("pathway_name" = path_name, "pathway_code" = path_code,
                "metabolites" = out_met, "PCA" = res_pca, "model" = res_lmm,
                "test_pathway" = res_test)
    return(res)
  })
  
  names(res_pathway) <- unique(pathways$pathway_name)
  class(res_pathway) <- "pathwayRes"
  
  return(res_pathway)
}
