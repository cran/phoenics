# Query KEGGREST package

use_KEGGREST <- function(expr) {
  kegg_out <- try({ expr }, silent = TRUE)
  if (inherits(kegg_out, "try-error")) {
    message("'KEGGREST' query is not available at the moment. Are you connected",
            " to internet and if so, please, try again later.")
    return(NULL)
  }
  return(kegg_out)
}
