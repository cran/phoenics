## ----loadLib, message=FALSE---------------------------------------------------
library("phoenics")

## ----loadData-----------------------------------------------------------------
data("MTBLS422")
ls()
head(quantif)

## ----quantif------------------------------------------------------------------
quantif_f <- from_ASICS_to_PHOENICS(quantif)
head(quantif_f)

## ----pathways-----------------------------------------------------------------
head(pathways)

## ----searchPaths, eval=FALSE--------------------------------------------------
# pathways <- pathway_search(metab = colnames(quantif), organism = "mmu")

## ----design-------------------------------------------------------------------
head(design)

## ----test, message=FALSE------------------------------------------------------
out_test <- test_pathway(quantif_f, design, pathways, 
                         fixed = c("Age", "Treatment"), random = "Mouse", 
                         npc = 2, model = "blmer")
out_test

## ----pathwayNames-------------------------------------------------------------
names(out_test)

## ----extract------------------------------------------------------------------
res_1 <- extract(out_test, "Galactose metabolism")

## ----detailedRes--------------------------------------------------------------
res_1[[1]]$test_pathway
res_1[[1]]$model

## ----BH-----------------------------------------------------------------------
adjust_pval(out_test)

## -----------------------------------------------------------------------------
res_1[[1]]$PCA

## ----plotPCA------------------------------------------------------------------
plot(out_test, pathway_id = "Galactose metabolism", plot = "var")
plot(out_test, pathway_id = "Galactose metabolism", plot = "ind", 
     habillage = "Age")
plot(out_test, pathway_id = "Galactose metabolism", plot = "eig")
plot(out_test, pathway_id = "Galactose metabolism", plot = "var")
plot(out_test, pathway_id = "Galactose metabolism", plot = "group")

## ----testMFA------------------------------------------------------------------
out_test <- test_pathway(quantif_f, design, pathways, 
                         fixed = c("Age", "Treatment"), random = "Mouse", 
                         npc = 2, model = "blmer", analysis = "MFA")
out_test

## ----auto, eval=FALSE---------------------------------------------------------
# out_test3 <- test_pathway(quantif, design, pathways = "auto",
#                           fixed = c("Age", "Treatment"), random = "Mouse",
#                           npc = 2, model = "blmer", organism = "mmu")
# out_test3

## ----session------------------------------------------------------------------
sessionInfo()

