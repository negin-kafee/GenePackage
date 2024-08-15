## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GenePackage)

## ----create-gene-example------------------------------------------------------
gene <- createGene("protein_coding", 
                   ID = "ENSG000001", 
                   symbol = "TP53", 
                   name = "Tumor protein p53", 
                   description = "Plays a role in cell cycle regulation", 
                   chromosome = "17", 
                   start = 7565097, 
                   end = 7590856, 
                   strand = "+", 
                   proteinID = "NP_000537", 
                   proteinSequence = "MDM2")

## ----get-symbol-example-------------------------------------------------------
getSymbol(gene)

## ----length-product-example---------------------------------------------------
lengthProduct(gene)

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

