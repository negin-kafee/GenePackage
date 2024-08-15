library(testthat)
library(GenePackage)

test_that("getSymbol works correctly", {
  gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53")
  expect_equal(getSymbol(gene), "TP53")
})

test_that("lengthProduct works for ProteinCodingGene", {
  gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53", proteinSequence = "MDM2")
  expect_equal(lengthProduct(gene), 4)
})
