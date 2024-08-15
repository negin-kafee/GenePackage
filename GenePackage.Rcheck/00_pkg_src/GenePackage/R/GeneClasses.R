# Gene Class Definition
# This class defines a general gene structure.

#' Gene Class
#'
#' A class representing a general gene structure with various attributes.
#' @name Gene-class
#' @slot ID A character string representing the gene ID.
#' @slot symbol A character string representing the gene symbol.
#' @slot name A character string representing the gene name.
#' @slot description A character string representing the gene description.
#' @slot chromosome A character string representing the chromosome.
#' @slot start A numeric value representing the start position of the gene.
#' @slot end A numeric value representing the end position of the gene.
#' @slot strand A character string representing the DNA strand.
#' @export
setClass("Gene",
  slots = list(
    ID = "character",
    symbol = "character",
    name = "character",
    description = "character",
    chromosome = "character",
    start = "numeric",
    end = "numeric",
    strand = "character"
  )
)

#' Get the Symbol of a Gene Object
#'
#' This method returns the symbol of a gene object.
#' @param object An object of class Gene.
#' @return The symbol of the gene.
#' @aliases getSymbol getSymbol,Gene-method
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53")
#' getSymbol(gene)
#' @export
#' @usage \S4method{getSymbol}{Gene}(object)
#' @rdname getSymbol-method
setGeneric("getSymbol", function(object) standardGeneric("getSymbol"))

#' @export
setMethod("getSymbol", "Gene", function(object) object@symbol)

#' Set the Symbol of a Gene Object
#'
#' This method sets the symbol of a gene object.
#' @param object An object of class Gene.
#' @param value The new symbol value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53")
#' gene <- setSymbol(gene, "NEW_SYMBOL")
#' @export
#' @aliases setSymbol,Gene-method
#' @rdname setSymbol-method
setGeneric("setSymbol", function(object, value) standardGeneric("setSymbol"))

setMethod("setSymbol", "Gene", function(object, value) {
  object@symbol <- value
  return(object)
})

#' Set the Name of a Gene Object
#'
#' This method sets the name of a gene object.
#' @param object An object of class Gene.
#' @param value The new name value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", name = "Tumor protein p53")
#' gene <- setName(gene, "New Gene Name")
#' @export
#' @aliases setName,Gene-method
#' @rdname setName-method
setGeneric("setName", function(object, value) standardGeneric("setName"))

setMethod("setName", "Gene", function(object, value) {
  object@name <- value
  return(object)
})

#' Set the Description of a Gene Object
#'
#' This method sets the description of a gene object.
#' @param object An object of class Gene.
#' @param value The new description value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", description = "Tumor protein p53")
#' gene <- setDescription(gene, "New Description")
#' @export
#' @aliases setDescription,Gene-method
#' @rdname setDescription-method
setGeneric("setDescription", function(object, value) standardGeneric("setDescription"))

setMethod("setDescription", "Gene", function(object, value) {
  object@description <- value
  return(object)
})

#' Set the Chromosome of a Gene Object
#'
#' This method sets the chromosome of a gene object.
#' @param object An object of class Gene.
#' @param value The new chromosome value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", chromosome = "17")
#' gene <- setChromosome(gene, "18")
#' @export
#' @aliases setChromosome,Gene-method
#' @rdname setChromosome-method
setGeneric("setChromosome", function(object, value) standardGeneric("setChromosome"))

setMethod("setChromosome", "Gene", function(object, value) {
  object@chromosome <- value
  return(object)
})

#' Set the Start Position of a Gene Object
#'
#' This method sets the start position of a gene object.
#' @param object An object of class Gene.
#' @param value The new start position value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", start = 7565097)
#' gene <- setStart(gene, 7566000)
#' @export
#' @aliases setStart,Gene-method
#' @rdname setStart-method
setGeneric("setStart", function(object, value) standardGeneric("setStart"))

setMethod("setStart", "Gene", function(object, value) {
  object@start <- value
  return(object)
})

#' Set the End Position of a Gene Object
#'
#' This method sets the end position of a gene object.
#' @param object An object of class Gene.
#' @param value The new end position value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", end = 7590856)
#' gene <- setEnd(gene, 7591000)
#' @export
#' @aliases setEnd,Gene-method
#' @rdname setEnd-method
setGeneric("setEnd", function(object, value) standardGeneric("setEnd"))

setMethod("setEnd", "Gene", function(object, value) {
  object@end <- value
  return(object)
})

#' Set the Strand of a Gene Object
#'
#' This method sets the strand of a gene object.
#' @param object An object of class Gene.
#' @param value The new strand value to assign.
#' @return The modified gene object.
#' @examples
#' gene <- createGene("protein_coding", ID = "ENSG000001", strand = "+")
#' gene <- setStrand(gene, "-")
#' @export
#' @aliases setStrand,Gene-method
#' @rdname setStrand-method
setGeneric("setStrand", function(object, value) standardGeneric("setStrand"))

setMethod("setStrand", "Gene", function(object, value) {
  object@strand <- value
  return(object)
})

# Protein-Coding Gene Class Definition
# This class defines a protein-coding gene, inheriting from the Gene class.

#' Protein-Coding Gene Class
#'
#' A class representing a protein-coding gene with specific attributes.
#' @name ProteinCodingGene-class
#' @slot proteinID A character string representing the protein ID.
#' @slot proteinSequence A character string representing the protein sequence.
#' @slot exonCount An integer representing the number of exons.
#' @export
setClass("ProteinCodingGene",
  contains = "Gene",
  slots = list(
    proteinID = "character",
    proteinSequence = "character",
    exonCount = "integer"
  )
)

# Long Non-Coding RNA Gene Class Definition
# This class defines a long non-coding RNA gene, inheriting from the Gene class.

#' Long Non-Coding RNA Gene Class
#'
#' A class representing a long non-coding RNA gene with specific attributes.
#' @name LongNonCodingRNA-class
#' @slot lncRNAID A character string representing the lncRNA ID.
#' @slot RNASequence A character string representing the RNA sequence.
#' @export
setClass("LongNonCodingRNA",
  contains = "Gene",
  slots = list(
    lncRNAID = "character",
    RNASequence = "character"
  )
)

# MicroRNA Gene Class Definition
# This class defines a microRNA gene, inheriting from the Gene class.

#' MicroRNA Gene Class
#'
#' A class representing a microRNA gene with specific attributes.
#' @name MicroRNA-class
#' @slot miRNAID A character string representing the miRNA ID.
#' @slot seedSequence A character string representing the miRNA seed sequence.
#' @export
setClass("MicroRNA",
  contains = "Gene",
  slots = list(
    miRNAID = "character",
    seedSequence = "character"
  )
)
# Transcript Class Definition
# This class defines a transcript, inheriting from the Gene class.

#' Transcript Class
#'
#' A class representing a transcript, inheriting from the Gene class.
#' @name Transcript-class
#' @slot transcriptID A character string representing the transcript ID.
#' @slot sequence A character string representing the transcript sequence.
#' @export
setClass("Transcript",
  contains = "Gene",
  slots = list(
    transcriptID = "character",
    sequence = "character"
  )
)

#' Create a Gene Object
#'
#' This function creates a gene object of the specified type.
#' @param type A character string specifying the type of gene ('protein_coding', 'lncRNA', 'miRNA', 'transcript').
#' @param ... Additional parameters for the gene object.
#' @return An object of the specified gene class.
#' @examples
#' gene <- createGene("protein_coding",
#'   ID = "ENSG000001", symbol = "TP53",
#'   name = "Tumor protein p53", description = "Plays a role in cell cycle regulation",
#'   chromosome = "17", start = 7565097, end = 7590856, strand = "+",
#'   proteinID = "NP_000537", proteinSequence = "MDM2", exonCount = 11L
#' )
#' getSymbol(gene)
#' lengthProduct(gene)
#' gcContent(gene)
#' @export
createGene <- function(type, ...) {
  if (type == "protein_coding") {
    return(new("ProteinCodingGene", ...))
  } else if (type == "lncRNA") {
    return(new("LongNonCodingRNA", ...))
  } else if (type == "miRNA") {
    return(new("MicroRNA", ...))
  } else if (type == "transcript") {
    return(new("Transcript", ...))
  } else {
    stop("Unknown gene type")
  }
}

#' Get the Length of the Gene Product
#'
#' This method returns the length of the gene product, depending on the gene type.
#' @param object An object of a specific gene class.
#' @return The length of the gene product.
#' @examples
#' gene <- createGene("protein_coding",
#'   ID = "ENSG000001", symbol = "TP53",
#'   proteinSequence = "MDM2"
#' )
#' lengthProduct(gene)
#' @export
#' @aliases lengthProduct lengthProduct,ProteinCodingGene-method lengthProduct,LongNonCodingRNA-method lengthProduct,MicroRNA-method lengthProduct,Transcript-method
#' @rdname lengthProduct-method
setGeneric("lengthProduct", function(object) standardGeneric("lengthProduct"))

#' @export
setMethod("lengthProduct", "ProteinCodingGene", function(object) {
  nchar(object@proteinSequence)
})

#' @export
setMethod("lengthProduct", "LongNonCodingRNA", function(object) {
  nchar(object@RNASequence)
})

#' @export
setMethod("lengthProduct", "MicroRNA", function(object) {
  nchar(object@seedSequence)
})

#' @export
setMethod("lengthProduct", "Transcript", function(object) {
  nchar(object@sequence)
})

#' Calculate the GC Content of a Gene Sequence
#'
#' This method calculates the GC content (percentage of G and C nucleotides) of a gene's sequence.
#' @param object An object of class ProteinCodingGene or Transcript.
#' @return The GC content as a percentage.
#' @examples
#' gene <- createGene("protein_coding",
#'   ID = "ENSG000001", symbol = "TP53",
#'   proteinSequence = "MDM2"
#' )
#' gcContent(gene)
#' transcript <- createGene("transcript", ID = "ENST000002", sequence = "ATGCATGC")
#' gcContent(transcript)
#' @export
#' @aliases gcContent gcContent,ProteinCodingGene-method gcContent,Transcript-method
#' @rdname gcContent-method
setGeneric("gcContent", function(object) standardGeneric("gcContent"))

#' @export
setMethod("gcContent", "ProteinCodingGene", function(object) {
  sequence <- object@proteinSequence
  gc <- sum(toupper(strsplit(sequence, NULL)[[1]]) %in% c("G", "C"))
  return((gc / nchar(sequence)) * 100)
})

#' @export
setMethod("gcContent", "Transcript", function(object) {
  sequence <- object@sequence
  gc <- sum(toupper(strsplit(sequence, NULL)[[1]]) %in% c("G", "C"))
  return((gc / nchar(sequence)) * 100)
})
