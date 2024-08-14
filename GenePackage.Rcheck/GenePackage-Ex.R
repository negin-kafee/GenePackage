pkgname <- "GenePackage"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GenePackage')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("createGene")
### * createGene

flush(stderr()); flush(stdout())

### Name: createGene
### Title: Create a Gene Object
### Aliases: createGene

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53", 
                   name = "Tumor protein p53", description = "Plays a role in cell cycle regulation", 
                   chromosome = "17", start = 7565097, end = 7590856, strand = "+", 
                   proteinID = "NP_000537", proteinSequence = "MDM2", exonCount = 11L)
getSymbol(gene)
lengthProduct(gene)
gcContent(gene)



cleanEx()
nameEx("gcContent-method")
### * gcContent-method

flush(stderr()); flush(stdout())

### Name: gcContent
### Title: Calculate the GC Content of a Gene Sequence
### Aliases: gcContent gcContent,ProteinCodingGene-method
###   gcContent,Transcript-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53", 
                   proteinSequence = "MDM2")
gcContent(gene)
transcript <- createGene("transcript", ID = "ENST000002", sequence = "ATGCATGC")
gcContent(transcript)



cleanEx()
nameEx("getSymbol-method")
### * getSymbol-method

flush(stderr()); flush(stdout())

### Name: getSymbol
### Title: Get the Symbol of a Gene Object
### Aliases: getSymbol getSymbol,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53")
getSymbol(gene)



cleanEx()
nameEx("lengthProduct-method")
### * lengthProduct-method

flush(stderr()); flush(stdout())

### Name: lengthProduct
### Title: Get the Length of the Gene Product
### Aliases: lengthProduct lengthProduct,ProteinCodingGene-method
###   lengthProduct,LongNonCodingRNA-method lengthProduct,MicroRNA-method
###   lengthProduct,Transcript-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53", 
                   proteinSequence = "MDM2")
lengthProduct(gene)



cleanEx()
nameEx("setChromosome-method")
### * setChromosome-method

flush(stderr()); flush(stdout())

### Name: setChromosome
### Title: Set the Chromosome of a Gene Object
### Aliases: setChromosome setChromosome,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", chromosome = "17")
gene <- setChromosome(gene, "18")



cleanEx()
nameEx("setDescription-method")
### * setDescription-method

flush(stderr()); flush(stdout())

### Name: setDescription
### Title: Set the Description of a Gene Object
### Aliases: setDescription setDescription,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", description = "Tumor protein p53")
gene <- setDescription(gene, "New Description")



cleanEx()
nameEx("setEnd-method")
### * setEnd-method

flush(stderr()); flush(stdout())

### Name: setEnd
### Title: Set the End Position of a Gene Object
### Aliases: setEnd setEnd,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", end = 7590856)
gene <- setEnd(gene, 7591000)



cleanEx()
nameEx("setName-method")
### * setName-method

flush(stderr()); flush(stdout())

### Name: setName
### Title: Set the Name of a Gene Object
### Aliases: setName setName,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", name = "Tumor protein p53")
gene <- setName(gene, "New Gene Name")



cleanEx()
nameEx("setStart-method")
### * setStart-method

flush(stderr()); flush(stdout())

### Name: setStart
### Title: Set the Start Position of a Gene Object
### Aliases: setStart setStart,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", start = 7565097)
gene <- setStart(gene, 7566000)



cleanEx()
nameEx("setStrand-method")
### * setStrand-method

flush(stderr()); flush(stdout())

### Name: setStrand
### Title: Set the Strand of a Gene Object
### Aliases: setStrand setStrand,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", strand = "+")
gene <- setStrand(gene, "-")



cleanEx()
nameEx("setSymbol-method")
### * setSymbol-method

flush(stderr()); flush(stdout())

### Name: setSymbol
### Title: Set the Symbol of a Gene Object
### Aliases: setSymbol setSymbol,Gene-method

### ** Examples

gene <- createGene("protein_coding", ID = "ENSG000001", symbol = "TP53")
gene <- setSymbol(gene, "NEW_SYMBOL")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
