
# GenePackage

`GenePackage` is an R package designed to represent different gene types
using S4 classes. It provides functionality for creating, manipulating,
and analyzing gene objects, including protein-coding genes, long
non-coding RNAs, microRNAs, and transcripts.

## Features

- S4 classes to represent various gene types.
- Constructor functions to easily create gene objects.
- Accessor methods for retrieving gene attributes such as gene symbol,
  name, chromosome location, etc.
- Methods to calculate gene product lengths.
- Functions to compute GC content for gene sequences.

## Installation

You can install the development version of `GenePackage` from GitHub
using the following commands:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install GenePackage
devtools::install_github("negin-kafee/GenePackage")
```

## Usage

Here is an example of how to use the `GenePackage`:

``` r
library(GenePackage)

# Creating a protein-coding gene object
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

# Accessing gene attributes
getSymbol(gene)
#> [1] "TP53"

# Calculating the length of the gene product
lengthProduct(gene)
#> [1] 4

# Calculating the GC content of the gene product sequence
gcContent(gene)
#> [1] 0.5
```

## Vignettes

For a comprehensive guide on how to use `GenePackage`, please refer to
the vignette:

``` r
browseVignettes(package = "GenePackage")
```

## Contributing

If you would like to contribute to `GenePackage`, please fork the
repository and submit a pull request. We welcome contributions of all
kinds.

## License

`GenePackage` is licensed under the MIT License.

## Author

- **Negin Kafee** - *Author and Maintainer* - <negin.kafee@yahoo.com>

<!-- -->
