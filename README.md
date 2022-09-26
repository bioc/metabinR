
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metabinR

<!-- badges: start -->

[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/metabinR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/metabinR)
<!-- badges: end -->

Metagenomics holds great promises for deepening our knowledge of key
bacterial driven processes, but metagenome assembly remains problematic,
typically resulting in representation biases and discarding significant
amounts of non-redundant sequence information. In order to alleviate
constraints assembly can impose on downstream analyses, and/or to
increase the fraction of raw reads assembled via targeted assemblies
relying on pre-assembly binning steps, we developed a set of binning
modules and evaluated their combination in a new “assembly-free” binning
protocol.

The goal of metabinR is to provide functions for performing abundance
and compositional based binning on metagenomic samples, directly from
FASTA or FASTQ files.

## Requirements

A JDK, at least 8, is required and needs to be present before installing
`metabinR`.

## Installation

To install `metabinR` package:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("metabinR")
```

You can install the development version of `metabinR` like so:

``` r
devtools::install_github("gkanogiannis/metabinR")
```

## Sample data

Toy fasta metagenome sample data files are provided in `inst/extdata`.

### samples.fasta.gz

Sample metagenome FASTA file containing 1000 fasta sequences. Original
sequences obtained from the 2nd CAMI Challenge (Rhizosphere) available
at [https://data.cami-challenge.org/participate](CAMI%20Challenge),
[https://www.microbiome-cosi.org/cami/cami/cami2](CAMI2%20Challenge) and
[https://frl.publisso.de/data/frl:6425521/](CAMI2%20raw%20data)..

``` r
fastaFile <- system.file("extdata", "samples.metagenome.fasta.gz", package="metabinR")
```

## Memory requirements

At minimum, make sure to allocate for JVM at least 10 bytes per variant
per sample. If there are `n` samples and `m` variants allocate
`10 x n x m` bytes of RAM. For example, for processing a VCF file
containing data for 1 million variants and 1 thousand samples, allocate
at least : 10^6 x 10^3 x 10 = 10^10 bytes = 10GB of RAM. For optimal
execution, allocate more RAM than minimum. This will trigger less times
garbage collections and hence less pauses.

In order to allocate RAM, a special parameter needs to be passed while
JVM initializes. JVM parameters can be passed by setting
`java.parameters` option. The `-Xmx` parameter, followed (without space)
by an integer value and a letter, is used to tell JVM what is the
maximum amount of heap RAM that it can use. The letter in the parameter
(uppercase or lowercase), indicates RAM units. For example, parameters
`-Xmx1024m` or `-Xmx1024M` or `-Xmx1g` or `-Xmx1G`, all allocate 1
Gigabyte or 1024 Megabytes of maximum RAM for JVM.

In order to allocate 3GB of RAM for the JVM, through R code, use:

``` r
options(java.parameters="-Xmx3G")
```

A rough estimation for the required RAM, if sample and variant numbers
are not known, is half the size of the uncompressed VCF file. For
example for processing a VCF file, which uncompressed occupies 2GB of
disk space, allocate 1GB of RAM.

## Section

## Section

## Section

## Section
