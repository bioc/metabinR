
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metabinR

<!-- badges: start
[![BioC status]
(http://www.bioconductor.org/shields/build/release/bioc/metabinR.svg)]
(https://bioconductor.org/checkResults/release/bioc-LATEST/metabinR)
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
and composition based binning on metagenomic samples, directly from
FASTA or FASTQ files.

Abundance based binning is performed by analyzing sequences with long
kmers (k\>8), whereas composition based binning is performed by
utilizing short kmers (k\<8).

## Requirements

A JDK, at least 8, is required and needs to be present before installing
`metabinR`.

## Citation

If you find my work useful, please cite the following publication : A
scalable assembly-free variable selection algorithm for biomarker
discovery from metagenomes.
[10.1186/s12859-016-1186-3](https://doi.org/10.1186/s12859-016-1186-3)

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

Toy gzip compressed fasta simulated metagenome file is provided in
`inst/extdata`.

We use [CAMISIM](https://github.com/CAMI-challenge/CAMISIM) and simulate
a toy metagenome. See file `CAMISIM_config.ini` for parameters.

Paired-end 2x150bp Illumina reads are sampled from 10 bacterial genomes
(see file `genome_to_id.tsv`) with abundances drawn from a log-normal
distribution (see file `distribution_0.txt`).

In total 4Mbases are sampled, yielding 26664 reads (13332 pairs) (see
file `reads.metagenome.fasta.gz`).

Original mapping information of each simulated read (from which genome
each read is originating from) is available in `reads_mapping.tsv.gz`.

### CAMISIM_config.ini

Parameters used in `CAMISIM` to generate the simulated metagenome.

### genome_to_id.tsv

Bacterial genomes and their corresponding id that have been used to
generate the simulated metagenome.

### distribution_0.txt

Abundances of sampled genomes and the Abundance class each one belongs
to. For this toy simulated metagenome, we assume 2 Abundance classes.
Class 1 of most abundant taxa and class 2 of less abundant ones.

### reads.metagenome.fasta.gz

Gzip compressed fasta file of the simulated 26664 reads (13332 pairs,
2x150bp).

### reads_mapping.tsv.gz

Gzip compressed tabular file of the original mapping information of each
simulated read.

Read id from `reads.metagenome.fasta.gz` matches column
`anonymous_read_id` and read genome id matches column `genome_id`.

## Memory requirements

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

## Abundance based Binning

The Abundance based Binning module analyzes input fasta/fastq sequences
by long kmers (k\>8), in order to capture differences in the abundance
signal of organisms in the input data.

Initially a kmer count dictionary is created. An efficient kmer counter
is implemented in the java backend. It supports multi-threaded counting
and small memory footprint. This is the most expensive operation of the
module, due to the exponential size of kmer space (number of distinct
kmers is 4^k). For example the number of distinct 20mers is one
trillion.

Next an EM routine tries to identify Poisson mixtures in the kmer count
spectrum. A Poisson distribution is fitted for each abundance cluster.
The number of abundance clusters cannot be automatically determined and
therefore needs to be user defined. The characteristics of each fitted
Poisson distribution are linked to the abundance level of the taxa
belonging to this abundance class and the their total genome size.

For each abundance class and fitted Poisson distribution, a vector
structure in the kmer spectrum space is constructed. A vector in the
same kmer space is also constructed for each input fasta sequence,
utilizing a tf\*idf weighting scheme.

Finally, cosine distance is calculated between input sequence vectors
and abundance class vectors. Each input sequence is assigned to the
abundance class of least distance.

## Composition based Binning

The Composition based Binning module analyzes input sequences by short
kmers (k\<8), to capture differences in the composition signal and finer
details between similar organisms in the input data.

This module does not involve an expensive kmer counting step, since
short kmer space is relatively small. For example the number of distinct
8mers is around 60 thousands.

For every input sequence a vector is constructed in the kmer space. A
kmeans clustering algorithms is run on the input vectors. Distance
measure is calculated with Spearman footrule. An efficient
multi-threaded implementation of the kmeans algorithm is provided by the
java backend. The number of composition clusters cannot be automatically
determined and therefore needs to be user defined.

Finally, input sequences are assigned to the closest composition
cluster.

## 2step
