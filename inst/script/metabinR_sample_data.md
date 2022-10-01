## Sample data

Toy gzip compressed fasta simulated metagenome file 
is provided in `inst/extdata`.

We use [CAMISIM](https://github.com/CAMI-challenge/CAMISIM) and simulate a toy
metagenome. See file `CAMISIM_config.ini` for parameters.

Paired-end 2x150bp Illumina reads are sampled from 10 bacterial genomes 
(see file `genome_to_id.tsv`) with abundances drawn 
from a log-normal distribution (see file `distribution_0.txt`). 

In total 4Mbases are sampled, yielding 26664 reads (13332 pairs) 
(see file `reads.metagenome.fasta.gz`).

Original mapping information of each simulated read 
(from which genome each read is originating from) is available in 
`reads_mapping.tsv.gz`.

### CAMISIM_config.ini
Parameters used in `CAMISIM` to generate the simulated metagenome.

### genome_to_id.tsv
Bacterial genomes and their corresponding id that have been used to generate 
the simulated metagenome.

### distribution_0.txt
Abundances of sampled genomes and the Abundance class each one belongs to. 
For this toy simulated metagenome, we assume 2 Abundance classes. 
Class 1 of most abundant taxa and class 2 of less abundant ones. 

### reads.metagenome.fasta.gz
Gzip compressed fasta file of the simulated 26664 reads (13332 pairs, 2x150bp).

### reads_mapping.tsv.gz
Gzip compressed tabular file of the original mapping information of each 
simulated read.

Read id from `reads.metagenome.fasta.gz` matches column `anonymous_read_id` and
read genome id matches column `genome_id`.
