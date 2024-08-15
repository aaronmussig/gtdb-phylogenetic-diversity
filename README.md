# gtdb-phylogenetic-diversity

Snakemake pipeline to calculate the phylogenetic diversity of each phylum in GTDB R220.

## Dependencies

The only dependency is Snakemake (version >= 7.32.4) which is assumed to be on the command line.

## Installation

Clone the repository:

```bash
git clone git@github.com:aaronmussig/gtdb-phylogenetic-diversity.git
cd gtdb-phylogenetic-diversity
```

## Running the pipeline

To run the pipeline, simply execute the following command:

```bash
snakemake --use-conda --cores all
```

## Pipeline steps

1. Trees and taxonomy files for the Archaeal and Bacterial domains are downloaded.
2. Both trees are scaled by their relative evolutionary divergence (RED) using PhyloRank.
3. A script will then calculate the phylogenetic diversity of each phylum in GTDB R220 by the sum of their branch
   lengths.

Both the RED and non-scaled trees will produce values.

Output files are saved under: `./output/`
