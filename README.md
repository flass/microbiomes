Software for handling metagenomic whole-genome shotgun (WGS) data.  

The [scripts/](https://github.com/flass/microbiomes/tree/master/scripts) folder contains script resource (mostly Python and R) for running and analysing results from popular microbiome programs/pipelines:  
- Kraken, a k-mer-based tool for fast taxonomic classification of metagenomic reads (see Kraken [paper](http://doi.org/10.1186/gb-2014-15-3-r46)). The [parseKrona Python module](https://github.com/flass/microbiomes/blob/master/scripts/kraken/parseKrona.py) annd associated scripts notably allow one to filter reads based on their taxonomic classification by Kraken and efficiently extract them from FASTQ/FASTA sequence files.  
- Phylosift, a pipeline for core taxon abundance estimation in microbiome data using phylogenetic placement of reads matching conserved marker genes (see Phylosift [paper](https://peerj.com/articles/243/))  
- Interproscan (notably as part of the EBI Metagenomics pipeline), a tool for functional annotation of metagenomic reads (see Interproscan [paper](http://doi.org/10.1093/nar/gkv1195))  

This software suit was originally developped for the study published by Lassalle et al. in *Molecular Ecology* (2017) comparing human oral microbiomes from hunter-gatherers vs. traditional farmers  from the Philippines, for which metadata required for reproduing the study can be found in the [data/](https://github.com/flass/microbiomes/tree/master/data) folder.  
Related links:
- published *Molecular Ecology* paper: [http://dx.doi.org/10.1111/mec.14435](http://dx.doi.org/10.1111/mec.14435)
- underlying data and supplementary result files: [Figshare project](https://figshare.com/projects/Oral_microbiomes_from_hunter-gatherers_and_traditional_farmers_from_the_Philippines/23425)
