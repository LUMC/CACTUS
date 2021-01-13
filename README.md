# CACTUS

This is the directory containing source code necessary to reproduce analyses presented in the manuscript:  

> **CACTUS: integrating clonal architecture with genomic clustering and transcriptome profiling of single tumor cells**  
> Shadi Darvish Shafighi, Szymon M Kiełbasa, Julieta Sepúlveda Yáñez, Ramin Monajemi, Davy Cats, Hailiang Mei, Roberta Menafra, Susan Kloet, Hendrik Veelken, Cornelis A.M. van Bergen, Ewa Szczurek

CACTUS extends [cardelino](https://github.com/single-cell-genetics/cardelino) and maps single cells to their clones based on comparing the allele specific transcript counts on mutated positions to given clonal genotypes. It additionally takes advantage of having information about clustering of the cells (here, BCR clusters).  

# Input files

Input files are provided in a separate [git repository CACTUS-data](https://github.com/LUMC/CACTUS-data) referenced as a git submodule in directory `input`.  
See the `README.md` file there for the description of the input files.

# Running

CACTUS can be run from the `code` directory as follows:

```
bash cactus_parallel_samples.sh  'WES R-object file'  'allele transcript counts'  'result directory'  'genotype directory' 'clustering file'
```

The arguments of the bash script ('WES R-object file', 'allele transcript counts', 'genotype directory')  are explained in the input section. 'result directory' is the name for the directory which will be generated for the results of CACTUS.

# Manuscript run

For the data presented in the manuscript the following command was used:

```
bash cactus_parallel_samples.sh  "../input/WES/wes.rds"  "../input/scRNA/ac.rds"  "../results/"   "../tree/"  "../input/scBCR_GEX.rds"
```
