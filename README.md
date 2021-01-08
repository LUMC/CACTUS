# CACTUS

CACTUS: Clonal Architecture of genomic Clustering of TUmor Single cells

# Input files

## WES

Varscan somatic variant calling of whole exome sequencing of FL and stromal samples produced the following files with mutation calls:

- for subject [S144](input/WES/S144.vcf.gz)
- for subject [S12118](input/WES/S12118.vcf.gz)

Using Varscan output, we also generated the r object of the mutations for both subjects provided in [input/WES/wes.rds](input/WES/wes.rds)

## Copy number analysis

Falconx produced copy number regarding each chromosome which can be found here: 

- for subject [S144](input/CopyNumber/S144)
- for subject [S12118](input/CopyNumber/S12118)

We also produced the r object of the whole chromosome for both subjects:

- for subject [S144](input/CopyNumber/S144.rds)
- for subject [S12118](input/CopyNumber/S12118.rds)


## scRNA-seq

Allele transcript counts are provided in [input/scRNA/ac.rds](input/scRNA/ac.rds) file.
The file is provided in R rds format. The contents represent a list which can be loaded with the following command:

```{r}
l <- readRDS("input/scRNA/ac.rds")
```

The list contains two data frames corresponding to the subjects described in the manuscript (S144, S12118).
Here is the beginning of one of the data frames:

```
  refContig  refPos refAllele base source                   cell cnt  MUTATION_ID cluster
1      chr1 1218897         G    A    K1B K1B:AACTTTCTCTTCCTTC-1   1 chr1 1218897    S144
2      chr1 1218897         G    A    K1B K1B:AGCTCCTTCTGGTATG-1   1 chr1 1218897    S144
3      chr1 1218897         G    A    K1B K1B:CTACGTCGTTTGGCGC-1   1 chr1 1218897    S144
4      chr1 1218897         G    A    K1B K1B:CTCGAGGAGCGATTCT-1   1 chr1 1218897    S144
5      chr1 1218897         G    A    K1B K1B:TCGAGGCTCTGTCCGT-1   1 chr1 1218897    S144
6      chr1 1218897         G    A    K1B K1B:TGCACCTTCACCTCGT-1   1 chr1 1218897    S144
```

Each row describes observed non-zero allelic counts at a certain mutation position in a certain cell.
The following columns are important for the analysis:

- `refContig` and `refPos` describe position of a mutation (overlapping with a position from WES)
- `refAllele` is the reference genome base at the position
- there are `cnt` transcripts with `base` nucleotide at the mutation position
- `cell` is the identifier of a single cell


## BCR and GEX clusterings

BCR and GEX clusterings are provided in [input/scBCR_GEX.rds](input/scBCR_GEX.rds) file stored in R rds format.
The file contains a list with elements corresponding to the studied subjects.
Each list element is a single data frame of the following format:

```
                     cell source subject cTypeHC cTypeLC genExprCluster
2  K1B:AAACCTGAGCCAGTTT-1    K1B    S144      u1       b              0
4  K1B:AAACCTGAGGCATGTG-1    K1B    S144       C       s              3
11 K1B:AAACCTGGTATTAGCC-1    K1B    S144       H       d              2
19 K1B:AAACCTGTCTGGGCCA-1    K1B    S144      u2       c              1
22 K1B:AAACGGGAGGATGGAA-1    K1B    S144       E       a              0
24 K1B:AAACGGGCATAGACTC-1    K1B    S144       A       b              2
```

Each row describes assignment of a single cell to a BCR heavy chain cluster or to gene expression cluster.
The following columns are provided:

- `cell` is the identifier of a single cell
- `subject` describes from which subject the cell originates
- `cTypeHC` is an identifier shared by cells with the same BCR sequence
- `genExprCluster` is an identifier shared by cells with similar gene expression profile

## Genotype of the clones

Genotype of the clones are inferred using [canopy](https://github.com/yuchaojiang/Canopy) method and the result are provided in [input/tree](input/tree). The file contains a list of trees which canopy returns. Each list element is a single data frame of the following format:

```
                clone1 clone2 clone3 clone4 clone5
chr1 1218897         0      1      0      0      0
chr1 15929547        0      0      0      1      1
chr1 16395227        0      1      0      0      0
chr1 21233759        0      1      0      0      0
chr1 33013252        0      0      0      0      1
chr1 33013270        0      0      1      0      0
chr1 33013276        0      0      1      0      0

```
Genotype is a zero/one matrix and each row describes mapping a mutation position to one or more clones. 

# Running

CACTUS, modified [cardelino](https://github.com/single-cell-genetics/cardelino) to take advantage of having information about clustering of the cells (here, BCR clonotypes). CACTUS can be run for the existing data in the code directory using the following script: 

```
bash cactus_parallel_samples.sh  'WES R-object file'  'Allele transcript counts'  'result directory'  'genotype directory' 'Clustering file' 
```
The argumans of the bash script ('WES R-object file', 'Allele transcript counts', 'genotype directory')  are explained in the input section. 'result directory' is the name for the directory which will be generated for the results of CACTUS.


## Test

For the data provided in this repository, the following command can be used:

```
bash cactus_parallel_samples.sh  "../input/WES/wes.rds"  "../input/scRNA/ac.rds"  "../results/"   "../tree/"  "../input/scBCR_GEX.rds" 
```


