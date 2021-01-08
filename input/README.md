# CACTUS Input files

This is the directory containing input data necessary to reproduce analyses presented in the manuscript:  

> **CACTUS: integrating clonal architecture with genomic clustering and transcriptome profiling of single tumor cells**  
> Shadi Darvish Shafighi, Szymon M Kiełbasa, Julieta Sepúlveda Yáñez, Ramin Monajemi, Davy Cats, Hailiang Mei, Roberta Menafra, Susan Kloet, Hendrik Veelken, Cornelis A.M. van Bergen, Ewa Szczurek

Two patients (denoted `S144` and `S12118`) suffering from follicular lymphoma (FL) are studied.

## WES

[VarScan](http://varscan.sourceforge.net/) somatic variant calling of whole exome sequencing of FL and stromal samples produced the following mutation calls:

- for subject [S144](WES/S144.vcf.gz)
- for subject [S12118](WES/S12118.vcf.gz)

Using VarScan output, we also generated the R object of the selected mutations for both subjects provided in [WES/wes.rds](WES/wes.rds). We selected mutations with somatic p-value (SPV in Varscan) less than 0.1 which have seen in the single-cell data as well. 

## Copy number analysis

[Falconx](https://cran.r-project.org/web/packages/falconx/index.html) produced copy number regarding each chromosome which can be found here:

- for subject [S144](CopyNumber/S144)
- for subject [S12118](CopyNumber/S12118)

We also produced the R objects of the whole chromosomes for both subjects:

- for subject [S144](CopyNumber/S144.rds)
- for subject [S12118](CopyNumber/S12118.rds)

## scRNA-seq

Allele-specific transcript counts are provided in [scRNA/ac.rds](scRNA/ac.rds) file.
The file is provided in the R rds format. The contents represent a list which can be loaded with the following command:

```{r}
l <- readRDS("scRNA/ac.rds")
```

The list contains two data frames corresponding to each of the subjects.
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

Each row describes observed (non-zero) allelic counts at a certain mutation position in a certain cell.  
The following columns are important for the analysis:

- `refContig` and `refPos` describe position of the mutation (overlapping with a position from WES)
- `refAllele` is the reference genome base at the position
- there are `cnt` transcripts with `base` nucleotide at the mutation position (this is the count of *unique* transcripts, calculated as a number of unique UMIs at that position)
- `cell` is the identifier of a single cell


## BCR and GEX clusterings

BCR and GEX clusterings are provided in [scBCR/scBCR_GEX.rds](scBCR/scBCR_GEX.rds) file stored in R rds format.  
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

Each row describes an assignment of a single cell to a BCR heavy chain cluster or to gene expression cluster.  
The following columns are provided:

- `cell` is the identifier of a single cell
- `subject` describes from which subject the cell originates
- `cTypeHC` is an identifier shared by cells with the same BCR sequence
- `genExprCluster` is an identifier shared by cells with similar gene expression profile

## Genotype of the clones

Genotype of the clones are inferred using [canopy](https://github.com/yuchaojiang/Canopy) method and the results are provided in [tree](tree).  
The file contains a list of trees which canopy returns.  
Each list element is a single data frame of the following format:

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
