# CACTUS

CACTUS: Clonal Architecture of genomic Clustering of TUmor Single cells

# Input files

## WES

Varscan somatic variant calling of whole exome sequencing of FL and stromal samples produced the following files with mutation calls:

- for subject [S144](input/WES/S144.vcf.gz)
- for subject [S12118](input/WES/S12118.vcf.gz)

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
