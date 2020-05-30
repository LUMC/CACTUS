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
