# Biofilm-rnaseq-practice
RNA-seq differential expression practice on a bacterial biofilm dataset (DESeq2, R)
# Biofilm RNA-seq practice: P. aeruginosa

Practice differential expression analysis of a public RNA-seq dataset,
comparing biofilm and planktonic growth. This is a learning project using
published data. I did not generate the data.

## Dataset
- GEO accession: GSE302777
- Organism: Pseudomonas aeruginosa (PAO1; the count file name says PAO1DK, so confirm on the GEO page)
- Platform: Illumina NovaSeq 6000
- Samples: 12
- Design: 4 conditions x 3 biological replicates
  - planktonic, untreated control
  - planktonic, midazolam 10 ug/ml
  - biofilm, untreated control
  - biofilm, midazolam 10 ug/ml
- Main comparison for this project: biofilm vs planktonic (untreated controls)

## Source paper
Barcelo IM et al. Analyzing the effects of benzodiazepines on the virulence
and biofilm formation of Pseudomonas aeruginosa. Scientific Reports 2025;16:2923.
DOI: 10.1038/s41598-025-32848-4

## Files
- data/ : count table and sample metadata downloaded from GEO
- scripts/ : analysis code (R, DESeq2)
- results/ : result tables
- figures/ : plots

## Status
In progress. Updated as each analysis step is finished.