# Biofilm RNA-seq practice: P. aeruginosa

Practice differential expression analysis of a public RNA-seq dataset
(DESeq2, R), comparing biofilm and planktonic growth. This is a learning
project using published data. I did not generate the data.

## Dataset
- GEO accession: GSE302777
- Organism: Pseudomonas aeruginosa (PAO1; count file labelled PAO1DK)
- Platform: Illumina NovaSeq 6000
- Samples: 12
- Design: 4 conditions x 3 biological replicates
  - planktonic, untreated control
  - planktonic, midazolam 10 ug/ml
  - biofilm, untreated control
  - biofilm, midazolam 10 ug/ml
- Main comparison: biofilm vs planktonic (untreated controls, 3 vs 3)

## Source paper
Barcelo IM et al. Analyzing the effects of benzodiazepines on the virulence
and biofilm formation of Pseudomonas aeruginosa. Scientific Reports 2025;16:2923.
DOI: 10.1038/s41598-025-32848-4

## Method
- Raw counts filtered to genes with at least 10 reads in total (5,678 -> 5,642 genes).
- Differential expression with DESeq2 (design: ~ growth), planktonic as reference.
- Fold changes shrunk with apeglm for plots; counts of genes use unshrunk results.
- PCA and heatmap use variance-stabilised counts (vst).
- Top 30 genes annotated with KEGG product names (KEGG REST, retrieved 8 November 2026);
  categories were assigned by hand from the product names.
- KEGG pathway enrichment (clusterProfiler, organism pae) run separately for genes
  higher in biofilm and genes higher in planktonic, with all filtered genes as background.

## How to run
1. Open the .Rproj file in RStudio.
2. Run scripts/analysis.R from top to bottom (internet needed for KEGG).

Tested with R 4.6.2, DESeq2 1.52.0, clusterProfiler 4.20.3, KEGGREST 1.52.2,
EnhancedVolcano 1.30.0, ggplot2 4.0.3.

## Results
- Thresholds: adjusted p < 0.05 and |log2 fold change| > 1
- Higher in biofilm: 1,489 genes. Higher in planktonic: 1,254 genes.
- Total counts per sample: about 5.4 to 9.9 million, no low-depth outliers.
- PCA: PC1 explains 98% of the variance and separates the two groups.
- Heatmap of the 30 most significant genes: 22 higher in biofilm, 8 higher in planktonic.
- KEGG, higher in biofilm: two-component system, bacterial chemotaxis, flagellar assembly.
- KEGG, higher in planktonic: ribosome, oxidative phosphorylation, carbon metabolism,
  citrate cycle (TCA cycle).

![PCA](figures/pca.png)
![Volcano plot](figures/volcano.png)
![Heatmap of top 30 genes](figures/heatmap_top30.png)

## Interpretation
Untreated biofilm and planktonic cells of P. aeruginosa differed widely in gene
expression: 1,489 genes were higher in biofilm and 1,254 were higher in planktonic
(adjusted p < 0.05, |log2FC| > 1, n = 3 per group), and PC1 explained 98% of the
variance. Among the 30 most significant genes, 22 were higher in biofilm and 8 were
higher in planktonic, including genes annotated as the fucose-binding lectin LecB, a
chemotaxis transducer and the potassium-transporting ATPase subunits KdpA and KdpB
(biofilm), and glyceraldehyde-3-phosphate dehydrogenase and a cytochrome c oxidase
subunit (planktonic). KEGG enrichment of the genes higher in biofilm highlighted
two-component systems, bacterial chemotaxis and flagellar assembly, whereas genes
higher in planktonic were enriched in the ribosome, oxidative phosphorylation, carbon
metabolism and the TCA cycle. The data show association, not mechanism. This raises
the question of whether the lower expression of protein-synthesis and energy-metabolism
genes in biofilm cells reflects slower growth, and how that relates to stress tolerance.

## Track B: image analysis (practice)
Fiji analysis of one public BiofilmQ sample stack (22 slices): Otsu threshold,
thresholded area per slice, approximate biovolume 45.2 um^3. Practice only, not my own
data. See notes/trackB_imaging_note.md and notes/trackB_data_source.txt.

## Limitations
- Public dataset analysed for practice; I did not generate the data.
- Three biological replicates per group; bulk RNA-seq averages over cells.
- Biofilm and planktonic cultures also differ in growth conditions.
- KEGG enrichment shows where changes cluster, not the mechanism.
- Gene categories for the top 30 were assigned by hand from product names; 9 of 30 are
  hypothetical or of unknown function.
- No experimental validation.

## Files
- data/ : count table and metadata from GEO
- scripts/analysis.R : full RNA-seq analysis
- scripts/trackB_biovolume.R : Track B area and biovolume script
- results/ : result tables (top-30 annotation, KEGG tables, Track B measurements)
- figures/ : QC plot, PCA, volcano plot, heatmap, Track B figures
- notes/ : data source and imaging notes