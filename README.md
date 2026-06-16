The spz_scRNAseq repository contains scripts used to compare Plasmodium falciparum sporozoites obtained from Anopheles gambiae mosquitos raised under 2 conditions:

"1BF" mosquitos provided 1 infectious bloodmeal; and "2BF" mosquitos provided 1 infectious bloodmeal plus an additional non-infectious bloodmeal 3 days later

Sporozoites were purified from the salivary glands via density gradient centrifugation 11 days after the infections bloodmeal unless otherwise indicated

Single-cell transcriptomes were generated using 10x Genomics + Illumina NovaSeq S4

The first set of scripts describes initial processing from raw Counts (QC, normalization, and integration, dimensionality reduction)

The second set describes downstream analyses of the Seurat object (differential expression and abundance tests, pseudo-time analysis, etc.)

The scripts were run in R 4.2.2 on a ThinkPad P15v with 32 Gb RAM
