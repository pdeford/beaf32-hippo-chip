# beaf32-hippo-chip
This repository contains the code used to generate a figure showing the presence of BEAF binding to the genes encoding for members of the Hippo pathway in *Drosophila Melanogaster*.

## ChIP Analysis:
BEAF32/BEAF32A ChIP for *Drosophila Melanogaster* was obtained from [modENCODE](http://www.modencode.org). There were six available studies, after excluding those with a biological RNAi target. For the purpose of this analysis, the seven sets of preprocessed peak calls available in the ['dmel-interpreted-1'](ftp://data.modencode.org/all_files/dmel-interpreted-1/)  directory were used. The midpoint of each identified peak region was used to annotate the gene diagrams, based on the FlyBase dm3 gene annotation on the [UCSC Genome Browser](genome.ucsc.edu). The scripts used can be found [here](https://github.com/pdeford/beaf32-hippo-chip)

## Additional Details:
* The `hippo_genes.txt` file was created manually by looking up the ID for each gene in the pathway on [FlyBase](flybase.org).
* `beaf_links.txt` was produced by selecting the **BEAF-32** and **BEAF32A and BEAF42B** assay factor filters, and selecting the appropriate studies (those that did not do RNAi for a specific factor). The download URLs were attained and filtered for the interpreted results.
* `flyBase_genes_dm3.bed` was downloaded using the [UCSC Genome Browser's](genome.ucsc.edu) Table Browser
* `hippo_gene_exons.tsv` was obtained by uploading the Hippo gene identifiers to the [UCSC Genome Browser](genome.ucsc.edu), and specifying the Name, Chrom, Strand, ExonStarts, ExonEnds fields from the table.
