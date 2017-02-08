# hiv-crispr-review

This repository contains a bioinformatics analysis of HIV-1 targeted CRISRP/Cas9 gRNAs. HIV-1 sequence data was downloaded from the Los Alamos Database and then processed using the MIT CRISPR score proposed by Hsu et al. Figures and summary statistics are included in the repository.

This analysis should be run in following order:
  1. `SequenceProcess.ipynb` - This loads the sequence data and calculates the gRNA cleavage ability across patient samples.
  2. `ProcessMITResults.ipynb` - This loads the processed data from the previous script, merges the data with that from the MIT Webserver, and generates the Supplemental Table 1 results.
  3. `VisualizeResults.ipynb` - This loads the results from Supplemental Table 1 and generates the publication quality figures.
  
Repository Inventory:
  - `data/LANLData.tar.gz` - Compressed alignments from the LANL database. Downloaded 1/17/2017
  - `data/MITData.tar.gz` - Compressed Genbank files from the MIT Webserver when uploading HXB2
  - `data/gRNAList.xlsx` - List of gRNAs found in publications as of 1/17/2017
  - `results/summary_res_all_gRNAs.xlsx` - Aggregated cleavage results across all gRNAs.
  - `results/SupTable1.xlsx` - Processed results.
  - `results/Table1.xlsx` - Summary results of the gRNAs aggregated by region.
  - `results/Fig1ABCD.png` - Figure describing the results across the HXB2 genome.
