# Analysis of Adaptation Patterns in Chironomus riparius

This repository contains scripts and analysis pipelines used to investigate mainly nucelotide diversity and Fst of rapid adaptation patterns in *Chironomus riparius* using Pool-Seq data across multiple generations.

## Methods Overview

The analysis workflow implemented in these scripts consists of several key steps:

1. **Pool-Seq_Data_Processing_workflow**: General workflow: Quality filtering and mapping to the *C. riparius* reference genome, preparation for Popoolation2,  Multiple testing correction and outlier detection
2. **Pool-Seq_Data_Analysis.R**: R-script: Identification of genomic regions with significant differentiation using FST and FET. Genetic Drift Simulation
3. **Genetic_Diversity_Pipeline**: Calculation of genetic diversity measures (π, θ, Tajima's D) using PoPoolation
4. **Processing_InterProScan_Annotations.py**: Analysis of protein domains and gene functions in adaptation-related regions


## Requirements

- PoPoolation (v1.2.2)
- PoPoolation2 (v1.201)
- R (v4.0+) with packages:
  - dplyr
  - ggplot2
- Python (v3.6+) with packages:
  - pandas
  - matplotlib
- Samtools (v1.9+)
- BWA-MEM (v0.7.17+)
- InterProScan (v5.52+)

## Usage

Each script contains detailed documentation within its header. For the main workflow:

3. Process annotations with: `python Processing_InterProScan_Annotation.py`
4. Statistical analysis: `Rscript r-analysis.R`

## Citation

If you use these scripts or methods in your research, please cite:

https://github.com/CoCaliendo/Analysis_AdaptationPattern

## Contact

For questions or issues, please open an issue on this repository or contact the repository owner.
