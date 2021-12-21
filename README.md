# SARS-CoV-2-RBD-Abs-HTDMS
High-throughput deep mutational scanning processing procedure

This repository contains code for processing deep mutational scanning data (it is a modified version of [J. Bloom's pipeline](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS)), and code for reproducing figures and analysis in [our paper](https://biorxiv.org/cgi/content/short/2021.12.07.470392v1).

If you want to begin with raw sequencing data, please download our raw data from [our SRA BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA787091) or prepare your sequencing data in the same .fastq.gz format.

File `raw_processing/datasets/sample_to_ref.csv` gives the corresponding reference file(s) for each antibody replication. Some experiments use more than 1 reference fastq files, and we advice you to merge them for further analysis.

Processed mutation escape score data can also be found at `processed/high_antibody_results.csv`. `processed/plot_with_R.ipynb` contains code for reproducing some figures in the manuscript.

Thanks to [the Bloom lab](https://research.fhcrc.org/bloom/en.html) for their fantastic deep mutational scanning technology and data processing pipeline. Refer to [Sunney Xie lab](https://sunneyxielab.org) for more information.
