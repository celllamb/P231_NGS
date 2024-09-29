
# Python script for P231 Study

## Introduction

This Python script is designed for the P231 study to analyze the proportion of specific meat in the SM2 sample using next-generation sequencing data. It performs reference-based mapping of sequencing reads to multiple species' genomes, calculates the proportion of each species in the sample, and provides confidence intervals for these proportions using bootstrap analysis.

## Features

- Combines multiple reference genomes for simultaneous mapping
- Performs read mapping using BWA MEM with customizable settings
- Analyzes unique mappings, single-species multi-mappings, and cross-species multi-mappings
- Calculates species proportions and their confidence intervals using bootstrap analysis
- Generates visualizations of mapping results
- Utilizes multi-core processing for improved performance
- Provides detailed progress updates during the analysis
- Outputs results to both screen and log file simultaneously

## Prerequisites

- Python 3.6 or higher
- BWA (Burrows-Wheeler Aligner)
- Samtools
- Python libraries: pandas, matplotlib, numpy, scipy, pysam

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/celllamb/P231_NGS.git
   cd P231_NGS
   ```

2. Install required Python libraries:
   ```bash
   pip install pandas matplotlib numpy scipy pysam
   ```

3. Ensure BWA and Samtools are installed and accessible in your PATH.

## Usage

1. Prepare your input files:
   - Place your trimmed paired-end sequencing reads (SM2_trim_1.fastq.gz and SM2_trim_2.fastq.gz) in the script's directory.
   - Place reference genome FASTA files (beef.fa, chicken.fa, goat.fa, horse.fa, pork.fa, sheep.fa) in the script's directory.

2. Run the script:
   ```bash
   python map2ratio.py
   ```

## Workflow

1. **Reference Preparation**: 
   - Combines all specified reference genomes into a single FASTA file.
   - Creates a BWA index for the combined reference.

2. **Read Mapping**: 
   - Maps reads to the combined reference using BWA MEM.
   - Performs mapping with three different stringency settings: default, strict, and very strict.

3. **Mapping Analysis**: 
   - Categorizes mapped reads as unique, single-species multi-mapped, or cross-species multi-mapped.
   - Calculates the proportion of reads mapping to each species.

4. **Bootstrap Analysis**: 
   - Performs bootstrap sampling to calculate confidence intervals for species proportions.

5. **Result Generation**: 
   - Saves detailed mapping statistics and species proportions to CSV files.
   - Generates visualizations of mapping results and confidence intervals.

## Output

The script generates the following outputs in the `output` directory:

1. `mapped_{setting}.bam`: BAM file of mapped reads for each mapping setting.
2. `cross_species_multi_map_{setting}.csv`: Details of cross-species multi-mappings.
3. `single_species_multi_map_{setting}.csv`: Counts of single-species multi-mappings.
4. `mapping_proportions_{setting}.csv`: Proportions of different mapping types for each species.
5. `confidence_intervals_{setting}.csv`: Bootstrap confidence intervals for species proportions.
6. `enhanced_mapping_proportions.png`: Visualization of mapping results and confidence intervals.
7. `analysis_log.txt`: Detailed log of the analysis process with timestamps.

## Interpreting Results

- The proportion of unique mappings indicates the confidence in species identification.
- Wide confidence intervals suggest higher uncertainty in the proportion estimate.
- Comparing results across different stringency settings can provide insights into mapping quality and cross-species similarities.
- The unmapped read ratio is provided for each mapping setting.

## Notes

- This tool is designed specifically for the P231 study and may require modifications for use in other contexts.
- The accuracy of results depends on the quality of input sequencing data and reference genomes.
- For samples with very low abundance of certain species, consider increasing sequencing depth for more accurate detection.
- The analysis progress is displayed in real-time on the screen and also logged to a file for later review.

## Contributing

Contributions to improve the tool are welcome. Please feel free to submit pull requests or open issues for bugs and feature requests.

## License

You can use this script freely with proper citations.

## Contact

For inquiries, please contact Seil Kim at [stapler@kriss.re.kr](mailto:stapler@kriss.re.kr).
