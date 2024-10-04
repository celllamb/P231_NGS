# Multi-species Genomic Analysis Tool for P231 Study

## Introduction

This Python script is designed for the P231 study to analyze the proportion of specific meat in the SM2 sample using next-generation sequencing data. It performs reference-based mapping of sequencing reads to multiple species' genomes, calculates the proportion of each species in the sample, and provides detailed analysis of mapping results.

## Features

- Performs multi-species read mapping and analysis
- Calculates species proportions with various metrics (Raw and Normalized proportions)
- Supports multiple mapping stringency settings
- Provides both raw and genome-size normalized proportions
- Generates detailed results including primary, secondary, and supplementary alignments
- Calculates statistics for total mapped reads across all species
- Implements parallel processing for improved performance on large datasets
- Dynamic chunk size determination for optimal memory usage
- Generates visualizations of mapping results

## Prerequisites

- Python 3.6 or higher
- BWA (Burrows-Wheeler Aligner)
- Samtools
- Python libraries: pandas, matplotlib, numpy, scipy, pysam, biopython, psutil

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/celllamb/P231_NGS.git
   cd P231_NGS
   ```

2. Install required Python libraries:
   ```bash
   pip install pandas matplotlib numpy scipy pysam biopython psutil
   ```

3. Ensure BWA and Samtools are installed and accessible in your PATH.

## Input Files

This tool requires specific input files to be present in the script's directory:

1. **Sequencing Reads**:
   - `SM2_trim_1.fastq.gz`: Forward reads of the trimmed paired-end sequencing data
   - `SM2_trim_2.fastq.gz`: Reverse reads of the trimmed paired-end sequencing data

2. **Reference Genomes**:
   The following FASTA files should be present for each species:
   - `beef.fa`: Reference genome for beef (Bos taurus)
   - `chicken.fa`: Reference genome for chicken (Gallus gallus)
   - `goat.fa`: Reference genome for goat (Capra hircus)
   - `horse.fa`: Reference genome for horse (Equus caballus)
   - `pork.fa`: Reference genome for pork (Sus scrofa)
   - `sheep.fa`: Reference genome for sheep (Ovis aries)

Ensure that all these files are present in the same directory as the script before running the analysis.

## Usage

1. Prepare your input files as described in the "Input Files" section.

2. Run the script:
   ```bash
   python map2ratio.py [options]
   ```

   Options:
   - `-i` or `--skip-index`: Skip the BWA indexing step
   - `-m` or `--skip-mapping`: Skip the BWA mapping step (implies -i)
   - `--force-recalculate`: Force recalculation of genome sizes
   - `-r` or `--report-progress`: Report progress during analysis

   Example:
   ```bash
   python map2ratio.py -m -r
   ```
   This will skip the mapping step, start from the analysis of existing BAM files, and report progress during the analysis.

## Output

The script generates the following outputs in the `output` directory:

1. `mapped_{setting}.bam`: BAM file of mapped reads for each mapping setting.
2. `genome_sizes.txt`: Calculated genome sizes for each species.
3. `reference_to_species_mapping.tsv`: Mapping of reference sequences to species.
4. `detailed_results.tsv`: Comprehensive results including read counts, proportions, and various metrics.
5. `enhanced_mapping_proportions.png`: Visualization of mapping results.
6. `analysis_log.txt`: Detailed log of the analysis process with timestamps.

## Read Classification Method

This section explains how reads are classified into primary, secondary, supplementary, multi_within, and multi_across categories.

### Definitions

- **Primary**: A read alignment that the aligner considers to be the "best" alignment. Each read can have at most one primary alignment across all species.
- **Secondary**: Alternative alignments of a read that are not considered the best.
- **Supplementary**: Alignments of parts of a read that did not align as part of the primary or secondary alignments.
- **Total_Multi_Within**: The sum of all secondary and supplementary alignments within a single species' genome.
- **Unique_Multi_Within**: The count of unique reads that have more than one alignment (including secondary and supplementary) within a single species' genome.
- **Total_Multi_Across**: The sum of all alignments (primary, secondary, and supplementary) for reads that map to multiple species.
- **Unique_Multi_Across**: The count of unique reads that have alignments to multiple species.

### Classification Process

1. For each alignment in the BAM file:
   - Classify as Primary, Secondary, or Supplementary based on alignment flags.
   - Record the species to which the read is mapped.

2. After processing all alignments:
   - Calculate Multi_Within and Multi_Across for each species.
   - Compute Unique_Multi_Within and Unique_Multi_Across.

3. Calculate additional statistics:
   - Total mapped reads, reads with/without primary alignments, primary alignment rate, etc.

## Species Proportion Calculations and Interpretation

Our analysis focuses on calculating the relative proportions of different species within the mapped reads. These proportions provide insights into the composition of the sample based on the sequencing data.

### Calculation of Category Totals

- Total Primary Alignments: Sum of Primary alignments for all species
- Total Secondary Alignments: Sum of Secondary alignments for all species
- Total Supplementary Alignments: Sum of Supplementary alignments for all species
- Total Mapped Reads: Sum of (Primary + Secondary + Supplementary) alignments for all species
- Total Unique Multi Within: Sum of Unique Multi Within alignments for all species
- Total Unique Multi Across: Sum of Unique Multi Across alignments for all species

Note: Total Mapped Reads represents the total number of all alignments, including multiple alignments of the same read. This may count some reads multiple times if they have multiple alignments.

### Raw Proportions

For each species, we calculate the following raw proportions:

1. **Raw_Primary_Proportion**: Primary alignments / Total Primary Alignments
2. **Raw_Secondary_Proportion**: Secondary alignments / Total Secondary Alignments
3. **Raw_Supplementary_Proportion**: Supplementary alignments / Total Supplementary Alignments
4. **Raw_Unique_Multi_Within_Proportion**: Unique_Multi_Within / Total Unique Multi Within
5. **Raw_Unique_Multi_Across_Proportion**: Unique_Multi_Across / Total Unique Multi Across
6. **Raw_Total_Proportion**: (Primary + Secondary + Supplementary) / Total Mapped Reads

### Normalized Proportion

To account for differences in genome sizes between species:

**Normalized_Proportion**: ((Primary + Secondary + Supplementary) / Genome Size) / Σ((Primary + Secondary + Supplementary) / Genome Size) for all species

### Interpretation Guidelines

- Each Raw Proportion category sums to 1 (or 100%) across all species.
- Raw_Primary_Proportion indicates the distribution of best alignments among species.
- Raw_Secondary_Proportion and Raw_Supplementary_Proportion show the distribution of alternative alignments.
- Raw_Unique_Multi_Within_Proportion reveals how multi-mapping reads within a single species are distributed.
- Raw_Unique_Multi_Across_Proportion shows the distribution of reads that map to multiple species.
- Raw_Total_Proportion gives an overall view of how all alignments are distributed among species.
- Normalized_Proportion is useful for comparing species with significantly different genome sizes.

When interpreting results:
1. Start with Raw_Total_Proportion for an overall view of species distribution.
2. Compare Raw_Primary_Proportion to understand the distribution of best alignments.
3. Use Raw_Unique_Multi_Within_Proportion and Raw_Unique_Multi_Across_Proportion to assess the level of ambiguous mappings.
4. Refer to Normalized_Proportion for the most balanced comparison between species, especially if there are large differences in genome sizes.
5. Always consider the actual read counts alongside the proportions for a comprehensive understanding.

Remember that these proportions are based on all alignments, which may count some reads multiple times. The overall mapping rate (Total Mapped Reads / Total Input Reads) should be considered to assess the completeness of the analysis.

## Notes

- This tool is designed specifically for the P231 study and may require modifications for use in other contexts.
- The accuracy of results depends on the quality of input sequencing data and reference genomes.
- For samples with very low abundance of certain species, consider increasing sequencing depth for more accurate detection.
- The analysis progress is displayed in real-time on the screen and also logged to a file for later review.
- Performance may vary depending on the system specifications and input data size.

## Troubleshooting

If you encounter any issues:
1. Check that all required input files are present and correctly named.
2. Ensure all dependencies are installed and up to date.
3. Review the `analysis_log.txt` file for detailed error messages and execution steps.
4. If issues persist, please contact me with the log file and a description of the problem.

## License

You can use this script freely with proper citations.

## Contact

For inquiries, please contact Seil Kim at [stapler@kriss.re.kr](mailto:stapler@kriss.re.kr).