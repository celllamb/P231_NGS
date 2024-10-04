# Multi-species Genomic Analysis Tool for P231 Study

## Introduction

This Python script is designed for the P231 study to analyze the proportion of specific meat in the SM2 sample using next-generation sequencing data. It performs reference-based mapping of sequencing reads to multiple species' genomes, calculates the proportion of each species in the sample, and provides detailed analysis of mapping results.

## Features

- Multi-species read mapping and analysis
- Species proportion calculation with confidence intervals
- Multiple mapping stringency settings
- Genome-size normalized proportions

## Alignment Definitions

To better understand the analysis results, it's important to be familiar with the following alignment types:

### Definitions

- **Primary**: A read alignment that the aligner considers to be the "best" alignment. Each read can have at most one primary alignment across all species.
- **Secondary**: Alternative alignments of a read that are not considered the best.
- **Supplementary**: Alignments of parts of a read that did not align as part of the primary or secondary alignments.

These definitions are crucial for interpreting the results of the multi-species genomic analysis. The tool primarily focuses on primary alignments to calculate species proportions, ensuring that each read contributes to only one species in the final analysis.

## Prerequisites

- Python 3.6 or higher
- BWA (Burrows-Wheeler Aligner)
- Samtools
- Python libraries: pandas, matplotlib, numpy, scipy, biopython

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/celllamb/P231_NGS.git
   cd P231_NGS
   ```

2. Install required Python libraries:
   ```bash
   pip install pandas matplotlib numpy scipy biopython
   ```

3. Ensure BWA and Samtools are installed and accessible in your PATH.

## Hardcoded Parameters

The following parameters are hardcoded in the script. If you need to modify these, you'll need to edit the script directly.

```python
# Species analyzed
species = ['beef', 'chicken', 'goat', 'horse', 'pork', 'sheep']

# FASTA file names for each species
FASTA_FILES = {
    'beef': 'beef.fa',
    'chicken': 'chicken.fa',
    'goat': 'goat.fa',
    'horse': 'horse.fa',
    'pork': 'pork.fa',
    'sheep': 'sheep.fa'
}

# Output directory
OUTPUT_DIR = "output"

# Sequencing read file names
UNKNOWN_SAMPLE_R1 = "SM2_trim_1.fastq.gz"
UNKNOWN_SAMPLE_R2 = "SM2_trim_2.fastq.gz"

# Mapping settings
MAPPING_SETTINGS = [
    {"name": "default", "params": "-a"},
    {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3 -a"},
    {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4 -a"}
]
```

These parameters define the species being analyzed, the input file names, output directory, and the different mapping stringency settings used in the analysis.

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

1. `mapped_{setting}.bam`: BAM file of alignments for each mapping setting.
2. `genome_sizes.txt`: Calculated genome sizes for each species.
3. `reference_to_species_mapping.tsv`: Mapping of reference sequences to species.
4. `primary_alignment_results.tsv`: Comprehensive results including read counts, proportions, and various metrics.
5. `primary_alignment_proportions.png`: Visualization of mapping results.
6. `analysis_log.txt`: Detailed log of the analysis process with timestamps.

## Species Proportion Calculations

This tool calculates species proportions in two ways:

1. Raw Proportion: The number of primary alignments for a species divided by the total number of primary alignments across all species.

2. Normalized Proportion: The raw proportion normalized by the genome size of each species. This accounts for differences in genome sizes between species.

The calculations are performed as follows:

```python
raw_proportion = primary_reads / category_total
normalized_proportion = (primary_reads / genome_size) / sum((primary_reads / genome_size) for all species)
```

Where:
- `primary_reads` is the number of primary alignments for a species
- `category_total` is the total number of primary alignments across all species
- `genome_size` is the size of the genome for each species

## Confidence Interval Calculations

Confidence intervals (CI) for the proportions are calculated using a multinomial distribution approach with bootstrapping. This method considers all species simultaneously, providing a more accurate estimation of uncertainty.

The process is as follows:

1. Count the number of primary alignments for each species.
2. Use these counts to define a multinomial distribution.
3. Perform bootstrap sampling from this multinomial distribution.
4. Calculate proportions for each bootstrapped sample.
5. Use the 2.5th and 97.5th percentiles of these bootstrapped proportions as the lower and upper bounds of the 95% confidence interval for each species.

The code for this calculation is:

```python
def calculate_multinomial_ci(counts, num_bootstraps=10000, confidence_level=0.95):
    total = sum(counts)
    probs = np.array(counts) / total
    
    bootstrapped_props = []
    for _ in range(num_bootstraps):
        sample = multinomial.rvs(n=total, p=probs)
        bootstrapped_props.append(sample / total)
    
    bootstrapped_props = np.array(bootstrapped_props)
    
    ci_lower = np.percentile(bootstrapped_props, (1 - confidence_level) / 2 * 100, axis=0)
    ci_upper = np.percentile(bootstrapped_props, (1 + confidence_level) / 2 * 100, axis=0)
    
    return ci_lower, ci_upper
```

### Advantages of the Multinomial Approach

1. Simultaneous Consideration: All species are considered together, reflecting the interdependence of their proportions.
2. Consistency: Ensures that the sum of proportions always equals 1, maintaining logical consistency.
3. Accuracy: Provides more accurate confidence intervals, especially for species with low representation.
4. Interpretability: Better reflects the nature of the data, where each read is assigned to exactly one species.

These confidence intervals provide a range of plausible values for the true proportion of each species, with 95% confidence. The multinomial approach ensures that the uncertainty in all species proportions is considered collectively, providing a more robust statistical framework for interpretation.

## Performance Considerations

- The script utilizes samtools for efficient BAM file processing.
- Multi-threading is implemented for various stages of the analysis, including samtools operations:
  - `samtools view` and `samtools sort` use the `-@` option to specify the number of threads.
  - `samtools sort` memory usage:
    - Memory per thread is set to 2GB (`-m 2G`)
    - Total memory usage is approximately 2GB * (number of threads)
    - This setting is optimized for large datasets, but may be excessive for smaller files
- The calculation of primary alignment ratios is efficient, especially for large datasets.
- The analysis of primary alignments uses samtools idxstats, which provides rapid statistics generation without reading through the entire BAM file.

Note: If you're working with smaller datasets or have limited system memory, you may want to adjust the memory usage for `samtools sort` in the script.

## Notes

- This tool is designed specifically for the P231 study and may require modifications for use in other contexts.
- The accuracy of results depends on the quality of input sequencing data and reference genomes.
- For samples with very low abundance of certain species, consider increasing sequencing depth for more accurate detection.
- The analysis progress is displayed in real-time on the screen and also logged to a file for later review.
- Performance may vary depending on the system specifications and input data size.

## Troubleshooting

If you encounter any issues:
1. Check that all required input files are present and correctly named.
2. Ensure all dependencies (BWA, samtools, Python libraries) are installed and up to date.
3. Review the `analysis_log.txt` file for detailed error messages and execution steps.
4. If issues persist, please contact me with the log file and a description of the problem.

## License

You can use this script freely with proper citations.

## Contact

For inquiries, please contact Seil Kim at [stapler@kriss.re.kr](mailto:stapler@kriss.re.kr).