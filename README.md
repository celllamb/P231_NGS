# Multi-species Genomic Analysis Tool for P231 Study

## Introduction

This Python script is designed for the P231 study to analyze the proportion of specific meat in the SM2 sample using next-generation sequencing data. It performs reference-based mapping of sequencing reads to multiple species' genomes, calculates the proportion of each species in the sample, and provides detailed analysis of mapping results.

## Features

- Multi-species read mapping and analysis
- Species proportion calculation with confidence intervals
- Multiple mapping stringency settings
- Genome-size normalized proportions

## Workflow

The script follows these main steps:

0. **Reference Sequence Preparation** (Performed before running this script)
   - Map sequencing reads of individual meat samples to NCBI RefSeq
   - Generate new reference sequences using the mapping results
   - Example command for beef:
     ```bash
     bgzip BEEF_filtered.vcf
     tabix -p vcf BEEF_filtered.vcf.gz
     bcftools consensus -f GCF_002263795.3_ARS-UCD2.0_genomic.fna -o beef.fa BEEF_filtered.vcf.gz 2>&1 | tee log.txt
     ```
   - Repeat this process for each species (chicken, goat, horse, pork, sheep)
   - Note: BEEF_filtered.vcf and similar files for other species are provided by the sequencing service provider

1. **Initialization and Parameter Setting**
   - Parse command-line arguments
   - Set up logging
   - Define species and file paths

2. **Genome Size Calculation**
   - Calculate or load genome sizes for each species
   - Save genome sizes for future use

3. **Reference Preparation**
   - Create a combined reference genome from all species
   - Generate BWA index for the combined reference (if not skipped)

4. **Read Counting and Mapping**
   - Calculate total number of read pairs from input FASTQ files
   - For each mapping setting (default, strict, very_strict):
     - Map reads to the combined reference using BWA MEM
     - Convert SAM to BAM, sort, and index the BAM file

5. **Alignment Analysis**
   - For each mapping setting:
     - Calculate primary alignment ratio
     - Analyze primary alignments and count reads per species

6. **Data Processing and Calculations**
   - Calculate raw and normalized proportions for each species
   - Compute confidence intervals using bootstrap method

7. **Results Output**
   - Save detailed results to TSV file

8. **Logging**
   - Log analysis completion

This workflow is designed to efficiently process large genomic datasets and provide accurate species proportion estimates with associated confidence intervals. The initial step (Step 0) ensures that the reference sequences used in the analysis are tailored to the specific samples being studied.

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
- seqkit
- Python libraries: scipy, biopython, numpy, pandas

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/celllamb/P231_NGS.git
   cd P231_NGS
   ```

2. Install required Python libraries:
   ```bash
   pip install scipy biopython numpy pandas
   ```

3. Ensure BWA, Samtools, and seqkit are installed and accessible in your PATH.

   **Important Note**: Do not install seqkit using conda. The script may not function correctly if seqkit is installed via conda. Please install seqkit directly from its official source or using a package manager that adds it to your system PATH.

## Hardcoded Parameters

The following parameters are hardcoded in the script. If you need to modify these, you'll need to edit the script directly.

```python
# Species analyzed
SPECIES = ['beef', 'chicken', 'goat', 'horse', 'pork', 'sheep']

# FASTA file names for each species
FASTA_FILES = {sp: f"{sp}.fa" for sp in SPECIES}

# Output directory
OUTPUT_DIR = "output"

# Sequencing read file names
UNKNOWN_SAMPLE_R1 = "SM2_trim_1.fastq.gz"
UNKNOWN_SAMPLE_R2 = "SM2_trim_2.fastq.gz"

# Mapping settings
MAPPING_SETTINGS = [
    {"name": "default", "params": ""},
    {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3"},
    {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4"}
]
```

These parameters define the species being analyzed, the input file names, output directory, and the different mapping stringency settings used in the analysis.

Note: The BWA MEM command uses the parameters specified in the `MAPPING_SETTINGS`. If you need to modify these settings, you can edit the `MAPPING_SETTINGS` in the script.

## Input Files

This tool requires specific input files to be present in the script's directory:

1. **Sequencing Reads**:
   - `SM2_trim_1.fastq.gz`: Forward reads of the trimmed paired-end sequencing data
   - `SM2_trim_2.fastq.gz`: Reverse reads of the trimmed paired-end sequencing data

2. **Reference Genomes**:
   The following FASTA files should be present for each species, generated using the process described in Step 0 of the Workflow:
   - `beef.fa`: Reference genome for beef (Bos taurus)
   - `chicken.fa`: Reference genome for chicken (Gallus gallus)
   - `goat.fa`: Reference genome for goat (Capra hircus)
   - `horse.fa`: Reference genome for horse (Equus caballus)
   - `pork.fa`: Reference genome for pork (Sus scrofa)
   - `sheep.fa`: Reference genome for sheep (Ovis aries)

   Note: These reference genomes are created using filtered VCF files (e.g., BEEF_filtered.vcf) provided by the sequencing service provider and the corresponding NCBI RefSeq genomes.

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

   Example:
   ```bash
   python map2ratio.py -m
   ```
   This will skip the mapping step and start from the analysis of existing BAM files.

## Output

The script generates the following outputs in the `output` directory:

1. `mapped_{setting}.bam`: BAM file of alignments for each mapping setting.
2. `genome_sizes.txt`: Calculated genome sizes for each species.
3. `detailed_results.tsv`: Comprehensive results including read counts, proportions, and various metrics.
4. `analysis_log.txt`: Detailed log of the analysis process with timestamps.

## Species Proportion Calculations

This tool calculates species proportions based on the number of mapped reads relative to the total number of reads. It's important to note that only primary alignments are considered as mapped reads in these calculations.

The tool calculates species proportions in two ways:

1. Raw Proportion: The number of primary alignments (mapped reads) for a species divided by the total number of primary alignments across all species.

2. Normalized Proportion: The raw proportion normalized by the genome size of each species. This accounts for differences in genome sizes between species.

The calculations are performed as follows:

```python
raw_proportion = mapped_reads / total_mapped_reads
normalized_proportion = (mapped_reads / genome_size) / sum((mapped_reads / genome_size) for all species)
```

Where:
- `mapped_reads` is the number of primary alignments for a species
- `total_mapped_reads` is the total number of primary alignments across all species
- `genome_size` is the size of the genome for each species

Note: Only primary alignments are considered as mapped reads. Secondary and supplementary alignments are excluded from these calculations to ensure each read contributes to only one species in the final analysis.

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