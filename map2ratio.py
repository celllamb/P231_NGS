#!/usr/bin/env python3

import argparse
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import multiprocessing
import sys
import time
import traceback
from Bio import SeqIO
from collections import defaultdict

# Analysis Configuration Constants
SPECIES = ['beef', 'chicken', 'goat', 'horse', 'pork', 'sheep']
FASTA_FILES = {sp: f"{sp}.fa" for sp in SPECIES}
MAPPING_SETTINGS = [
    {"name": "default", "params": "-a"},
    {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3 -a"},
    {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4 -a"}
]

# Input/Output Constants
OUTPUT_DIR = "output"
UNKNOWN_SAMPLE_R1 = "SM2_trim_1.fastq.gz"
UNKNOWN_SAMPLE_R2 = "SM2_trim_2.fastq.gz"
def print_progress(message, file=None):
    output = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    print(output)
    if file:
        file.write(output + "\n")
        file.flush()
    sys.stdout.flush()

def run_command(command, log=None):
    print_progress(f"Running command: {command}", log)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    while True:
        output = process.stdout.readline()
        if output == b'' and process.poll() is not None:
            break
        if output:
            print_progress(output.strip().decode(), log)
    
    rc = process.poll()
    if rc != 0:
        error = process.stderr.read().decode()
        print_progress(f"Error executing command: {command}", log)
        print_progress(f"Error message: {error}", log)
        raise subprocess.CalledProcessError(rc, command, error)
    
    return rc

def get_optimal_thread_count():
    total_cores = multiprocessing.cpu_count()
    physical_cores = total_cores // 2
    return max(1, physical_cores)

def calculate_genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            total_length += len(record.seq)
    return total_length

def get_genome_sizes(SPECIES, fasta_files, output_dir, force_recalculate=False):
    genome_sizes_file = os.path.join(output_dir, "genome_sizes.txt")
    
    if not force_recalculate and os.path.exists(genome_sizes_file):
        try:
            genome_sizes = {}
            with open(genome_sizes_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) != 2:
                        raise ValueError(f"Invalid line format in genome sizes file: {line}")
                    sp, size = parts
                    if sp not in SPECIES:
                        raise ValueError(f"Unknown species in genome sizes file: {sp}")
                    try:
                        genome_sizes[sp] = int(size)
                    except ValueError:
                        raise ValueError(f"Invalid genome size for species {sp}: {size}")
            return genome_sizes
        except (IOError, ValueError) as e:
            print_progress(f"Error reading genome sizes file: {str(e)}. Recalculating genome sizes.")
    
    genome_sizes = {}
    for sp in SPECIES:
        if sp not in fasta_files:
            raise ValueError(f"No FASTA file specified for species: {sp}")
        if not os.path.exists(fasta_files[sp]):
            raise FileNotFoundError(f"FASTA file not found for species {sp}: {fasta_files[sp]}")
        try:
            genome_sizes[sp] = calculate_genome_size(fasta_files[sp])
        except Exception as e:
            raise RuntimeError(f"Error calculating genome size for species {sp}: {str(e)}")
    
    try:
        with open(genome_sizes_file, 'w') as f:
            for sp, size in genome_sizes.items():
                f.write(f"{sp}\t{size}\n")
    except IOError as e:
        print_progress(f"Warning: Unable to save genome sizes file: {str(e)}")
    
    return genome_sizes

def calculate_multinomial_ci(counts, num_bootstraps=10000, confidence_level=0.95):
    total = sum(counts)
    probs = np.array(counts) / total
    
    bootstrapped_props = []
    for _ in range(num_bootstraps):
        sample = stats.multinomial.rvs(n=total, p=probs)
        bootstrapped_props.append(sample / total)
    
    bootstrapped_props = np.array(bootstrapped_props)
    
    ci_lower = np.percentile(bootstrapped_props, (1 - confidence_level) / 2 * 100, axis=0)
    ci_upper = np.percentile(bootstrapped_props, (1 + confidence_level) / 2 * 100, axis=0)
    
    return ci_lower, ci_upper

def load_reference_to_species_mapping(filename):
    mapping = {}
    with open(filename, 'r') as f:
        for line in f:
            reference, species = line.strip().split('\t')
            mapping[reference] = species
    return mapping

def save_reference_to_species_mapping(mapping, filename):
    with open(filename, 'w') as f:
        for reference, species in mapping.items():
            f.write(f"{reference}\t{species}\n")

def create_reference_to_species_mapping_from_fasta(fasta_files, SPECIES, log=None):
    reference_to_species = {}
    for sp, fasta_file in zip(SPECIES, fasta_files.values()):
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                reference_to_species[record.id] = sp
    return reference_to_species

def create_combined_reference(references, output_dir, log=None):
    print_progress("Creating combined reference...", log)
    combined_ref = os.path.join(output_dir, "combined_reference.fa")
    with open(combined_ref, 'w') as outfile:
        for ref in references:
            with open(ref, 'r') as infile:
                outfile.write(infile.read())
    print_progress("Combined reference created.", log)
    return combined_ref

def generate_bwa_index(combined_ref, log=None):
    print_progress("Generating BWA index...", log)
    command = f"bwa index -a bwtsw {combined_ref}"
    
    start_time = time.time()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    previous_size = 0
    while process.poll() is None:
        time.sleep(10)  # Check every 10 seconds
        current_size = os.path.getsize(combined_ref)
        if current_size != previous_size:
            elapsed_time = time.time() - start_time
            print_progress(f"Indexing in progress... Elapsed time: {elapsed_time:.0f} seconds. Reference file size: {current_size/1e6:.2f} MB", log)
            previous_size = current_size
    
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print_progress(f"Error executing command: {command}", log)
        print_progress(f"Error message: {stderr.decode()}", log)
        raise subprocess.CalledProcessError(process.returncode, command, stderr)
    
    print_progress("BWA index generated.", log)

def run_bwa_mem(ref_file, read1, read2, output_dir, settings, cpu_count, log=None):
    print_progress(f"Running BWA MEM with {settings['name']} settings...", log)
    
    all_bam_file = os.path.join(output_dir, f"mapped_{settings['name']}.bam")
    command = (f"bwa mem -t {cpu_count} {settings['params']} {ref_file} {read1} {read2} | "
               f"samtools view -b - | "
               f"samtools sort -@ {cpu_count} -m 2G - > {all_bam_file}")
    run_command(command, log)
    run_command(f"samtools index {all_bam_file}", log)
    
    print_progress(f"BWA MEM completed for {settings['name']} settings.", log)
    return all_bam_file

def calculate_primary_alignment_ratio(all_bam_file, cpu_count, log=None):
    print_progress("Calculating primary alignment ratio...", log)
    
    cmd = f"samtools view -c -F 4 {all_bam_file} && samtools view -c -F 260 {all_bam_file}"
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        total_mapped, primary_mapped = map(int, result.stdout.strip().split('\n'))
    except subprocess.CalledProcessError as e:
        print_progress(f"Error executing samtools: {e}", log)
        print_progress(f"Error output: {e.stderr}", log)
        raise
    
    ratio = primary_mapped / total_mapped if total_mapped > 0 else 0
    print_progress(f"Total mapped reads: {total_mapped}, Primary alignments: {primary_mapped}, Ratio: {ratio:.4f}", log)
    
    return total_mapped, primary_mapped, ratio

def analyze_primary_alignments(primary_bam_file, SPECIES, output_dir, settings, reference_to_species, genome_sizes, log=None, report_progress=False):
    print_progress(f"Analyzing primary alignments for {settings['name']} settings...", log)
    
    results = {sp: {'Primary': 0} for sp in SPECIES}
    
    cmd = f"samtools idxstats {primary_bam_file}"
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in result.stdout.strip().split('\n'):
            ref, length, mapped, unmapped = line.split()
            if ref in reference_to_species:
                species = reference_to_species[ref]
                results[species]['Primary'] += int(mapped)
    except subprocess.CalledProcessError as e:
        print_progress(f"Error executing samtools: {e}", log)
        print_progress(f"Error output: {e.stderr}", log)
        raise

    total_reads = sum(sp_data['Primary'] for sp_data in results.values())
    print_progress(f"Analysis completed. Total primary alignments: {total_reads:,}", log)
    
    return results, total_reads

def calculate_proportions_and_ci(results, SPECIES, genome_sizes, total_reads):
    proportions = {}
    category_total = sum(results[sp]['Primary'] for sp in SPECIES)
    
    species_counts = [results[sp]['Primary'] for sp in SPECIES]
    ci_lower, ci_upper = calculate_multinomial_ci(species_counts)
    
    for i, sp in enumerate(SPECIES):
        primary_reads = results[sp]['Primary']
        raw_proportion = primary_reads / category_total if category_total > 0 else 0
        normalized_proportion = (primary_reads / genome_sizes[sp]) / sum((results[s]['Primary'] / genome_sizes[s]) for s in SPECIES)
        
        proportions[sp] = {
            'Primary_Reads': primary_reads,
            'Raw_Proportion': raw_proportion,
            'Normalized_Proportion': normalized_proportion,
            'CI_Lower': ci_lower[i],
            'CI_Upper': ci_upper[i]
        }
    
    return proportions, category_total

def save_primary_alignment_results(results, total_reads, total_mapped_reads, primary_alignment_ratios, output_dir, genome_sizes, log=None):
    print_progress("Saving primary alignment results...", log)
    output_file = os.path.join(output_dir, "primary_alignment_results.tsv")
    
    with open(output_file, 'w') as f:
        f.write("Mapping_Setting\tSpecies\tPrimary_Reads\tRaw_Proportion\tNormalized_Proportion\tCI_Lower\tCI_Upper\tGenome_Size(Mb)\n")
        
        for setting, data in results.items():
            proportions, _ = calculate_proportions_and_ci(data, SPECIES, genome_sizes, total_reads[setting])
            for species, values in proportions.items():
                f.write(f"{setting}\t{species}\t{values['Primary_Reads']}\t{values['Raw_Proportion']:.6f}\t")
                f.write(f"{values['Normalized_Proportion']:.6f}\t{values['CI_Lower']:.6f}\t{values['CI_Upper']:.6f}\t")
                f.write(f"{genome_sizes.get(species, 'N/A') / 1e6:.2f}\n")
        
        # Add overall mapping statistics
        f.write("\nOverall Mapping Statistics\n")
        f.write("Mapping_Setting\tTotal_Mapped_Reads\tPrimary_Alignments\tPrimary_Alignment_Ratio\n")
        for setting in results.keys():
            f.write(f"{setting}\t{total_mapped_reads[setting]}\t{total_reads[setting]}\t{primary_alignment_ratios[setting]:.6f}\n")
    
    print_progress(f"Primary alignment results saved to {output_file}", log)

def plot_primary_alignment_results(results, output_dir, genome_sizes, total_reads, log=None):
    print_progress("Plotting primary alignment results...", log)
    fig, axes = plt.subplots(len(results), 1, figsize=(12, 6*len(results)), squeeze=False)
    
    for idx, (setting, data) in enumerate(results.items()):
        proportions, _ = calculate_proportions_and_ci(data, SPECIES, genome_sizes, total_reads[setting])
        df = pd.DataFrame(proportions).T
        df = df.sort_values('Raw_Proportion', ascending=False)
        
        ax = axes[idx, 0]
        df['Raw_Proportion'].plot(kind='bar', ax=ax, yerr=[df['Raw_Proportion'] - df['CI_Lower'], df['CI_Upper'] - df['Raw_Proportion']])
        ax.set_title(f'{setting.capitalize()} Primary Alignment Proportions')
        ax.set_ylabel('Proportion')
        ax.set_ylim(0, 1)
        
        for i, (_, values) in enumerate(df.iterrows()):
            ax.text(i, values['Raw_Proportion'], f'{values["Raw_Proportion"]:.2%}', ha='center', va='bottom')
        
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'primary_alignment_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print_progress("Primary alignment results plot saved.", log)

def main():
    parser = argparse.ArgumentParser(description="Multi-species Genomic Analysis Tool")
    parser.add_argument('-i', '--skip-index', action='store_true', help="Skip BWA indexing step")
    parser.add_argument('-m', '--skip-mapping', action='store_true', help="Skip BWA mapping step (implies -i)")
    parser.add_argument('--force-recalculate', action='store_true', help="Force recalculation of genome sizes")
    parser.add_argument('-r', '--report-progress', action='store_true', help="Report progress every 1 million reads")
    args = parser.parse_args()

    if args.skip_mapping:
        args.skip_index = True

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    log_file = os.path.join(OUTPUT_DIR, "analysis_log.txt")
    
    with open(log_file, 'a') as log:
        try:
            print_progress("Starting analysis...", log)
            
            cpu_count = get_optimal_thread_count()
            total_cores = multiprocessing.cpu_count()
            print_progress(f"Detected {total_cores} total CPU cores. Using {cpu_count} threads for optimal performance.", log)
            
            genome_sizes = get_genome_sizes(SPECIES, FASTA_FILES, OUTPUT_DIR, args.force_recalculate)

            print_progress("Genome sizes:", log)
            for sp, size in genome_sizes.items():
                print_progress(f"{sp}: {size:,} bp", log)
            print_progress("", log)
            
            combined_ref = os.path.join(OUTPUT_DIR, "combined_reference.fa")
            
            if not args.skip_index:
                combined_ref = create_combined_reference(list(FASTA_FILES.values()), OUTPUT_DIR, log)
                generate_bwa_index(combined_ref, log)
            else:
                print_progress("Skipping BWA indexing step...", log)
                if not os.path.exists(combined_ref):
                    raise FileNotFoundError(f"Combined reference file not found: {combined_ref}")
                if not os.path.exists(combined_ref + ".bwt"):
                    raise FileNotFoundError(f"BWA index files not found for: {combined_ref}")
            
            mapping_file = os.path.join(OUTPUT_DIR, "reference_to_species_mapping.tsv")
            if os.path.exists(mapping_file):
                print_progress("Loading existing reference to species mapping...", log)
                reference_to_species = load_reference_to_species_mapping(mapping_file)
            else:
                print_progress("Creating reference to species mapping from FASTA files...", log)
                reference_to_species = create_reference_to_species_mapping_from_fasta(FASTA_FILES, SPECIES, log)
                save_reference_to_species_mapping(reference_to_species, mapping_file)

            results = {}
            total_reads = {}
            total_mapped_reads = {}
            primary_alignment_ratios = {}

            if not args.skip_mapping:
                for settings in MAPPING_SETTINGS:
                    print_progress(f"\nProcessing {settings['name']} mapping...", log)
                    all_bam_file = run_bwa_mem(combined_ref, UNKNOWN_SAMPLE_R1, UNKNOWN_SAMPLE_R2, OUTPUT_DIR, settings, cpu_count, log)
                    
                    total_mapped, primary_mapped, ratio = calculate_primary_alignment_ratio(all_bam_file, cpu_count, log)
                    total_mapped_reads[settings['name']] = total_mapped
                    primary_alignment_ratios[settings['name']] = ratio
                    
                    analysis_results, primary_count = analyze_primary_alignments(
                        all_bam_file, SPECIES, OUTPUT_DIR, settings, reference_to_species, genome_sizes, log, args.report_progress
                    )
                    results[settings['name']] = analysis_results
                    total_reads[settings['name']] = primary_count
            else:
                print_progress("Skipping BWA mapping step...", log)
                for settings in MAPPING_SETTINGS:
                    all_bam_file = os.path.join(OUTPUT_DIR, f"mapped_{settings['name']}.bam")
                    if not os.path.exists(all_bam_file):
                        raise FileNotFoundError(f"BAM file not found: {all_bam_file}")
                    
                    total_mapped, primary_mapped, ratio = calculate_primary_alignment_ratio(all_bam_file, cpu_count, log)
                    total_mapped_reads[settings['name']] = total_mapped
                    primary_alignment_ratios[settings['name']] = ratio
                    
                    analysis_results, primary_count = analyze_primary_alignments(
                        all_bam_file, SPECIES, OUTPUT_DIR, settings, reference_to_species, genome_sizes, log, args.report_progress
                    )
                    results[settings['name']] = analysis_results
                    total_reads[settings['name']] = primary_count

            save_primary_alignment_results(results, total_reads, total_mapped_reads, primary_alignment_ratios, OUTPUT_DIR, genome_sizes, log)
            plot_primary_alignment_results(results, OUTPUT_DIR, genome_sizes, total_reads, log)
            
            print_progress(f"Analysis complete. Results are in {OUTPUT_DIR}", log)

        except Exception as e:
            print_progress(f"An error occurred: {str(e)}", log)
            print_progress("Traceback:", log)
            traceback.print_exc(file=log)
            print_progress("Script execution failed. Check the log file for details.", log)
            sys.exit(1)

if __name__ == "__main__":
    main()