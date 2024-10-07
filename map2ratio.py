#!/usr/bin/env python3

import argparse
import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import multinomial
import multiprocessing
import sys
import time
import traceback
from Bio import SeqIO
import warnings

# Analysis Configuration Constants
SPECIES = ['beef', 'chicken', 'goat', 'horse', 'pork', 'sheep']
FASTA_FILES = {sp: f"{sp}.fa" for sp in SPECIES}
MAPPING_SETTINGS = [
    {"name": "default", "params": ""},
    {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3"},
    {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4"}
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

def get_genome_sizes(SPECIES, fasta_files, output_dir, force_recalculate=False, log=None):
    genome_sizes_file = os.path.join(output_dir, "genome_sizes.txt")
    
    if not force_recalculate and os.path.exists(genome_sizes_file):
        try:
            genome_sizes = {}
            with open(genome_sizes_file, 'r') as f:
                for line in f:
                    sp, size = line.strip().split('\t')
                    genome_sizes[sp] = int(size)
            return genome_sizes
        except IOError as e:
            print_progress(f"Error reading genome sizes file: {str(e)}. Recalculating genome sizes.", log)
    
    genome_sizes = {}
    for sp in SPECIES:
        if sp not in fasta_files:
            raise ValueError(f"No FASTA file specified for species: {sp}")
        if not os.path.exists(fasta_files[sp]):
            raise FileNotFoundError(f"FASTA file not found for species {sp}: {fasta_files[sp]}")
        try:
            total_length = 0
            with open(fasta_files[sp], 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    total_length += len(record.seq)
            genome_sizes[sp] = total_length
        except Exception as e:
            raise RuntimeError(f"Error calculating genome size for species {sp}: {str(e)}")
    
    try:
        with open(genome_sizes_file, 'w') as f:
            for sp, size in genome_sizes.items():
                f.write(f"{sp}\t{size}\n")
    except IOError as e:
        print_progress(f"Warning: Unable to save genome sizes file: {str(e)}", log)
    
    return genome_sizes

def get_total_reads(fastq_file, cache_dir="cache"):
    os.makedirs(cache_dir, exist_ok=True)
    
    file_stats = os.stat(fastq_file)
    cache_file = os.path.join(cache_dir, f"{os.path.basename(fastq_file)}.{file_stats.st_mtime}.txt")
    
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            total_reads = int(f.read().strip())
        print(f"Using cached total read count for {fastq_file}")
        return total_reads
    
    cpu_count = multiprocessing.cpu_count()
    cmd_total_reads = f"seqkit stats -j {cpu_count} -T {fastq_file} | awk 'NR==2 {{print $4}}'"
    total_reads = int(subprocess.check_output(cmd_total_reads, shell=True))
    
    with open(cache_file, 'w') as f:
        f.write(str(total_reads))
    
    print(f"Calculated and cached total read count for {fastq_file}")
    return total_reads

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
    last_update_time = start_time
    while process.poll() is None:
        time.sleep(30)
        current_size = os.path.getsize(combined_ref)
        current_time = time.time()
        
        if current_size != previous_size:
            elapsed_time = current_time - start_time
            print_progress(f"Indexing in progress... Elapsed time: {elapsed_time:.0f} seconds. Reference file size: {current_size/1e6:.2f} MB", log)
            previous_size = current_size
            last_update_time = current_time
        elif current_time - last_update_time > 600:
            elapsed_time = current_time - start_time
            print_progress(f"Indexing still in progress... Elapsed time: {elapsed_time:.0f} seconds. No file size change in the last 10 minutes, but process is still running.", log)
            last_update_time = current_time
    
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
               f"samtools sort -@ {cpu_count} - > {all_bam_file}")
    run_command(command, log)
    run_command(f"samtools index {all_bam_file}", log)
    
    # Log mapped read counts
    mapped_count = int(subprocess.check_output(f"samtools view -c -F 4 {all_bam_file}", shell=True))
    print_progress(f"BWA MEM completed for {settings['name']} settings. Mapped reads: {mapped_count}", log)
    return all_bam_file

def calculate_alignment_statistics(bam_file, total_read_pairs, log=None):
    print_progress(f"Calculating alignment statistics for {bam_file}...", log)
    
    try:
        # Count properly paired reads (only primary alignments, exclude secondary, supplementary, and unmapped)
        cmd_mapped = f"samtools view -c -f 3 -F 2316 {bam_file}"
        mapped_read_pairs = max(0, int(subprocess.check_output(cmd_mapped, shell=True)) // 2)
        
        # Calculate unmapped read pairs
        unmapped_read_pairs = max(0, total_read_pairs - mapped_read_pairs)
        
        # Calculate mapping ratio
        mapping_ratio = mapped_read_pairs / total_read_pairs if total_read_pairs > 0 else 0
        
        print_progress(f"Total read pairs: {total_read_pairs}", log)
        print_progress(f"Mapped read pairs (primary alignments): {mapped_read_pairs}", log)
        print_progress(f"Unmapped read pairs: {unmapped_read_pairs}", log)
        print_progress(f"Mapping ratio: {mapping_ratio:.4f}", log)
        
        if mapped_read_pairs > total_read_pairs:
            print_progress(f"Warning: Number of mapped read pairs ({mapped_read_pairs}) exceeds total read pairs count ({total_read_pairs})", log)
        
        return mapped_read_pairs, unmapped_read_pairs, mapping_ratio
    
    except subprocess.CalledProcessError as e:
        print_progress(f"Error executing samtools command: {e}", log)
        raise
    except ValueError as e:
        print_progress(f"Error parsing samtools output: {e}", log)
        raise

def analyze_alignments(bam_file, SPECIES, reference_to_species, genome_sizes, log=None):
    print_progress(f"Analyzing alignments for {bam_file}...", log)
    
    results = {sp: {'Mapped': 0} for sp in SPECIES}
    
    try:
        # Analyze properly paired primary alignments by species (exclude secondary, supplementary, and unmapped)
        cmd_idxstats = f"samtools view -f 3 -F 2316 {bam_file} | cut -f 3 | sort | uniq -c"
        result = subprocess.run(cmd_idxstats, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        for line in result.stdout.strip().split('\n'):
            count, ref = line.strip().split()
            count = int(count) // 2  # Convert to read pairs
            if ref in reference_to_species:
                species = reference_to_species[ref]
                results[species]['Mapped'] += count
        
        total_mapped = sum(results[sp]['Mapped'] for sp in SPECIES)
        
        for sp in SPECIES:
            results[sp]['Raw_Proportion'] = results[sp]['Mapped'] / total_mapped if total_mapped > 0 else 0
            results[sp]['Normalized_Proportion'] = (results[sp]['Mapped'] / genome_sizes[sp]) / sum((results[s]['Mapped'] / genome_sizes[s]) for s in SPECIES) if total_mapped > 0 else 0
        
        print_progress(f"Total properly paired primary alignments: {total_mapped}", log)
        
        return results, total_mapped
    
    except subprocess.CalledProcessError as e:
        print_progress(f"Error executing samtools command: {e}", log)
        print_progress(f"Error output: {e.stderr}", log)
        raise
    except ValueError as e:
        print_progress(f"Error parsing samtools output: {e}", log)
        raise

def calculate_multinomial_ci(counts, num_bootstraps=10000, confidence_level=0.95):
    counts = np.array(counts)
    if not np.all(counts >= 0) or not np.issubdtype(counts.dtype, np.integer):
        print_progress(f"Warning: Invalid counts detected: {counts}. Setting negative values to 0 and rounding to integers.")
        counts = np.round(np.maximum(counts, 0)).astype(int)
    
    total = sum(counts)
    if total == 0:
        warnings.warn("All counts are zero, returning zero confidence intervals")
        return np.zeros_like(counts), np.zeros_like(counts)
    
    probs = counts / total
    
    bootstrapped_props = multinomial.rvs(n=total, p=probs, size=num_bootstraps) / total
    
    ci_lower = np.percentile(bootstrapped_props, (1 - confidence_level) / 2 * 100, axis=0)
    ci_upper = np.percentile(bootstrapped_props, (1 + confidence_level) / 2 * 100, axis=0)
    
    # Ensure CIs are within [0, 1] and lower <= upper
    ci_lower = np.maximum(0, ci_lower)
    ci_upper = np.minimum(1, ci_upper)
    ci_lower = np.minimum(ci_lower, probs)
    ci_upper = np.maximum(ci_upper, probs)
    
    return ci_lower, ci_upper

def save_primary_alignment_results(results, total_read_pairs, primary_mapped_read_pairs, unmapped_read_pairs, mapping_ratio, output_dir, genome_sizes, log=None):
    print_progress("Saving primary alignment results...", log)
    output_file = os.path.join(output_dir, "detailed_results.tsv")
    
    with open(output_file, 'w') as f:
        f.write("Mapping_Setting\tSpecies\tMapped_Read_Pairs\tRaw_Proportion\tRaw_CI_Lower\tRaw_CI_Upper\tNormalized_Proportion\tNormalized_CI_Lower\tNormalized_CI_Upper\tGenome_Size(Mb)\n")
        
        for setting, data in results.items():
            counts = [data[species]['Mapped'] for species in SPECIES]
            raw_ci_lower, raw_ci_upper = calculate_multinomial_ci(counts)
            
            normalized_counts = [int(data[species]['Mapped'] / genome_sizes[species] * 1e6) for species in SPECIES]
            norm_ci_lower, norm_ci_upper = calculate_multinomial_ci(normalized_counts)
            
            for i, species in enumerate(SPECIES):
                mapped_read_pairs = data[species]['Mapped']
                raw_proportion = data[species]['Raw_Proportion']
                normalized_proportion = data[species]['Normalized_Proportion']
                f.write(f"{setting}\t{species}\t{mapped_read_pairs}\t{raw_proportion:.6f}\t{raw_ci_lower[i]:.6f}\t{raw_ci_upper[i]:.6f}\t{normalized_proportion:.6f}\t{norm_ci_lower[i]:.6f}\t{norm_ci_upper[i]:.6f}\t{genome_sizes[species]/1e6:.2f}\n")
        
        # Add overall mapping statistics
        f.write("\nOverall Mapping Statistics\n")
        f.write("Mapping_Setting\tTotal_Read_Pairs\tPrimary_Mapped_Read_Pairs\tUnmapped_Read_Pairs\tMapping_Ratio\n")
        for setting in results.keys():
            f.write(f"{setting}\t{total_read_pairs}\t{primary_mapped_read_pairs[setting]}\t{unmapped_read_pairs[setting]}\t{mapping_ratio[setting]:.6f}\n")
    
    print_progress(f"Primary alignment results saved to {output_file}", log)

def plot_primary_alignment_results(results, output_dir, genome_sizes, log=None):
    print_progress("Plotting primary alignment results...", log)
    fig, axes = plt.subplots(len(results), 2, figsize=(20, 10*len(results)), squeeze=False)
    
    for idx, (setting, data) in enumerate(results.items()):
        df = pd.DataFrame({
            'Species': SPECIES,
            'Mapped_Read_Pairs': [data[sp]['Mapped'] for sp in SPECIES],
            'Raw_Proportion': [data[sp]['Raw_Proportion'] for sp in SPECIES],
            'Normalized_Proportion': [data[sp]['Normalized_Proportion'] for sp in SPECIES]
        })
        df = df.sort_values('Raw_Proportion', ascending=False)
        
        counts = df['Mapped_Read_Pairs'].values
        raw_ci_lower, raw_ci_upper = calculate_multinomial_ci(counts)
        normalized_counts = np.round((df['Mapped_Read_Pairs'] / [genome_sizes[sp] for sp in df['Species']]).values * 1e6).astype(int)
        norm_ci_lower, norm_ci_upper = calculate_multinomial_ci(normalized_counts)
        
        # Raw proportions plot
        ax = axes[idx, 0]
        yerr_lower = df['Raw_Proportion'] - raw_ci_lower
        yerr_upper = raw_ci_upper - df['Raw_Proportion']
        yerr = np.array([yerr_lower, yerr_upper])
        df.plot(kind='bar', x='Species', y='Raw_Proportion', ax=ax, yerr=yerr, capsize=5)
        ax.set_title(f'{setting.capitalize()} Raw Alignment Proportions')
        ax.set_ylabel('Proportion of Mapped Read Pairs')
        ax.set_ylim(0, 1)
        
        for i, row in df.iterrows():
            ax.text(i, row['Raw_Proportion'], f'{row["Raw_Proportion"]:.2%}', ha='center', va='bottom')
        
        # Normalized proportions plot
        ax = axes[idx, 1]
        yerr_lower = df['Normalized_Proportion'] - norm_ci_lower
        yerr_upper = norm_ci_upper - df['Normalized_Proportion']
        yerr = np.array([yerr_lower, yerr_upper])
        df.plot(kind='bar', x='Species', y='Normalized_Proportion', ax=ax, yerr=yerr, capsize=5)
        ax.set_title(f'{setting.capitalize()} Normalized Alignment Proportions')
        ax.set_ylabel('Normalized Proportion of Mapped Read Pairs')
        ax.set_ylim(0, 1)
        
        for i, row in df.iterrows():
            ax.text(i, row['Normalized_Proportion'], f'{row["Normalized_Proportion"]:.2%}', ha='center', va='bottom')
        
        for ax in axes[idx]:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'alignment_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print_progress("Alignment results plot saved.", log)

def main():
    parser = argparse.ArgumentParser(description="Multi-species Genomic Analysis Tool")
    parser.add_argument('-i', '--skip-index', action='store_true', help="Skip BWA indexing step")
    parser.add_argument('-m', '--skip-mapping', action='store_true', help="Skip BWA mapping step (implies -i)")
    parser.add_argument('--force-recalculate', action='store_true', help="Force recalculation of genome sizes")
    args = parser.parse_args()

    if args.skip_mapping:
        args.skip_index = True

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    log_file = os.path.join(OUTPUT_DIR, "analysis_log.txt")
    
    with open(log_file, 'a') as log:
        try:
            print_progress("Starting analysis...", log)
            
            cpu_count = max(1, multiprocessing.cpu_count() // 2)
            print_progress(f"Detected {multiprocessing.cpu_count()} CPU cores. Using {cpu_count} cores for processing.", log)
            
            genome_sizes = get_genome_sizes(SPECIES, FASTA_FILES, OUTPUT_DIR, args.force_recalculate, log)

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
            
            reference_to_species = {record.id: sp for sp, file in FASTA_FILES.items() for record in SeqIO.parse(file, 'fasta')}

            results = {}
            primary_mapped_read_pairs = {}
            unmapped_read_pairs = {}
            mapping_ratio = {}

            # Calculate total read pairs from FASTQ files
            total_reads = get_total_reads(UNKNOWN_SAMPLE_R1)
            total_read_pairs = total_reads // 2

            if not args.skip_mapping:
                for settings in MAPPING_SETTINGS:
                    print_progress(f"\nProcessing {settings['name']} mapping...", log)
                    all_bam_file = run_bwa_mem(combined_ref, UNKNOWN_SAMPLE_R1, UNKNOWN_SAMPLE_R2, OUTPUT_DIR, settings, cpu_count, log)
                    
                    mapped_read_pairs, unmapped, ratio = calculate_alignment_statistics(all_bam_file, total_read_pairs, log)
                    
                    try:
                        analysis_results, _ = analyze_alignments(
                            all_bam_file, SPECIES, reference_to_species, genome_sizes, log
                        )
                        results[settings['name']] = analysis_results
                        primary_mapped_read_pairs[settings['name']] = mapped_read_pairs
                        unmapped_read_pairs[settings['name']] = unmapped
                        mapping_ratio[settings['name']] = ratio
                    except Exception as e:
                        print_progress(f"Error in analyze_alignments for {settings['name']}: {str(e)}", log)
                        raise
            else:
                print_progress("Skipping BWA mapping step...", log)
                for settings in MAPPING_SETTINGS:
                    all_bam_file = os.path.join(OUTPUT_DIR, f"mapped_{settings['name']}.bam")
                    if not os.path.exists(all_bam_file):
                        raise FileNotFoundError(f"BAM file not found: {all_bam_file}")
                    
                    mapped_read_pairs, unmapped, ratio = calculate_alignment_statistics(all_bam_file, total_read_pairs, log)
                    
                    try:
                        analysis_results, _ = analyze_alignments(
                            all_bam_file, SPECIES, reference_to_species, genome_sizes, log
                        )
                        results[settings['name']] = analysis_results
                        primary_mapped_read_pairs[settings['name']] = mapped_read_pairs
                        unmapped_read_pairs[settings['name']] = unmapped
                        mapping_ratio[settings['name']] = ratio
                    except Exception as e:
                        print_progress(f"Error in analyze_alignments for {settings['name']}: {str(e)}", log)
                        raise

            save_primary_alignment_results(results, total_read_pairs, primary_mapped_read_pairs, unmapped_read_pairs, mapping_ratio, OUTPUT_DIR, genome_sizes, log)
            plot_primary_alignment_results(results, OUTPUT_DIR, genome_sizes, log)
            
            print_progress(f"Analysis complete. Results are in {OUTPUT_DIR}", log)

        except Exception as e:
            print_progress(f"An error occurred: {str(e)}", log)
            print_progress("Traceback:", log)
            traceback.print_exc(file=log)
            print_progress("Script execution failed. Check the log file for details.", log)
            sys.exit(1)

if __name__ == "__main__":
    main()