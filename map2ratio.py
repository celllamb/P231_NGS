#!/usr/bin/env python3

import argparse
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import multiprocessing
import pysam
import sys
import time
import traceback
from Bio import SeqIO
from collections import defaultdict

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

def get_genome_sizes(species, fasta_files, output_dir, force_recalculate=False):
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
                    if sp not in species:
                        raise ValueError(f"Unknown species in genome sizes file: {sp}")
                    try:
                        genome_sizes[sp] = int(size)
                    except ValueError:
                        raise ValueError(f"Invalid genome size for species {sp}: {size}")
            return genome_sizes
        except (IOError, ValueError) as e:
            print_progress(f"Error reading genome sizes file: {str(e)}. Recalculating genome sizes.")
    
    genome_sizes = {}
    for sp in species:
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

def calculate_confidence_intervals(data, num_bootstraps=10000, confidence_level=0.95):
    bootstrapped_proportions = []
    n = len(data)
    for _ in range(num_bootstraps):
        bootstrap_sample = np.random.choice(data, size=n, replace=True)
        bootstrapped_proportions.append(np.mean(bootstrap_sample))
    
    ci_lower = np.percentile(bootstrapped_proportions, (1 - confidence_level) / 2 * 100)
    ci_upper = np.percentile(bootstrapped_proportions, (1 + confidence_level) / 2 * 100)
    
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

def create_reference_to_species_mapping_from_fasta(fasta_files, species, log=None):
    reference_to_species = {}
    for sp, fasta_file in zip(species, fasta_files.values()):
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
        time.sleep(1)  # Check every second
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
    bam_file = os.path.join(output_dir, f"mapped_{settings['name']}.bam")
    command = f"bwa mem -t {cpu_count} {settings['params']} {ref_file} {read1} {read2} | samtools view -bS - | samtools sort - > {bam_file}"
    run_command(command, log)
    run_command(f"samtools index {bam_file}", log)
    print_progress(f"BWA MEM completed for {settings['name']} settings.", log)
    return bam_file

def analyze_multi_mapped_reads(bam_file, species, output_dir, settings, cpu_count, reference_to_species, genome_sizes, log=None, report_progress=False):
    print_progress(f"Analyzing multi-mapped reads for {settings['name']} settings...", log)
    
    results = {sp: {'Primary': 0, 'Secondary': 0, 'Supplementary': 0, 'Total_Multi_Within': 0, 'Unique_Multi_Within': 0, 'Total_Multi_Across': 0, 'Unique_Multi_Across': 0} for sp in species}
    read_mappings = defaultdict(lambda: defaultdict(lambda: {'primary': 0, 'secondary': 0, 'supplementary': 0}))
    
    total_reads = 0
    mapped_reads = 0
    
    start_time = time.time()
    last_progress_time = start_time

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            total_reads += 1
            
            if read.is_unmapped:
                continue
            
            mapped_reads += 1
            read_name = read.query_name
            mapped_species = reference_to_species[read.reference_name]
            
            if read.is_secondary:
                results[mapped_species]['Secondary'] += 1
                read_mappings[read_name][mapped_species]['secondary'] += 1
            elif read.is_supplementary:
                results[mapped_species]['Supplementary'] += 1
                read_mappings[read_name][mapped_species]['supplementary'] += 1
            else:  # Primary alignment
                results[mapped_species]['Primary'] += 1
                read_mappings[read_name][mapped_species]['primary'] += 1
            
            if report_progress and (mapped_reads % 10000000 == 0 or time.time() - last_progress_time > 300):
                current_time = time.time()
                print_progress(f"Processed {mapped_reads:,} mapped reads out of {total_reads:,} total reads. Elapsed time: {current_time - start_time:.2f} seconds", log)
                last_progress_time = current_time

    # Calculate Multi_Within and Multi_Across
    for sp in species:
        results[sp]['Total_Multi_Within'] = results[sp]['Secondary'] + results[sp]['Supplementary']
        results[sp]['Unique_Multi_Within'] = sum(1 for read in read_mappings.values() if read[sp]['secondary'] + read[sp]['supplementary'] > 0)
        results[sp]['Total_Multi_Across'] = sum(read[sp]['primary'] + read[sp]['secondary'] + read[sp]['supplementary'] for read in read_mappings.values() if len(read) > 1)
        results[sp]['Unique_Multi_Across'] = sum(1 for read in read_mappings.values() if sp in read and len(read) > 1)

    unmapped_reads = total_reads - mapped_reads
    print_progress(f"Analysis completed. Total reads: {total_reads:,}, Mapped reads: {mapped_reads:,}, Unmapped reads: {unmapped_reads:,}", log)
    print_progress(f"Total processing time: {time.time() - start_time:.2f} seconds", log)
    
    return results, total_reads, mapped_reads, unmapped_reads

def calculate_species_proportions(results, species, genome_sizes, total_reads):
    category_totals = {
        'Primary': sum(results[sp]['Primary'] for sp in species),
        'Secondary': sum(results[sp]['Secondary'] for sp in species),
        'Supplementary': sum(results[sp]['Supplementary'] for sp in species),
        'Unique_Multi_Within': sum(results[sp]['Unique_Multi_Within'] for sp in species),
        'Unique_Multi_Across': sum(results[sp]['Unique_Multi_Across'] for sp in species)
    }
    category_totals['Total'] = category_totals['Primary'] + category_totals['Secondary'] + category_totals['Supplementary']
    
    proportions = {}
    for sp in species:
        total_for_species = results[sp]['Primary'] + results[sp]['Secondary'] + results[sp]['Supplementary']
        proportions[sp] = {
            'Primary_Reads': results[sp]['Primary'],
            'Secondary_Reads': results[sp]['Secondary'],
            'Supplementary_Reads': results[sp]['Supplementary'],
            'Unique_Multi_Within': results[sp]['Unique_Multi_Within'],
            'Unique_Multi_Across': results[sp]['Unique_Multi_Across'],
            'Total_Reads': total_for_species,
            'Raw_Primary_Proportion': results[sp]['Primary'] / category_totals['Primary'],
            'Raw_Secondary_Proportion': results[sp]['Secondary'] / category_totals['Secondary'] if category_totals['Secondary'] > 0 else 0,
            'Raw_Supplementary_Proportion': results[sp]['Supplementary'] / category_totals['Supplementary'] if category_totals['Supplementary'] > 0 else 0,
            'Raw_Unique_Multi_Within_Proportion': results[sp]['Unique_Multi_Within'] / category_totals['Unique_Multi_Within'] if category_totals['Unique_Multi_Within'] > 0 else 0,
            'Raw_Unique_Multi_Across_Proportion': results[sp]['Unique_Multi_Across'] / category_totals['Unique_Multi_Across'] if category_totals['Unique_Multi_Across'] > 0 else 0,
            'Raw_Total_Proportion': total_for_species / category_totals['Total']
        }
    
    # Calculate Normalized Proportion
    total_normalized = sum((results[s]['Primary'] + results[s]['Secondary'] + results[s]['Supplementary']) / genome_sizes[s] for s in species)
    for sp in species:
        proportions[sp]['Normalized_Proportion'] = ((results[sp]['Primary'] + results[sp]['Secondary'] + results[sp]['Supplementary']) / genome_sizes[sp]) / total_normalized

    # Add totals
    proportions['Total'] = {
        'Primary_Reads': category_totals['Primary'],
        'Secondary_Reads': category_totals['Secondary'],
        'Supplementary_Reads': category_totals['Supplementary'],
        'Unique_Multi_Within': category_totals['Unique_Multi_Within'],
        'Unique_Multi_Across': category_totals['Unique_Multi_Across'],
        'Total_Reads': category_totals['Total'],
        'Raw_Primary_Proportion': 1.0,
        'Raw_Secondary_Proportion': 1.0,
        'Raw_Supplementary_Proportion': 1.0,
        'Raw_Unique_Multi_Within_Proportion': 1.0,
        'Raw_Unique_Multi_Across_Proportion': 1.0,
        'Raw_Total_Proportion': 1.0,
        'Normalized_Proportion': 1.0
    }

    return proportions, category_totals
	
def save_detailed_results(results, total_reads, mapped_reads, unmapped_reads, output_dir, genome_sizes, log=None):
    print_progress("Saving detailed results...", log)
    output_file = os.path.join(output_dir, "detailed_results.tsv")
    
    with open(output_file, 'w') as f:
        f.write("Mapping_Setting\tSpecies\tPrimary_Reads\tSecondary_Reads\tSupplementary_Reads\t")
        f.write("Unique_Multi_Within\tUnique_Multi_Across\tTotal_Reads\t")
        f.write("Raw_Primary_Proportion\tRaw_Secondary_Proportion\tRaw_Supplementary_Proportion\t")
        f.write("Raw_Unique_Multi_Within_Proportion\tRaw_Unique_Multi_Across_Proportion\tRaw_Total_Proportion\t")
        f.write("Normalized_Proportion\tGenome_Size(Mb)\n")
        
        for setting in results.keys():
            for species, data in results[setting].items():
                f.write(f"{setting}\t{species}\t{data['Primary_Reads']}\t{data['Secondary_Reads']}\t{data['Supplementary_Reads']}\t")
                f.write(f"{data['Unique_Multi_Within']}\t{data['Unique_Multi_Across']}\t{data['Total_Reads']}\t")
                f.write(f"{data['Raw_Primary_Proportion']:.6f}\t{data['Raw_Secondary_Proportion']:.6f}\t{data['Raw_Supplementary_Proportion']:.6f}\t")
                f.write(f"{data['Raw_Unique_Multi_Within_Proportion']:.6f}\t{data['Raw_Unique_Multi_Across_Proportion']:.6f}\t{data['Raw_Total_Proportion']:.6f}\t")
                f.write(f"{data['Normalized_Proportion']:.6f}\t")
                f.write(f"{genome_sizes.get(species, 'N/A') / 1e6:.2f}\n")
        
        # Add overall mapping statistics
        f.write("\nOverall Mapping Statistics\n")
        f.write("Mapping_Setting\tTotal_Reads\tMapped_Reads\tUnmapped_Reads\tMapping_Rate\n")
        for setting in results.keys():
            total = total_reads[setting]
            mapped = mapped_reads[setting]
            unmapped = unmapped_reads[setting]
            mapping_rate = mapped / total if total > 0 else 0
            f.write(f"{setting}\t{total}\t{mapped}\t{unmapped}\t{mapping_rate:.6f}\n")
    
    print_progress(f"Detailed results saved to {output_file}", log)

def plot_enhanced_results(results, output_dir, log=None):
    print_progress("Plotting enhanced results...", log)
    fig, axes = plt.subplots(len(results), 1, figsize=(12, 6*len(results)), squeeze=False)
    
    for idx, (setting, data) in enumerate(results.items()):
        df = pd.DataFrame(data).T
        df = df.sort_values('Raw_Total_Proportion', ascending=False)
        
        ax = axes[idx, 0]
        df[['Raw_Primary_Proportion', 'Raw_Secondary_Proportion', 'Raw_Supplementary_Proportion']].plot(kind='bar', stacked=True, ax=ax)
        ax.set_title(f'{setting.capitalize()} Mapping Proportions')
        ax.set_ylabel('Proportion')
        ax.set_ylim(0, 1)
        ax.legend(title='Mapping Type')
        
        # Add total percentage text
        for i, (species, values) in enumerate(df.iterrows()):
            total = values['Raw_Total_Proportion']
            ax.text(i, total, f'{total:.2%}', ha='center', va='bottom')
        
        # Rotate x-axis labels for better readability
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'enhanced_mapping_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print_progress("Enhanced results plot saved.", log)

def main():
    parser = argparse.ArgumentParser(description="Multi-species Genomic Analysis Tool")
    parser.add_argument('-i', '--skip-index', action='store_true', help="Skip BWA indexing step")
    parser.add_argument('-m', '--skip-mapping', action='store_true', help="Skip BWA mapping step (implies -i)")
    parser.add_argument('--force-recalculate', action='store_true', help="Force recalculation of genome sizes")
    parser.add_argument('-r', '--report-progress', action='store_true', help="Report progress every 10 million reads")
    args = parser.parse_args()

    if args.skip_mapping:
        args.skip_index = True

    OUTPUT_DIR = "output"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    log_file = os.path.join(OUTPUT_DIR, "analysis_log.txt")
    
    with open(log_file, 'a') as log:
        try:
            print_progress("Starting analysis...", log)
            
            cpu_count = get_optimal_thread_count()
            total_cores = multiprocessing.cpu_count()
            print_progress(f"Detected {total_cores} total CPU cores. Using {cpu_count} threads for optimal performance.", log)

            UNKNOWN_SAMPLE_R1 = "SM2_trim_1.fastq.gz"
            UNKNOWN_SAMPLE_R2 = "SM2_trim_2.fastq.gz"
            SPECIES = ['beef', 'chicken', 'goat', 'horse', 'pork', 'sheep']
            FASTA_FILES = {sp: f"{sp}.fa" for sp in SPECIES}
            
            MAPPING_SETTINGS = [
                {"name": "default", "params": "-a"},
                {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3 -a"},
                {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4 -a"}
            ]
            
            # Get genome sizes
            genome_sizes = get_genome_sizes(SPECIES, FASTA_FILES, OUTPUT_DIR, args.force_recalculate)

            # Print genome sizes
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
            
            if not args.skip_mapping:
                bam_files = {}
                for settings in MAPPING_SETTINGS:
                    bam_files[settings['name']] = run_bwa_mem(combined_ref, UNKNOWN_SAMPLE_R1, UNKNOWN_SAMPLE_R2, OUTPUT_DIR, settings, cpu_count, log)
            else:
                print_progress("Skipping BWA mapping step...", log)
                bam_files = {settings['name']: os.path.join(OUTPUT_DIR, f"mapped_{settings['name']}.bam") for settings in MAPPING_SETTINGS}
                for bam_file in bam_files.values():
                    if not os.path.exists(bam_file):
                        raise FileNotFoundError(f"BAM file not found: {bam_file}")
            
            # Create or load reference_to_species mapping
            mapping_file = os.path.join(OUTPUT_DIR, "reference_to_species_mapping.tsv")
            if os.path.exists(mapping_file) and os.path.getsize(mapping_file) > 0:
                print_progress("Loading existing reference to species mapping...", log)
                reference_to_species = load_reference_to_species_mapping(mapping_file)
            else:
                print_progress("Creating reference to species mapping from FASTA files...", log)
                reference_to_species = create_reference_to_species_mapping_from_fasta(FASTA_FILES, SPECIES, log)
                if reference_to_species:
                    save_reference_to_species_mapping(reference_to_species, mapping_file)
                    print_progress(f"Reference to species mapping saved to {mapping_file}", log)
                else:
                    print_progress("Error: Failed to create reference to species mapping. Exiting program.", log)
                    sys.exit(1)
            
            results = {}
            total_reads = {}
            mapped_reads = {}
            unmapped_reads = {}
            for settings in MAPPING_SETTINGS:
                print_progress(f"\nAnalyzing results for {settings['name']} mapping...", log)
                analysis_results, total, mapped, unmapped = analyze_multi_mapped_reads(
                    bam_files[settings['name']], SPECIES, OUTPUT_DIR, settings, cpu_count, reference_to_species, genome_sizes, log, args.report_progress
                )
                results[settings['name']], _ = calculate_species_proportions(analysis_results, SPECIES, genome_sizes, total)
                total_reads[settings['name']] = total
                mapped_reads[settings['name']] = mapped
                unmapped_reads[settings['name']] = unmapped
            
            print_progress("\nSaving detailed results...", log)
            save_detailed_results(results, total_reads, mapped_reads, unmapped_reads, OUTPUT_DIR, genome_sizes, log)
            
            print_progress("\nPlotting results...", log)
            plot_enhanced_results(results, OUTPUT_DIR, log)
            
            print_progress(f"Analysis complete. Results are in {OUTPUT_DIR}", log)

        except Exception as e:
            print_progress(f"An error occurred: {str(e)}", log)
            print_progress("Traceback:", log)
            traceback.print_exc(file=log)
            print_progress("Script execution failed. Check the log file for details.", log)
            sys.exit(1)

if __name__ == "__main__":
    main()