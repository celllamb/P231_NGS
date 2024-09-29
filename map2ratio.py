#!/usr/bin/python3
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

def print_progress(message, file=None):
    output = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}"
    print(output)
    if file:
        print(output, file=file)
    sys.stdout.flush()

def get_cpu_count():
    return multiprocessing.cpu_count()

def run_command(command, log=None):
    print_progress(f"Running command: {command}", log)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    if process.returncode != 0:
        print_progress(f"Error executing command: {command}", log)
        print_progress(error.decode(), log)
        return None
    return output.decode()

def create_combined_reference(references, output_dir, log=None):
    print_progress("Creating combined reference...", log)
    combined_ref = os.path.join(output_dir, "combined_reference.fa")
    with open(combined_ref, 'w') as outfile:
        for ref in references:
            with open(ref, 'r') as infile:
                outfile.write(infile.read())
    print_progress("Combined reference created.", log)
    return combined_ref

def run_bwa_mem(ref_file, read1, read2, output_dir, settings, cpu_count, log=None):
    print_progress(f"Running BWA MEM with {settings['name']} settings...", log)
    bam_file = os.path.join(output_dir, f"mapped_{settings['name']}.bam")
    command = f"bwa mem -t {cpu_count} {settings['params']} {ref_file} {read1} {read2} | samtools view -bS - | samtools sort - > {bam_file}"
    run_command(command, log)
    run_command(f"samtools index {bam_file}", log)
    print_progress(f"BWA MEM completed for {settings['name']} settings.", log)
    return bam_file

def bootstrap_confidence_interval(data, num_bootstraps=1000, confidence_level=0.95):
    bootstrapped_means = []
    for _ in range(num_bootstraps):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrapped_means.append(np.mean(sample))
    ci_lower, ci_upper = np.percentile(bootstrapped_means, 
                                       [(1-confidence_level)/2 * 100, 
                                        (1+confidence_level)/2 * 100])
    return ci_lower, ci_upper

def analyze_bam_file(bam_file, log=None):
    print_progress(f"Analyzing BAM file: {bam_file}", log)
    total_reads = 0
    unmapped_reads = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            total_reads += 1
            if read.is_unmapped:
                unmapped_reads += 1
    
    mapped_reads = total_reads - unmapped_reads
    unmapped_ratio = (unmapped_reads / total_reads) * 100 if total_reads > 0 else 0
    
    print_progress(f"Total reads: {total_reads}", log)
    print_progress(f"Mapped reads: {mapped_reads}", log)
    print_progress(f"Unmapped reads: {unmapped_reads}", log)
    print_progress(f"Unmapped ratio: {unmapped_ratio:.2f}%", log)
    
    return {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'unmapped_reads': unmapped_reads,
        'unmapped_ratio': unmapped_ratio
    }

def analyze_multi_mapped_reads(bam_file, species, output_dir, settings, cpu_count, log=None):
    print_progress(f"Analyzing multi-mapped reads for {settings['name']} settings...", log)
    multi_mapped_file = os.path.join(output_dir, f"multi_mapped_details_{settings['name']}.tsv")
    command = f"samtools view -@ {cpu_count} -F 4 {bam_file} | awk '{{if ($5 < 30) print $0}}' > {multi_mapped_file}"
    run_command(command, log)
    
    multi_map_counts = {sp: {other_sp: 0 for other_sp in species} for sp in species}
    single_species_multi_map = {sp: 0 for sp in species}
    unique_map_counts = {sp: 0 for sp in species}
    
    read_mappings = {}
    
    with open(multi_mapped_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            read_name = fields[0]
            mapped_species = fields[2].split('.')[0]
            
            if read_name not in read_mappings:
                read_mappings[read_name] = set()
            read_mappings[read_name].add(mapped_species)
    
    for read_name, mapped_species_set in read_mappings.items():
        if len(mapped_species_set) == 1:
            species_name = list(mapped_species_set)[0]
            single_species_multi_map[species_name] += 1
        else:
            for sp1 in mapped_species_set:
                for sp2 in mapped_species_set:
                    if sp1 != sp2:
                        multi_map_counts[sp1][sp2] += 1
    
    # Count unique mappings
    unique_command = f"samtools view -@ {cpu_count} -F 4 -q 30 {bam_file} | cut -f 3 | cut -d '.' -f 1 | sort | uniq -c"
    unique_output = run_command(unique_command, log)
    for line in unique_output.strip().split('\n'):
        count, species_name = line.strip().split()
        unique_map_counts[species_name] = int(count)
    
    # Calculate total reads and proportions
    total_reads = sum(unique_map_counts.values()) + sum(single_species_multi_map.values())
    for cross_maps in multi_map_counts.values():
        total_reads += sum(cross_maps.values()) // 2  # Divide by 2 to avoid double counting
    
    proportions = {}
    confidence_intervals = {}
    for sp in species:
        unique_prop = unique_map_counts[sp] / total_reads * 100
        single_multi_prop = single_species_multi_map[sp] / total_reads * 100
        cross_multi_prop = sum(multi_map_counts[sp].values()) / total_reads * 100
        total_prop = unique_prop + single_multi_prop + cross_multi_prop
        
        proportions[sp] = {
            'Unique': unique_prop,
            'Single-species Multi': single_multi_prop,
            'Cross-species Multi': cross_multi_prop,
            'Total': total_prop
        }
        
        # Bootstrap for confidence interval
        sp_data = [sp] * (unique_map_counts[sp] + single_species_multi_map[sp]) + \
                  [other_sp for other_sp in species for _ in range(multi_map_counts[sp][other_sp])]
        ci_lower, ci_upper = bootstrap_confidence_interval(sp_data)
        confidence_intervals[sp] = (ci_lower * 100, ci_upper * 100)  # Convert to percentage
    
    # Save results
    pd.DataFrame(multi_map_counts).fillna(0).to_csv(os.path.join(output_dir, f"cross_species_multi_map_{settings['name']}.csv"))
    pd.DataFrame(single_species_multi_map, index=['Count']).T.to_csv(os.path.join(output_dir, f"single_species_multi_map_{settings['name']}.csv"))
    pd.DataFrame(proportions).T.to_csv(os.path.join(output_dir, f"mapping_proportions_{settings['name']}.csv"))
    pd.DataFrame(confidence_intervals, index=['CI_Lower', 'CI_Upper']).T.to_csv(os.path.join(output_dir, f"confidence_intervals_{settings['name']}.csv"))
    
    print_progress(f"Multi-mapped read analysis completed for {settings['name']} settings.", log)
    return proportions, confidence_intervals

def plot_enhanced_results(results_dict, ci_dict, unmapped_dict, output_dir, log=None):
    print_progress("Plotting enhanced results...", log)
    fig, axes = plt.subplots(len(results_dict), 1, figsize=(12, 6*len(results_dict)), squeeze=False)
    
    for idx, (setting, data) in enumerate(results_dict.items()):
        df = pd.DataFrame(data).T
        df.plot(kind='bar', stacked=True, ax=axes[idx, 0])
        axes[idx, 0].set_title(f'{setting.capitalize()} Mapping Proportions')
        axes[idx, 0].set_ylabel('Percentage')
        axes[idx, 0].legend(title='Mapping Type')
        
        # Add confidence intervals
        ci_data = ci_dict[setting]
        x = range(len(ci_data))
        for i, (species, (ci_lower, ci_upper)) in enumerate(ci_data.items()):
            axes[idx, 0].plot([i, i], [ci_lower, ci_upper], color='black', linewidth=2)
        
        # Add unmapped ratio text
        unmapped_ratio = unmapped_dict[setting]['unmapped_ratio']
        axes[idx, 0].text(0.95, 0.95, f"Unmapped: {unmapped_ratio:.2f}%", 
                          transform=axes[idx, 0].transAxes, ha='right', va='top')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'enhanced_mapping_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print_progress("Enhanced results plot saved.", log)

def main():
    # Configuration
    OUTPUT_DIR = "output"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Set up log file
    log_file = os.path.join(OUTPUT_DIR, "analysis_log.txt")
    
    with open(log_file, 'a') as log:
        print_progress("Starting analysis...", log)
        
        # Detect CPU cores
        cpu_count = get_cpu_count()
        print_progress(f"Detected {cpu_count} CPU cores. Will use all available cores.", log)

        UNKNOWN_SAMPLE_R1 = "SM2_trim_1.fastq.gz"
        UNKNOWN_SAMPLE_R2 = "SM2_trim_2.fastq.gz"
        REFERENCES = ["beef.fa", "chicken.fa", "goat.fa", "horse.fa", "pork.fa", "sheep.fa"]
        SPECIES = ["beef", "chicken", "goat", "horse", "pork", "sheep"]
        
        MAPPING_SETTINGS = [
            {"name": "default", "params": "-a"},
            {"name": "strict", "params": "-B 12 -O 18,18 -E 3,3 -a"},
            {"name": "very_strict", "params": "-B 16 -O 24,24 -E 4,4 -a"}
        ]
        
        # Create combined reference
        combined_ref = create_combined_reference(REFERENCES, OUTPUT_DIR, log)
        
        # Generate BWA index
        print_progress("Generating BWA index...", log)
        run_command(f"bwa index {combined_ref}", log)
        
        results = {}
        ci_results = {}
        unmapped_results = {}
        for settings in MAPPING_SETTINGS:
            print_progress(f"\nStarting analysis for {settings['name']} mapping...", log)
            bam_file = run_bwa_mem(combined_ref, UNKNOWN_SAMPLE_R1, UNKNOWN_SAMPLE_R2, OUTPUT_DIR, settings, cpu_count, log)
            
            print_progress(f"Analyzing BAM file for {settings['name']} mapping...", log)
            unmapped_results[settings['name']] = analyze_bam_file(bam_file, log)
            
            print_progress(f"Performing multi-mapped analysis for {settings['name']} mapping...", log)
            results[settings['name']], ci_results[settings['name']] = analyze_multi_mapped_reads(bam_file, SPECIES, OUTPUT_DIR, settings, cpu_count, log)
            
            print_progress(f"Results for {settings['name']} mapping:", log)
            print(pd.DataFrame(results[settings['name']]))
            print(pd.DataFrame(results[settings['name']]), file=log)
            print_progress("\nConfidence Intervals:", log)
            print(pd.DataFrame(ci_results[settings['name']], index=['CI_Lower', 'CI_Upper']).T)
            print(pd.DataFrame(ci_results[settings['name']], index=['CI_Lower', 'CI_Upper']).T, file=log)
        
        print_progress("Plotting results...", log)
        plot_enhanced_results(results, ci_results, unmapped_results, OUTPUT_DIR, log)
        
        print_progress(f"Analysis complete. Results are in {OUTPUT_DIR}", log)

    print(f"\nAnalysis complete. Results and log are in {OUTPUT_DIR}")

if __name__ == "__main__":
    main()