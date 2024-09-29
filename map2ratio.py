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

def create_combined_reference(references, output_dir, log=None):
    print_progress("Creating combined reference...", log)
    combined_ref = os.path.join(output_dir, "combined_reference.fa")
    with open(combined_ref, 'w') as outfile:
        for ref in references:
            with open(ref, 'r') as infile:
                outfile.write(infile.read())
    print_progress("Combined reference created.", log)
    return combined_ref

def generate_bwa_index(combined_ref, cpu_count, log=None):
    print_progress("Generating BWA index...", log)
    command = f"bwa index -a bwtsw -t {cpu_count} {combined_ref}"
    run_command(command, log)
    print_progress("BWA index generated.", log)

def run_bwa_mem(ref_file, read1, read2, output_dir, settings, cpu_count, log=None):
    print_progress(f"Running BWA MEM with {settings['name']} settings...", log)
    bam_file = os.path.join(output_dir, f"mapped_{settings['name']}.bam")
    command = f"bwa mem -t {cpu_count} {settings['params']} {ref_file} {read1} {read2} | samtools view -bS - | samtools sort - > {bam_file}"
    run_command(command, log)
    run_command(f"samtools index {bam_file}", log)
    print_progress(f"BWA MEM completed for {settings['name']} settings.", log)
    return bam_file

def analyze_multi_mapped_reads(bam_file, species, output_dir, settings, cpu_count, log=None):
    print_progress(f"Analyzing multi-mapped reads for {settings['name']} settings...", log)
    multi_mapped_file = os.path.join(output_dir, f"multi_mapped_details_{settings['name']}.tsv")
    
    # Use FLAG field to identify multi-mapped reads
    # FLAG & 256 : not primary alignment
    # FLAG & 2048 : supplementary alignment
    command = f"samtools view -@ {cpu_count} -F 4 {bam_file} | awk '{{if (and($2, 256) || and($2, 2048)) print $0}}' > {multi_mapped_file}"
    run_command(command, log)
    
    multi_map_counts = {sp: {other_sp: 0 for other_sp in species} for sp in species}
    single_species_multi_map = {sp: 0 for sp in species}
    unique_map_counts = {sp: 0 for sp in species}
    
    read_mappings = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            
            read_name = read.query_name
            mapped_species = read.reference_name.split('.')[0]
            
            if read_name not in read_mappings:
                read_mappings[read_name] = set()
            read_mappings[read_name].add(mapped_species)
            
            if read.is_secondary or read.is_supplementary:
                if len(read_mappings[read_name]) == 1:
                    single_species_multi_map[mapped_species] += 1
                else:
                    for sp in read_mappings[read_name]:
                        if sp != mapped_species:
                            multi_map_counts[mapped_species][sp] += 1
            elif len(read_mappings[read_name]) == 1:
                unique_map_counts[mapped_species] += 1
    
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

def bootstrap_confidence_interval(data, num_bootstraps=1000, confidence_level=0.95):
    bootstrapped_means = []
    for _ in range(num_bootstraps):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrapped_means.append(np.mean(sample))
    ci_lower, ci_upper = np.percentile(bootstrapped_means, 
                                       [(1-confidence_level)/2 * 100, 
                                        (1+confidence_level)/2 * 100])
    return ci_lower, ci_upper

def plot_enhanced_results(results_dict, ci_dict, output_dir, log=None):
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
        
        # Add total percentage text
        for i, (species, values) in enumerate(data.items()):
            total = values['Total']
            axes[idx, 0].text(i, total, f'{total:.2f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'enhanced_mapping_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print_progress("Enhanced results plot saved.", log)

def main():
    parser = argparse.ArgumentParser(description="Multi-species Genomic Analysis Tool")
    parser.add_argument('-i', '--skip-index', action='store_true', help="Skip BWA indexing step")
    parser.add_argument('-m', '--skip-mapping', action='store_true', help="Skip BWA mapping step (implies -i)")
    args = parser.parse_args()

    # If -m is used, automatically set -i as well
    if args.skip_mapping:
        args.skip_index = True

    OUTPUT_DIR = "output"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    log_file = os.path.join(OUTPUT_DIR, "analysis_log.txt")
    
    with open(log_file, 'a') as log:
        try:
            print_progress("Starting analysis...", log)
            
            cpu_count = multiprocessing.cpu_count()
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
            
            combined_ref = os.path.join(OUTPUT_DIR, "combined_reference.fa")
            
            if not args.skip_index:
                combined_ref = create_combined_reference(REFERENCES, OUTPUT_DIR, log)
                generate_bwa_index(combined_ref, cpu_count, log)
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
            
            results = {}
            ci_results = {}
            for settings in MAPPING_SETTINGS:
                print_progress(f"\nAnalyzing results for {settings['name']} mapping...", log)
                results[settings['name']], ci_results[settings['name']] = analyze_multi_mapped_reads(bam_files[settings['name']], SPECIES, OUTPUT_DIR, settings, cpu_count, log)
                
                print_progress(f"Results for {settings['name']} mapping:", log)
                print(pd.DataFrame(results[settings['name']]))
                print(pd.DataFrame(results[settings['name']]), file=log)
                print_progress("\nConfidence Intervals:", log)
                print(pd.DataFrame(ci_results[settings['name']], index=['CI_Lower', 'CI_Upper']).T)
                print(pd.DataFrame(ci_results[settings['name']], index=['CI_Lower', 'CI_Upper']).T, file=log)
            
            print_progress("Plotting results...", log)
            plot_enhanced_results(results, ci_results, OUTPUT_DIR, log)
            
            print_progress(f"Analysis complete. Results are in {OUTPUT_DIR}", log)

        except Exception as e:
            print_progress(f"An error occurred: {str(e)}", log)
            print_progress("Traceback:", log)
            traceback.print_exc(file=log)
            print_progress("Script execution failed. Check the log file for details.", log)
            sys.exit(1)

if __name__ == "__main__":
    main()