"""
FASTA Sequence Analysis and Summary Tool
--------------------------------------

This script analyses a log file containing a list of alignment FASTA files, .out files
containing 'raw' summary stats from MGE, and cleaning statistics from a CSV file.
It generates comprehensive summary statistics in CSV format. 

Usage:
    python mge_stats.py -a/--alignment_log <log_file> -o/--output <output_csv> -od/--out_file_dir <out_file_dir> -c/--cleaning_csv <cleaning_csv>

Arguments:
    -a, --alignment_log : str
        Path to a text file containing a list of FASTA file paths (one per line)
    -o, --output : str
        Name of the output CSV file
    -od, --out_file_dir : str
        Directory containing .out files with additional sequence statistics
    -c, --cleaning_csv : str
        Path to CSV file containing cleaning statistics

Outputs:
    - <output_file>.csv: Main summary file containing all statistics
    - mge_stats.log: Log file with processing information

The script generates the following metrics for each sample:
    - ID: Identifier extracted from the filename
    - mge_params: Parameters used for MGE (e.g., r_1.3_s_50)
    - n_reads_in: Number of input sequences (from .out file)
    - n_reads_aligned: Number of aligned sequences in the FASTA file
    - n_reads_skipped: Number of sequences that were successfully aligned but not in the FASTA file
    - ref_length: Length of alignment (from .out file)
    - cov_min: Minimum coverage depth across alignment
    - cov_max: Maximum coverage depth across alignment
    - cov_avg: Average coverage depth across alignment
    - cov_med: Median coverage depth across alignment
    - cleaning_input_reads: Number of input sequences for cleaning
    - cleaning_kept_reads: Number of sequences kept after cleaning
    - cleaning_removed_human: Number of sequences removed due to human similarity
    - cleaning_removed_at: Number of sequences removed due to AT content
    - cleaning_removed_outlier: Number of sequences removed as statistical outliers
    - cleaning_ambig_bases: Number of ambiguous bases after cleaning
    - cleaning_cov_percent: Coverage percentage after cleaning
    - cleaning_cov_avg: Average coverage after cleaning
    - cleaning_cov_max: Maximum coverage after cleaning
    - cleaning_cov_min: Minimum coverage after cleaning
"""



import os
import csv
import argparse
import re
import logging
from Bio import SeqIO
import numpy as np
import glob


# Set up logging
logger = logging.getLogger('mge_stats')
logger.setLevel(logging.INFO)
# Use mode='w' to overwrite the log file each time
file_handler = logging.FileHandler('mge_stats.log', mode='w')
file_handler.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

def extract_sample_info(filename):
    """
    Extract the process ID and parameters from the filename.
    E.g., from "BGSNL096-23_r_1.3_s_50_align_BGSNL096-23.fas"
    Returns:
        - base_id: "BGSNL096-23"
        - full_id: "BGSNL096-23_r_1.3_s_50"
        - params: "r_1.3_s_50"
    """
    filename_no_ext = os.path.splitext(filename)[0]
    
    # Extract base ID (everything before first underscore)
    base_id_match = re.match(r'^([^_]+)', filename_no_ext)
    base_id = base_id_match.group(1) if base_id_match else None
    
    # Extract parameters (r_X_s_Y pattern)
    params_match = re.search(r'(r_[0-9.]+_s_[0-9.]+)', filename_no_ext)
    params = params_match.group(1) if params_match else ""
    
    # Combine to form the full ID used for matching files
    full_id = f"{base_id}_{params}" if params else base_id
    
    return base_id, full_id, params

def extract_process_id(filename):
    """Extract the process ID from the filename (for backward compatibility)."""
    base_id, _, _ = extract_sample_info(filename)
    return base_id

def parse_out_file(out_file_path):
    """Parse the .out file for additional statistics."""
    try:
        with open(out_file_path, 'r') as file:
            content = file.read()
    except Exception as e:
        logger.error(f"Error reading .out file {out_file_path}: {e}")
        return None

    process_id_data = {
        'n_reads_in': None,
        'ref_length': None,
        'successful_aligned': None
    }

    # Extract number of input sequences (updated pattern)
    reads_in_match = re.search(r'Number of input sequences:\s*(\d+)', content)
    if reads_in_match:
        process_id_data['n_reads_in'] = int(reads_in_match.group(1))
    
    # Simplified pattern - just target the line with the number
    aligned_match = re.search(r'to the amino acid sequence, see vulgar file:\s*(\d+)', content)
    if aligned_match:
        process_id_data['successful_aligned'] = int(aligned_match.group(1))
    
    # Reference length pattern
    ref_length_match = re.search(r'Length of alignment:\s*(\d+)', content)
    if ref_length_match:
        process_id_data['ref_length'] = int(ref_length_match.group(1))
 
    return process_id_data

def parse_cleaning_csv(file_path):
    """Parse the cleaning CSV file for filtering statistics."""
    cleaning_data = {}
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Extract sample name which should match with full_id
                sample_name = row.get('sample_name', '')
                if sample_name:
                    # Map the CSV columns to our output structure
                    cleaning_stats = {
                        'input_reads': int(row.get('input_reads', 0)),
                        'removed_human': int(row.get('removed_human', 0)),
                        'removed_at': int(row.get('removed_at_distance', 0)),
                        'removed_outlier': int(row.get('removed_outliers', 0)),
                        'kept_reads': int(row.get('cleaned_reads', 0)),
                        'ambig_bases': int(row.get('final_ambig_bases', 0)),
                        'cov_percent': float(row.get('cov_percent', 0)),
                        'cov_avg': float(row.get('cov_mean', 0)),
                        'cov_max': float(row.get('cov_max', 0)),
                        'cov_min': float(row.get('cov_min', 0))
                    }
                    
                    # Store under original sample name
                    cleaning_data[sample_name] = cleaning_stats
                    
                    # Store under versions without suffixes
                    # Handle _fcleaner suffix
                    if '_fcleaner' in sample_name:
                        clean_name = sample_name.replace('_fcleaner', '')
                        cleaning_data[clean_name] = cleaning_stats
                    
                    # Handle _fcleaner_merge suffix
                    if '_fcleaner_merge' in sample_name:
                        clean_name = sample_name.replace('_fcleaner_merge', '')
                        cleaning_data[clean_name] = cleaning_stats
                        
                        # Also store without _merge suffix in case there's inconsistency
                        clean_name_alt = sample_name.replace('_merge', '')
                        cleaning_data[clean_name_alt] = cleaning_stats
            
        logger.info(f"Parsed cleaning data for {len(cleaning_data)} samples from CSV")
        return cleaning_data
            
    except Exception as e:
        logger.error(f"Error reading cleaning CSV file {file_path}: {e}")
        return {}

def process_fasta_file(file_path):
    """Process a FASTA file and extract sequence statistics."""
    logger.info(f"Processing file: {file_path}")  # Log which file is being processed
    
    if os.path.getsize(file_path) == 0:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        logger.warning(f"Empty file: {file_path}")
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    if not sequences:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        logger.warning(f"No sequences found in file: {file_path}")
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    unique_sequences = {}
    for seq in sequences:
        if seq.id not in unique_sequences:
            unique_sequences[seq.id] = seq
    
    # Detect sequence length mismatches
    sequence_lengths = [len(seq.seq) for seq in unique_sequences.values()]
    if len(set(sequence_lengths)) > 1:
        logger.error(f"Error in file {file_path}: Found sequences with different lengths: {set(sequence_lengths)}")
        # Log some example sequence IDs with their lengths
        for seq_id, seq in list(unique_sequences.items())[:5]:  # Log up to 5 examples
            logger.error(f"  Sequence {seq_id}: length {len(seq.seq)}")
        
        # Use the most common length
        from collections import Counter
        most_common_length = Counter(sequence_lengths).most_common(1)[0][0]
        logger.info(f"Using most common length: {most_common_length} for file {file_path}")
        
        # Filter sequences to only include those with the most common length
        filtered_sequences = {seq_id: seq for seq_id, seq in unique_sequences.items() 
                             if len(seq.seq) == most_common_length}
        
        if not filtered_sequences:
            logger.error(f"No sequences with common length in {file_path}")
            base_id, _, _ = extract_sample_info(os.path.basename(file_path))
            return {
                'ID': base_id,
                'n_reads_aligned': len(unique_sequences),
                'cov_min': 0,
                'cov_max': 0,
                'cov_avg': 0,
                'cov_med': 0,
            }
        
        logger.info(f"Filtered from {len(unique_sequences)} to {len(filtered_sequences)} sequences")
        unique_sequences = filtered_sequences

    sequence_count = len(unique_sequences)
    
    try:
        first_seq = next(iter(unique_sequences.values()))
        coverage = np.zeros(len(first_seq.seq))
        
        for seq_id, seq in unique_sequences.items():
            try:
                seq_array = np.array([1 if base != '-' else 0 for base in seq.seq])
                if coverage.shape != seq_array.shape:
                    logger.error(f"Shape mismatch in {file_path}: Expected {coverage.shape}, got {seq_array.shape} for sequence {seq_id}")
                    # Skip this sequence
                    continue
                coverage += seq_array
            except Exception as e:
                logger.error(f"Error processing sequence {seq_id} in file {file_path}: {e}")
                # Skip this sequence and continue with others
                continue

        # Calculate coverage statistics
        if len(coverage) > 0:
            min_coverage = np.min(coverage)
            max_coverage = np.max(coverage)
            mean_coverage = np.mean(coverage)
            median_coverage = np.median(coverage)
        else:
            min_coverage = max_coverage = mean_coverage = median_coverage = 0
            
    except Exception as e:
        logger.error(f"Error calculating coverage for file {file_path}: {e}")
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': sequence_count,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    base_id, _, _ = extract_sample_info(os.path.basename(file_path))

    return {
        'ID': base_id,
        'n_reads_aligned': sequence_count,
        'cov_min': min_coverage,
        'cov_max': max_coverage,
        'cov_avg': mean_coverage,
        'cov_med': median_coverage,
    }

def summarise_fasta(log_file, output_file, out_file_dir, cleaning_csv=None):
    """Summarise FASTA files based on log file with cleaning statistics from CSV."""
    try:
        with open(log_file, 'r') as f:
            file_paths = [line.strip() for line in f if line.strip().endswith(('.fasta', '.fas'))]
    except Exception as e:
        logger.error(f"Error reading log file {log_file}: {e}")
        return

    # Log file information
    logger.info(f"The input log_file contains paths to {len(file_paths)} files for processing")
    
    if not file_paths:
        logger.warning(f"No valid FASTA files found in log file: {log_file}")
        return

    # Extract and log sample names
    sample_info = [extract_sample_info(os.path.basename(path)) for path in file_paths]
    base_ids = [info[0] for info in sample_info if info[0]]
    
    logger.info(f"Sample base IDs found: {', '.join(set(base_ids))}")
    logger.info(f"Total number of samples being processed: {len(file_paths)}")

    # Get alignment stats from .out files
    out_files = [f for f in os.listdir(out_file_dir) if f.endswith('.out')]
    logger.info(f"The out_file_dir contains {len(out_files)} .out files for processing")
    
    out_file_data = {}
    for file in out_files:
        _, full_id, _ = extract_sample_info(file)
        if full_id:
            out_data = parse_out_file(os.path.join(out_file_dir, file))
            if out_data:
                out_file_data[full_id] = out_data
    
    # Get cleaning stats from CSV file
    cleaning_data = {}
    if cleaning_csv and os.path.exists(cleaning_csv):
        cleaning_data = parse_cleaning_csv(cleaning_csv)
    else:
        logger.warning(f"Cleaning CSV file not provided or does not exist")

    # Define fieldnames with updated column headings
    fieldnames = [
        'Filename', 'ID', 'mge_params', 'n_reads_in', 'n_reads_aligned', 'n_reads_skipped', 'ref_length', 
        'cov_min', 'cov_max', 'cov_avg', 'cov_med',
        'cleaning_input_reads', 'cleaning_kept_reads', 'cleaning_removed_human', 'cleaning_removed_at', 
        'cleaning_removed_outlier', 'cleaning_ambig_bases', 'cleaning_cov_percent', 'cleaning_cov_avg', 
        'cleaning_cov_max', 'cleaning_cov_min'
    ]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file_path in file_paths:
            base_id, full_id, params = extract_sample_info(os.path.basename(file_path))
            if not base_id:
                continue

            result = process_fasta_file(file_path)
            if result:
                result['Filename'] = os.path.basename(file_path).replace('.fasta', '').replace('.fas', '').replace('_align_', '_')
                result['ID'] = base_id
                result['mge_params'] = params

                # Add alignment stats with updated keys - use full_id for matching
                out_data = out_file_data.get(full_id, {})
                result['n_reads_in'] = out_data.get('n_reads_in', '')
                result['ref_length'] = out_data.get('ref_length', '')
                
                # Calculate n_reads_skipped (successful_aligned - n_reads_aligned)
                successful_aligned = out_data.get('successful_aligned', None)
                n_reads_aligned = result.get('n_reads_aligned', 0)
                
                if successful_aligned is not None and n_reads_aligned is not None:
                    result['n_reads_skipped'] = max(0, successful_aligned - n_reads_aligned)
                else:
                    result['n_reads_skipped'] = ''

                # Look for the sample in the cleaning data - try different format variants
                # Try multiple formats including different suffixes
                sample_key = full_id
                alt_key = f"{full_id}_{base_id}"
                fcleaner_key = f"{full_id}_{base_id}_fcleaner"
                fcleaner_merge_key = f"{full_id}_{base_id}_fcleaner_merge"
                
                # Try each key in sequence
                cleaning_stats = None
                for key in [sample_key, alt_key, fcleaner_key, fcleaner_merge_key]:
                    if key in cleaning_data:
                        cleaning_stats = cleaning_data[key]
                        logger.info(f"Found cleaning stats using key: {key}")
                        break
                
                if not cleaning_stats:
                    # Additional debug: list some of the keys actually available in cleaning_data
                    if cleaning_data:
                        sample_keys = list(cleaning_data.keys())[:10]  # Show first 10 keys
                        logger.warning(f"No cleaning stats found for {base_id}. Available keys include: {sample_keys}")
                    else:
                        logger.warning(f"No cleaning stats found for {base_id} and cleaning_data is empty")
                    cleaning_stats = {}
                
                # Add cleaning stats with renamed columns
                result['cleaning_input_reads'] = cleaning_stats.get('input_reads', '')
                result['cleaning_kept_reads'] = cleaning_stats.get('kept_reads', '')
                result['cleaning_removed_human'] = cleaning_stats.get('removed_human', '')
                result['cleaning_removed_at'] = cleaning_stats.get('removed_at', '')
                result['cleaning_removed_outlier'] = cleaning_stats.get('removed_outlier', '')
                result['cleaning_ambig_bases'] = cleaning_stats.get('ambig_bases', '')
                result['cleaning_cov_percent'] = cleaning_stats.get('cov_percent', '')
                result['cleaning_cov_avg'] = cleaning_stats.get('cov_avg', '')
                result['cleaning_cov_max'] = cleaning_stats.get('cov_max', '')
                result['cleaning_cov_min'] = cleaning_stats.get('cov_min', '')

                # Write all results to the CSV file
                writer.writerow(result)

    output_file_abs = os.path.abspath(output_file)
    logger.info(f"CSV summary created: {output_file_abs}")
    print(f"CSV summary created: '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise FASTA file statistics including sequence counts, coverage, and cleaning results.')
    parser.add_argument('-a', '--alignment_log', required=True, type=str, help='The log file containing the paths to FASTA files.')
    parser.add_argument('-o', '--output', required=True, type=str, help='The output CSV file name.')
    parser.add_argument('-od', '--out_file_dir', required=True, type=str, help='The directory containing .out files with additional statistics.')
    parser.add_argument('-c', '--cleaning_csv', type=str, help='Path to CSV file containing cleaning statistics')

    args = parser.parse_args()

    if not os.path.isfile(args.alignment_log):
        parser.error(f"The log file '{args.alignment_log}' does not exist.")
    if not os.path.isdir(args.out_file_dir):
        parser.error(f"The directory '{args.out_file_dir}' does not exist or is not a directory.")
    if args.cleaning_csv and not os.path.isfile(args.cleaning_csv):
        parser.error(f"The cleaning CSV file '{args.cleaning_csv}' does not exist.")

    summarise_fasta(args.alignment_log, args.output, args.out_file_dir, args.cleaning_csv)
