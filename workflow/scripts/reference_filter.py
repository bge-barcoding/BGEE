#!/usr/bin/env python3
"""
Reference Filter
Removes sequences that are outliers compared to reference sequences.

1. Finds corresponding reference sequences for each alignment file
2. Compares each sequence to the reference (not to a consensus generated from the alignment)
3. Calculates deviation scores for each sequence compared to the reference
4. Uses statistical thresholds to identify sequences that are outliers relative to the reference
5. Removes outlier sequences that deviate too much from the reference pattern
6. Outputs cleaned alignments with sequences that are similar to the reference
"""

import os
import sys
import csv
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Dict, Optional
import gc
from datetime import datetime
from collections import Counter

def log_message(message: str, log_file=None, stdout=False):
    """Log message to file and optionally stdout"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {message}"
    
    if log_file:
        log_file.write(formatted_msg + '\n')
        log_file.flush()
    
    if stdout:
        print(formatted_msg, flush=True)
        
def get_base_filename(filepath: str) -> str:
    """Extract base filename without extension and filter suffixes"""
    basename = os.path.basename(filepath)
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Remove MGE-specific patterns first
    import re
    basename = re.sub(r'_r_\d+_s_\d+_align_.*', '', basename)
    
    # Remove filter suffixes in order  
    for suffix in ['_outlier_filtered', '_at_filtered', '_human_filtered', '_align']:
        basename = basename.replace(suffix, '')
    
    return basename

def get_reference_sequence(reference_file: str) -> str:
    """Extract reference sequence from a FASTA file"""
    try:
        references = list(SeqIO.parse(reference_file, "fasta"))
        if not references:
            raise ValueError("No sequences found in reference file")
        if len(references) > 1:
            print(f"Warning: Multiple sequences found in {reference_file}, using first one")
        return str(references[0].seq)
    except Exception as e:
        raise ValueError(f"Error reading reference file: {str(e)}")

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """Calculate residue frequencies at each position"""
    if not sequences:
        return []
    
    max_len = max(len(seq) for seq in sequences)
    position_freqs = []
    
    for i in range(max_len):
        residues = []
        for seq in sequences:
            if i < len(seq) and seq[i] != '-':
                residues.append(seq[i])
        
        if residues:
            counter = Counter(residues)
            total = len(residues)
            frequencies = {res: count/total for res, count in counter.items()}
        else:
            frequencies = {}
            
        position_freqs.append(frequencies)
    
    return position_freqs

def generate_consensus_sequence(sequences: List[str], threshold: float = 0.5) -> Tuple[str, List[Dict[str, float]]]:
    """Generate consensus sequence and position frequencies"""
    if not sequences:
        return "", []
    
    frequencies = calculate_position_frequencies(sequences)
    consensus = []
    
    for pos_freqs in frequencies:
        if pos_freqs:
            most_common = max(pos_freqs.items(), key=lambda x: x[1])
            if most_common[1] >= threshold:
                consensus.append(most_common[0])
            else:
                consensus.append('-')
        else:
            consensus.append('-')
    
    return ''.join(consensus), frequencies

def calculate_unweighted_deviation(sequence: str, reference: str) -> float:
    """Calculate unweighted deviation score"""
    if len(sequence) != len(reference):
        min_len = min(len(sequence), len(reference))
        sequence = sequence[:min_len]
        reference = reference[:min_len]
    
    differences = 0
    valid_positions = 0
    
    for seq_res, ref_res in zip(sequence, reference):
        if seq_res != '-' and ref_res != '-':
            valid_positions += 1
            if seq_res != ref_res:
                differences += 1
    
    return differences / valid_positions if valid_positions > 0 else 0.0

def calculate_weighted_deviation(sequence: str, reference: str, frequencies: List[Dict[str, float]]) -> float:
    """Calculate weighted deviation score based on conservation"""
    if len(sequence) != len(reference):
        min_len = min(len(sequence), len(reference))
        sequence = sequence[:min_len]
        reference = reference[:min_len]
    
    if len(frequencies) < len(sequence):
        frequencies.extend([{}] * (len(sequence) - len(frequencies)))
    
    total_score = 0.0
    total_weight = 0.0
    
    for i, (seq_res, ref_res) in enumerate(zip(sequence, reference)):
        if seq_res != '-' and ref_res != '-' and i < len(frequencies):
            conservation_weight = frequencies[i].get(ref_res, 0)
            total_weight += conservation_weight
            
            if seq_res != ref_res:
                total_score += conservation_weight
    
    return total_score / total_weight if total_weight > 0 else 0.0


def process_single_file(file_path: str, reference_dir: str, outlier_percentile: float, 
                       consensus_threshold: float, output_dir: str) -> Dict:
    """Process a single FASTA file for reference-based filtering"""
    try:
        base_name = get_base_filename(file_path)
        
        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'empty_file',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Find reference file
        reference_file = os.path.join(reference_dir, f"{base_name}_reference.fasta")
        if not os.path.exists(reference_file):
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'no_reference_file',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Read reference sequence
        try:
            reference_seq = get_reference_sequence(reference_file)
        except Exception as e:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'error',
                'reason': f'reference_read_error: {str(e)}',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        # Read input sequences
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
        except Exception as e:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'error',
                'reason': f'parse_error: {str(e)}',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        if not records:
            return {
                'file_path': file_path,
                'base_name': base_name,
                'status': 'skipped',
                'reason': 'no_sequences',
                'input_count': 0,
                'kept_count': 0,
                'removed_count': 0,
                'removed_sequences': []
            }
        
        input_count = len(records)
        
        # Generate consensus sequence from input for frequency calculation
        sequences = [str(record.seq).upper() for record in records]
        consensus_seq, frequencies = generate_consensus_sequence(sequences, consensus_threshold)
        
        # Pad reference sequence to match alignment length
        max_len = max(len(seq) for seq in sequences) if sequences else len(reference_seq)
        padded_reference = reference_seq.ljust(max_len, '-').upper()
        
        # Calculate deviation scores for all sequences compared to reference
        deviation_scores = []
        for record in records:
            sequence = str(record.seq).upper()
            
            # Pad sequence to match reference length
            padded_sequence = sequence.ljust(len(padded_reference), '-')
            
            unweighted_dev = calculate_unweighted_deviation(padded_sequence, padded_reference)
            weighted_dev = calculate_weighted_deviation(padded_sequence, padded_reference, frequencies)
            
            deviation_scores.append({
                'record': record,
                'unweighted_deviation': unweighted_dev,
                'weighted_deviation': weighted_dev
            })
        
        # Calculate outlier thresholds
        unweighted_scores = [s['unweighted_deviation'] for s in deviation_scores if s['unweighted_deviation'] > 0]
        weighted_scores = [s['weighted_deviation'] for s in deviation_scores if s['weighted_deviation'] > 0]
        
        unweighted_threshold = np.percentile(unweighted_scores, outlier_percentile) if unweighted_scores else float('inf')
        weighted_threshold = np.percentile(weighted_scores, outlier_percentile) if weighted_scores else float('inf')
        
        # Filter sequences
        kept_records = []
        removed_sequences = []
        
        for score_info in deviation_scores:
            record = score_info['record']
            unweighted_dev = score_info['unweighted_deviation']
            weighted_dev = score_info['weighted_deviation']
            
            # Check if sequence is an outlier
            is_outlier = (unweighted_dev > unweighted_threshold or weighted_dev > weighted_threshold)
            
            if is_outlier:
                removal_reason = 'reference_outlier'
                if unweighted_dev > unweighted_threshold and weighted_dev > weighted_threshold:
                    removal_reason = 'reference_outlier_both'
                elif unweighted_dev > unweighted_threshold:
                    removal_reason = 'reference_outlier_unweighted'
                elif weighted_dev > weighted_threshold:
                    removal_reason = 'reference_outlier_weighted'
                
                removed_sequences.append({
                    'sequence_id': record.id,
                    'removal_reason': removal_reason,
                    'unweighted_deviation': unweighted_dev,
                    'weighted_deviation': weighted_dev,
                    'unweighted_threshold': unweighted_threshold,
                    'weighted_threshold': weighted_threshold
                })
            else:
                kept_records.append(record)
        
        # Write filtered sequences
        if kept_records:
            output_file = os.path.join(output_dir, f"{base_name}_reference_filtered.fasta")
            with open(output_file, 'w') as handle:
                SeqIO.write(kept_records, handle, "fasta")
        
        return {
            'file_path': file_path,
            'base_name': base_name,
            'status': 'success',
            'reason': 'processed',
            'input_count': input_count,
            'kept_count': len(kept_records),
            'removed_count': len(removed_sequences),
            'removed_sequences': removed_sequences,
            'output_file': os.path.join(output_dir, f"{base_name}_reference_filtered.fasta") if kept_records else None,
            'reference_file': reference_file,
            'thresholds': {
                'unweighted_threshold': unweighted_threshold,
                'weighted_threshold': weighted_threshold,
                'outlier_percentile': outlier_percentile
            }
        }
        
    except Exception as e:
        return {
            'file_path': file_path,
            'base_name': get_base_filename(file_path),
            'status': 'error',
            'reason': f'processing_error: {str(e)}',
            'input_count': 0,
            'kept_count': 0,
            'removed_count': 0,
            'removed_sequences': []
        }

def main():
    parser = argparse.ArgumentParser(description='Filter sequences based on reference comparison')
    parser.add_argument('--input-files-list', required=True, help='File containing list of FASTA files to process')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered files')
    parser.add_argument('--filtered-files-list', required=True, help='Output file listing successfully filtered files')
    parser.add_argument('--metrics-csv', required=True, help='Output CSV file with filtering metrics')
    parser.add_argument('--reference-dir', required=True, help='Directory containing reference sequences')
    parser.add_argument('--outlier-percentile', type=float, default=90.0, help='Percentile threshold for outlier detection')
    parser.add_argument('--consensus-threshold', type=float, default=0.5, help='Consensus generation threshold')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel processing')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate reference directory
    if not os.path.exists(args.reference_dir):
        print(f"Error: Reference directory does not exist: {args.reference_dir}")
        sys.exit(1)
    
    # Read input files
    with open(args.input_files_list, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]
    
    print(f"Processing {len(input_files)} files with reference comparison (outlier percentile: {args.outlier_percentile})")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        future_to_file = {
            executor.submit(process_single_file, file_path, args.reference_dir, 
                          args.outlier_percentile, args.consensus_threshold, args.output_dir): file_path
            for file_path in input_files
        }
        
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'success':
                    thresholds = result.get('thresholds', {})
                    print(f"✓ {result['base_name']}: {result['kept_count']}/{result['input_count']} sequences kept")
                    print(f"  Reference: {os.path.basename(result.get('reference_file', 'N/A'))}")
                    print(f"  Thresholds: unweighted={thresholds.get('unweighted_threshold', 'N/A'):.4f}, weighted={thresholds.get('weighted_threshold', 'N/A'):.4f}")
                elif result['status'] == 'skipped':
                    print(f"⚠ {result['base_name']}: skipped ({result['reason']})")
                else:
                    print(f"✗ {result['base_name']}: error ({result['reason']})")
                    
            except Exception as e:
                print(f"✗ {file_path}: processing failed - {str(e)}")
                results.append({
                    'file_path': file_path,
                    'base_name': os.path.basename(file_path),
                    'status': 'error',
                    'reason': f'executor_error: {str(e)}',
                    'input_count': 0,
                    'kept_count': 0,
                    'removed_count': 0,
                    'removed_sequences': []
                })
    
    # Write filtered files list
    successful_files = []
    with open(args.filtered_files_list, 'w') as f:
        for result in results:
            if result['status'] == 'success' and result.get('output_file') and os.path.exists(result['output_file']):
                f.write(result['output_file'] + '\n')
                successful_files.append(result['output_file'])
    
    # Write detailed metrics CSV
    with open(args.metrics_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'file_path', 'base_name', 'sequence_id', 'removal_reason', 
            'unweighted_deviation', 'weighted_deviation', 'unweighted_threshold', 'weighted_threshold',
            'reference_file', 'step_name', 'input_count', 'kept_count', 'removed_count'
        ])
        
        for result in results:
            # Write file-level summary
            thresholds = result.get('thresholds', {})
            writer.writerow([
                result['file_path'],
                result['base_name'],
                'FILE_SUMMARY',
                result['reason'],
                '',
                '',
                thresholds.get('unweighted_threshold', ''),
                thresholds.get('weighted_threshold', ''),
                result.get('reference_file', ''),
                'reference_filter',
                result['input_count'],
                result['kept_count'],
                result['removed_count']
            ])
            
            # Write individual removed sequences
            for seq_info in result['removed_sequences']:
                writer.writerow([
                    result['file_path'],
                    result['base_name'],
                    seq_info['sequence_id'],
                    seq_info['removal_reason'],
                    f"{seq_info['unweighted_deviation']:.6f}",
                    f"{seq_info['weighted_deviation']:.6f}",
                    f"{seq_info['unweighted_threshold']:.6f}",
                    f"{seq_info['weighted_threshold']:.6f}",
                    result.get('reference_file', ''),
                    'reference_filter',
                    '',
                    '',
                    ''
                ])
    
    # Summary statistics
    total_input = sum(r['input_count'] for r in results)
    total_kept = sum(r['kept_count'] for r in results)
    total_removed = sum(r['removed_count'] for r in results)
    successful_count = len([r for r in results if r['status'] == 'success'])
    skipped_no_ref = len([r for r in results if r['status'] == 'skipped' and r['reason'] == 'no_reference_file'])
    
    print(f"\nReference Filtering Summary:")
    print(f"Files processed successfully: {successful_count}/{len(input_files)}")
    print(f"Files skipped (no reference): {skipped_no_ref}")
    print(f"Total input sequences: {total_input}")
    print(f"Sequences kept: {total_kept}")
    print(f"Sequences removed (reference outliers): {total_removed}")
    print(f"Filtered files written: {len(successful_files)}")

if __name__ == "__main__":
    main()
