#!/usr/bin/env python3

import sys
import os
import argparse
from subprocess import run, CalledProcessError
import pandas as pd
from tqdm import tqdm
from collections import defaultdict, Counter
from pyfaidx import Fasta

# ================== Step 0: Filter PAF Reads that Span Gaps ==================

def load_gaps(bed_file):
    """
    Load gaps from a BED file.

    Args:
    bed_file (str): Path to the BED file containing gap information.

    Returns:
    dict: A dictionary where keys are chromosome names and values are lists of (start, end) tuples.
    """
    gaps = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[:3]
            if chrom not in gaps:
                gaps[chrom] = []
            gaps[chrom].append((int(start), int(end)))
    return gaps

def is_spanning_gap(read_start, read_end, gaps):
    """
    Check if a read spans any gap.

    Args:
    read_start (int): Start position of the read.
    read_end (int): End position of the read.
    gaps (list): List of (start, end) tuples representing gaps.

    Returns:
    bool: True if the read spans any gap, False otherwise.
    """
    for gap_start, gap_end in gaps:
        if read_start <= gap_start and read_end >= gap_end:
            return True
    return False

def filter_paf_reads(input_paf, gaps, output_paf):
    """
    Filter PAF reads that span gaps and write to output PAF file.

    Args:
    input_paf (str): Path to input PAF file.
    gaps (dict): Dictionary of gaps by chromosome.
    output_paf (str): Path to output filtered PAF file.
    """
    with open(input_paf, 'r') as infile, open(output_paf, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            chrom = fields[5]
            try:
                read_start = int(fields[7])
                read_end = int(fields[8])
            except ValueError:
                continue
            if chrom in gaps and is_spanning_gap(read_start, read_end, gaps[chrom]):
                outfile.write(line)
    print(f"Filtered PAF written to {output_paf}")

# ================== Step 1: Split PAF by Each Gap ==================

def load_gaps_df(gaps_bed):
    """
    Load gaps from BED file into a DataFrame.

    Args:
    gaps_bed (str): Path to gaps BED file.

    Returns:
    DataFrame: Pandas DataFrame with columns ['chrom', 'start', 'end'].
    """
    gaps = pd.read_csv(gaps_bed, sep='\t', header=None, names=['chrom', 'start', 'end'])
    return gaps

def split_paf_by_gap(paf_file, gaps, output_dir, min_reads=1):
    """
    Split PAF file into gap-specific PAF files with a minimum number of reads.

    Args:
    paf_file (str): Path to input PAF file.
    gaps (DataFrame): DataFrame containing gaps.
    output_dir (str): Directory to save split PAF files.
    min_reads (int): Minimum number of reads required to create a gap-specific PAF file.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Convert gaps DataFrame to list of dictionaries for faster access
    gaps_list = gaps.to_dict(orient='records')

    # Create a dictionary to hold gap_id and corresponding lines
    gap_paf_dict = defaultdict(list)

    print("Splitting PAF by gaps...")
    total_lines = sum(1 for _ in open(paf_file))
    with open(paf_file, 'r') as paf, tqdm(total=total_lines, desc="Splitting PAF by Gap") as pbar:
        for line in paf:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                pbar.update(1)
                continue
            ref_name, ref_start, ref_end = cols[5], int(cols[7]), int(cols[8])
            for gap in gaps_list:
                if ref_name == gap['chrom'] and ref_start <= gap['start'] and ref_end >= gap['end']:
                    gap_id = f"{gap['chrom']}_{gap['start']}_{gap['end']}"
                    gap_paf_dict[gap_id].append(line)
            pbar.update(1)

    # Write to individual PAF files if the number of reads meets min_reads
    for gap_id, lines in gap_paf_dict.items():
        if len(lines) >= min_reads:
            with open(os.path.join(output_dir, f"{gap_id}.paf"), 'w') as gap_paf:
                gap_paf.writelines(lines)
    print(f"PAF files split into directory {output_dir} with min_reads={min_reads}")

# ================== Step 2: Extract and Split Reads from FASTQ ==================

def parse_all_paf_files(paf_dir):
    """
    Parse all PAF files to build a dictionary mapping each gap to its reads.

    Args:
    paf_dir (str): Directory containing split PAF files.

    Returns:
    dict: {gap_id: [read_name1, read_name2, ...]}
    """
    gap_reads = defaultdict(list)
    for paf_file in os.listdir(paf_dir):
        if paf_file.endswith(".paf"):
            gap_id = paf_file.replace(".paf", "")
            with open(os.path.join(paf_dir, paf_file), 'r') as f:
                for line in f:
                    read_name = line.split('\t')[0]
                    gap_reads[gap_id].append(read_name)
    return gap_reads

def write_total_reads_list(gap_reads, output_file):
    """
    Write all unique read names to a list file.

    Args:
    gap_reads (dict): {gap_id: [read_name1, read_name2, ...]}
    output_file (str): Path to output read list file.
    """
    all_reads = set()
    for reads in gap_reads.values():
        all_reads.update(reads)
    with open(output_file, 'w') as f:
        f.write("\n".join(all_reads))
    print(f"Total {len(all_reads)} unique reads written to {output_file}")

def extract_reads_with_seqkit(reads_list, fastq_files, output_fastq):
    """
    Use seqkit to extract reads based on a list.

    Args:
    reads_list (str): Path to file containing read names.
    fastq_files (list): List of FASTQ file paths.
    output_fastq (str): Path to output extracted FASTQ file.
    """
    reads_list_quoted = f"'{reads_list}'"
    fastq_files_quoted = " ".join([f"'{f}'" for f in fastq_files])
    cmd = f"seqkit grep -f {reads_list_quoted} {fastq_files_quoted} -o {output_fastq}"
    try:
        run(cmd, shell=True, check=True)
        print(f"Extracted reads written to {output_fastq}")
    except CalledProcessError as e:
        print(f"Error running seqkit: {e}")
        sys.exit(1)

def split_fastq_by_gap(extracted_fastq, gap_reads, output_dir):
    """
    Split the extracted FASTQ into gap-specific FASTQ files.

    Args:
    extracted_fastq (str): Path to extracted FASTQ file.
    gap_reads (dict): {gap_id: [read_name1, read_name2, ...]}
    output_dir (str): Directory to save gap-specific FASTQ files.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Map read names to gaps
    read_to_gaps = defaultdict(list)
    for gap_id, reads in gap_reads.items():
        for read in reads:
            read_to_gaps[read].append(gap_id)

    # Open file handles for each gap FASTQ
    gap_fastq_handles = {gap_id: open(os.path.join(output_dir, f"{gap_id}.fastq"), 'a') for gap_id in gap_reads.keys()}

    print("Splitting extracted FASTQ into gap-specific FASTQ files...")
    # Read and split FASTQ
    with open(extracted_fastq, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            read_name = header[1:].split()[0]  # Remove "@" and any suffix
            if read_name in read_to_gaps:
                for gap_id in read_to_gaps[read_name]:
                    gap_fastq_handles[gap_id].write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    # Close all file handles
    for handle in gap_fastq_handles.values():
        handle.close()
    print(f"FASTQ files split into directory {output_dir}")

def step2_extract_and_split(args):
    """
    Execute Step 2: Extract and split reads from FASTQ files.

    Args:
    args: Parsed command-line arguments.
    """
    # Parse PAF files to get reads per gap
    print("Parsing PAF files to get reads per gap...")
    gap_reads = parse_all_paf_files(args.paf_dir)

    # Write all unique reads to a temporary list
    print("Writing unique reads to temporary list...")
    write_total_reads_list(gap_reads, args.temp_reads_list)

    # Extract reads from FASTQ files using seqkit
    print("Extracting reads from FASTQ files using seqkit...")
    extract_reads_with_seqkit(args.temp_reads_list, args.fastq_files, args.extracted_fastq)

    # Split the extracted FASTQ into gap-specific FASTQ files
    print("Splitting extracted FASTQ into gap-specific FASTQ files...")
    split_fastq_by_gap(args.extracted_fastq, gap_reads, args.output_dir)

    # Optionally, remove the temporary reads list and extracted FASTQ
    if os.path.exists(args.temp_reads_list):
        os.remove(args.temp_reads_list)
        print(f"Temporary reads list {args.temp_reads_list} removed.")
    if os.path.exists(args.extracted_fastq):
        os.remove(args.extracted_fastq)
        print(f"Temporary extracted FASTQ {args.extracted_fastq} removed.")

# ================== Step 3: Generate Consensus Sequences ==================

def parse_fastq(file_path):
    """
    Parse FASTQ file and return list of sequences.

    Args:
    file_path (str): Path to FASTQ file.

    Returns:
    list: List of sequences.
    """
    sequences = []
    with open(file_path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            sequences.append(seq)
    return sequences

def generate_consensus_from_reads(reads):
    """
    Generate consensus sequence from multiple reads using majority voting.

    Args:
    reads (list): List of sequences.

    Returns:
    str: Consensus sequence.
    """
    if len(reads) == 1:
        return reads[0]

    consensus = []
    for i in range(len(reads[0])):
        bases = [read[i] for read in reads if i < len(read)]
        if not bases:
            consensus.append('N')  # Unknown base
            continue
        most_common_base = Counter(bases).most_common(1)[0][0]
        consensus.append(most_common_base)
    return "".join(consensus)

def run_hifiasm_assembly(input_fastq, output_prefix, threads=4):
    """
    Run hifiasm to assemble reads.

    Args:
    input_fastq (str): Path to input FASTQ file.
    output_prefix (str): Prefix for hifiasm output files.
    threads (int): Number of threads.

    Returns:
    str or None: Path to GFA file if successful, else None.
    """
    gfa_file = f"{output_prefix}.bp.p_ctg.gfa"
    cmd = f"hifiasm -o {output_prefix} -t {threads} -f0 {input_fastq}"
    try:
        run(cmd, shell=True, check=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
        if os.path.exists(gfa_file) and os.path.getsize(gfa_file) > 0:
            return gfa_file
    except CalledProcessError as e:
        print(f"Hifiasm assembly failed for {input_fastq}: {e}")
    return None

def extract_sequence_from_gfa(gfa_file):
    """
    Extract assembled sequence from GFA file.

    Args:
    gfa_file (str): Path to GFA file.

    Returns:
    str or None: Assembled sequence or None if not found.
    """
    assembled_sequences = []
    with open(gfa_file) as f:
        for line in f:
            if line.startswith('S'):  # Sequence line in GFA
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    assembled_sequences.append(parts[2])
    return "".join(assembled_sequences) if assembled_sequences else None

def select_best_read(reads):
    """
    Select the longest read as the best read.

    Args:
    reads (list): List of sequences.

    Returns:
    str or None: The longest read or None if no reads.
    """
    return max(reads, key=len) if reads else None

def process_gap(fastq_file, output_file, assembly_prefix, threads=4):
    """
    Process a single gap to generate consensus sequence.

    Args:
    fastq_file (str): Path to gap-specific FASTQ file.
    output_file (str): Path to output consensus FASTA file.
    assembly_prefix (str): Prefix for assembly output.
    threads (int): Number of threads.
    """
    reads = parse_fastq(fastq_file)

    if not reads:
        print(f"No reads found in {fastq_file}. Skipping.")
        return

    if len(reads) == 1:
        consensus_sequence = reads[0]
    else:
        # Attempt to assemble with hifiasm
        gfa_file = run_hifiasm_assembly(fastq_file, assembly_prefix, threads)
        if gfa_file:
            consensus_sequence = extract_sequence_from_gfa(gfa_file)
            if not consensus_sequence:
                print(f"Assembly failed to produce a consensus for {fastq_file}.")
                consensus_sequence = select_best_read(reads) or generate_consensus_from_reads(reads)
        else:
            # Assembly failed, select best read or generate consensus
            print(f"Assembly failed for {fastq_file}. Selecting best read or generating consensus.")
            best_read = select_best_read(reads)
            if best_read:
                consensus_sequence = best_read
            else:
                consensus_sequence = generate_consensus_from_reads(reads)

    # Write consensus to output file
    with open(output_file, 'w') as out:
        out.write(f">{os.path.basename(output_file).replace('_consensus.fasta', '')}\n")
        out.write(consensus_sequence + "\n")
    print(f"Consensus sequence written to {output_file}")

def process_all_gaps(input_dir, output_dir, threads=4):
    """
    Process all gaps to generate consensus sequences.

    Args:
    input_dir (str): Directory containing gap-specific FASTQ files.
    output_dir (str): Directory to save consensus FASTA files.
    threads (int): Number of threads.
    """
    os.makedirs(output_dir, exist_ok=True)
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq")]
    for fastq_file in tqdm(fastq_files, desc="Generating Consensus Sequences"):
        gap_id = fastq_file.replace(".fastq", "")
        input_fastq = os.path.join(input_dir, fastq_file)
        output_file = os.path.join(output_dir, f"{gap_id}_consensus.fasta")
        assembly_prefix = os.path.join(output_dir, f"{gap_id}_assembly")
        process_gap(input_fastq, output_file, assembly_prefix, threads)

# ================== Step 4: Replace Gaps in Reference Genome ==================

def combine_consensus_sequences(consensus_dir, combined_fasta):
    """
    Combine all consensus sequences into a single FASTA file.

    Args:
    consensus_dir (str): Directory containing consensus FASTA files.
    combined_fasta (str): Path to output combined FASTA file.

    Returns:
    dict: {gap_id: sequence}
    """
    consensus_dict = {}
    with open(combined_fasta, 'w') as out_f:
        for consensus_file in os.listdir(consensus_dir):
            if consensus_file.endswith("_consensus.fasta"):
                gap_id = consensus_file.replace("_consensus.fasta", "")
                with open(os.path.join(consensus_dir, consensus_file), 'r') as f:
                    header = f.readline().strip()
                    if not header.startswith('>'):
                        print(f"Invalid FASTA format in {consensus_file}. Skipping.")
                        continue
                    sequence = f.readline().strip()
                    consensus_dict[gap_id] = sequence
                    out_f.write(f">{gap_id}\n{sequence}\n")
    print(f"Combined all consensus sequences into {combined_fasta}")
    return consensus_dict

def run_minimap2(reference_fasta, combined_fasta, paf_output, threads=4):
    """
    Run minimap2 to align consensus sequences to the reference genome.

    Args:
    reference_fasta (str): Path to reference genome FASTA.
    combined_fasta (str): Path to combined consensus FASTA.
    paf_output (str): Path to output PAF file.
    threads (int): Number of threads.
    """
    cmd = f"minimap2 -x asm5 -t {threads} {reference_fasta} {combined_fasta} > {paf_output}"
    try:
        run(cmd, shell=True, check=True)
        print(f"Minimap2 mapping completed. Results saved to {paf_output}")
    except CalledProcessError as e:
        print(f"Error running minimap2: {e}")
        sys.exit(1)

def parse_paf_and_select_best(paf_file, gaps, consensus_dict):
    """
    Parse PAF file and select the best alignment for each consensus sequence.

    Args:
    paf_file (str): Path to PAF file.
    gaps (DataFrame): DataFrame containing gaps.
    consensus_dict (dict): {gap_id: sequence}

    Returns:
    dict: {gap_id: (ref_start, ref_end, coverage, sequence)}
    """
    gap_info = {f"{row['chrom']}_{row['start']}_{row['end']}": (row['chrom'], row['start'], row['end'])
                for _, row in gaps.iterrows()}

    best_matches = {}

    print("Parsing PAF and selecting best matches...")
    with open(paf_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            query_name, ref_name = cols[0], cols[5]
            try:
                ref_start, ref_end = int(cols[7]), int(cols[8])
            except ValueError:
                continue
            overlap = ref_end - ref_start

            if query_name in gap_info:
                chrom, gap_start, gap_end = gap_info[query_name]

                # Ensure alignment covers the gap
                if ref_name == chrom and ref_start <= gap_end and ref_end >= gap_start:
                    replace_start = max(ref_start, gap_start)
                    replace_end = min(ref_end, gap_end)
                    coverage = replace_end - replace_start

                    # Extract the corresponding sequence from consensus
                    consensus_seq = consensus_dict.get(query_name, "")
                    if not consensus_seq:
                        continue
                    # Calculate slice indices relative to consensus sequence
                    # Assuming consensus sequence corresponds exactly to the gap region
                    # This may need adjustment based on actual alignment
                    best_seq = consensus_seq[replace_start - gap_start:replace_end - gap_start]

                    # Select the best coverage
                    if query_name not in best_matches or coverage > best_matches[query_name][2]:
                        best_matches[query_name] = (replace_start, replace_end, coverage, best_seq)

    return best_matches

def replace_gaps(reference_fasta, best_matches, output_fasta):
    """
    Replace gaps in the reference genome with consensus sequences.

    Args:
    reference_fasta (str): Path to reference genome FASTA.
    best_matches (dict): {gap_id: (ref_start, ref_end, coverage, sequence)}
    output_fasta (str): Path to output updated FASTA.
    """
    ref_genome = Fasta(reference_fasta)
    new_sequences = {chrom: list(str(ref_genome[chrom][:].seq)) for chrom in ref_genome.keys()}

    print("Replacing gaps in the reference genome with consensus sequences...")
    for gap_id, (start, end, _, seq) in best_matches.items():
        chrom = gap_id.split('_')[0]
        if chrom not in new_sequences:
            print(f"Chromosome {chrom} not found in reference genome. Skipping.")
            continue
        # Replace the sequence
        new_sequences[chrom][start:end] = list(seq)
        print(f"Replaced {chrom}:{start}-{end} with consensus sequence for {gap_id}")

    # Write the new genome to output FASTA
    with open(output_fasta, 'w') as f_out:
        for chrom, seq in new_sequences.items():
            f_out.write(f">{chrom}\n")
            for i in range(0, len(seq), 80):
                f_out.write("".join(seq[i:i + 80]) + "\n")
    print(f"New genome with replaced gaps saved to {output_fasta}")

def step4_replace_gaps(args):
    """
    Execute Step 4: Replace gaps in reference genome with consensus sequences.

    Args:
    args: Parsed command-line arguments.
    """
    # Combine all consensus sequences
    print("Combining consensus sequences...")
    combined_fasta = "combined_consensus.fasta"
    consensus_dict = combine_consensus_sequences(args.consensus_dir, combined_fasta)

    # Run minimap2 to map consensus sequences to the reference
    print("Running minimap2 to map consensus sequences to the reference genome...")
    unified_paf = "combined_mapping.paf"
    run_minimap2(args.reference, combined_fasta, unified_paf, args.threads)

    # Parse PAF and select best matches
    print("Parsing PAF and selecting best matches...")
    gaps = load_gaps_df(args.gaps)
    best_matches = parse_paf_and_select_best(unified_paf, gaps, consensus_dict)

    # Replace gaps in the reference genome
    print("Replacing gaps in the reference genome with consensus sequences...")
    replace_gaps(args.reference, best_matches, args.output)

    # Optionally, remove temporary files
    if os.path.exists(combined_fasta):
        os.remove(combined_fasta)
        print(f"Temporary combined FASTA {combined_fasta} removed.")
    if os.path.exists(unified_paf):
        os.remove(unified_paf)
        print(f"Temporary PAF file {unified_paf} removed.")

# ================== Main Pipeline Controller ==================

def main():
    parser = argparse.ArgumentParser(description="Genome Gap Filling Pipeline")
    subparsers = parser.add_subparsers(dest='step', required=True, help='Pipeline steps')

    # Step 0 Parser
    parser_step0 = subparsers.add_parser('step0_filter_paf', help='Filter PAF reads that span gaps')
    parser_step0.add_argument('--input_paf', required=True, help='Input PAF file')
    parser_step0.add_argument('--gaps_bed', required=True, help='Gaps BED file')
    parser_step0.add_argument('--output_paf', required=True, help='Output filtered PAF file')

    # Step 1 Parser
    parser_step1 = subparsers.add_parser('step1_split_paf', help='Split PAF by each gap')
    parser_step1.add_argument('--filtered_paf', required=True, help='Filtered PAF file from step0')
    parser_step1.add_argument('--gaps_bed', required=True, help='Gaps BED file')
    parser_step1.add_argument('--output_dir', required=True, help='Output directory for split PAF files')
    parser_step1.add_argument('--min_reads', type=int, default=1, help='Minimum number of reads required per gap (default: 1)')

    # Step 2 Parser
    parser_step2 = subparsers.add_parser('step2_extract_split_fastq', help='Extract and split reads from FASTQ files')
    parser_step2.add_argument('--paf_dir', required=True, help='Directory containing split PAF files from step1')
    parser_step2.add_argument('--fastq_files', nargs='+', required=True, help='List of input FASTQ files')
    parser_step2.add_argument('--output_dir', required=True, help='Output directory for gap-specific FASTQ files')
    parser_step2.add_argument('--temp_reads_list', default="temp_reads_list.txt", help='Temporary reads list file')
    parser_step2.add_argument('--extracted_fastq', default="extracted_reads.fastq", help='Temporary extracted FASTQ file')

    # Step 3 Parser
    parser_step3 = subparsers.add_parser('step3_generate_consensus', help='Generate consensus sequences for gaps')
    parser_step3.add_argument('--input_dir', required=True, help='Directory containing gap-specific FASTQ files from step2')
    parser_step3.add_argument('--output_dir', required=True, help='Output directory for consensus sequences')
    parser_step3.add_argument('--threads', type=int, default=4, help='Number of threads for hifiasm')

    # Step 4 Parser
    parser_step4 = subparsers.add_parser('step4_replace_gaps', help='Replace gaps in reference genome with consensus sequences')
    parser_step4.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser_step4.add_argument('--gaps', required=True, help='Gaps BED file')
    parser_step4.add_argument('--consensus_dir', required=True, help='Directory containing consensus FASTA files from step3')
    parser_step4.add_argument('--output', required=True, help='Output reference genome FASTA with gaps replaced')
    parser_step4.add_argument('--threads', type=int, default=4, help='Number of threads for minimap2')

    # Full Pipeline Parser
    parser_full = subparsers.add_parser('full_pipeline', help='Run the full pipeline from step0 to step4')
    parser_full.add_argument('--input_paf', required=True, help='Input PAF file for step0')
    parser_full.add_argument('--gaps_bed', required=True, help='Gaps BED file for step0 and step1')
    parser_full.add_argument('--fastq_files', nargs='+', required=True, help='List of input FASTQ files for step2')
    parser_full.add_argument('--consensus_output_dir', required=True, help='Output directory for consensus sequences in step3')
    parser_full.add_argument('--final_output', required=True, help='Final reference genome FASTA with gaps replaced')
    parser_full.add_argument('--threads', type=int, default=4, help='Number of threads for hifiasm and minimap2')
    parser_full.add_argument('--min_reads', type=int, default=1, help='Minimum number of reads required per gap (default: 1)')

    args = parser.parse_args()

    if args.step == 'step0_filter_paf':
        gaps = load_gaps(args.gaps_bed)
        filter_paf_reads(args.input_paf, gaps, args.output_paf)

    elif args.step == 'step1_split_paf':
        gaps_df = load_gaps_df(args.gaps_bed)
        split_paf_by_gap(args.filtered_paf, gaps_df, args.output_dir, args.min_reads)

    elif args.step == 'step2_extract_split_fastq':
        step2_extract_and_split(args)

    elif args.step == 'step3_generate_consensus':
        process_all_gaps(args.input_dir, args.output_dir, args.threads)

    elif args.step == 'step4_replace_gaps':
        step4_replace_gaps(args)

    elif args.step == 'full_pipeline':
        # Step 0
        print("=== Step 0: Filter PAF Reads that Span Gaps ===")
        gaps = load_gaps(args.gaps_bed)
        filtered_paf = "filtered.paf"
        filter_paf_reads(args.input_paf, gaps, filtered_paf)

        # Step 1
        print("=== Step 1: Split PAF by Each Gap ===")
        gaps_df = load_gaps_df(args.gaps_bed)
        split_paf_dir = "split_paf"
        split_paf_by_gap(filtered_paf, gaps_df, split_paf_dir, args.min_reads)

        # Step 2
        print("=== Step 2: Extract and Split Reads from FASTQ Files ===")
        temp_reads_list = "temp_reads_list.txt"
        extracted_fastq = "extracted_reads.fastq"
        split_fastq_dir = "gap_fastq"
        # Create an object with the necessary attributes
        class ArgsStep2:
            def __init__(self, paf_dir, fastq_files, output_dir, temp_reads_list, extracted_fastq):
                self.paf_dir = paf_dir
                self.fastq_files = fastq_files
                self.output_dir = output_dir
                self.temp_reads_list = temp_reads_list
                self.extracted_fastq = extracted_fastq

        step2_args = ArgsStep2(
            paf_dir=split_paf_dir,
            fastq_files=args.fastq_files,
            output_dir=split_fastq_dir,
            temp_reads_list=temp_reads_list,
            extracted_fastq=extracted_fastq
        )
        step2_extract_and_split(step2_args)

        # Step 3
        print("=== Step 3: Generate Consensus Sequences for Each Gap ===")
        consensus_dir = args.consensus_output_dir
        process_all_gaps(split_fastq_dir, consensus_dir, args.threads)

        # Step 4
        print("=== Step 4: Replace Gaps in Reference Genome ===")
        output_fasta = args.final_output
        # Create an object with the necessary attributes
        class ArgsStep4:
            def __init__(self, reference, gaps, consensus_dir, output, threads):
                self.reference = reference
                self.gaps = gaps
                self.consensus_dir = consensus_dir
                self.output = output
                self.threads = threads

        step4_args = ArgsStep4(
            reference=args.reference,
            gaps=args.gaps_bed,
            consensus_dir=consensus_dir,
            output=output_fasta,
            threads=args.threads
        )
        step4_replace_gaps(step4_args)

        # Clean up intermediate files
        intermediate_files = [filtered_paf, "combined_consensus.fasta", "combined_mapping.paf"]
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
                print(f"Temporary file {file} removed.")

        print("=== Full Pipeline Completed Successfully ===")

if __name__ == "__main__":
    main()
