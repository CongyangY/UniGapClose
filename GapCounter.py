import argparse
import concurrent.futures

def count_n_and_gaps(sequence):
    n_count = sequence.count('N')
    gap_count = 0
    i = 0
    while i < len(sequence):
        if sequence[i] == 'N':
            gap_count += 1
            while i < len(sequence) and sequence[i] == 'N':
                i += 1
        else:
            i += 1
    return n_count, gap_count

def process_sequence(seq_name, sequence):
    n_count, gap_count = count_n_and_gaps(sequence)
    return seq_name, n_count, gap_count

def process_fasta(file_path, output_file):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    sequences = []
    current_seq = ''
    seq_name = ''

    for line in lines:
        if line.startswith('>'):
            if current_seq:
                sequences.append((seq_name, current_seq))
            seq_name = line.strip().split('>')[1]
            current_seq = ''
        else:
            current_seq += line.strip()

    if current_seq:
        sequences.append((seq_name, current_seq))

    results = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_sequence, seq[0], seq[1]) for seq in sequences]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    total_n = sum([result[1] for result in results])
    total_gaps = sum([result[2] for result in results])

    with open(output_file, 'w') as out_file:
        out_file.write("File, Total Ns, Total Gaps\n")
        out_file.write(f"{output_file}, {total_n}, {total_gaps}\n")
        out_file.write("-" * 30 + "\n")
        out_file.write("Chromosome, Ns, Gaps\n")
        for result in results:
            out_file.write(f"{result[0]}, {result[1]}, {result[2]}\n")

def main():
    parser = argparse.ArgumentParser(description='Count the number of Ns and continuous N gaps in genome sequences from a FASTA file using parallel processing.')
    parser.add_argument('fasta_file', help='Path to the FASTA file')
    parser.add_argument('output_file', help='Path to the output file')

    args = parser.parse_args()
    process_fasta(args.fasta_file, args.output_file)
    print(f"Results written to {args.output_file}")

if __name__ == "__main__":
    main()

