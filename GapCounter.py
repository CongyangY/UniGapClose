import argparse
import concurrent.futures
import re

def count_n_and_gaps(sequence):
    # 使用正则表达式查找连续的 'N' 序列
    gaps = re.findall('N+', sequence)
    return sequence.count('N'), len(gaps)

def process_sequences(data):
    seq_name, sequence = data
    return seq_name, *count_n_and_gaps(sequence)

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        return file.read().split('>')[1:]  # 分割每个序列，跳过空的第一部分

def split_sequences(raw_sequences):
    sequences = []
    for raw_sequence in raw_sequences:
        lines = raw_sequence.split('\n', 1)
        seq_name = lines[0].strip()
        sequence = lines[1].replace('\n', '')  # 移除换行符
        sequences.append((seq_name, sequence))
    return sequences

def process_fasta(file_path, output_file):
    raw_sequences = read_fasta(file_path)
    sequences = split_sequences(raw_sequences)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(process_sequences, sequences))

    total_n = sum(n for _, n, _ in results)
    total_gaps = sum(gaps for _, _, gaps in results)

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

