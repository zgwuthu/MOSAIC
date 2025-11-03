# This code was completed with the assistance of Deepseek
# The file Max_cds_ALL.txt contains reference sequences for all genes. 
# In this code, the paired-end sequencing files are initially assembled by identifying overlapping regions, 
# generating a file of concatenated sequences (joined_sequences.fasta). 
# Subsequently, the reference sequences are searched against the assembly results, 
# with occurrence frequencies being quantified and statistically analyzed. 


port gzip
import os
import time
from collections import defaultdict, Counter, deque
import multiprocessing as mp
from multiprocessing import Pool, Manager
import itertools
import heapq
import numpy as np
import psutil


class AhoCorasick:
    """Aho-Corasick automaton implementation - memory optimized version"""

    def __init__(self, patterns):
        self.patterns = patterns
        self.build()

    def build(self):
        """Construct the automaton"""
        self.goto = [{}]
        self.fail = [0]
        self.output = [[]]

        # Build goto table
        for pattern in self.patterns:
            current = 0
            for char in pattern:
                if char in self.goto[current]:
                    current = self.goto[current][char]
                else:
                    new_state = len(self.goto)
                    self.goto.append({})
                    self.fail.append(0)
                    self.output.append([])
                    self.goto[current][char] = new_state
                    current = new_state
            self.output[current].append(pattern)

        # Build failure table
        queue = deque()
        for char, state in self.goto[0].items():
            self.fail[state] = 0
            queue.append(state)

        while queue:
            current = queue.popleft()
            for char, next_state in self.goto[current].items():
                queue.append(next_state)
                fail_state = self.fail[current]
                while fail_state != 0 and char not in self.goto[fail_state]:
                    fail_state = self.fail[fail_state]
                self.fail[next_state] = self.goto[fail_state].get(char, 0)
                self.output[next_state].extend(self.output[self.fail[next_state]])

    def search(self, text):
        """Search for all patterns in the input text"""
        current = 0
        results = set()

        for char in text:
            while current != 0 and char not in self.goto[current]:
                current = self.fail[current]
            if char in self.goto[current]:
                current = self.goto[current][char]
            else:
                current = 0

            if self.output[current]:
                results.update(self.output[current])

        return results


def reverse_complement(seq):
    """Generate reverse complement sequence - optimized version"""
    complement = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(complement)[::-1]


def average_quality_fast(qualities, threshold=30):
    """Rapid quality assessment - utilizing numpy"""
    return np.mean(qualities) >= threshold


def find_overlap_fast(seq1, seq2, min_overlap=30, max_mismatch=0.05):
    """
    Rapid overlap detection - employing sliding window with early termination
    """
    len1, len2 = len(seq1), len(seq2)
    max_possible_overlap = min(len1, len2)

    # Search from maximum possible overlap downwards
    for overlap in range(max_possible_overlap, min_overlap - 1, -1):
        s1_end = seq1[-overlap:]
        s2_start = seq2[:overlap]

        # Rapid mismatch calculation
        mismatches = sum(1 for a, b in zip(s1_end, s2_start) if a != b)

        if mismatches <= overlap * max_mismatch:
            # Acceptable overlap identified
            # Select base from first sequence at mismatched positions
            merged_overlap = ''.join(s1_end[i] if s1_end[i] == s2_start[i] else s1_end[i]
                                     for i in range(overlap))
            return overlap, merged_overlap

    return 0, ""


def process_fastq_chunk(args):
    """
    Process FASTQ file chunk - multiprocessing version
    Returns read pairs passing quality filtration
    """
    chunk, quality_threshold = args
    filtered_pairs = []

    for i in range(0, len(chunk), 8):  # Each 8 lines = 2 reads (4 lines each)
        if i + 7 >= len(chunk):
            break

        # Parse R1
        r1_header = chunk[i].strip()
        r1_seq = chunk[i + 1].strip()
        r1_plus = chunk[i + 2].strip()
        r1_qual = chunk[i + 3].strip()

        # Parse R2
        r2_header = chunk[i + 4].strip()
        r2_seq = chunk[i + 5].strip()
        r2_plus = chunk[i + 6].strip()
        r2_qual = chunk[i + 7].strip()

        # Quality assessment
        r1_qual_values = [ord(c) - 33 for c in r1_qual]
        r2_qual_values = [ord(c) - 33 for c in r2_qual]

        if (average_quality_fast(r1_qual_values, quality_threshold) and
                average_quality_fast(r2_qual_values, quality_threshold)):
            filtered_pairs.append((r1_seq, r2_seq))

    return filtered_pairs


def process_read_pair_fast(args):
    """Rapid processing of read pair"""
    r1_seq, r2_seq = args

    # Generate reverse complement of r2
    r2_rc = reverse_complement(r2_seq)

    # Identify overlap region and perform assembly
    overlap, merged_overlap = find_overlap_fast(r1_seq, r2_rc, min_overlap=30, max_mismatch=0.05)

    if overlap > 0:
        # Sequence assembly
        joined_seq = r1_seq[:-overlap] + merged_overlap + r2_rc[overlap:]
        return joined_seq

    return None


def read_fastq_in_chunks(file_path, chunk_size=100000):
    """
    Stream FASTQ file reading to prevent memory overload
    Yields chunks (each containing chunk_size*4 lines)
    """
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'

    with open_func(file_path, mode) as f:
        chunk = []
        for i, line in enumerate(f):
            chunk.append(line)
            if len(chunk) >= chunk_size * 8:  # Each read 4 lines, each pair 8 lines
                yield chunk
                chunk = []

        if chunk:  # Process final segment
            yield chunk


def filter_and_join_reads_parallel(r1_file, r2_file, output_file, quality_threshold=30, num_processes=30):
    """
    Parallel read filtering and assembly - optimized for high memory and multi-core systems
    """
    print(f"Initiating processing of {r1_file} and {r2_file}")
    print(f"Utilizing {num_processes} processes")

    total_pairs = 0
    passed_pairs = 0
    successfully_joined = 0
    start_time = time.time()

    # Create output file
    with open(output_file, 'w') as out_f:
        # Parallel processing of each chunk - employing larger process pool
        with Pool(processes=num_processes) as pool:
            # Concurrent reading of R1 and R2 chunks
            r1_chunks = read_fastq_in_chunks(r1_file, chunk_size=50000)  # Larger chunks
            r2_chunks = read_fastq_in_chunks(r2_file, chunk_size=50000)

            chunk_count = 0
            for r1_chunk, r2_chunk in zip(r1_chunks, r2_chunks):
                chunk_count += 1

                # Merge chunks from both files
                combined_chunk = []
                for i in range(0, len(r1_chunk), 4):
                    if i + 3 < len(r1_chunk) and i + 3 < len(r2_chunk):
                        # R1 4 lines + R2 4 lines
                        combined_chunk.extend(r1_chunk[i:i + 4])
                        combined_chunk.extend(r2_chunk[i:i + 4])

                # Quality filtration - larger batch sizes
                filtered_results = pool.map(process_fastq_chunk,
                                            [(combined_chunk[i:i + 40000], quality_threshold)
                                             for i in range(0, len(combined_chunk), 40000)])

                # Collect filtered read pairs
                filtered_pairs = []
                for result in filtered_results:
                    filtered_pairs.extend(result)

                # Sequence assembly - larger batch sizes
                if filtered_pairs:
                    join_results = pool.map(process_read_pair_fast,
                                            [filtered_pairs[i:i + 5000]
                                             for i in range(0, len(filtered_pairs), 5000)])

                    # Write results
                    for result_batch in join_results:
                        for joined_seq in result_batch:
                            if joined_seq is not None:
                                successfully_joined += 1
                                out_f.write(f">read_{successfully_joined}\n{joined_seq}\n")

                # Update statistics
                total_pairs += len(combined_chunk) // 8
                passed_pairs += len(filtered_pairs)

                # Progress reporting
                if chunk_count % 5 == 0:
                    elapsed = time.time() - start_time
                    memory_usage = psutil.virtual_memory().percent
                    print(f"Processed {total_pairs:,} read pairs, passed {passed_pairs:,}, assembled {successfully_joined:,}, "
                          f"Memory utilization: {memory_usage}%, Elapsed: {elapsed:.2f} seconds")

    total_time = time.time() - start_time
    print(f"\nProcessing completed!")
    print(f"Total read pairs: {total_pairs:,}")
    print(f"Passed quality filtration: {passed_pairs:,}")
    print(f"Successfully assembled: {successfully_joined:,}")
    print(f"Total duration: {total_time:.2f} seconds")
    print(f"Processing rate: {total_pairs / total_time:.0f} reads/second")

    return successfully_joined


def build_gene_automaton(ref_file):
    """Construct Aho-Corasick automaton from reference file"""
    print("Constructing gene automaton...")
    start_time = time.time()

    # Load reference sequences
    ref_genes = {}
    with open(ref_file, 'r') as f:
        current_gene = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene and current_seq:
                    ref_genes[current_gene] = ''.join(current_seq)
                current_gene = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_gene and current_seq:
            ref_genes[current_gene] = ''.join(current_seq)

    # Construct automaton
    gene_sequences = list(ref_genes.values())
    automaton = AhoCorasick(gene_sequences)

    build_time = time.time() - start_time
    print(f"Automaton construction completed, duration: {build_time:.2f} seconds")
    print(f"Loaded {len(ref_genes)} genes")

    return automaton, ref_genes


def search_genes_in_sequence_batch(args):
    """Search for genes in a batch of sequences (multiprocessing)"""
    sequences_batch, automaton, ref_genes = args
    gene_counts = defaultdict(int)

    # Create reverse mapping from sequence to gene identifier
    seq_to_gene = {seq: gene for gene, seq in ref_genes.items()}

    for seq in sequences_batch:
        found_genes = automaton.search(seq)
        for gene_seq in found_genes:
            gene_name = seq_to_gene[gene_seq]
            gene_counts[gene_name] += 1

    return gene_counts


def read_joined_sequences_in_batches(joined_file, batch_size=50000):
    """Batch reading of assembled sequences, employing larger batch sizes"""
    sequences_batch = []
    with open(joined_file, 'r') as f:
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences_batch.append(''.join(current_seq))
                    if len(sequences_batch) >= batch_size:
                        yield sequences_batch
                        sequences_batch = []
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_seq:
            sequences_batch.append(''.join(current_seq))

    if sequences_batch:
        yield sequences_batch


def match_genes_parallel(joined_file, ref_file, num_processes=30):
    """Parallel gene matching - optimized for high memory and multi-core systems"""
    print("Initiating gene matching...")
    start_time = time.time()

    # Construct automaton
    automaton, ref_genes = build_gene_automaton(ref_file)

    # Multiprocess gene searching - larger process pool
    gene_counts = defaultdict(int)
    total_seqs = 0

    with Pool(processes=num_processes) as pool:
        # Batch processing of assembled sequences - larger batch sizes
        for batch_num, sequences_batch in enumerate(read_joined_sequences_in_batches(joined_file, batch_size=50000)):
            total_seqs += len(sequences_batch)

            # Prepare arguments
            search_args = [(sequences_batch, automaton, ref_genes)]

            # Multiprocess searching
            batch_results = pool.map(search_genes_in_sequence_batch, search_args)

            # Aggregate results
            for result in batch_results:
                for gene, count in result.items():
                    gene_counts[gene] += count

            # Progress reporting
            if batch_num % 5 == 0:
                elapsed = time.time() - start_time
                memory_usage = psutil.virtual_memory().percent
                print(f"Processed {total_seqs:,} assembled sequences, Memory utilization: {memory_usage}%, Elapsed: {elapsed:.2f} seconds")

    match_time = time.time() - start_time
    print(f"Gene matching completed, duration: {match_time:.2f} seconds")
    print(f"Matching rate: {total_seqs / match_time:.0f} sequences/second")

    return gene_counts


def main():
    """Main function - optimized for high memory and multi-core systems"""
    # File paths
    r1_file = "xxx_1.fq.gz"
    r2_file = "xxx_2.fq.gz"
    ref_file = "Max_cds_ALL.txt"
    output_file = "joined_sequences.fasta"

    # Verify file existence
    for f in [r1_file, r2_file, ref_file]:
        if not os.path.exists(f):
            print(f"Error: File {f} not found")
            return

    # Set process count - utilizing 40-core CPU efficiently
    num_processes = 35  # Reserve cores for system operations

    print("=" * 60)
    print("NGS Data Processing Pipeline - High Memory Multi-Core Optimized Version")
    print(f"CPU Cores: {mp.cpu_count()}, Allocated Processes: {num_processes}")
    print("=" * 60)

    # Step 1: Filtration and assembly
    print("\nStep 1: Quality Filtration and Sequence Assembly")
    start_total = time.time()

    successfully_joined = filter_and_join_reads_parallel(
        r1_file, r2_file, output_file,
        quality_threshold=30,
        num_processes=num_processes
    )

    # Step 2: Gene matching
    if successfully_joined > 0:
        print("\nStep 2: Gene Matching")
        gene_counts = match_genes_parallel(output_file, ref_file, num_processes=num_processes)

        # Output results
        print("\n" + "=" * 60)
        print("Final Results")
        print("=" * 60)
        print(f"Successfully assembled sequences: {successfully_joined:,}")
        print(f"Detected gene types: {len(gene_counts)}")

        # Display matching results
        if gene_counts:
            print("\nGene Matching Statistics (Top 20):")
            sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
            for gene, count in sorted_genes[:20]:
                print(f"  {gene}: {count:,} occurrences")

            # Statistics for unmatched genes
            ref_genes = set()
            with open(ref_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        ref_genes.add(line[1:].strip())

            unmatched = ref_genes - set(gene_counts.keys())
            if unmatched:
                print(f"\nUnmatched gene count: {len(unmatched)}")
                print(f"Examples: {list(unmatched)[:5]}")
        else:
            print("No genes detected")

    total_time = time.time() - start_total
    print(f"\nTotal processing time: {total_time:.2f} seconds ({total_time / 60:.2f} minutes)")
    print(f"Average processing rate: {successfully_joined / total_time:.0f} sequences/second")


if __name__ == "__main__":
    # Set multiprocessing start method
    mp.set_start_method('spawn', force=True)
    main()
