


import os
import time
from collections import defaultdict
from multiprocessing import Pool, cpu_count, Manager
import psutil
import mmap
from typing import Dict, List, Tuple, Set


def read_cds_file_fast(cds_file: str) -> Dict[str, str]:
    """Fast read CDS reference file, properly handle duplicate gene names"""
    print("Reading CDS reference file...")
    start_time = time.time()
    cds_genes = {}
    duplicate_count = 0
    
    with open(cds_file, 'r') as f:
        current_gene = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene and current_seq:
                    if current_gene in cds_genes:
                        duplicate_count += 1
                        new_gene_name = f"{current_gene}_dup{duplicate_count}"
                        cds_genes[new_gene_name] = ''.join(current_seq)
                    else:
                        cds_genes[current_gene] = ''.join(current_seq)
                
                current_gene = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_gene and current_seq:
            if current_gene in cds_genes:
                duplicate_count += 1
                new_gene_name = f"{current_gene}_dup{duplicate_count}"
                cds_genes[new_gene_name] = ''.join(current_seq)
            else:
                cds_genes[current_gene] = ''.join(current_seq)

    elapsed = time.time() - start_time
    print(f"Read {len(cds_genes)} genes from CDS file, time: {elapsed:.2f} seconds")
    if duplicate_count > 0:
        print(f"Found and processed {duplicate_count} duplicate gene names")
    
    return cds_genes


def read_joined_sequences_mmap(joined_file: str) -> List[str]:
    """Fast read joined sequences using memory mapping"""
    print("Reading joined sequences file...")
    start_time = time.time()
    
    sequences = []
    with open(joined_file, 'r') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            current_seq = []
            for line in iter(mm.readline, b""):
                line = line.decode('utf-8').strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq))
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_seq:
                sequences.append(''.join(current_seq))
    
    elapsed = time.time() - start_time
    print(f"Read {len(sequences)} joined sequences, time: {elapsed:.2f} seconds")
    
    return sequences


def exact_match_worker(args):
    """Exact match worker process - stop searching gene once match found"""
    batch_seqs, genes_dict, found_genes = args
    newly_found = set()
    
    for gene, seq in genes_dict.items():
        if gene in found_genes or gene in newly_found:
            continue
            
        for read_seq in batch_seqs:
            if seq in read_seq:
                newly_found.add(gene)
                break
    
    return newly_found


def exact_match_all_genes_fast(genes_dict: Dict[str, str], joined_seqs: List[str], 
                              num_processes: int = 35) -> Set[str]:
    """Fast exact match all genes - stop searching gene once match found"""
    print("Starting fast exact match analysis...")
    start_time = time.time()
    
    manager = Manager()
    found_genes = manager.list()
    
    batch_size = 100000
    batches = [joined_seqs[i:i + batch_size] for i in range(0, len(joined_seqs), batch_size)]
    
    print(f"Using {num_processes} processes to handle {len(batches)} batches...")
    
    total_found = 0
    processed_batches = 0
    
    with Pool(processes=num_processes) as pool:
        for i in range(0, len(batches), 10):
            batch_group = batches[i:i+10]
            
            args_list = [(batch, genes_dict, found_genes) for batch in batch_group]
            results = pool.map(exact_match_worker, args_list)
            
            for result in results:
                for gene in result:
                    if gene not in found_genes:
                        found_genes.append(gene)
                        total_found += 1
            
            processed_batches += len(batch_group)
            
            elapsed = time.time() - start_time
            memory_usage = psutil.virtual_memory().percent
            print(f"Processed {processed_batches}/{len(batches)} batches, found {total_found} genes,"
                  f" memory usage: {memory_usage}%, time: {elapsed:.2f} seconds")
            
            if len(found_genes) >= len(genes_dict):
                break
    
    elapsed = time.time() - start_time
    print(f"Exact match completed, time: {elapsed:.2f} seconds")
    print(f"Found {len(found_genes)} genes through exact matching")
    
    return set(found_genes)


def extract_feature_sequences(gene_seq: str) -> Tuple[str, str, str]:
    """Extract three feature sequences from gene"""
    feature1 = gene_seq[:20]
    
    if len(gene_seq) >= 120:
        feature2 = gene_seq[100:120]
    else:
        feature2 = ""
    
    feature3 = gene_seq[-20:] if len(gene_seq) >= 20 else ""
    
    return feature1, feature2, feature3


def feature_match_worker(args):
    """Feature match worker process - collect matching reads"""
    batch_seqs, unmatched_genes_dict = args
    gene_reads = defaultdict(list)
    
    for gene, seq in unmatched_genes_dict.items():
        feature1, feature2, feature3 = extract_feature_sequences(seq)
        if not (feature1 and feature2 and feature3):
            continue
            
        for read_seq in batch_seqs:
            if (feature1 in read_seq and 
                feature2 in read_seq and 
                feature3 in read_seq):
                gene_reads[gene].append(read_seq)
    
    return gene_reads


def feature_match_all_genes(unmatched_genes_dict: Dict[str, str], joined_seqs: List[str],
                           num_processes: int = 35) -> Tuple[Dict[str, List[str]], Set[str]]:
    """Use feature sequences to match unfound genes and collect matching reads"""
    print("Starting feature sequence match analysis...")
    start_time = time.time()
    
    batch_size = 50000
    batches = [joined_seqs[i:i + batch_size] for i in range(0, len(joined_seqs), batch_size)]
    
    print(f"Using {num_processes} processes to handle {len(batches)} batches...")
    
    all_gene_reads = defaultdict(list)
    list_a = set()
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(feature_match_worker, [(batch, unmatched_genes_dict) for batch in batches])
    
    # Merge results
    for result in results:
        for gene, reads in result.items():
            all_gene_reads[gene].extend(reads)
            list_a.add(gene)
    
    elapsed = time.time() - start_time
    print(f"Feature sequence match completed, time: {elapsed:.2f} seconds")
    print(f"Found {len(list_a)} genes through feature sequence matching")
    
    return all_gene_reads, list_a


def save_feature_matched_reads(gene_reads: Dict[str, List[str]], output_file: str = "feature_matched_reads.fasta"):
    """Save feature-matched reads to file"""
    print("Saving feature-matched reads...")
    
    with open(output_file, 'w') as f:
        read_count = 0
        for gene, reads in gene_reads.items():
            for read in reads:
                f.write(f">{gene}_read{read_count}\n{read}\n")
                read_count += 1
    
    print(f"Feature-matched reads saved to {output_file}, total {read_count} reads")
    return read_count


def save_final_results(exact_matched: Set[str], feature_matched: Set[str], 
                      all_genes: Dict[str, str], output_file: str = "gene_detection_results.tsv"):
    """Save final results in table format"""
    print("Saving final results...")
    
    # Statistics for different match types
    only_exact = exact_matched
    only_feature = feature_matched - exact_matched
    not_found = set(all_genes.keys()) - exact_matched - feature_matched
    
    with open(output_file, 'w') as f:
        f.write("Gene\tDetection_Type\n")
        
        for gene in sorted(only_exact):
            f.write(f"{gene}\tExact_Match\n")
        
        for gene in sorted(only_feature):
            f.write(f"{gene}\tFeature_Match\n")
        
        for gene in sorted(not_found):
            f.write(f"{gene}\tNot_Found\n")
    
    print(f"Final results saved to {output_file}")
    
    print(f"\nStatistics:")
    print(f"Total genes: {len(all_genes)}")
    print(f"Exact matched genes: {len(only_exact)}")
    print(f"Feature matched genes: {len(only_feature)}")
    print(f"Unmatched genes: {len(not_found)}")


def main_optimized():
    """Highly optimized main function - only up to feature matching and output reads"""
    joined_file = "joined_sequences.fasta"
    cds_file = "Max_cds_ALL.txt"
    
    if not os.path.exists(joined_file):
        print(f"Error: Joined sequence file {joined_file} does not exist")
        return
    
    if not os.path.exists(cds_file):
        print(f"Error: CDS reference file {cds_file} does not exist")
        return
    
    num_processes = min(38, cpu_count())
    
    print("=" * 70)
    print("Gene Detection Analysis Tool - Feature Match Version")
    print(f"CPU cores: {cpu_count()}, Assigned processes: {num_processes}")
    print(f"Memory: {psutil.virtual_memory().total / (1024**3):.1f} GB")
    print("=" * 70)
    
    start_total = time.time()
    
    try:
        # Step 1: Read files
        cds_genes = read_cds_file_fast(cds_file)
        joined_seqs = read_joined_sequences_mmap(joined_file)
        
        # Step 2: Fast exact matching
        exact_matched = exact_match_all_genes_fast(cds_genes, joined_seqs, num_processes)
        
        # Step 3: Feature sequence matching for unmatched genes
        unmatched_genes = set(cds_genes.keys()) - exact_matched
        if unmatched_genes:
            print(f"\nStarting feature sequence matching for {len(unmatched_genes)} unmatched genes...")
            unmatched_genes_dict = {gene: cds_genes[gene] for gene in unmatched_genes}
            gene_reads, list_a = feature_match_all_genes(unmatched_genes_dict, joined_seqs, num_processes)
            
            # Save feature-matched reads to file
            total_reads = save_feature_matched_reads(gene_reads)
            print(f"Feature-matched reads file contains {total_reads} reads")
        else:
            list_a = set()
            gene_reads = {}
        
        # Step 4: Save final results
        save_final_results(exact_matched, list_a, cds_genes)
        
        print(f"\nAnalysis completed!")
        print(f"- Exact matched genes: {len(exact_matched)}")
        print(f"- Feature matched genes: {len(list_a)}")
        if gene_reads:
            print(f"- Feature matched reads saved to file, total {total_reads} reads")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
    
    total_time = time.time() - start_total
    print(f"\nTotal processing time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")


if __name__ == "__main__":
    main_optimized()
