

import os
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import psutil
import glob


def read_cds_file(cds_file: str) -> Dict[str, str]:
    """Read CDS reference file"""
    print("Reading CDS reference file...")
    cds_genes = {}
    
    with open(cds_file, 'r') as f:
        current_gene = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene and current_seq:
                    cds_genes[current_gene] = ''.join(current_seq)
                
                current_gene = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_gene and current_seq:
            cds_genes[current_gene] = ''.join(current_seq)

    print(f"Read {len(cds_genes)} genes from CDS file")
    return cds_genes


def read_feature_matched_reads(feature_reads_file: str) -> Dict[str, List[str]]:
    """Read feature-matched reads file, extract gene names and sequences"""
    print("Reading feature-matched reads file...")
    gene_reads = defaultdict(list)
    
    with open(feature_reads_file, 'r') as f:
        current_gene = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene and current_seq:
                    gene_reads[current_gene].append(''.join(current_seq))
                
                # Extract gene name from header (remove _readX suffix)
                header = line[1:]
                gene_name = header.split('_read')[0]
                current_gene = gene_name
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_gene and current_seq:
            gene_reads[current_gene].append(''.join(current_seq))
    
    print(f"Read reads for {len(gene_reads)} genes from file")
    total_reads = sum(len(reads) for reads in gene_reads.values())
    print(f"Total reads: {total_reads}")
    
    return gene_reads


def find_best_alignment(gene_seq: str, read_seq: str) -> Tuple[int, int]:
    """Find best alignment position for gene sequence in read sequence"""
    gene_len = len(gene_seq)
    read_len = len(read_seq)
    
    # If read is shorter than gene, cannot align
    if read_len < gene_len:
        return 0, 0
    
    best_score = 0
    best_start = 0
    
    # Slide window to find best match position
    for start in range(0, read_len - gene_len + 1):
        end = start + gene_len
        window = read_seq[start:end]
        
        # Calculate match score
        score = sum(1 for a, b in zip(gene_seq, window) if a == b)
        
        if score > best_score:
            best_score = score
            best_start = start
    
    return best_start, best_score


def count_mutations_with_consolidation(gene_seq: str, aligned_read_seq: str) -> Tuple[int, List[Tuple[str, int, str, str]]]:
    """
    Count mutation numbers, consolidate consecutive insertions/deletions as one mutation
    Returns: (total mutations, mutation details list)
    """
    mutations = []
    i, j = 0, 0
    gene_len = len(gene_seq)
    read_len = len(aligned_read_seq)
    
    while i < gene_len or j < read_len:
        if i < gene_len and j < read_len and gene_seq[i] == aligned_read_seq[j]:
            # Match, move forward
            i += 1
            j += 1
        else:
            # Found difference
            mutation_pos = i  # Position based on gene sequence
            
            # Check if substitution, deletion or insertion
            if i < gene_len and j < read_len:
                # Substitution
                mutations.append(('substitution', mutation_pos, gene_seq[i], aligned_read_seq[j]))
                i += 1
                j += 1
            elif j < read_len:
                # Insertion (read has more bases than gene)
                inserted_bases = []
                while j < read_len and (i >= gene_len or gene_seq[i] != aligned_read_seq[j]):
                    inserted_bases.append(aligned_read_seq[j])
                    j += 1
                # Consecutive insertions count as one mutation
                mutations.append(('insertion', mutation_pos, '', ''.join(inserted_bases)))
            else:
                # Deletion (gene has more bases than read)
                deleted_bases = []
                while i < gene_len and (j >= read_len or gene_seq[i] != aligned_read_seq[j]):
                    deleted_bases.append(gene_seq[i])
                    i += 1
                # Consecutive deletions count as one mutation
                mutations.append(('deletion', mutation_pos, ''.join(deleted_bases), ''))
    
    # Consolidate consecutive insertions/deletions
    consolidated_mutations = []
    i = 0
    while i < len(mutations):
        current_mut = mutations[i]
        
        if current_mut[0] in ['insertion', 'deletion']:
            # Check if subsequent mutations are consecutive same type
            j = i + 1
            while (j < len(mutations) and 
                   mutations[j][0] == current_mut[0] and 
                   mutations[j][1] == current_mut[1] + (j - i)):
                j += 1
            
            if j > i + 1:
                # Consolidate consecutive insertions/deletions
                if current_mut[0] == 'insertion':
                    all_inserted = current_mut[3]
                    for k in range(i + 1, j):
                        all_inserted += mutations[k][3]
                    consolidated_mutations.append(('insertion', current_mut[1], '', all_inserted))
                else:  # deletion
                    all_deleted = current_mut[2]
                    for k in range(i + 1, j):
                        all_deleted += mutations[k][2]
                    consolidated_mutations.append(('deletion', current_mut[1], all_deleted, ''))
                i = j
            else:
                consolidated_mutations.append(current_mut)
                i += 1
        else:
            consolidated_mutations.append(current_mut)
            i += 1
    
    return len(consolidated_mutations), consolidated_mutations


def analyze_mutations_for_gene(gene_name: str, gene_seq: str, read_seqs: List[str]) -> Dict[str, any]:
    """Analyze mutation situation for all reads of one gene"""
    results = {
        'gene_name': gene_name,
        'gene_length': len(gene_seq),
        'total_reads': len(read_seqs),
        'reads_analysis': [],
        'min_mutation_count': float('inf'),
        'min_mutation_reads': []  # Store read information with minimum mutation count
    }
    
    for read_idx, read_seq in enumerate(read_seqs):
        # Find best alignment position
        best_start, match_score = find_best_alignment(gene_seq, read_seq)
        
        # Extract aligned read sequence portion
        aligned_read = read_seq[best_start:best_start + len(gene_seq)]
        
        # If extracted sequence is not long enough, pad with N
        if len(aligned_read) < len(gene_seq):
            aligned_read += 'N' * (len(gene_seq) - len(aligned_read))
        
        # Count mutations
        mutation_count, mutation_details = count_mutations_with_consolidation(gene_seq, aligned_read)
        
        # Calculate match rate
        match_rate = match_score / len(gene_seq) if len(gene_seq) > 0 else 0
        
        read_analysis = {
            'read_index': read_idx,
            'read_length': len(read_seq),
            'alignment_start': best_start,
            'match_score': match_score,
            'match_rate': match_rate,
            'mutation_count': mutation_count,
            'mutation_details': mutation_details,
            'aligned_read_seq': aligned_read
        }
        
        results['reads_analysis'].append(read_analysis)
        
        # Update minimum mutation count
        if mutation_count < results['min_mutation_count']:
            results['min_mutation_count'] = mutation_count
            results['min_mutation_reads'] = [read_analysis]
        elif mutation_count == results['min_mutation_count']:
            results['min_mutation_reads'].append(read_analysis)
    
    return results


def save_mutation_analysis(results: Dict[str, any], output_dir: str = "mutation_analysis"):
    """Save mutation analysis results"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    gene_name = results['gene_name']
    
    # Save detailed mutation information
    detail_file = os.path.join(output_dir, f"{gene_name}_mutations_detail.tsv")
    with open(detail_file, 'w') as f:
        f.write("Read_Index\tRead_Length\tAlignment_Start\tMatch_Score\tMatch_Rate\tMutation_Count\tMutation_Type\tPosition\tReference\tVariant\n")
        
        for read_analysis in results['reads_analysis']:
            base_info = f"{read_analysis['read_index']}\t{read_analysis['read_length']}\t{read_analysis['alignment_start']}\t{read_analysis['match_score']}\t{read_analysis['match_rate']:.4f}\t{read_analysis['mutation_count']}"
            
            if read_analysis['mutation_count'] == 0:
                f.write(f"{base_info}\tNone\tNone\tNone\tNone\n")
            else:
                for mut_type, pos, ref, var in read_analysis['mutation_details']:
                    f.write(f"{base_info}\t{mut_type}\t{pos}\t{ref}\t{var}\n")
    
    # Save summary information
    summary_file = os.path.join(output_dir, f"{gene_name}_summary.tsv")
    with open(summary_file, 'w') as f:
        f.write("Gene_Name\tGene_Length\tTotal_Reads\tMin_Mutation_Count\tReads_With_Min_Mutation\tAvg_Mutations\tAvg_Match_Rate\tPerfect_Reads\n")
        
        total_mutations = sum(ra['mutation_count'] for ra in results['reads_analysis'])
        avg_mutations = total_mutations / len(results['reads_analysis']) if results['reads_analysis'] else 0
        
        total_match_rate = sum(ra['match_rate'] for ra in results['reads_analysis'])
        avg_match_rate = total_match_rate / len(results['reads_analysis']) if results['reads_analysis'] else 0
        
        perfect_reads = sum(1 for ra in results['reads_analysis'] if ra['mutation_count'] == 0)
        
        f.write(f"{gene_name}\t{results['gene_length']}\t{results['total_reads']}\t{results['min_mutation_count']}\t{len(results['min_mutation_reads'])}\t{avg_mutations:.2f}\t{avg_match_rate:.4f}\t{perfect_reads}\n")
    
    # Save detailed information for reads with minimum mutations
    if results['min_mutation_reads']:
        min_mutation_file = os.path.join(output_dir, f"{gene_name}_min_mutation_reads.tsv")
        with open(min_mutation_file, 'w') as f:
            f.write("Gene_Name\tMin_Mutation_Count\tRead_Index\tRead_Length\tAlignment_Start\tMatch_Score\tMatch_Rate\tMutation_Details\n")
            
            for read_analysis in results['min_mutation_reads']:
                mutation_details_str = "; ".join([f"{mut_type}@{pos}:{ref}->{var}" 
                                                for mut_type, pos, ref, var in read_analysis['mutation_details']])
                if not mutation_details_str:
                    mutation_details_str = "None"
                
                f.write(f"{gene_name}\t{results['min_mutation_count']}\t{read_analysis['read_index']}\t{read_analysis['read_length']}\t{read_analysis['alignment_start']}\t{read_analysis['match_score']}\t{read_analysis['match_rate']:.4f}\t{mutation_details_str}\n")
    
    return detail_file, summary_file


def save_min_mutations_summary(all_results: Dict[str, Dict], output_file: str = "min_mutations_summary.tsv"):
    """Save summary of minimum mutation counts for all genes"""
    print("Generating minimum mutation count summary file...")
    
    with open(output_file, 'w') as f:
        f.write("Gene_Name\tGene_Length\tTotal_Reads\tMin_Mutation_Count\tReads_With_Min_Mutation\tAvg_Mutations\tPerfect_Reads\n")
        
        for gene_name, results in all_results.items():
            total_mutations = sum(ra['mutation_count'] for ra in results['reads_analysis'])
            avg_mutations = total_mutations / len(results['reads_analysis']) if results['reads_analysis'] else 0
            perfect_reads = sum(1 for ra in results['reads_analysis'] if ra['mutation_count'] == 0)
            
            f.write(f"{gene_name}\t{results['gene_length']}\t{results['total_reads']}\t{results['min_mutation_count']}\t{len(results['min_mutation_reads'])}\t{avg_mutations:.2f}\t{perfect_reads}\n")
    
    print(f"Minimum mutation count summary saved to: {output_file}")
    return output_file


def main_mutation_analysis():
    """Main function: analyze mutations in feature-matched reads"""
    cds_file = "Max_cds_ALL.txt"
    feature_reads_file = "feature_matched_reads.fasta"
    
    if not os.path.exists(cds_file):
        print(f"Error: CDS reference file {cds_file} does not exist")
        return
    
    if not os.path.exists(feature_reads_file):
        print(f"Error: Feature-matched reads file {feature_reads_file} does not exist")
        return
    
    print("=" * 70)
    print("Mutation Analysis Tool - Includes Minimum Mutation Count Statistics")
    print(f"Memory: {psutil.virtual_memory().total / (1024**3):.1f} GB")
    print("=" * 70)
    
    start_total = time.time()
    
    try:
        # Step 1: Read CDS reference sequences
        cds_genes = read_cds_file(cds_file)
        
        # Step 2: Read feature-matched reads
        gene_reads = read_feature_matched_reads(feature_reads_file)
        
        # Step 3: Perform mutation analysis for each gene
        print("\nStarting mutation analysis...")
        
        analyzed_genes = 0
        total_reads_analyzed = 0
        all_results = {}
        
        for gene_name, read_seqs in gene_reads.items():
            if gene_name not in cds_genes:
                print(f"Warning: Gene {gene_name} not found in CDS reference file, skipping")
                continue
            
            print(f"Analyzing gene {gene_name}, total {len(read_seqs)} reads...")
            
            gene_seq = cds_genes[gene_name]
            results = analyze_mutations_for_gene(gene_name, gene_seq, read_seqs)
            all_results[gene_name] = results
            
            # Save results
            detail_file, summary_file = save_mutation_analysis(results)
            
            # Statistics
            total_mutations = sum(ra['mutation_count'] for ra in results['reads_analysis'])
            perfect_reads = sum(1 for ra in results['reads_analysis'] if ra['mutation_count'] == 0)
            
            print(f"  - Minimum mutation count: {results['min_mutation_count']}")
            print(f"  - Number of reads with minimum mutations: {len(results['min_mutation_reads'])}")
            print(f"  - Total mutations: {total_mutations}")
            print(f"  - Average mutations per read: {total_mutations/len(read_seqs):.2f}")
            print(f"  - Perfect match reads: {perfect_reads}/{len(read_seqs)}")
            print(f"  - Results saved to: {detail_file}")
            
            analyzed_genes += 1
            total_reads_analyzed += len(read_seqs)
        
        # Step 4: Save minimum mutation count summary
        summary_file = save_min_mutations_summary(all_results)
        
        # Step 5: Statistics of minimum mutation count distribution
        min_mutation_counts = [results['min_mutation_count'] for results in all_results.values()]
        zero_min_count = sum(1 for count in min_mutation_counts if count == 0)
        one_min_count = sum(1 for count in min_mutation_counts if count == 1)
        two_min_count = sum(1 for count in min_mutation_counts if count == 2)
        three_plus_min_count = sum(1 for count in min_mutation_counts if count >= 3)
        
        print(f"\nMinimum mutation count distribution statistics:")
        print(f"- Genes with minimum mutation count 0: {zero_min_count}")
        print(f"- Genes with minimum mutation count 1: {one_min_count}")
        print(f"- Genes with minimum mutation count 2: {two_min_count}")
        print(f"- Genes with minimum mutation count ≥3: {three_plus_min_count}")
        
        print(f"\nAnalysis completed!")
        print(f"- Total analyzed genes: {analyzed_genes}")
        print(f"- Total analyzed reads: {total_reads_analyzed}")
        print(f"- Detailed results saved in mutation_analysis directory")
        print(f"- Minimum mutation count summary: {summary_file}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
    
    total_time = time.time() - start_total
    print(f"\nTotal processing time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")


if __name__ == "__main__":
    main_mutation_analysis()
