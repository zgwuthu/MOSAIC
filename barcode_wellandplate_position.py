# This code was completed with the assistance of Deepseek
# For each sequencing read, the identity was determined based on three sequential barcodes: 
# the gene barcode, located between the sequences TGCAAGAGACTTCCATCCAG and GCTTCCGGTCTGGTTCGCTTTGAAGCTCGA; 
# the well barcode, situated between GCTTCCGGTCTGGTTCGCTTTGAAGCTCGA and CCGAGCCCACGAGACTATCG; 
# and the plate barcode, identified by the five bases immediately following CCGAGCCCACGAGACTATCG.
# Barcode assignment was achieved using three reference files: 
# one listing 100 gene barcodes (formatted as Gene Name [Tab] Barcode Sequence), one containing 96 well barcodes, and one containing 50 plate barcodes.
# The primary objective was to assign each read to its specific plate, well, and gene.
# Acknowledging the high-throughput nature of NGS and the potential for experimental cross-contamination to introduce multiple genes per well, 
# the final output quantified the top two most frequently observed genes within every unique well, 
# reporting both the gene identity and the specific read counts in a single, clearly formatted Excel spreadsheet for subsequent analysis.




import gzip
import pandas as pd
from collections import defaultdict, Counter
import re
from Bio.Seq import Seq
import os


def load_barcode_references(gene_ref_file, well_ref_file, plate_ref_file):
    """
    Load three reference files
    """
    # Load gene barcode reference file (each line: gene_name\tbarcode_sequence)
    gene_barcodes = {}
    with open(gene_ref_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            # Split by tab
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_name = parts[0].strip()
                barcode = parts[1].strip().upper()  # Convert to uppercase
                gene_barcodes[barcode] = gene_name
                if line_num <= 5:  # Print first 5 as examples
                    print(f"Gene barcode: {gene_name} -> {barcode}")
            else:
                print(f"Warning: Gene barcode file line {line_num} has incorrect format: {line}")

    # Load well barcode reference file (each line: one barcode)
    well_barcodes = {}
    with open(well_ref_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            barcode = line.strip().upper()  # Convert to uppercase
            well_name = f"Well_{line_num:02d}"
            well_barcodes[barcode] = well_name
            if line_num <= 5:  # Print first 5 as examples
                print(f"Well barcode: {well_name} -> {barcode}")

    # Load plate barcode reference file (each line: one barcode)
    plate_barcodes = {}
    with open(plate_ref_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            barcode = line.strip().upper()  # Convert to uppercase
            plate_name = f"Plate_{line_num:02d}"
            plate_barcodes[barcode] = plate_name
            if line_num <= 5:  # Print first 5 as examples
                print(f"Plate barcode: {plate_name} -> {barcode}")

    print(f"Loaded {len(gene_barcodes)} gene barcodes")
    print(f"Loaded {len(well_barcodes)} well barcodes")
    print(f"Loaded {len(plate_barcodes)} plate barcodes")

    return gene_barcodes, well_barcodes, plate_barcodes


def extract_barcodes(sequence):
    """
    Extract three types of barcodes from sequence
    """
    sequence_upper = sequence.upper()

    # Define fixed sequences
    gene_start = "TGCAAGAGACTTCCATCCAG"
    gene_end = "GCTTCCGGTCTGGTTCGCTTTGAAGCTCGA"  # Also the start of well barcode
    well_end = "CCGAGCCCACGAGACTATCG"  # Also the start of plate barcode

    gene_barcode = None
    well_barcode = None
    plate_barcode = None

    # Find gene barcode (between gene_start and gene_end)
    start_pos = sequence_upper.find(gene_start)
    end_pos = sequence_upper.find(gene_end)

    if start_pos != -1 and end_pos != -1 and start_pos < end_pos:
        gene_barcode_start = start_pos + len(gene_start)
        gene_barcode_end = end_pos
        gene_barcode = sequence_upper[gene_barcode_start:gene_barcode_end]

    # Find well barcode (between gene_end and well_end)
    well_end_pos = sequence_upper.find(well_end)

    if end_pos != -1 and well_end_pos != -1 and end_pos < well_end_pos:
        # Well barcode
        well_barcode_start = end_pos + len(gene_end)
        well_barcode_end = well_end_pos
        well_barcode = sequence_upper[well_barcode_start:well_barcode_end]

    # Find plate barcode (5 bases after well_end)
    if well_end_pos != -1:
        plate_barcode_start = well_end_pos + len(well_end)
        plate_barcode_end = plate_barcode_start + 5
        if plate_barcode_end <= len(sequence_upper):
            plate_barcode = sequence_upper[plate_barcode_start:plate_barcode_end]

    return gene_barcode, well_barcode, plate_barcode


def find_best_match(barcode, reference_dict, max_mismatches=2):
    """
    Find best matching barcode in reference dictionary, allowing few mismatches
    """
    if not barcode:
        return None

    if barcode in reference_dict:
        return reference_dict[barcode]

    # If no exact match, try allowing mismatches
    best_match = None
    best_mismatches = max_mismatches + 1

    for ref_barcode, name in reference_dict.items():
        if len(barcode) != len(ref_barcode):
            continue

        # Calculate number of mismatches
        mismatches = sum(1 for a, b in zip(barcode, ref_barcode) if a != b)
        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_match = name

    if best_mismatches <= max_mismatches:
        return best_match

    return None


def analyze_ngs_barcodes(input_file, gene_ref_file, well_ref_file, plate_ref_file, output_excel):
    """
    Analyze barcode information in NGS file
    """
    # Load reference files
    gene_ref, well_ref, plate_ref = load_barcode_references(gene_ref_file, well_ref_file, plate_ref_file)

    # Statistics data structure: {(plate, well): {gene: count}}
    stats = defaultdict(lambda: defaultdict(int))
    total_reads = 0
    valid_reads = 0

    # Auto-detect file format
    open_func = gzip.open if input_file.endswith('.gz') else open
    file_mode = 'rt' if input_file.endswith('.gz') else 'r'

    print("Start analyzing barcodes in NGS file...")

    try:
        with open_func(input_file, file_mode) as f:
            read_count = 0
            while True:
                header = f.readline().strip()
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()

                if not header:
                    break

                total_reads += 1
                read_count += 1

                # Extract barcodes
                gene_barcode, well_barcode, plate_barcode = extract_barcodes(sequence)

                # Find corresponding names
                gene_name = find_best_match(gene_barcode, gene_ref) if gene_barcode else None
                well_name = find_best_match(well_barcode, well_ref) if well_barcode else None
                plate_name = find_best_match(plate_barcode, plate_ref) if plate_barcode else None

                # If all barcodes successfully identified
                if gene_name and well_name and plate_name:
                    valid_reads += 1
                    stats[(plate_name, well_name)][gene_name] += 1

                # Show progress
                if read_count % 100000 == 0:
                    print(f"Processed {read_count} reads, valid reads: {valid_reads}")

        print(f"\nProcessing completed!")
        print(f"Total reads: {total_reads}")
        print(f"Valid reads: {valid_reads}")
        print(f"Efficiency: {valid_reads / total_reads * 100:.2f}%")

        # Generate result table
        result_data = []

        for (plate, well), gene_counts in stats.items():
            # Get top two genes by frequency
            top_genes = Counter(gene_counts).most_common(2)
            total_count = sum(gene_counts.values())

            row = {
                'Plate': plate,
                'Well': well,
                'Total_Reads': total_count,
                'Gene_Diversity': len(gene_counts)  # Number of different genes detected in this well
            }

            # Add top two genes and their counts
            for i, (gene, count) in enumerate(top_genes, 1):
                row[f'Gene_{i}_Name'] = gene
                row[f'Gene_{i}_Count'] = count
                row[f'Gene_{i}_Percentage'] = f"{count / total_count * 100:.2f}%"

            # If only one gene, set second gene as empty
            if len(top_genes) < 2:
                row['Gene_2_Name'] = 'N/A'
                row['Gene_2_Count'] = 0
                row['Gene_2_Percentage'] = '0%'

            result_data.append(row)

        # Create DataFrame and sort
        df = pd.DataFrame(result_data)

        # Sort by plate and well
        df = df.sort_values(['Plate', 'Well'])

        # Save to Excel
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            # Main analysis sheet
            df.to_excel(writer, sheet_name='Well_Gene_Analysis', index=False)

            # Create summary sheet
            summary_data = {
                'Metric': ['Total Reads', 'Valid Reads', 'Efficiency', 'Detected Plates', 'Detected Wells',
                           'Detected Genes'],
                'Value': [
                    total_reads,
                    valid_reads,
                    f"{valid_reads / total_reads * 100:.2f}%",
                    len(set([plate for plate, well in stats.keys()])),
                    len(stats),
                    len(set(gene for gene_counts in stats.values() for gene in gene_counts))
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

            # Create gene distribution sheet
            gene_dist = defaultdict(int)
            for (plate, well), gene_counts in stats.items():
                for gene, count in gene_counts.items():
                    gene_dist[gene] += count

            gene_dist_df = pd.DataFrame({
                'Gene': list(gene_dist.keys()),
                'Total_Count': list(gene_dist.values()),
                'Percentage': [f"{count / sum(gene_dist.values()) * 100:.2f}%" for count in gene_dist.values()]
            }).sort_values('Total_Count', ascending=False)
            gene_dist_df.to_excel(writer, sheet_name='Gene_Distribution', index=False)

            # Create plate statistics sheet
            plate_stats = defaultdict(lambda: {'wells': 0, 'reads': 0})
            for (plate, well), gene_counts in stats.items():
                plate_stats[plate]['wells'] += 1
                plate_stats[plate]['reads'] += sum(gene_counts.values())

            plate_df = pd.DataFrame({
                'Plate': list(plate_stats.keys()),
                'Detected_Wells': [stats['wells'] for stats in plate_stats.values()],
                'Total_Reads': [stats['reads'] for stats in plate_stats.values()]
            }).sort_values('Plate')
            plate_df.to_excel(writer, sheet_name='Plate_Statistics', index=False)

        print(f"Results saved to: {output_excel}")
        print(f"Excel file contains the following sheets:")
        print("- Well_Gene_Analysis: Top two gene statistics for each well")
        print("- Summary: Analysis summary information")
        print("- Gene_Distribution: Distribution of all genes")
        print("- Plate_Statistics: Detection status of each plate")

        return df, valid_reads, total_reads

    except Exception as e:
        print(f"Error processing file: {e}")
        import traceback
        traceback.print_exc()
        return None, 0, 0


# Usage example
if __name__ == "__main__":
    # File path configuration
    input_file = "xxx.fq.gz"  # Replace with your NGS file
    gene_ref_file = "gene_barcodes.txt"  # Gene barcode reference file (each line: gene_name\tbarcode_sequence)
    well_ref_file = "well_barcodes.txt"  # Well barcode reference file (each line: one barcode)
    plate_ref_file = "plate_barcodes.txt"  # Plate barcode reference file (each line: one barcode)
    output_excel = "barcode_analysis_results.xlsx"

    print("Start analyzing NGS barcode...")

    # Run analysis
    df, valid, total = analyze_ngs_barcodes(
        input_file, gene_ref_file, well_ref_file, plate_ref_file, output_excel
    )

    if df is not None:
        print("\nAnalysis completed!")
        # Display some statistics
        print(f"Detected plates: {df['Plate'].nunique()}")
        print(f"Detected wells: {len(df)}")
        print(f"Detected gene types: {len(set(df['Gene_1_Name'].unique()) | set(df['Gene_2_Name'].unique()) - {'N/A'})}")
