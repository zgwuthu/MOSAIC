GCwindow.py is used for calculating the GC content of GC rich fragment

sequence_split.py is used for generating the oligos for 100 fragments


100genebarcode.R: Following the ligation of barcodes to the termini of the 100-gene pool, high-throughput sequencing was performed on the barcode regions of the pooled plasmid library. This script processes the resulting sequencing data to decode the barcode associations.

FP_library_code.R and PETase.R: These scripts are employed for the analysis and quantification of mutation spectra, specifically to identify the varieties and respective frequencies of observed mutations within the library.

barcode_wellandplate_position.py: This script is designed to process the sequencing data derived from the PCR amplicons collected across 50 multi-well plates. Its function is to map each sequencing result back to its specific physical location, thereby enabling the identification of the gene contained within any given well of a specific plate.

High_throughput_error_free.py: This script analyzes the sequencing results from a high-throughput gene synthesis campaign. It is utilized to identify which target gene sequences are represented by perfectly matching, error-free reads and to determine the abundance (frequency) of such reads for each gene.

mutation_less_than_10_part1.py and mutation_less_than_10_part2.py: These scripts perform a complementary analysis on the subset of genes for which no error-free reads were identified. They determine the prevalence of reads containing a low number of mutations (specifically, between 1 and 10 mutations) for these genes, providing a metric for near-correct synthesis
