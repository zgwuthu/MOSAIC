GCwindow.py is used for calculating the GC content of GC rich fragment

sequence_split.py is used for generating the oligos for 100 fragments


100genebarcode.R: Following the ligation of barcodes to the termini of the 100-gene pool, high-throughput sequencing was performed on the barcode regions of the pooled plasmid library. This script processes the resulting sequencing data to decode the barcode associations.
Barcode regions were recognized and extracted by surrounded conservative region (TGCAAGAGACTTCCATCCAG and GCTTCCGGTCTGGTTCGCTTTGAAGCTCGA). Occurrences of each barcode would be counted in the extracted results. 


FP_library_code.R and PETase.R: These scripts are employed for the analysis and quantification of mutation spectra, specifically to identify the varieties and respective frequencies of observed mutations within the library.
12 variant regions for FP and 8 variant regions for PETase were strictly recognized (without mismatch allowance) and extracted by sequences of conservative region. 
For FP library: 
Position		Conservative Region						Recognized Position
1, 2, 3		GGCGTGCAGTGCTTC					-9 to -1
4, 5			CGCTACCCCGACCACATG				-3 to -1, +1 to +3
6, 7, 8		ATCCTGGGGCACAAGCTGGAGTACAAC	+1 to +6, +13 to +15
9, 10, 11		GCCGACAAGCAGAAGAACGGCATCAAG	-3 to -1, +1 to +3, +13 to +15
12			CCGACAACCACTACCTGAGC				+1 to +3
For PETase library: 
Position		Conservative Region						Recognized Position
1			GTTCTTATGGTGTAGTAAGT				+1 to +3
2			CCAATGCTACGGAGCGTTTT				+1 to +3
3			GTGCTCTCGATTGGGCGTCA				+1 to +3
4			TAGCACCTTGGCATACGACT				-3 to -1
5			TTTTGGGTGGTCAAAATGAT				-3 to -1
6			ATTGCCCCTGTCTCT						-3 to -1
7			CACATGCAATTCCAATGTAT				-3 to -1
8			GCCATAATTTTCCGAATTCG				-3 to -1
Reads would be ignored if there was no matching for some of conservative regions. Variant regions were translated to amino acid sequences. And occurrences of each theoretical combination would be counted in the results. 


barcode_wellandplate_position.py: This script is designed to process the sequencing data derived from the PCR amplicons collected across 50 multi-well plates. Its function is to map each sequencing result back to its specific physical location, thereby enabling the identification of the gene contained within any given well of a specific plate.

High_throughput_error_free.py: This script analyzes the sequencing results from a high-throughput gene synthesis campaign. It is utilized to identify which target gene sequences are represented by perfectly matching, error-free reads and to determine the abundance (frequency) of such reads for each gene.

mutation_less_than_10_part1.py and mutation_less_than_10_part2.py: These scripts perform a complementary analysis on the subset of genes for which no error-free reads were identified. They determine the prevalence of reads containing a low number of mutations (specifically, between 1 and 10 mutations) for these genes, providing a metric for near-correct synthesis
