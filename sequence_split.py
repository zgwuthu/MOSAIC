# This script was utilized to standardize the length of the gene sequence via padding 
# and subsequent cleavage into uniform-length oligonucleotides.
# sequence.fasta is the sequence file of target genes

# In this specific demonstration, we showcase the processing of 100 genes, 
# all initially shorter than 1590 bp. These genes were first padded to a uniform length of 1590 bp. 
# Subsequently, they were flanked by two homologous ends (left and right) and then segmented into component oligonucleotides 
# at a fixed length of 60 bp per oligo.

import random

# define a function to read the fasta file
def read_fasta_file(file_path):
    sequences = {}
    current_sequence_name = ''
    current_sequence_content = ''

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence_name != '':
                    sequences[current_sequence_name] = current_sequence_content
                current_sequence_name = line[1:]
                current_sequence_content = ''
            else:
                current_sequence_content += line

    if current_sequence_name != '':
        sequences[current_sequence_name] = current_sequence_content

    return sequences

# read the fasta file
file_path = 'sequence.fasta'
sequences = read_fasta_file(file_path)
output = open('output.txt', mode='a', encoding='utf-8')

# the length here can be changed as demanded 
for name, content in sequences.items():
    if len(content) > 1590:
        print(name, 'is longer than 1590bp')
        continue
    elif len(content) == 1590:
        content = content
    # if this fragment is shorter than 1590 bp, it would be padded by nonsense sequence
    while len(content) < 1590:
        nucleotide = random.choice(['A', 'G', 'C', 'T'])
        content += nucleotide
    
    # these two sequences are two homologous regions appended two ends
    content = 'GCAATGCAGACTCAGAGAGAACCCGCCACC' + content.upper() + 'GCTTCCGGTCTGGTTCGCTTTGAAGCTCGA'
    complementary_sequence = ''
    ture_complementary_sequence = ''
    for base in content:
        if base == 'A':
            complementary_sequence += 'T'
        elif base == 'T':
            complementary_sequence += 'A'
        elif base == 'C':
            complementary_sequence += 'G'
        elif base == 'G':
            complementary_sequence += 'C'
    ture_complementary_sequence = complementary_sequence[::-1]

    # the 1650-bp sequence is segmented to 54 60-mers and 2 30-mers
    a = 30
    for i in range(27):
        b = a + 60
        in_put = content[a:b]
        k = i + 1,
        print(name, "_F", "%s" % k, '\t', "%s" % in_put, sep='', file=output)
        a += 60

    a = 30
    for i in range(27):
        b = a + 60
        in_put = ture_complementary_sequence[a:b]
        k = i + 1,
        print(name, "_R", "%s" % k, '\t', "%s" % in_put, sep='', file=output)
        a += 60
output.close()

