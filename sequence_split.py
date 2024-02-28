import random


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


file_path = 'sequence.fasta'
sequences = read_fasta_file(file_path)
output = open('output.txt', mode='a', encoding='utf-8')

for name, content in sequences.items():
    if len(content) > 1590:
        print(name, 'is longer than 1590bp')
    elif len(content) == 1590:
        content = content
    while len(content) < 1590:
        nucleotide = random.choice(['A', 'G', 'C', 'T'])
        content += nucleotide

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

