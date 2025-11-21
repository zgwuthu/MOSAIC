import random


content = ''
while len(content) < 1001:
        nucleotide = random.choice(['A', 'T'])
        content += nucleotide

output = open('output.txt', mode='a', encoding='utf-8')
print(content, file=output)
output.close()