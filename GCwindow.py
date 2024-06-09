import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

# the GC content is defined as the ratio of G and C in a DNA sequence 
def gc_content(seq):
    return (seq.count('g') + seq.count('c')) / len(seq)

def sliding_window_gc(seq, window_size):
    gc_contents = []
    for i in range(0, len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        gc_contents.append(gc_content(window))
    return gc_contents

def plot_gc_content(gc_contents,):
    plt.plot(gc_contents, color = 'black')
    plt.xlabel('Position in Sequence', font = 'Times New Roman')
    plt.ylabel('GC Content',font = 'Times New Roman')
    plt.title('GC Content Sliding Window',font = 'Times New Roman')
    plt.show()

# read DNA sequence file
with open('sequence.txt', 'r') as file:
    sequence = file.read().replace('\n', '')

# calculate the GC content
gc_contents = sliding_window_gc(sequence, window_size=200)

# plot the GC content figure
plot_gc_content(gc_contents)
