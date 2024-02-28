import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')


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

# 读取DNA序列文件
with open('sequence.txt', 'r') as file:
    sequence = file.read().replace('\n', '')

# 计算滑动窗口GC含量
gc_contents = sliding_window_gc(sequence, window_size=200)

# 生成GC含量滑窗图
plot_gc_content(gc_contents)
