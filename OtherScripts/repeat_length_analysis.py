#Use this script to create a both a distribution of repeat lengths in the human genome, as well as create an output file containing their values.

from pylab import hist, show, xticks
import numpy as np
import matplotlib.pyplot as plt
import math

with open('hg19_Repeatmasker_UCSC.gff') as f:
    content = f.readlines()
    
repeat_lengths = []
acceptable=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

for i in range(1, len(content)):
    chr_num = content[i].split()[5]
    if chr_num in acceptable:
        genome_start = int(content[i].split()[6])
        genome_end = int(content[i].split()[7])
        length = genome_end-genome_start
        repeat_lengths.append(length)

print('min: ' + str(min(repeat_lengths)))
print('max: ' + str(max(repeat_lengths)))
print('median: ' + str(np.median(repeat_lengths)))
print('mean: ' + str(np.mean(repeat_lengths)))
print('std dev: ' + str(np.std(repeat_lengths)))

bins = []
start = 0
bin_value = 5
num_bins = 200
for i in range(0, num_bins):
    bins.append(start)
    start = start + bin_value
    
plt.hist(repeat_lengths, bins, density = True, stacked = True, alpha = 0.5, ec = 'black')
plt.gca().set(xlabel = 'Repeat lengths', ylabel = 'Probability (divided by bin width)');
plt.axis([0, 1000, 0, 0.007])
plt.savefig('C:/Users/kelle/Downloads/repeat_lengths_distribution.pdf')
show()

repeat_lengths_file = open('C:/Users/kelle/OneDrive/slim_folder/raw_repeat_lengths.txt', 'w+')
for i in repeat_lengths:
    repeat_lengths_file.write(str(i) + '\n')
repeat_lengths_file.close()
