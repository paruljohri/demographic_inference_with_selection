#This script creates an SFS from a single chromosome segment using all available SNPs

import sys

input_chromosome = sys.argv[1] #path to the input .ms file
size = int(sys.argv[2]) #size of the input genome in base pairs
all_SNPs_ratio_file_name = sys.argv[3] #path and name of the ratio file
SFS_file_name = sys.argv[4] #path and name of the output SFS file

num_chr = 1 #number of input genomes
num_haplotypes = 100 #number of haploytypes in the .ms file

total_SFS = [0] * (num_haplotypes + 1)
with open(input_chromosome) as f:
    content = f.readlines()

segsites = int(content[1].strip('segsites: ').strip('\n'))

del content[0:3]

all_haplotypes = [haplotype.strip('\n') for haplotype in content]

ntons = [0] * segsites
for individual in all_haplotypes:
    for d in range(0, len(individual)):
        site = int(individual[d])
        if site == 1:
            ntons[d] += 1

SFS = [0]
max_nton = max(ntons)
for e in range(1, max_nton + 1):
    count = ntons.count(e)
    SFS.append(count)
total = sum(SFS)
SFS[0] = int(size) - total
SFS[-1] = 0
print(SFS)

for w in range(0, len(SFS)):
    total_SFS[w] = total_SFS[w] + SFS[w]
print(total_SFS)
        
all_SNPs_ratio = total_SFS[0] / sum(total_SFS)
ratio = open(all_SNPs_ratio_file_name, 'w+')
ratio.write(str(all_SNPs_ratio))
ratio.close()
        
build_d0=''
build_info=''
for f in range(0, len(total_SFS)):
    d0 = 'd0_' + str(f)
    build_d0 = build_d0 + d0 + '\t'
    info = str(total_SFS[f])
    build_info = build_info + info + '\t'

edit = open(SFS_file_name, 'w+')
edit.write('1 observations\n')
edit.write(build_d0 + '\n')
edit.write(build_info)
edit.close()
