#This script creates an SFS from a single chromosome segment using only SNPs a minimum distance away from one another

import sys

input_chromosome = sys.argv[1] #path to the input .ms file
size = int(sys.argv[2]) #size of the input genome in base pairs
all_SNPs_ratio_file_name = sys.argv[3] #path and name of the ratio file
SFS_file_name = sys.argv[4] #path and name of the output SFS file
distance_apart = sys.argv[5] #minimum distance between SNPs counted towards the SFS in base pairs

num_chr = 1 #number of input genomes
num_haplotypes = 100 #number of haploytypes in the .ms file

total_SFS = [0] * (num_haplotypes + 1)
with open(input_chromosome) as f:
    content = f.readlines()

segsites = int(content[1].strip('segsites: ').strip('\n'))

SNP_positions = [round(size * float(i)) for i in content[2].strip('positions: ').split()]

start = SNP_positions[0]
positions_to_take=[0]
for i in range(1, segsites):
    difference = SNP_positions[i] - start
    if difference >= distance_apart:
        start = SNP_positions[i]
        positions_to_take.append(i)    
positions_to_take = [(g + 1) for g in positions_to_take]

del content[0:3]
all_haplotypes = [haplotype.strip('\n') for haplotype in content]

sampled_haplotypes = []
for individual in all_haplotypes:
    build_individual = ''
    for q in positions_to_take:
        sampled_SNP = individual[q - 1]
        build_individual = build_individual + sampled_SNP
    sampled_haplotypes.append(build_individual)

ntons = [0] * (len(positions_to_take))
for individual in sampled_haplotypes:
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
SFS[-1] = 0
print(SFS)

for w in range(0, len(SFS)):
    total_SFS[w] = total_SFS[w] + SFS[w]
print(total_SFS)
        
with open(all_SNPs_ratio_file_name) as g:
        ratio = g.readlines()
all_SNPs_ratio = float(ratio[0])

total_non_zero_type = sum(total_SFS[1:])
total_zero_type = total_non_zero_type * all_SNPs_ratio / (1 - all_SNPs_ratio)
total_SFS[0] = int(total_zero_type)    
        
build_d0 = ''
build_info = ''
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
