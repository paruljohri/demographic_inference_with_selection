#This script is designed to mask the centromere region of genomes simulated in accordance with characteristics of human chromosome 6.
#The centromere region of human chromosome 6 ranges roughly from base pair positions 58500000 - 62500000. However, human chromosome 6 is about 171Mb long. 
#Given that the genomes simulated in this project are only 150Mb long, the centromere region for the purposes of this masking procedure was defined to range from base pair positions 48500000 - 52500000 instead. 

import numpy as np
import sys

size = int(sys.argv[1]) #Size of chromosome in base pairs as integer
inputt = sys.argv[2] #Input path to .ms file to have the centromere masked
outputt = sys.argv[3] #Output path to newly masked .ms file

#chromosome_6_size = 170805979
#chromosome_6_centromere_start_position=58500000
#chromosome_6_centromere_end_position=62500000

masking_position_start = 48500000 #Specifies the starting position of the centromere region
masking_position_end = 52500000 #Specifies the ending position of the centromere region

with open(inputt) as f:
    content = f.readlines()

segsites = int(content[1].strip('segsites: ').strip('\n'))
print('segsites: ' + str(segsites))
SNP_positions = [round(size * float(i)) for i in content[2].strip('positions: ').split()]
del content[0:3]
all_haplotypes = [haplotype.strip('\n') for haplotype in content]


SNP_positions_to_take = []
for i in range(0, len(SNP_positions)):
    if masking_position_start < SNP_positions[i] < masking_position_end:
        pass
    else:
        SNP_positions_to_take.append(i)

new_positions = []
for q in SNP_positions_to_take:
    position = SNP_positions[q] / size
    position_formatted = '{:.8f}'.format(position)
    new_positions.append(position_formatted)

masked_haplotypes = []
for haplotype in all_haplotypes:
    masked = ''
    for q in SNP_positions_to_take:
        value = haplotype[q]
        masked = masked + str(value)
    masked_haplotypes.append(masked)

a = len(all_haplotypes[0])
b = len(masked_haplotypes[0])
print('____________________________')
print(new_positions[0:100])
print(SNP_positions_to_take[0:100])
print(a)
print(b)

num_sites_masked = a-b
percentage_masked = num_sites_masked / segsites * 100
print('An additional ' + str(num_sites_masked) + ' sites in the centromere were masked which is ' + str(round(percentage_masked, 2)) + '% of all sites')

masked_ms_file = open(outputt, 'w+')
masked_ms_file.write('//' + '\n')
masked_ms_file.write('segsites: ' + str(b) + '\n')
masked_ms_file.write('positions: ' + ' '.join(map(str, new_positions)) + '\n')
for sequence in masked_haplotypes:
    masked_ms_file.write(sequence + '\n')
masked_ms_file.close()
