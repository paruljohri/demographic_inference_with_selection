#This script is designed to randomly mask an additional percentage of a genome.
#For the purposes of this project, an additional 10% of the genome was masked, though this is customizable.
#The length of any particular new masked region is randomly drawn from a distribution of repeat lengths.
#This length is then applied to some random region of the genome.
#This script accounts for the fact that the exons in the input genome have already been masked,
#and ensures that any newly masked areas do not overlap with them.
#Because the centromere region is also masked in another script,
#this script does not allow for any of the random masking to be applied to this pre-defined centromere region.

import numpy as np
import random
import sys



### Recieve values and file pathways
size = int(sys.argv[1]) #Size of chromosome in base pairs as integer
percentage_to_mask = float(sys.argv[2]) #Percentage of genome to be masked as a decimal from 0 to 1 (i.e. write in '0.25' if you desire 25% of the genome to be masked)
inputt = str(sys.argv[3]) #Input path to .ms file to have random masking performed on
outputt = str(sys.argv[4]) #Output path to newly masked .ms file
genomeX_positions_file = str(sys.argv[5]) #File containing all exon positions in the genome
repeat_lengths_file = str(sys.argv[6]) #File containing all repeat region lengths



### Read in the chromosome to mask
with open(inputt) as f:
    content = f.readlines()
segsites = int(content[1].strip('segsites: ').strip('\n'))
SNP_positions = [round(size * float(i)) for i in content[2].strip('positions: ').split()]
SNP_positions_to_take = list(SNP_positions)
del content[0 : 3]
all_haplotypes = [haplotype.strip('\n') for haplotype in content]
content.clear()



### Mask windows from previous exon masking
with open(genomeX_positions_file, 'r+') as f:
    genomeX_positions_content = f.readlines()

masking_windows = []
for line in genomeX_positions_content:
    if 'g3' in line:
        windows = line[28 : -3]
        comma_position = windows.find(',')
        window_lower_bound = int(windows[0 : comma_position]) + 1
        window_upper_bound = int(windows[comma_position + 1 : ]) + 1
        true_window=[window_lower_bound, window_upper_bound]
        masking_windows.append(true_window)

        for position in SNP_positions_to_take:
            if window_lower_bound <= position < window_upper_bound:
                SNP_positions_to_take.remove(position)
genomeX_positions_content.clear()



### Read in all raw repeat lengths
repeat_lengths = []
with open(repeat_lengths_file, 'r+') as f:
    repeat_lengths_content = f.readlines()
for length in repeat_lengths_content:
    repeat_lengths.append(int(length))
repeat_lengths_content.clear()



### Define centromere boundaries
chromosome_6_size = 170805979
ratio = size / chromosome_6_size
# centromere_start_position = 58500000
# centromere_end_position = 62500000
centromere_window_start = 48500000
centromere_window_end = 52500000



### Find the new windows to mask
total_length_to_mask = size * percentage_to_mask
lower_end_tolerance = 0.9999
upper_end_tolarance = 1.0001
lower_end_length = round(total_length_to_mask * lower_end_tolerance)
upper_end_length = round(total_length_to_mask * upper_end_tolarance)

individual_masking_lengths = []
criteria_fulfilled = False
while criteria_fulfilled == False:
    current_masked_length = sum(individual_masking_lengths)
    if lower_end_length <= current_masked_length <= upper_end_length:
        criteria_fulfilled = True
    else:
        suitable_position_found = False
        while suitable_position_found == False:
            current_position = np.random.randint(1, size)
            length = random.choice(repeat_lengths)
            print('____________________________________')
            print('current position: ' + str(current_position))
            print('length: ' + str(length))

            if (current_position + length) > size:
                print('new window goes past chromosome tip')
                print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                continue

            if (current_masked_length + length) > upper_end_length:
                print('new length makes masking % too much')
                print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                continue

            if (centromere_window_start <= current_position < centromere_window_end):
                print('starting window position is inside centromere')
                print('interfering with centromere ' + str(centromere_window_start) + ' - ' + str(centromere_window_end))
                print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                continue
            elif (centromere_window_start <= current_position + length < centromere_window_end):
                print('new window cuts into centromere')
                print('interfering with centromere ' + str(centromere_window_start) + ' - ' + str(centromere_window_end))
                print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                continue
            elif (current_position < centromere_window_start < current_position + length) and (current_position < centromere_window_end < current_position + length):
                print('new window overlaps entirety of centromere')
                print('interfering with centromere ' + str(centromere_window_start) + ' - ' + str(centromere_window_end))
                print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                continue

            window_interference = False
            for window in masking_windows:
                window_lower_bound = int(window[0])
                window_upper_bound = int(window[1])
                if (window_lower_bound <= current_position < window_upper_bound):
                    print('starting window position is inside a previous boundary')
                    print('interfering with window ' + str(window_lower_bound) + ' - ' + str(window_upper_bound))
                    print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                    window_interference = True
                    break
                elif (window_lower_bound <= current_position + length < window_upper_bound):
                    print('new window cuts into a previous boundary')
                    print('interfering with window ' + str(window_lower_bound) + ' - ' + str(window_upper_bound))
                    print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                    window_interference = True
                    break
                elif (current_position < window_lower_bound < current_position + length) and (current_position < window_upper_bound < current_position + length):
                    print('new window overlaps entirety of a previous boundary')
                    print('interfering with window ' + str(window_lower_bound) + ' - ' + str(window_upper_bound))
                    print('would\'ve made new window ' + str(current_position) + ' - ' + str(current_position + length))
                    window_interference = True
                    break
            if window_interference == False:
                suitable_position_found = True
                print('current masked length: ' + str(current_masked_length))

        window_lower_bound = current_position
        window_upper_bound = current_position + length
        window = [window_lower_bound, window_upper_bound]
        masking_windows.append(window)

        difference = window_upper_bound - window_lower_bound
        individual_masking_lengths.append(difference)

        for position in SNP_positions_to_take:
            if window_lower_bound <= position < window_upper_bound:
                SNP_positions_to_take.remove(position)

        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')

### Take new segsite positions
new_positions = []
for q in SNP_positions_to_take:
    new_position = q / size
    new_position_formatted = '{:.8f}'.format(new_position)
    new_positions.append(new_position_formatted)



### Mask the haplotypes
SNP_indexes_to_take = []
for i in SNP_positions_to_take:
    index = SNP_positions.index(i)
    SNP_indexes_to_take.append(index)

masked_haplotypes = []
for haplotype in all_haplotypes:
    masked = ''
    for r in SNP_indexes_to_take:
        value = haplotype[r]
        masked = masked + str(value)
    masked_haplotypes.append(masked)



### Results of how many sites were masked
print('\n')
print('RESULTS')
print('total masking length: ' + str(current_masked_length))
# print('windows: ')
# print(masking_windows)
initial_num_segsites = len(all_haplotypes[0])
final_num_segsites = len(masked_haplotypes[0])
num_sites_masked = initial_num_segsites - final_num_segsites
percentage_masked = num_sites_masked / segsites * 100
print('initial # of segsites: ' + str(initial_num_segsites))
print('final # of segsites: ' + str(final_num_segsites))
print(str(num_sites_masked) + ' sites were masked which is ' + str(round(percentage_masked, 2)) + '% of all sites')



### Create new masked .ms file
masked_ms_file=open(outputt, 'w+')
masked_ms_file.write('//' + '\n')
masked_ms_file.write('segsites: ' + str(final_num_segsites) + '\n')
masked_ms_file.write('positions: ' + ' '.join(map(str, new_positions)) + '\n')
for sequence in masked_haplotypes:
    masked_ms_file.write(sequence + '\n')
masked_ms_file.close()


### Data printout
# print(SNP_positions)
# print(SNP_positions_to_take)
# print(new_positions)
# print(len(SNP_positions))
# print(len(all_haplotypes[0]))
# print(all_haplotypes[0])
# print(len(SNP_positions_to_take))
# print(len(masked_haplotypes[0]))
# print(masked_haplotypes[0])
