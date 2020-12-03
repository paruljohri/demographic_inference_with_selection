import math
import itertools
import numpy as np

def analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size): 
    def find_dashed_line_start_index(max_plotted_generation, bin_size):
        all_replicate_max_generation_boundaries=[]
        for k in range(1,number_of_msmc_replicates+1):
            msmc_results_file = '$$$/SupplementalInformation/MSMC_OutputFiles/transformed//TRANSFORMED_rep' + str(k) + '_' + demographic_model + '_DFE' + str(DFE_number) + '_' + genome_percentage + '_' + exon_state + '.txt'
            with open(msmc_results_file) as f:
                content = f.readlines()[-1]
            final_line=content.split()
            all_replicate_max_generation_boundaries.append(float(final_line[1]))
        max_replicate_generation = max(all_replicate_max_generation_boundaries)
        boom=(max_plotted_generation/bin_size)+1
        bop=int(math.floor(max_replicate_generation / 100)) * 100
        bing=(max_plotted_generation-bop) / bin_size
        dashed_line_start_index=int(boom-bing)+1
        return dashed_line_start_index
        
    def find_indeces(array, bin_size, bins):
        array_start_and_end = []
        proportions_list=[]
        start = array[0]
        end = array[1]
        for x in range(0, len(bins)-1):
            if bins[x] <= start and start < bins[x+1]:
                start_index = int(start / bin_size)
                end_index = int(end / bin_size)
                index_range = end_index - start_index;
                array_start_and_end.append(start_index)
                if start_index == end_index:
                    proportions_list=[[start_index, 1.0]]
                    return proportions_list
                if index_range > 1:
                    start_point = start_index + 1
                    while start_point < end_index:
                        array_start_and_end.append(start_point)
                        start_point = start_point + 1
                array_start_and_end.append(end_index)
                for r in range(0, len(array_start_and_end)):
                    proportions=[]
                    proportions.append(array_start_and_end[r])
                    lower_proportion=(array_start_and_end[1] * bin_size - start) / bin_size
                    upper_proportion=(end - array_start_and_end[-1] * bin_size) / bin_size
                    if r == 0:
                        proportions.append(lower_proportion)
                    elif r == len(array_start_and_end) - 1:
                        proportions.append(upper_proportion)
                    else:
                        proportions.append(1.0)
                    proportions_list.append(proportions)
                return proportions_list
                
    bins=[a for a in range(0, max_plotted_generation + 2 * bin_size, bin_size)]
    msmc_average_population_sizes = []
    msmc_replicate_population_sizes = []
    dashed_line_start_index = find_dashed_line_start_index(max_plotted_generation, bin_size) 
    
    for b in range(1, number_of_msmc_replicates + 1):
        msmc_results_file = '$$$/SupplementalInformation/MSMC_OutputFiles/transformed//TRANSFORMED_rep' + str(b) + '_' + demographic_model + '_DFE' + str(DFE_number) + '_' + genome_percentage + '_' + exon_state + '.txt'
        with open(msmc_results_file) as f:
                content = f.readlines()[1 : ]
        data = []
        for line in content:
            line_values = line.split()
            left_generation_boundary = int(line_values[1])
            if line_values[2] == 'inf':
                right_generation_boundary = bins[-1] - bin_size
            else:
                right_generation_boundary = int(line_values[2])
            population_size = int(line_values[3])
            values = [left_generation_boundary, right_generation_boundary, population_size]
            data.append(values)
        bin_totals = [0] * (len(bins) - 1)
        for d in data:
            data_indeces = find_indeces(d, bin_size, bins)
            for index in data_indeces:
                bin_totals[index[0]] = bin_totals[index[0]] + d[2] * index[1]
        msmc_replicate_population_sizes.append(bin_totals)
    del bins[-1]
    
    for a in range(0,len(msmc_replicate_population_sizes[0])):
        average_per_generation=float(np.mean([i[:][a] for i in msmc_replicate_population_sizes]))
        msmc_average_population_sizes.append(average_per_generation)
    
    msmc_average_population_sizes[-1] = msmc_average_population_sizes[-2]
    
    return bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index

graph_demographic_model = ['eqm', 'decline', 'growth']
graph_exon_state = ['unmasked', 'masked']
graph_genome_precentage = ['genome05', 'genome10', 'genome20', 'genome05_2fold', 'genome10_2fold', 'genome20_2fold']
graph_DFE = [0, 1, 2, 3, 4, 5, 6]

number_of_individuals_sampled = 2
bin_size = 100

number_of_msmc_replicates = 10

parameter_combinations = list(itertools.product(graph_demographic_model, graph_exon_state, graph_genome_precentage, graph_DFE))
for combo in parameter_combinations:
    demographic_model = combo[0]
    exon_state = combo[1]
    genome_percentage = combo[2]
    DFE_number = combo[3]
    if demographic_model == 'eqm' and ('2fold' in genome_percentage):
        continue
    
    if '2_fold' not in genome_percentage:
        if demographic_model == 'eqm':
            max_plotted_generation = 100000       
        if demographic_model == 'growth':
            max_plotted_generation = 25000
        if demographic_model == 'decline':
            max_plotted_generation = 100000 
    else:
        if demographic_model == 'growth':
            max_plotted_generation = 25000
        if demographic_model == 'decline':
            max_plotted_generation = 75000
            
    bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index = analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size)
    
    average_file_location = 'C:/Users/kelle/Downloads/final/new_average_msmc_outputs/AVERAGE_'+demographic_model+'_DFE'+str(DFE_number)+'_'+genome_percentage+'_'+exon_state+'.txt'
    average_file = open(average_file_location, "w+")
    average_file.write('time_index\tleft_generation_boundary\tright_generation_boundary\tN\n')
    
    x = 0
    for i in range(0, dashed_line_start_index):
        line = str(x) + '\t' + str(bins[i]) + '\t' + str(bins[i + 1]) + '\t' + str(int(msmc_average_population_sizes[i])) + '\n'
        average_file.write(line)
        x += 1
    line = str(x) + '\t' + str(dashed_line_start_index * bin_size) + '\t' + 'inf' + '\t' + str(int(msmc_average_population_sizes[-1])) + '\n'
    average_file.write(line)    
    average_file.close()