#This script allows the user to graph the results of MSMC runs on chromosomal segments.
#The scenario parameters this script will be run for can be altered towards the bottom by modifiyng the parameters in the the lists with a beginning with "GRAPH".
#The current parameters already specified showcase all possible parameters that can be entered for each list. 
#Running it as such would produce all possible graphs. 
#Areas requiring a file path alteration are marked with '$$$'

import math
import itertools
import numpy as np
import matplotlib.pyplot as plt

#This function builds the true population size line (in black), specifies the maximum generation that replicate lines are plotted out to backward in time ('max_plotted_generation'), and specifies the y axis range in terms of population size ('min_viewable_population_size' and 'max_viewable_population_size'). 
def demographic_model_parameters(chr_size):
    true_generations = []
    true_population_sizes = []
        
    min_viewable_population_size = 500
    if chr_size == '1Mb':
        max_viewable_population_size = 100000
    if chr_size == '10Mb':
        max_viewable_population_size = 100000
    if chr_size == '50Mb':
        max_viewable_population_size = 100000
    if chr_size == '200Mb':
        max_viewable_population_size = 100000
    if chr_size == '1Gb':
        max_viewable_population_size = 100000
        
    max_plotted_generation = 75000
    for generation in range(0, max_plotted_generation):
        true_generations.append(generation)
        size = 5000
        true_population_sizes.append(size)
                
    return true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size

#This function analyzes the MSMC results and builds output lines for each replicate, as well as an average output line
def analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size, chr_size, num_of_individuals, mutation_rate, simulation_tool):         
    #This sub-function determines what values are added to each bin when creating the average MSMC output line
    def find_indeces(array, bin_size, bins):
        array_start_and_end = []
        proportions_list = []
        start = array[0]
        end = array[1]
        for x in range(0, len(bins) - 1):
            if bins[x] <= start and start < bins[x + 1]:
                start_index = int(start / bin_size)
                end_index = int(end / bin_size)
                index_range = end_index - start_index;
                array_start_and_end.append(start_index)
                if start_index == end_index:
                    proportions_list = [[start_index, 1.0]]
                    return proportions_list
                if index_range > 1:
                    start_point = start_index + 1
                    while start_point < end_index:
                        array_start_and_end.append(start_point)
                        start_point = start_point + 1
                array_start_and_end.append(end_index)
                for r in range(0, len(array_start_and_end)):
                    proportions = []
                    proportions.append(array_start_and_end[r])
                    lower_proportion = (array_start_and_end[1] * bin_size - start) / bin_size
                    upper_proportion = (end - array_start_and_end[-1] * bin_size) / bin_size
                    if r == 0:
                        proportions.append(lower_proportion)
                    elif r == len(array_start_and_end) - 1:
                        proportions.append(upper_proportion)
                    else:
                        proportions.append(1.0)
                    proportions_list.append(proportions)
                return proportions_list
                
    bins = [a for a in range(0, max_plotted_generation + 2 * bin_size, bin_size)]
    msmc_average_population_sizes = []
    msmc_replicate_population_sizes = []
    
    all_replicate_max_generation_boundaries = []
    for b in range(1, number_of_msmc_replicates + 1):
        if simulation_tool == 'slim':
            msmc_results_file = '$$$/slim/raw_output_single_chr_' + chr_size + '_' + num_of_individuals + 'indv_rep' + str(b)+'.final.txt'
        if simulation_tool == 'msprime':
            msmc_results_file = '$$$/msprime/msmc1_run_' + chr_size + '_' + num_of_individuals + 'indv_' + str(b) + '.final.txt'
        with open(msmc_results_file) as f:
            content = f.readlines()[-1]
        final_line = content.split()
        all_replicate_max_generation_boundaries.append(int(float(final_line[1]) / mutation_rate ))
    max_replicate_generation = max(all_replicate_max_generation_boundaries)
    boom = (max_plotted_generation / bin_size) + 1
    bop = int(math.floor(max_replicate_generation / 100)) * 100
    bing = (max_plotted_generation - bop) / bin_size
    dashed_line_start_index = int(boom - bing) + 1
    
    for b in range(1, number_of_msmc_replicates + 1):
        if simulation_tool == 'slim':
            msmc_results_file = '$$$/slim/raw_output_single_chr_' + chr_size + '_' + num_of_individuals + 'indv_rep' + str(b)+'.final.txt'
        if simulation_tool == 'msprime':
            msmc_results_file = '$$$/msprime/msmc1_run_' + chr_size + '_' + num_of_individuals + 'indv_' + str(b) + '.final.txt'
        with open(msmc_results_file) as f:
            content = f.readlines()[1 : ]
        data = []
        for line in content:
            line_values = line.split()
            left_generation_boundary = int(float(line_values[1]) / mutation_rate)
            if line_values[2] == 'inf':
                right_generation_boundary = bins[-1] - bin_size
            else:
                right_generation_boundary = int(float(line_values[2]) / mutation_rate)
            population_size = int(1 / (float(line_values[3]) * 2 * mutation_rate))
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
        average_per_generation = float(np.mean([i[:][a] for i in msmc_replicate_population_sizes]))
        msmc_average_population_sizes.append(average_per_generation)
    
    msmc_average_population_sizes[-1] = msmc_average_population_sizes[-2]
    for j in range(0, len(msmc_replicate_population_sizes)):
        msmc_replicate_population_sizes[j][-1] = msmc_replicate_population_sizes[j][-2]
    
    return bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index

#This function graphs the data for MSMC depending on one's specifications
def graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, plot_msmc_average_line, plot_msmc_replicate_lines, chr_size, num_of_individuals, simulation_tool, save_graph_as_pdf, display_graph):
    fig = plt.figure(figsize = (20, 11.3))
    ax = fig.add_subplot(111)
    
#     The x and y axis parameters can be customized here
    plt.rc('xtick', labelsize = 50)
    ax.xaxis.set_tick_params(length = 8, width = 3)
    ax.xaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_xlabel('Number of Generations', fontsize = 47)
    ax.set_xscale('log')
    ax.set_xlim([50, max_plotted_generation])
    
    plt.rc('ytick', labelsize = 50)
    ax.yaxis.set_tick_params(length = 8, width = 3)
    ax.yaxis.set_tick_params(length = 5, width = 2, which = 'minor')
    ax.set_ylabel('${N}$', fontsize = 47)
    ax.set_yscale('log')
    ax.set_ylim([min_viewable_population_size, max_viewable_population_size])
    
    
    
    ax.plot(true_generations, true_population_sizes, color = (0, 0, 0, 1), linewidth = 3.0)
    
    if plot_msmc_average_line == True:
        if dashed_line_start_index > 0:
            ax.plot(bins[ : dashed_line_start_index], msmc_average_population_sizes[ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 1), linewidth = 2.0)
            ax.plot(bins[dashed_line_start_index - 1 : ], msmc_average_population_sizes[dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 1), linewidth = 2.0)
        else:
            ax.plot(bins, msmc_average_population_sizes, ls = 'steps', color = (1, 0, 0, 1), linewidth = 3.0)
    if plot_msmc_replicate_lines == True:
        for c in range(0, number_of_msmc_replicates):
            if dashed_line_start_index > 0:
                ax.plot(bins[ : dashed_line_start_index], msmc_replicate_population_sizes[c][ : dashed_line_start_index], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
                ax.plot(bins[dashed_line_start_index - 1 : ], msmc_replicate_population_sizes[c][dashed_line_start_index - 1 : ], ls = 'dashed', color = (1, 0, 0, 0.2), linewidth = 2.0)
            else:
                ax.plot(bins, msmc_replicate_population_sizes[c], ls = 'steps', color = (1, 0, 0, 0.2), linewidth = 2.0)
    
    if save_graph_as_pdf == True:
        fig.savefig('$$$/graph_msmc_' + chr_size + '_' + num_of_individuals + 'indv.pdf')
    
    if display_graph != True:
        plt.close()

#This section is used to specify which scenarios to create results figures for. All options currently listed encompass the full range of possible options
GRAPH_chr_size = ['1Mb','10Mb','50Mb','200Mb','1Gb'] #Specify the chromosome segment size as a string
GRAPH_num_of_individuals = [2,4,8] #Specify the number of indivudals MSMC was run with

simulation_tool = 'slim' #Specify either 'slim' or 'msprime' for which results will be analyzed

bin_size = 100 #Specify the binsize (in generations) that will be used in the creation of the MSMC average line
mutation_rate = 10**-8 #Specify the mutation rate the genomes were simulated under for the purposes of converting the MSMC values 

plot_msmc_replicate_lines = True #Specify either 'True' or 'False' for whether or not to include MSMC replicate results
plot_msmc_average_line = True #Specify either 'True' or 'False' for whether or not to include an average line of MSMC replicate results in the plot
display_graph = True #Specify either 'True' or 'False' for whther or not to display the plots
save_graph_as_pdf = True #Specify either 'True' or 'False' for whther or not to save the plots

parameter_combinations = list(itertools.product(GRAPH_chr_size, GRAPH_num_of_individuals))
for combo in parameter_combinations:
    chr_size = combo[0]
    num_of_individuals = str(combo[1])
    
    if simulation_tool == 'msprime':
        number_of_msmc_replicates = 100
        if num_of_individuals == '4':
            if chr_size == '50Mb' or chr_size == '200Mb':
                number_of_msmc_replicates = 99
            elif chr_size == '1Gb':
                number_of_msmc_replicates = 95
        elif num_of_individuals == '8':
            if chr_size == '1Gb':
                print('Script has stopped as there are no results for 1Gb - 8 individuals')
                break
            
    elif simulation_tool == 'slim':
        number_of_msmc_replicates = 100
        if chr_size == '1Gb':
            if num_of_individuals == '2':
                number_of_msmc_replicates = 8
            if num_of_individuals == '4':
                number_of_msmc_replicates = 9
            if num_of_individuals == '8':
                number_of_msmc_replicates = 2
    
    true_generations, true_population_sizes, max_plotted_generation, min_viewable_population_size, max_viewable_population_size = demographic_model_parameters(chr_size)
    bins, msmc_replicate_population_sizes, msmc_average_population_sizes, dashed_line_start_index = analyze_msmc_results(number_of_msmc_replicates, max_plotted_generation, bin_size, chr_size, num_of_individuals, mutation_rate, simulation_tool)
    graph_results(max_plotted_generation, min_viewable_population_size, max_viewable_population_size, true_generations, true_population_sizes, bins, msmc_average_population_sizes, msmc_replicate_population_sizes, dashed_line_start_index, number_of_msmc_replicates, plot_msmc_average_line, plot_msmc_replicate_lines, chr_size, num_of_individuals, simulation_tool, save_graph_as_pdf, display_graph)
