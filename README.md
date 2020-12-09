# demographic_inference_with_selection
This is the list of files and folders that are included in the Supplemental information for the manuscript titled "The impact of purifying and background selection on the inference of population history: problems and prospects" by Johri et al. Contact pjohri1@asu.edu for questions.

1. General workflow and command lines for performing simulations and demographic inference:
>> Pipeline_MSMC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with MSMC.\
>> Pipeline_FSC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with fastsimcoal2.\
>> Pipeline_MSMC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection with MSMC.\
>> Pipeline_FSC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection with fastsimcoal2.

2. Scripts to perform simulations in SliM.
Folder: SlimScripts
>> eqm_single_chr_1Mb.slim //Slim script for 1Mb chromosome under neutral equilibrium\
>> eqm_single_chr_10Mb.slim //Slim script for 10Mb chromosome under neutral equilibrium\
>> eqm_single_chr_50Mb.slim //Slim script for 50Mb chromosome under neutral equilibrium\
>> eqm_single_chr_200Mb.slim //Slim script for 200Mb chromosome under neutral equilibrium\
>> eqm_single_chr_1Gb.slim //Slim script for 1Gb chromosome under neutral equilibrium\
>> eqm_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for equilibrium (20% of the genome is exonic)\
>> eqm_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for equilibrium (10% of the genome is exonic)\
>> eqm_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for equilibrium (5% of the genome is exonic)\
>> growth_dfe_150Mb_22chr_genome20.slim //Slim script for 30-fold exponential growth (20% of the genome is exonic)\
>> growth_dfe_150Mb_22chr_genome10.slim //Slim script for 30-fold exponential growth (10% of the genome is exonic)\
>> growth_dfe_150Mb_22chr_genome05.slim //Slim script for 30-fold exponential growth (05% of the genome is exonic)\
>> growth_2fold_dfe_150Mb_22chr_genome20.slim //Slim script for 2-fold exponential growth (20% of the genome is exonic)\
>> growth_2fold_dfe_150Mb_22chr_genome10.slim //Slim script for 2-fold exponential growth (10% of the genome is exonic)\
>> growth_2fold_dfe_150Mb_22chr_genome05.slim //Slim script for 2-fold exponential growth (5% of the genome is exonic)\
>> decline_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 6-fold instantaneous decline (20% of the genome is exonic)\
>> decline_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 6-fold instantaneous decline (10% of the genome is exonic)\
>> decline_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 6-fold instantaneous decline (5% of the genome is exonic)\
>> decline_2fold_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 2-fold instantaneous decline (20% of the genome is exonic)\
>> decline_2fold_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 2-fold instantaneous decline (10% of the genome is exonic)\
>> decline_2fold_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 2-fold instantaneous decline (5% of the genome is exonic)

3. Script to perform simulations in msprime
Folder: MsprimeScripts
>> simulate_msprime.py //Python script to simulate neutral equilibrium for varying chromosome segments and providing an output in .ms format

4. Other scripts used to process the input and output files generated for/by MSMC and fastsimcoal2.
Folder: OtherScripts
>> AIC_model_choice.py //Model selection in fastsimcoal2\
>> genome05_positions.txt, genome10_positions.txt, genome20_positions.txt //genomic position of the exons, introns and intergenic regions simulated.\
>> mask_coding_regions_human.py //Remove or mask the exonic regions from a .ms file.\
>> msmc_average_file_maker.py //Average MSMC outputs from replicates over fixed time intervals.\
>> plot_msmc_and_fsc_results.py //Plot the population size history estimates of either fastsimcoal2, MSMC, or both togethor.\ 
>> SFS_maker_human_genome_all_SNPs.py //Create an SFS for a human-like genome from the 22 .ms chromosome files.\
>> SFS_maker_human_genome_thinned_SNPs.py //Create an SFS for a human-like genome from the 22 .ms chromosome files only using SNPs that are at least 5kb apart from each other.\
>> SFS_maker_single_segment_all_SNPs.py //Create an SFS for a single chromosome segment from the .ms file.\
>> SFS_maker_single_segment_thinned_SNPs.py //Create an SFS for a single chromosome segment from the .ms file only using SNPs that are at least 5kb apart from each other.

5. Input files of MSMC
Folder: MSMC_InputFiles
>> WORK IN PROGRESS\

6. Output files of MSMC
Folder: MSMC_OutputFiles
>> "raw" files //The direct output from an MSMC run\
>> "transformed" files //An altered version of the "raw" file that has changed the values to be in terms of generations and population size based on the mutation rate\
>> "average" files //Average MSMC outputs from replicates over fixed time intervals created from the "transformed" files using the msmc_average_file_maker.py script\

7. Input files for fastsimcoal2
Folder: Fastsimcoal2_InputFiles
>> .tpl and .est files for equilibrium\
>> .tpl and .est files for instantaneous change\
>> .tpl and .est files for exponential change\
>> .tpl and .est files for instant bottleneck

8. Input SFS for fastsimcoal2
Folder: Fastsimcoal2_InputFiles_SFS
>> SFS input files needed to run fastsimcoal2

9. Output files for fastsimcoal2
Folder: Fastsimcoal2_OutputFiles
>> "bestlhood" files //The ".bestlhood" output files from fastsimcoal2 containing the estimated parameters and maximum estimated likelhood value for any particular run.\

10. Scripts to perform analytical calculations
Folder: AnalyticalCalculations
>> Analytical_size_change_B.nb //Mathematica notebook to calculate the theoretical expectations of B with population size change\
>> Expgrow-prog.txt //An example program that calculates analytical expectations of B and the SFS for exponential growth model.
