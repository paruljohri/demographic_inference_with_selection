# demographic_inference_with_selection
#Updated on Dec 17th, 2020

This is the list of files and folders that are included in the Supplemental information for the manuscript titled "The impact of purifying and background selection on the inference of population history: problems and prospects" by Johri et al., which can be found here - https://www.biorxiv.org/content/10.1101/2020.04.28.066365v3 . Contact pjohri1@asu.edu for questions.

1. General workflow and command lines for performing simulations and demographic inference:
- Pipeline_MSMC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with MSMC.\
- Pipeline_FSC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with fastsimcoal2.\
- Pipeline_MSMC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection with MSMC.\
- Pipeline_FSC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection with fastsimcoal2.

2. Commandlines used to generate every Figure and table (including supplementary ones):\
Folder: CommandLines

3. Scripts to perform simulations in SliM.\
Folder: SlimScripts
- eqm_single_chr_1Mb.slim //Slim script for 1Mb chromosome under neutral equilibrium\
- eqm_single_chr_10Mb.slim //Slim script for 10Mb chromosome under neutral equilibrium\
- eqm_single_chr_50Mb.slim //Slim script for 50Mb chromosome under neutral equilibrium\
- eqm_single_chr_200Mb.slim //Slim script for 200Mb chromosome under neutral equilibrium\
- eqm_single_chr_1Gb.slim //Slim script for 1Gb chromosome under neutral equilibrium\
- eqm_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for equilibrium (20% of the genome is exonic)\
- eqm_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for equilibrium (10% of the genome is exonic)\
- eqm_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for equilibrium (5% of the genome is exonic)\
- growth_dfe_150Mb_22chr_genome20.slim //Slim script for 30-fold exponential growth (20% of the genome is exonic)\
- growth_dfe_150Mb_22chr_genome10.slim //Slim script for 30-fold exponential growth (10% of the genome is exonic)\
- growth_dfe_150Mb_22chr_genome05.slim //Slim script for 30-fold exponential growth (05% of the genome is exonic)\
- growth_2fold_dfe_150Mb_22chr_genome20.slim //Slim script for 2-fold exponential growth (20% of the genome is exonic)\
- growth_2fold_dfe_150Mb_22chr_genome10.slim //Slim script for 2-fold exponential growth (10% of the genome is exonic)\
- growth_2fold_dfe_150Mb_22chr_genome05.slim //Slim script for 2-fold exponential growth (5% of the genome is exonic)\
- decline_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 6-fold instantaneous decline (20% of the genome is exonic)\
- decline_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 6-fold instantaneous decline (10% of the genome is exonic)\
- decline_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 6-fold instantaneous decline (5% of the genome is exonic)\
- decline_2fold_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 2-fold instantaneous decline (20% of the genome is exonic)\
- decline_2fold_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 2-fold instantaneous decline (10% of the genome is exonic)\
- decline_2fold_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 2-fold instantaneous decline (5% of the genome is exonic)

4. Script to perform simulations in msprime\
Folder: MsprimeScripts
- simulate_msprime.py //Python script to simulate neutral equilibrium for varying chromosome segments and providing an output in .ms format

5. Other scripts used to process the input and output files generated for/by MSMC and fastsimcoal2.\
Folder: OtherScripts
- AIC_model_choice.py //Model selection in fastsimcoal2\
- genome05_positions.txt, genome10_positions.txt, genome20_positions.txt //genomic position of the exons, introns and intergenic regions simulated.\
- mask_coding_regions_human.py //Remove or mask the exonic regions from a .ms file.\
- msmc_average_file_maker.py //Average MSMC outputs from replicates over fixed time intervals.\
- plot_msmc_and_fsc_results.py //Plot the population size history estimates of either fastsimcoal2, MSMC, or both togethor.\ 
- get_sfs_for_fsc_chr22.py //Create an SFS for a human-like genome from the 22 .ms chromosome files.\
- get_sfs_for_fsc_chr22_thin.py //Create an SFS for a human-like genome from the 22 .ms chromosome files using thinned SNPs.\
- get_sfs_for_fsc.py //Create an SFS for a single chromosome segment from the .ms file.\
- get_sfs_for_fsc_thin.py //Create an SFS for a single chromosome segment from the .ms file using thinned SNPs.

6. Scripts used to perform Approximate Bayesian Computation (ABC)\
Folder: ABCScripts //Includes scripts for simulations, calculating statistics and performing ABC + plotting results

7. Scripts used for simulating variable mutation and recombination rates
Folder: RecMutnRateMapScripts

8. Output files of MSMC\
Folder: MSMC_OutputFiles //Output files of MSMC from all evolutionary scenarios.

9. Input files for fastsimcoal2\
Folder: Fastsimcoal2_InputFiles
- .tpl and .est files for equilibrium\
- .tpl and .est files for instantaneous change\
- .tpl and .est files for exponential change\
- .tpl and .est files for instant bottleneck

10. Input SFS for fastsimcoal2\
Folder: Fastsimcoal2_InputFiles_SFS \\SFS input files needed to run fastsimcoal2. View the Readme in the folder to  access the correct file.\

11. Output files from fastsimcoal2\
Folder: Fastsimcoal2_OutputFiles // *.bestlhood files //The ".bestlhood" output files from fastsimcoal2 containing the estimated parameters and maximum estimated likelhood value for any particular run.\

12. Scripts to perform analytical calculations\
Folder: AnalyticalCalculations
- Analytical_size_change_B.nb //Mathematica notebook to calculate the theoretical expectations of B with population size change\
- Expgrow-prog.txt //An example program that calculates analytical expectations of B and the SFS for exponential growth model.
