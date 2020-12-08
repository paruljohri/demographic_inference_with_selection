# demographic_inference_with_selection
This is the list of files and folders that are included in the Supplemental information for the manuscript titled "The impact of purifying and background selection on the inference of population history: problems and prospects" by Johri et al. Contact pjohri1@asu.edu for questions.\

1. General workflow and command lines for performing simulations and demographic inference:\
>> Pipeline_MSMC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with MSMC.\
>> Pipeline_FSC_chr_segments.txt //workflow for chromosome segments of varying sizes under neutral equilibrium with fastsimcoal2.\
>> Pipeline_MSMC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection (with MSMC).\
>> Pipeline_FSC_human_chr.txt //workflow for human-like genomes under different demographic histories and selection (with fastsimcoal2).\

2. Scripts to perform simulations in SliM.
Folder: SlimScripts
>> eqm_single_chr_10Mb.slim //Slim script for 10Mb chromosome under neutral equilibrium
>> eqm_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for equilibrium (20% of the genome is exonic)
>> eqm_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for equilibrium (10% of the genome is exonic)
>> eqm_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for equilibrium (5% of the genome is exonic)
>> growth_dfe_150Mb_22chr_genome20.slim //Slim script for 30-fold exponential growth (20% of the genome is exonic)
>> growth_dfe_150Mb_22chr_genome10.slim //Slim script for 30-fold exponential growth (10% of the genome is exonic)
>> growth_dfe_150Mb_22chr_genome05.slim //Slim script for 30-fold exponential growth (05% of the genome is exonic)
>> growth_2fold_dfe_150Mb_22chr_genome20.slim //Slim script for 2-fold exponential growth (20% of the genome is exonic)
>> growth_2fold_dfe_150Mb_22chr_genome10.slim //Slim script for 2-fold exponential growth (10% of the genome is exonic)
>> growth_2fold_dfe_150Mb_22chr_genome05.slim //Slim script for 2-fold exponential growth (5% of the genome is exonic)
>> decline_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 6-fold instantaneous decline (20% of the genome is exonic)
>> decline_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 6-fold instantaneous decline (10% of the genome is exonic)
>> decline_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 6-fold instantaneous decline (5% of the genome is exonic)
>> decline_2fold_dfe_150Mb_22chr_scaled_genome20.slim //Slim script for 2-fold instantaneous decline (20% of the genome is exonic)
>> decline_2fold_dfe_150Mb_22chr_scaled_genome10.slim //Slim script for 2-fold instantaneous decline (10% of the genome is exonic)
>> decline_2fold_dfe_150Mb_22chr_scaled_genome05.slim //Slim script for 2-fold instantaneous decline (5% of the genome is exonic)

3. Script to perform simulations in msprime
Folder: MsprimeScripts
>> simulate_msprime.py //Python script to simulate neutral equilibrium for varying chromosome segments and providing an output in .ms format

4. Other scripts used to process the input and output files generated for/by MSMC and fastsimcoal2.
Folder: OtherScripts
>> AIC_model_choice.py //Model selection in fastsimcoal2
>> genome05_positions.txt, genome10_positions.txt, genome20_positions.txt //genomic position of the exons, introns and intergenic regions simulated.
>> mask_coding_regions_human.py //Remove or mask the exonic regions from a .ms file.
>> msmc_average_file_maker.py //average the results from 10 replicates of MSMC
>> plot_msmc_and_fsc_results.py //
>> SFS_maker_all.py //
>> SFS_maker_single_chr_all.py //
>> SFS_maker_single_chr_thinned.py //
>> SFS_maker_thinned.py //
>> plotting_script.py //Python script to plot either only MSMC results, or only fastsimcoal2 results, or to plot them together.
>> msmc_average_file_maker.py //Python script to average MSMC outputs from replicates over fixed time intervals.

5. Output files of MSMC

6. Input files for fastsimcoal2
Folder: Fastsimcoal2_InputFiles
>> .tpl and .est files for equilibrium
>> .tpl and .est files for instantaneous change
>> .tpl and .est files for exponential change
>> .tpl and .est files for instant bottleneck

7. Input SFS for fastsimcoal2

8. Scripts to perform analytical calculations
Folder: AnalyticalCalculations
>> Analytical_size_change_B.nb //Mathematica notebook to calculate the theoretical expectations of B with population size change
>> Expgrow-prog.txt //An example program that calculates analytical expectations of B and the SFS for exponential growth model.
