#Key to access raw output files from MSMC

#For varying chromosome sizes simulated under neutral equilibrium, files are in folder - MSMC_Chr_Segments_OutputFiles.rar

#Files, for both SLiM and msprime results, are named as follows: raw_output_single_chr_${chr_size}_${num_of_individuals}indv_rep${replicate_num}.final.txt
#where
#${chr_size} = (1Mb, 10Mb, 50Mb, 200Mb, 1Gb) #size of the chromosomes simulated
#${num_of_individuals} = (2, 4, 8) #number of indiviudals MSMC was run using
#${replicate_num} = (1, 2, 3, 4, 5, ..., 99, 100) #For each chromosome size 100 replicates were simulated.



#For human-like genomes simulated, the files are in folder - MSMC_Human_Genome_OutputFiles.rar

#Files for results from genomes made with fixed recombination and mutation rates, found in the "fixed_rec_and_mut" folder, are named as follows: RAW_rep${replicate_num}_${demographic_model}_sim${DFE}_${exon_density}_${masking_state}.final.txt
#where
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, unmasked)
#${exon_density} = (genome05, genome10, genome20)
#${DFE} = (0, 1, 2, 3, 4, 5, 6) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3,4,5,6,7,8,9,10)

#Files for results from genomes made with variable recombination and mutation rates, found in the "variable_rec_and_mut" folder, are named as follows: RAW_rep${replicate_num}_${demographic_model}_sim${DFE}_${exon_density}_${masking_state}.final.txt
#where
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, masked_random_regions, masked_centromere, masked_centromere+random_regions)
#${exon_density} = (genome20)
#${DFE} = (0, 4) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3)
