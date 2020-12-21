#Key to access output files of fastsimcoal2.
#Currently only the .bestlhoods are provided.
#Fastsimcoal2_Chr_Segments_OutputFiles.rar //for output of simulated genomes of varying chromosome sizes under neutral equilibrium
#Fastsimcoal2_Human_Genome_OutputFiles_genome05.rar //for output of simulated human-like genomes with exons comprising 5% of the genome
#Fastsimcoal2_Human_Genome_OutputFiles_genome10.rar //for output of simulated human-like genomes with exons comprising 10% of the genome
#Fastsimcoal2_Human_Genome_OutputFiles_genome20.rar //for output of simulated human-like genomes with exons comprising 20% of the genome
#Fastsimcoal2_Human_Genome_OutputFiles_variable_rec_and_mut.rar //for output of simulated human-like genomes with variable recombination and mutation rate maps with various types of additional masking

#For varying chromosome sizes simulated under neutral equilibrium, files are in folder - Fastsimcoal2_Chr_Segments_OutputFiles.rar
#Files, for both SLiM and msprime results, are named as follows: ${model_type}_${SNP_density}_${chr_size}_rep${replicate_num}.bestlhoods
#where
#${model_type} = (eqm, size_change_exponential, size_change_inst, inst_bot)
#${SNP_density} = (all, thinned_5kb, thinned_50kb, thinned_100kb) #is the density at which SNPs were sampled from the simulated genomes
#${chr_size} = (1Mb, 10Mb, 50Mb, 200Mb, 1Gb) #size of the chromosomes simulated
#${replicate_num} = (1, 2, 3, 4, 5, ..., 99, 100) #For each chromosome size 100 replicates were simulated.


#For human-like genomes simulated, the files are in folder - Fastsimcoal2_Human_Genome_OutputFiles_${exon_density}.rar
#Files for results from genomes made with fixed recombination and mutation rates are named as follows: ${model_type}_${SNP_density}_${demographic_model}_${masking_state}_${exon_density}_sim${DFE}_rep${replicate_num}_trial${trial_num}.bestlhoods
#where
#${model_type} = (eqm, size_change_exponential, size_change_inst, inst_bot)
#${SNP_density} = (all, thinned_5kb)
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, unmasked)
#${exon_density} = (genome05, genome10, genome20)
#${DFE} = (0, 1, 2, 3, 4, 5, 6) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3,4,5,6,7,8,9,10)
#${trial_num} = (1,2,3,4,5,6,7,8,9,10)


#For human-like genomes simulated under variable recombination and mutation rate maps with various types of additional masking, the files are in folder - Fastsimcoal2_Human_Genome_OutputFiles_variable_rec_and_mut.rar
#Files for results from genomes made with fixed recombination and mutation rates are named as follows: ${model_type}_${SNP_density}_recombination+mutation_${demographic_model}_${masking_state}_${exon_density}_sim${DFE}_rep${replicate_num}_trial${trial_num}.bestlhoods
#where
#${model_type} = (eqm, size_change_exponential, size_change_inst, inst_bot)
#${SNP_density} = (all, thinned_5kb)
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, masked_random_regions, masked_centromere, masked_centromere+random_regions)
#${exon_density} = (genome05, genome10, genome20)
#${DFE} = (0, 4) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3)
#${trial_num} = (1,2,3,4,5,6,7,8,9,10)
