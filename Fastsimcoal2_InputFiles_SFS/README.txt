#This is the key to access the SFS files used as input files for fastsimcoal2 

#For varying chromosome sizes simulated under neutral equilibrium, files are in folder - Fastsimcoal2_Chr_Segments_SFS_InputFiles.rar
#Files, for both SLiM and msprime results, are named as follows, : MASTER_${SNP_density}_${chr_size}_rep${replicate_num}_DAFpop0.obs
#where
#${SNP_density} = (all, thinned_5kb, thinned_50kb, thinned_100kb) #is the density at which SNPs were sampled from the simulated genomes
#${chr_size} = (1Mb, 10Mb, 50Mb, 200Mb, 1Gb) #size of the chromosomes simulated
#${replicate_num} = (1, 2, 3, 4, 5, ..., 99, 100) #For each chromosome size 100 replicates were simulated.



#For human-like genomes simulated - Fastsimcoal2_Human_Genome_SFS_InputFiles.rar

#Files for results from genomes made with fixed recombination and mutation rates are named as follows: MASTER_${SNP_density}_${demographic_model}_${masking_state}_${exon_density}_sim${DFE}_rep${replicate_num}_DAFpop0.obs
#where
#${SNP_density} = (all, thinned_5kb)
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, unmasked)
#${exon_density} = (genome05, genome10, genome20)
#${DFE} = (0, 1, 2, 3, 4, 5, 6) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3,4,5,6,7,8,9,10)

#Files for results from genomes made with variable recombination and mutation rates are named as follows: MASTER_${SNP_density}_recombination+mutation_${demographic_model}_${masking_state}_${exon_density}_sim${DFE}_rep${replicate_num}_DAFpop0.obs
#where
#${SNP_density} = (all, thinned_5kb)
#${demographic_model} = (eqm, growth, decline)
#${masking_state} = (masked, masked_random_regions, masked_centromere, masked_centromere+random_regions)
#${exon_density} = (genome05, genome10, genome20)
#${DFE} = (0, 4) where ${DFE}=0 corresponds to neutrality
#${replicate_num} = (1,2,3)
