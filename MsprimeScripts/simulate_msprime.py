#This is to simulate a chromosome of size 1Gb using msprime and output data in ms format:
#python simulate_msprime.py eqm_neutral_50Mb_msp 50000000 num_chrom repID seed
import sys
import os
import msprime

folder = sys.argv[1]
chrom_len = int(sys.argv[2])
num_chrom = int(sys.argv[3])
repID = sys.argv[4]
s_seed = int(sys.argv[5])*1000
#num_reps = 100

#define some constants:
n_mutn_rate = 1e-8
n_rec_rate = 1e-8
n_Ne = 5000

num_rep = 1
while num_rep <= num_chrom:
    try:
        f_check = open("/scratch/pjohri1/MSMC_SIMS/" + folder + "/output_genome" + str(repID) + "_chr" + str(num_rep) + ".ms", 'r')
        f_check.close()
    except:
        tree = msprime.simulate(length=chrom_len, recombination_rate=n_rec_rate, mutation_rate=n_mutn_rate, Ne=n_Ne, random_seed=s_seed+num_rep, sample_size=100)
        result = open("/scratch/pjohri1/MSMC_SIMS/" + folder + "/output_genome" + str(repID) + "_chr" + str(num_rep) + ".ms", 'w+')
        result.write("//" + '\n')
        d_posn, d_geno = {}, {}
        l_sites = []
        for variant in tree.variants():
            l_sites.append(variant.site.id)
            d_posn[variant.site.id] = round(float(variant.site.position)/chrom_len, 7)
            d_geno[variant.site.id] = variant.genotypes
            #print ('\t' + variant.site.position)

        result.write("segsites: " + str(len(l_sites)) + '\n')
        result.write("positions:")
        for site in l_sites:
            result.write(" " + str(d_posn[site]))
        result.write('\n')

        i = 0
        while i < 100:
            for site in l_sites:
                result.write(str(d_geno[site][i]))
            result.write('\n')
            i = i + 1
        result.close()
        ####loop ends for this replicate number
    num_rep = num_rep + 1

print ("done")

