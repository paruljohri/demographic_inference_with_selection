#This is to make a table that summarizes all SFSs of a genome with a fixed number fo individuals:

import sys

num_indv = sys.argv[1]

l_models = ["eqm", "growth", "decline"]
l_dfe = ["neutral", "dfe1", "dfe2", "dfe3", "dfe4", "dfe5", "dfe6"]
l_sims = ["sim0", "sim1", "sim2", "sim3", "sim4", "sim5", "sim6"]
l_reps = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7", "rep8", "rep9", "rep10"]

result = open("/scratch/pjohri1/MSMC_SIMS/genome10/sfs/genome10_numindv" + str(num_indv) + ".sfs", 'w+')
result.write("model" + '\t' + "DFE" + '\t' + "replicate" + '\t' + "type")
s_af = 1
while s_af < int(num_indv):
    result.write('\t' + str(s_af))
    s_af += 1
result.write('\n')
for model in l_models:
    i = 0
    for sim in l_sims:
        for rep in l_reps:
            f_sfs = open("/scratch/pjohri1/MSMC_SIMS/genome10/sfs/" + model + "_" + l_sims[i] + "_" + str(rep) + "_masked_" + str(num_indv) + ".sfs", 'r')
            for line in f_sfs:
                line1 = line.strip('\n')
                if "#" not in line1:
                    if line1 != "":
                        result.write(model + '\t' + l_dfe[i] + '\t' + rep)
                        result.write('\t' + line1 + '\n')
            f_sfs.close()
        i += 1

result.close()
print ("done")

