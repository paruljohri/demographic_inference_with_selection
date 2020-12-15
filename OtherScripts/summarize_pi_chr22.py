#This is to make a table that summarizes mean pi from all genomes and reps:

import sys



l_models = ["eqm", "growth", "decline"]
l_dfe = ["neutral", "dfe1", "dfe2", "dfe3", "dfe4", "dfe5", "dfe6"]
l_sims = ["sim0", "sim1", "sim2", "sim3", "sim4", "sim5", "sim6"]
l_reps = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7", "rep8", "rep9", "rep10"]

result = open("/scratch/pjohri1/MSMC_SIMS/genome10/pi/genome10.pi", 'w+')
result.write("model" + '\t' + "DFE" + '\t' + "replicate" + '\t' + "model2" + '\t' + "simnum" + '\t' + "repnum" + '\t' + "pi" + '\n')
for model in l_models:
    i = 0
    for sim in l_sims:
        for rep in l_reps:
            try:
                f_pi = open("/scratch/pjohri1/MSMC_SIMS/genome10/pi/" + model + "_" + l_sims[i] + "_" + str(rep) + "_masked.pi", 'r')
                for line in f_pi:
                    line1 = line.strip('\n')
                    if "#" not in line1:
                        if line1 != "":
                            result.write(model + '\t' + l_dfe[i] + '\t' + rep)
                            result.write('\t' + line1 + '\n')
                f_pi.close()
            except:
                print ("File not present")
                result.write(model + '\t' + l_dfe[i] + '\t' + rep)
                result.write('\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\n')
        i += 1

result.close()
print ("done")

