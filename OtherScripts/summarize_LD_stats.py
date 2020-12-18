#get rsq mean and sd
import sys
#import statistics
import numpy

chr_size = "200Mb"

l_rsq = []
l_S = []
repID = 1
while repID <= 100:
    #f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "/output_single_chr_" + chr_size + "_rep" + str(repID) + "_10000.stats", 'r')
    #f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "/output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_5kb_50000.stats", 'r')
    #f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "/output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_50kb_500000.stats", 'r')
    f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "/output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_100kb_1000000.stats", 'r')
    if chr_size=="1Mb":
        f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(repID) + "_100000.stats", 'r')
    elif chr_size=="10Mb":
        f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(repID) + "_1000000.stats", 'r')
    #elif chr_size=="50Mb":
    #    f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(repID) + "_5000000.stats", 'r')
    #elif chr_size=="200Mb":
    #    f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(repID) + "_20000000.stats", 'r')
    elif chr_size == "1Gb":
        f_tab = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(repID) + "_100000000.stats", 'r')
    d_col = {}
    for line in f_tab:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if line2[0] == "simID":
            col = 0
            for x in line2:
                d_col[x] = col
                col += 1
        else:
            #print (line2[d_col["rsq"]])
            if line2[d_col["rsq"]] != "NA" and line2[d_col["rsq"]]!="":
                l_rsq.append(float(line2[d_col["rsq"]]))
                l_S.append(float(line2[d_col["S"]]))
    f_tab.close()
    repID += 1

print ("mean rsq:")
print (numpy.mean(l_rsq))
print ("sd rsq")
print (numpy.std(l_rsq))
print ("mean S")
print (numpy.mean(l_S))
print ("sd of S")
print (numpy.std(l_S))
print ("done")