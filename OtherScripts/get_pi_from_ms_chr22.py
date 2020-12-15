#This is to get just the mean pi from an ms file:
#For the genomes this is meant to get pi from masked files.

from __future__ import print_function
import sys
import libsequence
import pandas                                                                   
import math                                                                     
import argparse

#parsing user given constants                                                   
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-demo', dest = 'demo', action='store', nargs = 1, type = str, help = 'demographic model; growth/ decline/ eqm ')
parser.add_argument('-dfe', dest = 'dfe', action='store', nargs = 1, type = str, help = 'DFE')
parser.add_argument('-rep', dest = 'rep', action='store', nargs = 1, type = str, help = 'repID')
parser.add_argument('-genome', dest = 'genome', action='store', nargs = 1, type = str, help = 'genome number like genome10')
parser.add_argument('-burnin', dest = 'burnin', action='store', nargs = 1, type = str, help = 'using values at burnin (T) or after (F)')
#read input parameters                                                          
args = parser.parse_args()
s_demo = args.demo[0]
s_dfe = args.dfe[0]
repID = args.rep[0]
s_genome = args.genome[0]
s_burnin = args.burnin[0]

#function to read ms file
def read_ms_for_pylibseq(f_ms):
    l_Pos = [] #list of positions of SNPs as given by ms file      
    d_gt = {} #hash table with SNP position -> genotype of each column
    l_data = []
    for line in f_ms:
        #print (line)
        line1 = line.strip('\n')
        if "positions" in line1:
            line2 = line1.split()
            i = 0
            for x in line2:
                if "position" not in x:
                    l_Pos.append(float(x))
                    d_gt[str(i)] = ""
                    i = i + 1
        elif "//" not in line and "segsites" not in line:
            i = 0
            while i < len(line1):
                d_gt[str(i)] = d_gt[str(i)] + line1[i]
                i = i + 1
    i = 0
    while i < len(l_Pos):
        t_tmp = (l_Pos[i], d_gt[str(i)])
        l_data.append(t_tmp)
        i = i + 1
    return (l_data)

if s_genome=="genome7" or s_genome=="genome20":
    s_tot = 118425300
elif s_genome=="genome03":
    s_tot = 146118000 #for genome03
elif s_genome=="genome10":
    s_tot = 135570750
elif s_genome=="genome05":
    s_tot = 142355000

s_pi_sum = 0.0
s_tot_sum = 0.0
chr_num = 1
while chr_num <= 22:
    print ("chr number is: " + str(chr_num))
    if s_burnin == "T":
        f_ms = open("/scratch/pjohri1/MSMC_SIMS/genome10/" + s_demo + "_burnin/" + s_dfe + "/output_genome" + str(repID) + "_chr" + str(chr_num) + "_masked.ms", 'r')
    elif s_burnin == "F":
        f_ms = open("/scratch/pjohri1/MSMC_SIMS/genome10/" + s_demo + "/" + s_dfe + "/output_genome" + str(repID) + "_chr" + str(chr_num) + "_masked.ms", 'r')
    l_ms = read_ms_for_pylibseq(f_ms)
    sd = libsequence.SimData(l_ms)
    ps = libsequence.PolySIM(sd)
    pi = ps.thetapi()
    s_pi_sum = s_pi_sum + pi
    s_tot_sum = s_tot_sum + s_tot
    chr_num += 1
if s_burnin == "T":
    result = open("/scratch/pjohri1/MSMC_SIMS/genome10/pi_burnin/" + s_demo + "_" + s_dfe + "_rep" + str(repID) + "_masked.pi", 'w+')
elif s_burnin == "F":
    result = open("/scratch/pjohri1/MSMC_SIMS/genome10/pi/" + s_demo + "_" + s_dfe + "_rep" + str(repID) + "_masked.pi", 'w+')
result.write("#demo" + '\t' + "dfe" + '\t' + "repID" + '\t' + "pi" + '\n')
result.write(s_demo + '\t' + s_dfe + '\t' + str(repID) + '\t' + str(float(s_pi_sum)/float(s_tot_sum)) + '\n')
#print (pi)
#print (pi/float(s_tot))
print ("done")


