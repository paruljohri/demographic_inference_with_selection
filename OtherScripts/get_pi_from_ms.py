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
parser.add_argument('-filename', dest = 'filename', action='store', nargs = 1, type = str, help = 'the name of the file with full path')

#read input parameters                                                          
args = parser.parse_args()
s_filename = args.filename[0]

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

if "genome7" in s_filename or "genome20" in s_filename:
    s_tot = 118425300
elif "genome03" in s_filename:
    s_tot = 146118000 #for genome03
elif "genome10" in s_filename:
    s_tot = 135570750
elif "genome05" in s_filename:
    s_tot = 142355000
f_ms = open(s_filename, 'r')
l_ms = read_ms_for_pylibseq(f_ms)
sd = libsequence.SimData(l_ms)
ps = libsequence.PolySIM(sd)
pi = ps.thetapi()
print (pi)
print (pi/float(s_tot))
print ("done")


