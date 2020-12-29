#This is to get the SFS (including the 0 class) from a .ms file in the format for fsc:
#It includes the option of thinning SNPs
#How to run:
#python get_sfs_for_fsc_thin.py -inputFolder /scratch/kriall/fsc/REDO_single_chr/chromosomes/${chr_size} -filename output_single_chr_${chr_size}_rep${repID}.ms -outputFolder /scratch/pjohri1/MSMC_SCRIPTS/mySFS/eqm_neutral_${chr_size}_slim -chrLength ${chr_size} -num_indv 100 -thinningInterval 5Kb

import sys
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-inputFolder', dest = 'inputFolder', action='store', nargs = 1, type = str, help = 'path to input folder')
parser.add_argument('-filename', dest = 'filename', action='store', nargs = 1, type = str, help = 'ms filename')
parser.add_argument('-outputFolder', dest = 'outputFolder', action='store', nargs = 1, type = str, help = 'path to output folder')
parser.add_argument('-chrLength', dest = 'chrLength', action='store', nargs = 1, type = str, help = '1Mb/10Mb/50Mb/200Mb/1Gb')
parser.add_argument('-num_indv', dest = 'num_indv', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')
parser.add_argument('-thinningInterval', dest = 'thinningInterval', action='store', nargs = 1, type = str, help = '5Kb/50Kb/100Kb')
#read input parameters
args = parser.parse_args()
in_folder = args.inputFolder[0]
filename = args.filename[0]
out_folder = args.outputFolder[0]
chr_len = args.chrLength[0]
num_indv = args.num_indv[0]
thinning = args.thinningInterval[0]

#Set some parameters:
if chr_len == "1Mb":
    chr_size = 1000000
elif chr_len == "10Mb":
    chr_size = 10000000
elif chr_len == "50Mb":
    chr_size = 50000000
elif chr_len == "200Mb":
    chr_size = 200000000
elif chr_len == "1Gb":
    chr_size = 1000000000
if thinning == "5Kb":
    thin_size = 5000
elif thinning == "50Kb":
    thin_size = 50000
elif thinning == "100Kb":
    thin_size = 100000

def get_sfs(l_af):
    d_sfs = {}
    s_seg = 0 #total number of truly segregating sites
    s_not_anc = 0 #required to know the d0_0 class
    for x in l_af:
        try:
            d_sfs[x] = d_sfs[x] + 1
        except:
            d_sfs[x] = 1
        if int(x) > 0 and int(x) < int(num_indv):
            s_seg += 1
        if int(x) > 0:
            s_not_anc += 1
    return(d_sfs, s_seg, s_not_anc)

def thin_sfs(d_af, d_posns, s_interval):
    d_af_thin = {}
    s_posn0 = 1
    d_af_thin[s_posn0] = d_af[s_posn0]
    s_posn = s_posn0 + 1
    while s_posn <= len(d_af):
        if d_posns[s_posn]-d_posns[s_posn0]-1 >= int(s_interval):
            d_af_thin[s_posn] = d_af[s_posn]
            s_posn0 = s_posn
        s_posn += 1
    return(d_af_thin)

#Read input file and create output file
f_ms = open(in_folder + "/" + filename, 'r')
result = open(out_folder + "/" + filename.split(".")[0] + "_" + str(num_indv) + "_thinned_" + thinning + "_DAFpop0.obs", 'w+')
result.write("1 observations" + '\n')
i = 0
result.write("d0_0")
while i <= int(num_indv):
    result.write('\t' + "d0_" + str(i))
    i = i + 1
result.write('\n')

#Read in input file and store derived allele frequencies at all positions
d_af = {} #column numbers starting at 1 -> genotype
d_posns = {} #column numbers starting at 1 -> chromosomal positions
linecount = 0
for line in f_ms:
    line1 = line.strip('\n')
    if "positions" in line1:
        line2 = line1.split()
        col = 1
        for x in line2:
            if "positions" not in x:
                d_posns[col] = round(float(x)*chr_size)
                col += 1
    elif "//" not in line1 and "segsites" not in line1 and "positions" not in line1:
        linecount += 1
        if linecount <= int(num_indv):
            col = 1
            for x in line1:
                try:
                    d_af[col] = int(d_af[col]) + int(x)
                except:
                    d_af[col] = int(x)
                col += 1
f_ms.close()

#Thin the dictionary with allele frequencies:
d_af_thinned = thin_sfs(d_af, d_posns, thin_size)

#Convert the thinned AF dictionary into a list
l_af = []
for posn in d_af_thinned.keys():
    l_af.append(d_af_thinned[posn])

#Calulate the thinned SFS
t_sfs = get_sfs(l_af)
d_sfs_all = t_sfs[0]
s_seg = t_sfs[1] #Number of polymorphic sites in the thinned sample
s_not_anc = t_sfs[2] #Number of polymorphic sites in the thinned sample + number of fixed sites for the derived allele in the thinned sample
chr_size_thinned = round(chr_size*(len(d_af_thinned)/float(len(d_af))))#Re-scaling of the total chromosome size proportionally to the thinning
result.write(str(chr_size_thinned-s_not_anc))#write the d0_0 class
i = 1
while (i <= int(num_indv)):
    result.write('\t' + str(d_sfs_all.get(i, 0)))
    i = i + 1
result.write('\n')

print ("done")

