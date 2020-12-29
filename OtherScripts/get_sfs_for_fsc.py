#This is to get the SFS (including class 0) from a .ms file:
#How to run:
#python get_sfs_for_fsc.py -inputFolder /scratch/kriall/fsc/REDO_single_chr/chromosomes/${chr_size} -filename output_single_chr_${chr_size}_rep${repID}.ms -outputFolder /scratch/pjohri1/MSMC_SCRIPTS/mySFS/eqm_neutral_${chr_size}_slim -chrLength ${chr_size} -num_indv 100

import sys
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-inputFolder', dest = 'inputFolder', action='store', nargs = 1, type = str, help = 'path to input folder')
parser.add_argument('-filename', dest = 'filename', action='store', nargs = 1, type = str, help = 'ms filename')
parser.add_argument('-outputFolder', dest = 'outputFolder', action='store', nargs = 1, type = str, help = 'path to output folder')
parser.add_argument('-chrLength', dest = 'chrLength', action='store', nargs = 1, type = str, help = '1Mb/10Mb/50Mb/200Mb/1Gb')
parser.add_argument('-num_indv', dest = 'num_indv', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')

#read input parameters
args = parser.parse_args()
in_folder = args.inputFolder[0]
filename = args.filename[0]
out_folder = args.outputFolder[0]
chr_len = args.chrLength[0]
num_indv = args.num_indv[0]

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

f_ms = open(in_folder + "/" + filename, 'r')
result = open(out_folder + "/" + filename.split(".")[0] + "_" + str(num_indv) + "_all_DAFpop0.obs", 'w+')
result.write("1 observations" + '\n')
i = 0
result.write("d0_0")
while i <= int(num_indv):
    result.write('\t' + "d0_" + str(i))
    i = i + 1
result.write('\n')

#Read input file and read in derived allele frequency at each position:
d_af = {}
linecount = 0
for line in f_ms:
    line1 = line.strip('\n')
    if "//" not in line1 and "segsites" not in line1 and "positions" not in line1:
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

#Convert dictionary into a list
l_af = []
for posn in d_af.keys():
    l_af.append(d_af[posn])

#Calculate the SFS
t_sfs = get_sfs(l_af)
d_sfs_all = t_sfs[0]
s_seg = t_sfs[1]#Number of polymorphic sites, excludes sites with AF=1
s_not_anc = t_sfs[2] #Number of polymorphic sites + number of fixed sites in the sample

result.write(str(chr_size-s_not_anc))#Write the d0_0 class
i = 1
while (i <= int(num_indv)):
    result.write('\t' + str(d_sfs_all.get(i, 0)))
    i = i + 1
result.write('\n')

print ("done")

