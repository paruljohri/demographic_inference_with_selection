#This is to get the SFS from .ms file for all 22 chromosomes:
#To run:
#python get_sfs_chr22_modified_for_kellen.py -demo ${demo} -dfe ${dfe} -rep ${rep_num} -num_indv 100 -number_of_chromosomes 1

import sys
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-demo', dest = 'demo', action='store', nargs = 1, type = str, help = 'demographic model; growth/ decline/ eqm ')
parser.add_argument('-dfe', dest = 'dfe', action='store', nargs = 1, type = str, help = 'DFE')
parser.add_argument('-rep', dest = 'rep', action='store', nargs = 1, type = str, help = 'repID')
parser.add_argument('-num_indv', dest = 'num_indv', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')
parser.add_argument('-number_of_chromosomes', dest = 'number_of_chromosomes', action='store', nargs = 1, type = int, help = 'number  of chromosomes')

#read input parameters
args = parser.parse_args()
demo = args.demo[0]
dfe = args.dfe[0] 
rep = args.rep[0]
num_indv = args.num_indv[0]
number_of_chromosomes = int(args.number_of_chromosomes[0])

def get_sfs(l_af):
    d_sfs = {}
    s_tot = 0
    for x in l_af:
        try:
            d_sfs[x] = d_sfs[x] + 1
        except:
            d_sfs[x] = 1
        if int(x) > 0 and int(x) < int(num_indv): #throws away all sites that are monomorphic
            s_tot += 1
    return(d_sfs, s_tot)

result = open("/scratch/pjohri1/MSMC_SIMS/genome10/sfs/" + demo + "_" + dfe + "_rep" + str(rep) + "_masked_" + str(num_indv) + ".sfs", 'w+')#modify path for output
l_af = []
chr_num = 1
while chr_num <= number_of_chromosomes:
    print ("chr number is: " + str(chr_num))
    d_af = {}
    f_ms = open("/scratch/pjohri1/MSMC_SIMS/genome10/" + demo + "/" + dfe + "/output_genome" + str(rep) + "_chr" + str(chr_num) + "_masked.ms", 'r') #modify path for input
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
    for posn in d_af.keys():
        l_af.append(d_af[posn])
    chr_num += 1

t_sfs = get_sfs(l_af)
d_sfs_all = t_sfs[0]
s_tot = t_sfs[1]
print (s_tot)
print (d_sfs_all)
result.write("#type")
i = 1
while (i < int(num_indv)):
    #print (i)
    result.write('\t' + str(i))
    i = i + 1
result.write('\n')
result.write("count")
i = 1
while (i < int(num_indv)):
    result.write('\t' + str(d_sfs_all.get(i, 0)))
    i = i + 1
result.write('\n')
result.write("freq")
i = 1
while (i < int(num_indv)):
    result.write('\t' + str(d_sfs_all[i]/float(s_tot)))
    i = i + 1
result.write('\n')

print ("done")
