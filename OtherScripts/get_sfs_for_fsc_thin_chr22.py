#This is to get the SFS (including class 0) from an .ms file in the format for fsc:
#It includes the option of thinning SNPs
#How to run:
#python get_sfs_for_fsc_thin_chr22.py -inputFolder /scratch/kriall/simulated_chromosomes -outputFolder /scratch/pjohri1/MSMC_SCRIPTS/mySFS/ -demography decline -genome genome20 -dfe 0 -repNum 1 -masking masked -thinning 5kb -num_indv 100

import sys
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-inputFolder', dest = 'inputFolder', action='store', nargs = 1, type = str, help = 'path to input folder')
parser.add_argument('-outputFolder', dest = 'outputFolder', action='store', nargs = 1, type = str, help = 'path to output folder')
parser.add_argument('-demography', dest = 'demography', action='store', nargs = 1, type = str, help = 'eqm/decline/growth')
parser.add_argument('-genome', dest = 'genome', action='store', nargs = 1, type = str, help = 'genome20/genome10/genome05')
parser.add_argument('-dfe', dest = 'dfe', action='store', nargs = 1, type = str, help = '0/1/2/3/4/5/6')
parser.add_argument('-repNum', dest = 'repNum', action='store', nargs = 1, type = str, help = '1-10')
parser.add_argument('-masking', dest = 'masking', action='store', nargs = 1, type = str, help = 'masked/unmasked')
parser.add_argument('-thinning', dest = 'thinning', action='store', nargs = 1, type = str, help = '5kb/100kb')
parser.add_argument('-num_indv', dest = 'num_indv', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')

#read input parameters
args = parser.parse_args()
in_folder = args.inputFolder[0]
out_folder = args.outputFolder[0]
s_demo = args.demography[0]
genome = args.genome[0]
s_dfe = args.dfe[0]
repID = args.repNum[0]
masking = args.masking[0]
thinning = args.thinning[0]
num_indv = args.num_indv[0]

if "genome20" in genome:
    if masking == "unmasked":
        chr_size = 150003700
    elif masking == "masked":
        chr_size = 118424753
elif "genome10" in genome:
    if masking == "unmasked":
        chr_size = 150029950
    elif masking == "masked":
        chr_size = 135572119
elif "genome05" in genome:
    if masking == "unmasked":
        chr_size = "150018600"
    elif masking=="masked":
        chr_size = 142354366

if thinning == "5kb":
    thin_size = 5000
elif thinning == "100kb":
    thin_size = 100000
else:
    print ("Error: check thinning parameters")

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

def thin_snps(d_af, d_posns, s_interval):
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

#f_ms = open(sys.argv[1], 'r')
result = open(out_folder + "/MASTER_thinned_" + thinning + "_" + s_demo + "_" + masking + "_" + genome + "_sim" + s_dfe + "_rep" + repID + "_DAFpop0.obs", 'w+')
result.write("1 observations" + '\n')
i = 0
result.write("d0_0")
while i <= int(num_indv):
    result.write('\t' + "d0_" + str(i))
    i = i + 1
result.write('\n')

#Make a list of all .ms files:
#os.system("ls " + in_folder + "/*.ms > " + out_folder + "/tmp.list")
#Going through .ms files for all 22 chromosomes:
d_sfs_chr22 = {}
af_bin = 0
while af_bin <= num_indv:
    d_sfs_chr22[af_bin] = 0
    af_bin += 1
chr_num = 1
#f_list = open(out_folder + "/tmp.list", 'r')
while chr_num <= 22:
    #Bline = Aline.strip('\n')
    print ("chromosome number: " + str(chr_num))
    if masking == "unmasked":
        f_ms = open(in_folder + "/" + s_demo + "_dfe_150Mb_22chr_" + genome + "/sim" + s_dfe + "/output_genome" + repID + "_chr" + str(chr_num) + ".ms", 'r')
    elif masking == "masked":
        f_ms = open(in_folder + "/" + s_demo + "_dfe_150Mb_22chr_" + genome + "/sim" + s_dfe + "/output_genome" + repID + "_chr" + str(chr_num) + "_masked.ms", 'r')
    
    #reading in genotypes from ms
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
    
    #Thin the ms file:
    d_af_thinned = thin_snps(d_af, d_posns, thin_size)
    #Make SFS for thinned SNPs
    l_af = []
    for y in d_af_thinned.keys():
        l_af.append(d_af_thinned[y])
    t_sfs = get_sfs(l_af)
    d_sfs_all = t_sfs[0]
    s_seg = t_sfs[1]
    s_not_anc = t_sfs[2]
    chr_size_thinned = round(chr_size*(len(d_af_thinned)/float(len(d_af))))
    s_bin0 = chr_size_thinned-s_not_anc
    #Add to the cumulative SFS:
    d_sfs_chr22[0] = d_sfs_chr22[0] + s_bin0
    for x in d_sfs_all.keys():
        d_sfs_chr22[x] = int(d_sfs_chr22[x]) + int(d_sfs_all[x])
    chr_num += 1

#Write the full result:
result.write(str(d_sfs_chr22[0]))#write the d0_0 class
i = 1
while (i <= int(num_indv)):
    result.write('\t' + str(d_sfs_chr22[i]))
    i = i + 1
result.write('\n')

print ("done")

