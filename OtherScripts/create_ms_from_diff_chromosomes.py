#This is to create an ms file that has SNPs from different chromosomes
#It is to meausre the baseline unlinked LD
#How to run:
#python create_ms_from_diff_chromosomes.py -regionLen 10000000 -chr_size 10Mb -thinning all

import sys
import os
import argparse
import random

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-chr_size', dest = 'chr_size', action='store', nargs = 1, type = str, help = 'foldername with .ms files')
parser.add_argument('-thinning', dest = 'thinning', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
#parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
chr_len =  args.regionLen[0]
chr_size =  args.chr_size[0]
thinning = args.thinning[0]
#outfolder = args.output_folder[0]
#prefix = args.output_prefix[0]
num_indv = 100

if chr_size == "1Gb":
    num_reps = 10
else:
    num_reps = 100

#read ms file:
def read_subset_ms(f_ms, start, end):
    l_Pos = [] #list of positions of SNPs
    d_gt = {}
    l_int_posn = []
    for line in f_ms:
        line1 = line.strip('\n')
        if "positions" in line1:
            line2 = line1.split()
            i = 0
            for x in line2:
                if "position" not in x:
                    if (float(x) >= float(start)) and (float(x) <= float(end)):
                        l_Pos.append(float(x))
                        d_gt[str(i)] = ""
                        l_int_posn.append(i)#1,2...n
                    i = i + 1
        elif "//" not in line and "segsites" not in line:
            for j in l_int_posn:
                d_gt[str(j)] = d_gt[str(j)] + line1[j]
    return (d_gt)

#create a dictionary for 100 random craws from a chromosome each
d_random_gts = {}
draw_num = 1
while draw_num <= 100:
    d_random_gts[draw_num] = []
    draw_num += 1

#go through all simulation replicates and read data into pylibseq format
repID = 1
while repID <= num_reps:
    print ("Reading file:" + str(repID))
    if thinning == "all":
        f_ms = open("/scratch/kriall/fsc/REDO_single_chr/chromosomes/" + chr_size + "/" + "output_single_chr_" + chr_size + "_rep" + str(repID) + ".ms", 'r') #for all SNPs
    elif thinning == "5kb":
        f_ms = open("/scratch/kriall/fsc/REDO_single_chr/chromosomes/" + chr_size + "/" + "output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_5kb.ms", 'r') #for 1 in 5kb
    elif thinning == "50kb":
        f_ms = open("/scratch/kriall/fsc/REDO_single_chr/chromosomes/" + chr_size + "/" + "output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_50kb.ms", 'r') #for 1 in 50kb
    elif thinning == "100kb":    
        f_ms = open("/scratch/kriall/fsc/REDO_single_chr/chromosomes/" + chr_size + "/" + "output_single_chr_" + chr_size + "_rep" + str(repID) + "_thinned_100kb.ms", 'r') #for 1 in 100kb

    d_ms = read_subset_ms(f_ms, 0.0, 1.0)
    f_ms.close()
    if len(d_ms["1"]) == num_indv:#checking that the .ms file has the number of individuals assumed
        draw_num = 1
        while draw_num <= 100:
            s_rand = random.randint(1, len(d_ms))
            d_random_gts[draw_num].append(d_ms[str(s_rand)])
            draw_num += 1
    repID += 1
print(len(d_random_gts))

#write the 100 draws into .ms files with arbitrary positions
draw_num = 1
while draw_num <= 100:
    num_SNPs = len(d_random_gts[draw_num])
    s_posn_interval = float(chr_len)/float(num_SNPs)
    #print ("draw num: " + str(draw_num))
    result = open("/home/pjohri1/MSMC/thinned_SNPs_LD/" + chr_size + "_unlinked/draw" + str(draw_num) + ".ms", 'w+')
    result.write("//" + '\n' + "segsites: " + str(num_SNPs) + '\n' + "positions:")
    j=1
    while j<=num_SNPs:
        result.write(" " + str((s_posn_interval*float(j))/float(chr_len)))
        j += 1
    #write the actual genotypes:
    row_num=0
    while row_num<num_indv:
        for l_gts in d_random_gts[draw_num]:
            #print (len(l_gts))
            #print ("row_num" + str(row_num))
            result.write(l_gts[row_num])
        result.write('\n')
        row_num += 1
    result.close()
    draw_num += 1        

print("done")
 
