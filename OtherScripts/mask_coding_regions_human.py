#This script is to mask coding regions from an ms file (mostly just remove them) and then create anohter ms file.
#python mask_coding_regions_human.py growth_dfe_150Mb_22chr_slim simID chrID repID

import sys
folder = sys.argv[1]
simID = sys.argv[2]
rep_num = sys.argv[3]
chr_num = sys.argv[4]

#Read the human genome coordinates and figure out which sites to mask:
if "genome4" in folder:
    f_genome = open("/scratch/pjohri1/MSMC_SCRIPTS/genome4_positions.txt", 'r')
    tot_len = 150086999
elif "genome20" in folder:
    f_genome = open("/scratch/pjohri1/MSMC_SCRIPTS/genome20_positions.txt", 'r')
    tot_len = 150003699
elif "genome10" in folder:
    f_genome = open("/scratch/pjohri1/MSMC_SCRIPTS/genome10_positions.txt", 'r')
    tot_len = 150029949
elif "genome03" in folder:
    f_genome = open("/scratch/pjohri1/MSMC_SCRIPTS/human_genome_positions.txt", 'r')
    tot_len = 150012799
elif "genome05" in folder:
    f_genome = open("/scratch/pjohri1/MSMC_SCRIPTS/genome05_positions.txt", 'r')
    tot_len = 150018599
d_sites_mask = {}
for line in f_genome:
    if "initializeGenomicElement(" in line:
        line1 = line.strip('\n').replace("initializeGenomicElement", "")
        line2 = line1.replace(" ", "")
        line3 = line2.replace("(", "")
        line4 = line3.replace(")", "")
        line5 = line4.replace(";", "")
        line6 = line5.split(',')
        s_region = line6[0]
        s_start = int(line6[1])
        s_end = int(line6[2])
        if s_region == "g3":#g3 is coding region
            posn = s_start
            while posn <= s_end:
                s_site = round(float(posn)/float(tot_len),7)
                d_sites_mask[s_site] = "mask"
                posn = posn + 1
f_genome.close()

#tot_chr = 22
#tot_reps = 100
#chr_num = 14
l_unmasked_posns = []
#while chr_num <= tot_chr:
#    rep_num = 1
#    while rep_num <= tot_reps:
l_unmasked_posns = []
print ("chr_num: " + str(chr_num) + "  rep_num:  " + str(rep_num))
f_ms = open("/scratch/kriall/" + folder + "/sim" + str(simID) + "/output_genome" + str(rep_num) + "_chr" + str(chr_num) + ".ms", 'r')
result = open("/scratch/kriall/" + folder + "/sim" + str(simID) + "/output_genome" + str(rep_num) + "_chr" + str(chr_num) + "_masked.ms", 'w+')
#identify unmasked positions
for line in f_ms:
    line1 = line.strip('\n')
    if "positions" in line1:
        line2 = line1.split()
        col = 0
        for posn in line2:
            if "position" not in posn:
                s_site = float(posn)
                if d_sites_mask.get(s_site, "NA") != "mask":#the site is not masked
                    l_unmasked_posns.append(col)
                col = col + 1
print(max(l_unmasked_posns))
#re-write ms file
f_ms.seek(0)
result.write("//" + '\n')
result.write("segsites: " + str(len(l_unmasked_posns)) + '\n')
#result.write("positions:")
for line in f_ms:
    line1 = line.strip('\n')
    line2 = line1.split()
    #if "//" in line1:
    #    result.write("//" + '\n')
    #elif "segsites" in line1:
    #    result.write("segsites: " + str(len(l_unmasked_posns)) + '\n')
    if "positions" in line1:
        result.write("positions:")
        for col_num in l_unmasked_posns:
            result.write(" " + str(line2[col_num+1]))
        result.write('\n')
    elif "//" not in line1 and "segsites" not in line1:
        #print (line2)
        for col_num in l_unmasked_posns:
            result.write(line1[col_num])
        result.write('\n')
result.close()
f_ms.close()
#        rep_num = rep_num + 1
#    chr_num = chr_num + 1
        
print("done")
