#This is to get the final set of stats for SingExon sims:
#In this script we first take means and variance of the 10 replicates and then use means of those values.
#python get_final_statistics_SingExon.py 3 numbp50 demo_disc_5_SingExon_osg_stats
#python get_final_statistics_SingExon.py 1 numbp50 demo_disc_5_SingExon

import sys
import statistics 

simID = sys.argv[1]
cutoffLink = sys.argv[2]
folder = sys.argv[3]
f_stats = open("/scratch/pjohri1/" + folder + "/sim" + simID + "_bigwindow_" + cutoffLink + ".stats", 'r')

def average(l_numbers):
    l_floats = []
    for x in l_numbers:
        if x != "NA" and x!="nan" and x != "inf":
            l_floats.append(float(x))
    avg = statistics.mean(l_floats)
    return avg

def sd(l_numbers):
    l_floats = []
    for x in l_numbers:
        if x != "NA" and x!="nan" and x != "inf":
            l_floats.append(float(x))
    stdev = statistics.stdev(l_floats)
    return stdev

def get_replicateNum(s_simnum):
    sim_num = int(s_simnum.replace("sim",""))
    if sim_num%10 == 0:
        repNum = "rep10"
    else:
        repNum = "rep" + str(sim_num%10)
    return repNum

#To calculate mean and variance of all statistics per window, for 94 genes:
l_windows = ["functional", "linked", "neutral"]
l_rep = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7", "rep8", "rep9", "rep10"]
d_func = {}
d_link = {}
d_neu = {}
l_stats = []
d_col = {}

for line in f_stats:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "simID":
        col = 0
        for x in line2:
            d_col[col] = x
            d_func[x], d_link[x], d_neu[x] = {}, {}, {}
            if x != "simID"  and x != "WinType" and x != "WinSize":
                l_stats.append(x)
            for rep in l_rep:
                d_func[x][rep] = []
                d_link[x][rep] = []
                d_neu[x][rep] = []
            col = col + 1
    else:
        if "gene" in line2[0]:
            rep_num = line2[0].split("_")[1]
        else:
            rep_num = get_replicateNum(line2[0])
        if line2[1] == "functional":
            col = 0
            for x in line2:
                if d_col[col] == "thetapi" or d_col[col] == "thetaw" or d_col[col]=="thetah" or d_col[col]=="numSing" or d_col[col]=="div":
                    if x != "NA":
                        d_func[d_col[col]][rep_num].append(float(x)/float(line2[2]))
                else:
                    d_func[d_col[col]][rep_num].append(x)
                col = col + 1
        elif line2[1] == "linked":
            col = 0
            for x in line2:
                if d_col[col] == "thetapi" or d_col[col] == "thetaw" or d_col[col]=="thetah" or d_col[col]=="numSing" or d_col[col]=="div":
                    if x != "NA":
                        d_link[d_col[col]][rep_num].append(float(x)/float(line2[2]))
                else:
                    d_link[d_col[col]][rep_num].append(x)
                col = col + 1
        elif line2[1] == "neutral":
            col = 0
            for x in line2:
                if d_col[col] == "thetapi" or d_col[col] == "thetaw" or d_col[col]=="thetah" or d_col[col]=="numSing" or d_col[col]=="div":
                    if x != "NA":
                        d_neu[d_col[col]][rep_num].append(float(x)/float(line2[2]))
                else:
                    d_neu[d_col[col]][rep_num].append(x)
                col = col + 1
f_stats.close()
#print (d_func)
result = open("/scratch/pjohri1/" + folder + "/sim" + str(simID) + "_" + str(cutoffLink) + ".bigwinsummary", 'w+')
result.write("functional" + '\t' + "linked" + '\t' + "neutral" + '\n')

#get mean across genes and then take means across reps:
for stat in l_stats:
    result.write(stat + "_m" + '\t')
    l_means_func, l_means_link, l_means_neu = [], [], []
    for rep in l_rep:
        #print (stat + '\t' + rep)
        #print (d_func[stat][rep])
        l_means_func.append(average(d_func[stat][rep]))#put mean across genes for each replicate in a list
        l_means_link.append(average(d_link[stat][rep]))
        l_means_neu.append(average(d_neu[stat][rep]))
    result.write(str("{:.15f}".format(average(l_means_func))) + '\t' + str("{:.15f}".format(average(l_means_link))) + '\t' + str("{:.15f}".format(average(l_means_neu))) + '\n')#write the mean of mean

#get variance across genes and then take means across reps:
for stat in l_stats:
    result.write(stat + "_sd" + '\t')
    l_sd_func, l_sd_link, l_sd_neu = [], [], []
    for rep in l_rep:
        l_sd_func.append(sd(d_func[stat][rep]))#put variance across genes for each replicate in a list
        l_sd_link.append(sd(d_link[stat][rep]))
        l_sd_neu.append(sd(d_neu[stat][rep]))
    #print (stat + '\t' + rep)
    #print (l_sd_neu)
    result.write(str("{:.15f}".format(average(l_sd_func))) + '\t' + str("{:.15f}".format(average(l_sd_link))) + '\t' + str("{:.15f}".format(average(l_sd_neu))) + '\n')#write the mean of mean

result.close()
print ("done")







