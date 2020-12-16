#Basic stats, sliding window, SLIm ms output
#Uses coding sizes for every gene separately
#python statistics_slidingwindow_pylibseq_SingExon_osg.py -winSize 200 -stepSize 200 -simID 1 -noncodingLen 4000
from __future__ import print_function
#from libsequence.polytable import SimData
#from libsequence.summstats import PolySIM
#from libsequence.windows import Windows
#from libsequence.summstats import ld
import libsequence
import sys
import pandas
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, choices = range(1,10000), help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, choices = range(1,10000), help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-folder', dest = 'FolderName', action='store', nargs = 1, type = str, help = 'the name of the folder or simulation to run')
parser.add_argument('-simID', dest = 'simulationID', action='store', nargs = 1, type = str, help = 'the name of the subfolder or simulation to run')
parser.add_argument('-noncodingLen', dest = 'noncodingLen', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'noncoding length in bp')

#read input parameters
args = parser.parse_args()
noncoding_size =  args.noncodingLen[0]
folder = args.FolderName[0]
simID = args.simulationID[0]
subfolder = "sim" + str(simID)

#defined constants:
num_genes = 94 #total number of replicates
num_reps = 10 #number of replicates of each gene

#get length of coding from file:
f_sing = open("/home/pjohri1/sing_exon_500_2000_4000_prop.txt", 'r')
d_codingsize = {}
geneID = 0
for line in f_sing:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "gene":
        geneID = geneID + 1
        #while i <= num_rep:
        #    repID = repID + 1
        d_codingsize["gene" + str(geneID)] = line2[3]
        #    i = i + 1
f_sing.close()

#get burn-in period from file:
if folder == "demo_disc_5_SingExon_osg":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_osg.txt", 'r')
elif folder == "demo_disc_5_SingExon_dfealpha":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_dfealpha.txt", 'r')
elif folder == "demo_disc_5_SingExon_posterior":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_posterior.txt", 'r')
elif folder == "demo_disc_5_SingExon_positive":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_positive.txt", 'r')
elif folder == "demo_disc_5_SingExon_negative":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_negative.txt", 'r')
elif folder == "demo_disc_5_SingExon_mutation":
    f_par = open("/home/pjohri1/demo_disc_parameters+5_SingExon_mutation.txt", 'r')
d_col = {}
for line in f_par:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "sim":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    elif line2[0] == "sim" + str(simID):
        gen_burnin = 10*int(line2[d_col["Na"]])
f_par.close()

#read fixed mutations: output of Slim:
def read_fixed_mutations(f_fixed):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        line2 = line1.split()
        if line1[0]!="#" and line2[0]!="Mutations:":
            posn = float(line2[3])/float(chr_len)
            num_gen = line2.pop()
            if int(num_gen) >= gen_burnin:
                d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) <= end and float(posn) > start:
			s_sum = s_sum + 1
	return s_sum

#result files:
result = open("/scratch/pjohri1/" + folder + "_stats/sim" + str(simID) + "_" + str(args.winSize[0]) + ".stats", 'w+')
result.write("simnum" + '\t' + "posn" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')

#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
numsim = 1
s_absent = 0
geneID = 1
while geneID <= num_genes:
    repID = 1
    while repID <= num_reps:
        #print ("Reading file:" + str(numsim))
        coding_size = int(d_codingsize["gene" + str(geneID)])
        chr_len = noncoding_size + coding_size
        win_size = args.winSize[0]/float(chr_len)
        step_size = args.stepSize[0]/float(chr_len)
        try:
        #if numsim > 0:
            #if folder == "demo_disc_5_SingExon_negative":
            #    f_ms = open("/scratch/pjohri1/" + folder + "/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep"+ str(repID) + "_masked.ms", 'r')
            #else:
            f_ms = open("/scratch/pjohri1/" + folder + "/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep"+ str(repID) + ".ms", 'r')
            f_subs = open("/scratch/pjohri1/" + folder + "/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".fixed", 'r')
            d_subs = read_fixed_mutations(f_subs)
            l_Pos = [] #list of positions of SNPs
            l_Genos = [] #list of alleles
            d_tmp = {}
            for line in f_ms:
                line1 = line.strip('\n')
                if "positions" in line1:
                    line2 = line1.split()
                    i = 0
                    for x in line2:
                        if "position" not in x:
                            l_Pos.append(float(x))
                            d_tmp[str(i)] = ""
                            i = i + 1
                elif "//" not in line and "segsites" not in line:
                    #print (d_tmp)
                    i = 0
                    while i < len(line1):
                        d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                        i = i + 1
		    #print (d_tmp)
            l_data = []
            i = 0
            while i < len(l_Pos):
                l_Genos.append(d_tmp[str(i)])
                t_tmp = (l_Pos[i], d_tmp[str(i)])
                l_data.append(t_tmp)
                i = i + 1


		    #assign object
            sd = libsequence.SimData(l_data)
		    #sd.assign(l_Pos[10:100],l_Genos[10:100])

		    #define sliding windows:
            w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
		    #chromosome length = 30kb, window size = 5 kb
            num_win = len(w)

		    #calculate summary statistic in sliding window:
            print ("calculating stats in windows")
            win_name = 1
            for i in range(len(w)):
                wi = w[i]
                #print (wi)
                pswi = libsequence.PolySIM(wi)
                result.write("gene" + str(geneID) + "_rep" + str(repID) + '\t' + str(win_name) + '\t' + str("{:.5f}".format(pswi.thetapi())) + '\t' + str("{:.5f}".format(pswi.thetaw())) + '\t' + str("{:.5f}".format(pswi.thetah())) + '\t' + str("{:.5f}".format(pswi.hprime())) + '\t' + str("{:.5f}".format(pswi.tajimasd())) + '\t' + str(pswi.numexternalmutations()) + '\t' + str("{:.5f}".format(pswi.hapdiv())) + '\t')
			    #read data to calculate LD based stats:
                
                if len(wi.pos()) >= 5: #These are pairwise stats. If only 1 site exists, it'll show an error.
                    #print (i)
                    LD_tmp = libsequence.ld(wi)
                    LDstats = pandas.DataFrame(LD_tmp)
                    meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
                    meanD = sum(LDstats['D'])/len(LDstats['D'])
                    meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
                    result.write(str("{:.5f}".format(meanrsq)) + '\t' + str("{:.5f}".format(meanD)) + '\t' + str("{:.5f}".format(meanDprime)) + '\t')
                else:
                    result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t') 
        	    #divergence:
                s_start = (i)*(1.0/float(num_win))
                s_end = s_start + win_size
                result.write(str("{:.5f}".format(avg_divergence_win(d_subs, s_start, s_end))) + '\n')
                win_name = win_name + 1
        except:
        #else:
            s_absent = s_absent + 1
            print ("This file does not exist or cannot be read or is empty")
        repID = repID + 1
        numsim = numsim + 1
    geneID = geneID + 1
 
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






