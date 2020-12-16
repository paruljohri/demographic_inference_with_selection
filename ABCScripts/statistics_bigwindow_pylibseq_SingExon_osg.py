#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_bigwindow_pylibseq_SingExon_osg.py -LinkCutoff numbp50 -folder eqm_disc_5 -simID 1
from __future__ import print_function
import libsequence
#from libsequence.polytable import SimData
#from libsequence.summstats import PolySIM
#from libsequence.windows import Windows
#from libsequence.summstats import ld
import sys
import pandas
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-LinkCutoff', dest = 'LinkCutoff', action='store', nargs = 1, type = str, help = 'which cutoff- 50, 75 or 90 are you using?')
parser.add_argument('-folder', dest = 'FolderName', action='store', nargs = 1, type = str, help = 'the name of the folder or simulation to run')
parser.add_argument('-simID', dest = 'simulationID', action='store', nargs = 1, type = str, help = 'the name of the subfolder or simulation to run')
parser.add_argument('-noncodingLen', dest = 'noncodingLen', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'noncoding length in bp')
parser.add_argument('-numRep', dest = 'numRep', action='store', nargs = 1, type = int, choices = range(1,10001), help = 'number of replicates of each gene')
args = parser.parse_args()
cutoff = args.LinkCutoff[0]
folder = args.FolderName[0]
simID = args.simulationID[0]
subfolder = "sim" + str(simID)
noncoding_size = args.noncodingLen[0]
num_rep = args.numRep[0]
#defined constants:
win_size = "200"
num_genes = 94
#get length of coding from file:
f_sing = open("/home/pjohri1/sing_exon_500_2000_4000_prop.txt", 'r')
d_codingsize = {}
geneID = 0
for line in f_sing:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "gene":
        geneID = geneID + 1
        d_codingsize["gene" + str(geneID)] = line2[3]
f_sing.close()


#read ms file:
def read_subset_ms(f_ms, start, end):
	l_Pos = [] #list of positions of SNPs
	#l_Genos = [] #list of alleles
	d_tmp = {}
	l_int_posn = []
	for line in f_ms:
		line1 = line.strip('\n')
		if "positions" in line1:
			line2 = line1.split()
			i = 0
			for x in line2:
				if "position" not in x:
					if (float(x) > float(start)) and (float(x) < float(end)):
						l_Pos.append(float(x))
						d_tmp[str(i)] = ""
						l_int_posn.append(i)
					i = i + 1
		elif "//" not in line and "segsites" not in line:
			for j in l_int_posn:
				d_tmp[str(j)] = d_tmp[str(j)] + line1[j]
	return (l_Pos, d_tmp, l_int_posn)

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
print ("gen_burnin: " + str(gen_burnin))

#read fixed mutations: output of Slim:
def read_subset_fixed_mutations(f_fixed, start, end):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        line2 = line1.split()
        if line1[0]!="#" and line2[0]!="Mutations:":
            posn = float(line2[3])/float(chr_len)
            num_gen = line2[8]
            #print (num_gen)
            if int(num_gen) >= gen_burnin:#include this mutation only if it fixed after burnin period- 10Na.
                print (line2)
                if posn > float(start) and posn < float(end):
                    d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) < float(end) and float(posn) > float(start):
			s_sum = s_sum + 1
	return s_sum

def print_5_dig(s_number):
    #print (s_number)
    if s_number == "NA" or s_number == "nan" or s_number == "inf":
        to_print = s_number
    else:
        #print(s_num)
        to_print = str("{:.5f}".format(s_number))
    return to_print

#get numbpLinked:
f_rec = open("/scratch/pjohri1/" + folder + "_stats/sim" + str(simID) + "_" + str(win_size) + ".recoverypi", 'r')
for line in f_rec:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] == "slope":
		i = 0
		for x in line2:
			if x == cutoff:
				col = i
			i = i + 1
	else:
		numbp = line2[col]
if float(numbp) <= noncoding_size:
    startLinked_size = noncoding_size - float(numbp)
else:
	startLinked_size = 0
f_rec.close()
#print ("bins are:" + '\t' + str(start_neu) + '\t' + str(startLinked) + '\t' + str(start_func) + '\t' + str(end_func))
#try:
#	result = open("/home/pjohri/" + folder + "/" + subfolder + "_bigwindow" + ".stats", 'a')
#except:
result = open("/scratch/pjohri1/" + folder + "_stats/" + subfolder + "_bigwindow_" + cutoff +  ".stats", 'w+')
result.write("simID" + '\t' + "WinType" + '\t' + "WinSize" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')


#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
geneID = 1
s_absent = 0
while geneID <= num_genes:
    #start_neu = 0.0
    coding_size = d_codingsize["gene" + str(geneID)]
    chr_len = noncoding_size + int(coding_size)
    start_neu = float(2000)/float(chr_len)
    startLinked = float(startLinked_size)/float(chr_len)
    start_func = float(noncoding_size)/float(chr_len)
    end_func = 1.0
    print ("Reading files in: gene" + str(geneID))
    print ("bins are:" + '\t' + str(start_neu) + '\t' + str(startLinked) + '\t' + str(start_func) + '\t' + str(end_func))
    repID = 1
    while repID <= num_rep:
        #if repID > 0:
        try:
            print ("Reading file: sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID))
            #f_ms = open("/home/pjohri/" + folder + "/" + subfolder + "/output" + str(numsim) + ".ms", 'r')
            f_subs = open("/scratch/pjohri1/" + folder + "/" + subfolder + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".fixed", 'r')
            d_subs_func = read_subset_fixed_mutations(f_subs, start_func, end_func)
            f_subs.seek(0)
            d_subs_link = read_subset_fixed_mutations(f_subs, startLinked, start_func)
            f_subs.seek(0)
            d_subs_neu = read_subset_fixed_mutations(f_subs, start_neu, startLinked)
            f_subs.close()
            #reading ms file:
            #if folder == "demo_disc_5_SingExon_negative":
            #    f_ms = open("/scratch/pjohri1/" + folder + "/" + subfolder + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + "_masked.ms", 'r')
            #else:
            f_ms = open("/scratch/pjohri1/" + folder + "/" + subfolder + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".ms", 'r')
            t_ms_func = read_subset_ms(f_ms, start_func, end_func)
            f_ms.seek(0)
            t_ms_link = read_subset_ms(f_ms, startLinked, start_func)
            f_ms.seek(0)
            t_ms_neu = read_subset_ms(f_ms, start_neu, startLinked)
            f_ms.close()

            l_Pos_func, l_Pos_link, l_Pos_neu = t_ms_func[0], t_ms_link[0], t_ms_neu[0]
		    #l_Genos_func, l_Genos_link, l_Genos_neu = t_ms_func[1], t_ms_link[1], t_ms_neu[1]
            d_tmp_func, d_tmp_link, d_tmp_neu = t_ms_func[1], t_ms_link[1], t_ms_neu[1]
            l_intPos_func, l_intPos_link, l_intPos_neu = t_ms_func[2], t_ms_link[2], t_ms_neu[2]
            print (len(d_subs_func))
            print (len(d_subs_link))
            print (len(d_subs_neu))
            #func:
            l_data = []
            l_Genos_func = []
            j = 0
            for i in l_intPos_func:
                l_Genos_func.append(d_tmp_func[str(i)])
                t_tmp = (l_Pos_func[j], d_tmp_func[str(i)])
                l_data.append(t_tmp)
                j = j + 1
            sd_func = libsequence.SimData(l_data)
		    #link:
            l_data = []
            l_Genos_link = []
            j = 0
            for i in l_intPos_link:
                l_Genos_link.append(d_tmp_link[str(i)])
                t_tmp = (l_Pos_link[j], d_tmp_link[str(i)])
                l_data.append(t_tmp)
                j = j + 1
            sd_link = libsequence.SimData(l_data)
            #neu:
            l_data = []
            l_Genos_neu = []
            j = 0
            for i in l_intPos_neu:
                l_Genos_neu.append(d_tmp_neu[str(i)])
                t_tmp = (l_Pos_neu[j], d_tmp_neu[str(i)])
                l_data.append(t_tmp)
                j = j + 1
            sd_neu = libsequence.SimData(l_data)
            #calculating stats:
            ps_func, ps_link, ps_neu = libsequence.PolySIM(sd_func), libsequence.PolySIM(sd_link), libsequence.PolySIM(sd_neu)
            if len(sd_func.pos()) >= 5:
                LD_func = libsequence.ld(sd_func)
                LDstats_func = pandas.DataFrame(LD_func)
                meanrsq_func = sum(LDstats_func['rsq'])/len(LDstats_func['rsq'])
                meanD_func = sum(LDstats_func['D'])/len(LDstats_func['D'])
                meanDprime_func = sum(LDstats_func['Dprime'])/len(LDstats_func['Dprime'])
            else:
                LD_func, meanrsq_func, meanD_func, meanDprime_func = "NA", "NA", "NA", "NA"
            if len(sd_link.pos()) >= 5:
                LD_link = libsequence.ld(sd_link)
                LDstats_link = pandas.DataFrame(LD_link)
                meanrsq_link = sum(LDstats_link['rsq'])/len(LDstats_link['rsq'])
                meanD_link = sum(LDstats_link['D'])/len(LDstats_link['D'])
                meanDprime_link = sum(LDstats_link['Dprime'])/len(LDstats_link['Dprime'])
            else:
                LD_link, meanrsq_link, meanD_link, meanDprime_link = "NA", "NA", "NA", "NA"
            if len(sd_neu.pos()) >= 5:
                LD_neu = libsequence.ld(sd_neu)
                LDstats_neu = pandas.DataFrame(LD_neu)
                meanrsq_neu = sum(LDstats_neu['rsq'])/len(LDstats_neu['rsq'])
                meanD_neu = sum(LDstats_neu['D'])/len(LDstats_neu['D'])
                meanDprime_neu = sum(LDstats_neu['Dprime'])/len(LDstats_neu['Dprime'])
            else:
                LD_neu, meanrsq_neu, meanD_neu, meanDprime_neu = "NA", "NA", "NA", "NA"
		    #LDstats_func, LDstats_link, LDstats_neu = pandas.DataFrame(LD_func), pandas.DataFrame(LD_link), pandas.DataFrame(LD_neu)
		    #meanrsq_func, meanrsq_link, meanrsq_neu = sum(LDstats_func['rsq'])/len(LDstats_func['rsq']), sum(LDstats_link['rsq'])/len(LDstats_link['rsq']), sum(LDstats_neu['rsq'])/len(LDstats_neu['rsq'])
		    #meanD_func, meanD_link, meanD_neu = sum(LDstats_func['D'])/len(LDstats_func['D']), sum(LDstats_link['D'])/len(LDstats_link['D']), sum(LDstats_neu['D'])/len(LDstats_neu['D'])
		    #meanDprime_func, meanDprime_link, meanDprime_neu = sum(LDstats_func['Dprime'])/len(LDstats_func['Dprime']), sum(LDstats_link['Dprime'])/len(LDstats_link['Dprime']), sum(LDstats_neu['Dprime'])/len(LDstats_neu['Dprime'])
            div_func, div_link, div_neu = avg_divergence_win(d_subs_func, start_func, end_func), avg_divergence_win(d_subs_link, startLinked, start_func), avg_divergence_win(d_subs_neu, start_neu, startLinked)
		    #write results:
            result.write("gene" + str(geneID) + "_rep" + str(repID) + '\t' + str("functional") + '\t' + str((end_func-start_func)*float(chr_len)) + '\t' + print_5_dig(ps_func.thetapi()) + '\t' + print_5_dig(ps_func.thetaw()) + '\t' + print_5_dig(ps_func.thetah()) + '\t' + print_5_dig(ps_func.hprime()) + '\t' + print_5_dig(ps_func.tajimasd()) + '\t' + str(ps_func.numexternalmutations()) + '\t' + print_5_dig(ps_func.hapdiv()) + '\t' + print_5_dig(meanrsq_func) + '\t' + print_5_dig(meanD_func) + '\t' + print_5_dig(meanDprime_func) + '\t' + print_5_dig(div_func) + '\n')
            result.write("gene" + str(geneID) + "_rep" + str(repID) + '\t' + str("linked") + '\t' + str((start_func-startLinked)*float(chr_len)) + '\t' + print_5_dig(ps_link.thetapi()) + '\t' + print_5_dig(ps_link.thetaw()) + '\t' + print_5_dig(ps_link.thetah()) + '\t' + print_5_dig(ps_link.hprime()) + '\t' + print_5_dig(ps_link.tajimasd()) + '\t' + print_5_dig(ps_link.numexternalmutations()) + '\t' + print_5_dig(ps_link.hapdiv()) + '\t' + print_5_dig(meanrsq_link) + '\t' + print_5_dig(meanD_link) + '\t' + print_5_dig(meanDprime_link) + '\t' + print_5_dig(div_link) + '\n')
            result.write("gene" + str(geneID) + "_rep" + str(repID) + '\t' + str("neutral") + '\t' + str((startLinked-start_neu)*float(chr_len)) + '\t' + print_5_dig(ps_neu.thetapi()) + '\t' + print_5_dig(ps_neu.thetaw()) + '\t' + print_5_dig(ps_neu.thetah()) + '\t' + print_5_dig(ps_neu.hprime()) + '\t' + print_5_dig(ps_neu.tajimasd()) + '\t' + str(ps_neu.numexternalmutations()) + '\t' + print_5_dig(ps_neu.hapdiv()) + '\t' + print_5_dig(meanrsq_neu) + '\t' + print_5_dig(meanD_neu) + '\t' + print_5_dig(meanDprime_neu) + '\t' + print_5_dig(div_neu) + '\n')
        #else:
        except:
            s_absent = s_absent + 1
            print ("This file does not exist or cannot be read or is empty")
        repID = repID + 1
    geneID = geneID + 1

result.close()
print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






