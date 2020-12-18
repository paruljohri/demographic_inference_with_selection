#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_slidingwindow_pylibseq_general.py -winSize 200 -stepSize 200 -regionLen 2341 -sim_num 21
#python statistics_slidingwindow_pylibseq_general_1rep.py -winSize 2341 -stepSize 2341 -regionLen 2341 -folder /scratch/pjohri1/VIRUS_DFE -filename foll_data.ms
#python statistics_slidingwindow_pylibseq_general_1rep.py -winSize 2314 -stepSize 2314 -regionLen 2314 -folder /scratch/pjohri1/VIRUS_DFE/ms_files_passages12_control_lines -filename Run1034_ms.txt
#python statistics_slidingwindow_pylibseq_general_1rep.py -winSize 23500 -stepSize 23500 -regionLen 23500 -folder /scratch/pjohri1/VIRUS_DFE/HCMV -filename qual0_475_folded_no5_ms.txt
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
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-filename', dest = 'filename', action='store', nargs = 1, type = str, help = 'full filename like sim1.ms')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
args = parser.parse_args()
chr_len =  args.regionLen[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
filename = args.filename[0]
in_folder = args.input_folder[0]
out_folder = args.output_folder[0]
#num_reps = int(args.numRep[0])


def read_fixed_mutations(f_fixed):
	d_subs = {}
	for line in f_fixed:
		line1 = line.strip('\n')
		line2 = line1.split()
		if line1[0]!="#" and line2[0]!="Mutations:":
			posn = float(line2[3])/float(chr_len)
			#if "grow10" in input_folder  or "red10" in input_folder:
			num_gen = line2.pop()
			#if int(num_gen) >= 100000:
			d_subs[posn] = d_subs.get(posn, 0) + 1
			#else:
			#	 d_subs[posn] = d_subs.get(posn, 0) + 1
	return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
	s_sum = 0
	for posn in d_subs.keys():
		if float(posn) <= end and float(posn) > start:
			s_sum = s_sum + 1
			
	return s_sum

def get_S(f_ms):
	for line in f_ms:
		line1 = line.strip('\n')
		if "segsites" in line1:
			S = line1.split()[1]
	f_ms.seek(0)
	return S
#result files:

result = open(out_folder + "/" + filename.split(".")[0] + "_" + str(args.winSize[0]) +	".stats", 'w+')
result.write("simID" + '\t' + "posn" + '\t' + "S" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\t' + "div" + '\n')

#go through the ms file and divergence file if it exists:
f_ms = open(in_folder + "/" + filename, 'r')
try:
	f_subs = open(f_name + ".txt", 'r')
	d_subs = read_fixed_mutations(f_subs)
except:
	print ("no divergence file")
	d_subs = {}
S = get_S(f_ms)
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
#print (l_Pos)
#print (l_Genos)


#assign object
sd = libsequence.SimData(l_data)
#sd.assign(l_Pos[10:100],l_Genos[10:100])

#define sliding windows:
w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
#chromosome length = 30kb, window size = 5 kb
num_win = len(w)

#calculate summary statistic in sliding window:
numsim=1
print ("calculating stats in windows")
win_name = 1
for i in range(len(w)):
	wi = w[i]
	#print (wi)
	pswi = libsequence.PolySIM(wi)
	result.write("sim" + str(numsim) + '\t' + str(win_name) + '\t' + str(len(wi.pos())) + '\t' + str(pswi.thetapi()) + '\t' + str(pswi.thetaw()) + '\t' + str(pswi.thetah()) + '\t' + str(pswi.hprime()) + '\t' + str(pswi.tajimasd()) + '\t' + str(pswi.numexternalmutations()) + '\t' + str(pswi.hapdiv()) + '\t')
	#read data to calculate LD based stats:
			
	if len(wi.pos()) >= 5: #These are pairwise stats. If only 1 site exists, it'll show an error.
		#print (i)
		LD_tmp = libsequence.ld(wi)
		LDstats = pandas.DataFrame(LD_tmp)
		meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
		meanD = sum(LDstats['D'])/len(LDstats['D'])
		meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
		result.write(str(meanrsq) + '\t' + str(meanD) + '\t' + str(meanDprime) + '\t')
	else:
		result.write("NA" + '\t' + "NA" + '\t' + "NA" + '\t') 
	#calculate and write divergence:
	s_start = (i)*(1.0/float(num_win))
	s_end = s_start + win_size
	result.write(str(avg_divergence_win(d_subs, s_start, s_end)) + '\n')
	win_name = win_name + 1

print ("Finished")






