#This is to perform 1000 replicates of each slim script:
#This program will create a temporary slim script and save it's output, and then delete the script.
#python write_slim_scripts_SingExon_osg.py -numRep 10 -simID 1
#go over starting ID to end ID. Check if a set has already been simulated. Then submit it if it hasn't been done.

from random import *
import os
import sys
import argparse
import math
parser = argparse.ArgumentParser(description='Information about number of replicates and cores used')
parser.add_argument('-numRep', dest = 'numReplicates', action='store', nargs = 1, default = 1, type = int, choices = range(1,10001), help = 'number of replicates of the same simulation')
parser.add_argument('-simID', dest = 'simID', action='store', nargs = 1, default = 1, type = int, help = 'the ID of the simulation to start with')

args = parser.parse_args()
num_rep = args.numReplicates[0]
simID = args.simID[0]

#Other constants:
num_gen_change = 5000
scaling_factor = 320

#make a directory for results of this simualtion:
os.system("mkdir /scratch/pjohri1/demo_disc_5_SingExon_osg/sim" + str(simID))

##copy and write the submit file with that simulation ID:
#f_submit = open("/home/pjohri/SingExon_slim_scratch.submit", 'r')
#result_submit = open("/local-scratch/pjohri/submitfiles/SingExon_slim_sim" + str(simID) + ".submit", 'w+')
#for line in f_submit:
#	line1 = line.replace("simID", "sim" + str(simID))
#	result_submit.write(line1)
#result_submit.close()
#f_submit.close()

#get parameters and ID of that simulation:
f_para = open("/home/pjohri1/demo_disc_parameters+5_SingExon_osg.txt", 'r')
d_col = {}
for line in f_para:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "sim":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    if line2[0] == "sim" + str(simID):
        print (line)
        f0 = line2[d_col["f0"]]
        f1 = line2[d_col["f1"]]
        f2 = line2[d_col["f2"]]
        f3 = line2[d_col["f3"]]
        Na = int(line2[d_col["Na"]])
        Nc = int(line2[d_col["Nc"]])
        model = line2[d_col["model"]]
        Nfold = float(line2[d_col["Nfold"]])
        g_factor = float(line2[d_col["g_factor"]])
f_para.close()


#re-write the main script and write out specific parameters in it with the simulation ID in it

f_exons = open("/home/pjohri1/sing_exon_500_2000_4000_prop.txt", 'r')

geneID = 0
for aline in f_exons:
    aline1 = aline.strip('\n')
    aline2 = aline1.split('\t')
    #print (line2)
    if aline2[0] != "gene":
        geneID = geneID + 1
        repID = 1
        counter = 0
        f_cmd = open("/scratch/pjohri1/demo_disc_5_SingExon_osg_slimfiles/sim" + str(simID) + "_gene" + str(geneID) + ".list", 'w+')
        while repID <= 10:
            #check if this sim has already been done:
            try:
                f_check = open("/scratch/pjohri1/demo_disc_5_SingExon_osg/sim" + str(simID) + "/sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".ms", 'r')
                f_check.close()
            except:
                counter = counter + 1
                f_cmd.write("slim script_sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".slim" + '\n')
                gene_name = aline2[0]
                exon_size = aline2[3]
                #print (aline2[4] + '\t' + aline2[5])
                rec_rate = str((float(aline2[4]) + float(aline2[5]))/2.0)
                #print (rec_rate)
			    
				    #print ("started reading slim file")
                f_script = open("/home/pjohri1/SlimPrograms/demo_disc_5_SingExon_osg.slim", 'r')
                result = open("/scratch/pjohri1/demo_disc_5_SingExon_osg_slimfiles/script_sim" + str(simID) + "_gene" + str(geneID) + "_rep" + str(repID) + ".slim", 'w+')
                for line in f_script:
                    line1 = line.replace("f0", f0)
                    line2 = line1.replace("f1", f1)
                    line3 = line2.replace("f2", f2)
                    line4 = line3.replace("f3", f3)
                    line5 = line4.replace("end_exon_size", str(4000 + int(exon_size)-1))
                    line6 = line5.replace("rec_rate", str(rec_rate))
                    line7 = line6.replace("N_anc", str(Na))
                    line8 = line7.replace("N_cur", str(Nc))
                    line9 = line8.replace("gen_burnin", str(10*Na))
                    line10 = line9.replace("g_factor", str(g_factor))
                    line11 = line10.replace("gen_eqm", str((10*Na) + 20000))
                    line12 = line11.replace("gen_stop", str((10*Na) + 20000 + 5000))
                    line13 = line12.replace("simID", "sim" + str(simID))
                    line14 = line13.replace("geneID", "gene" + str(geneID))
                    line15 = line14.replace("repID", "rep" + str(repID))
                    result.write(line15)
				    #print("end of reading")
                result.close()
                f_script.close()
            repID = repID + 1
        f_cmd.close()
        if counter==0:
            os.system("rm /scratch/pjohri1/demo_disc_5_SingExon_osg_slimfiles/sim" + str(simID) + "_gene" + str(geneID) + ".list")


print ("done")


