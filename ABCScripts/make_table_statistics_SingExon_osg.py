#This is to make a final table of all statistics:
#python make_table_statistics_SingExon_osg.py demo_disc_5_SingExon_osg_stats 3543
import sys
folder = sys.argv[1]
tot_numsims = int(sys.argv[2])
outputfile = folder.replace("/", "_")

#Going through stats files:
result_pi = open("/scratch/pjohri1/" + folder + "/" + outputfile + "_sumstats_pi.txt", 'w+')
result_pi.write("simID" + '\t' + "f0" + '\t' + "f1" + '\t' + "f2" + '\t' + "f3" + '\t' + "N_anc" + '\t' + "N_cur" + '\t' + "model" + '\t' + "N_fold" + '\t' + "g_factor" + '\t' + "pi_slope" + '\t' + "pi_intercept" + '\t' + "pi_max" + '\t' + "pi_predicted_max	pi_reduction	pi_predicted_reduction	pi_numbp50	pi_numbp75	pi_numbp90" + '\n')
result_50 = open("/scratch/pjohri1/" + folder + "/" + outputfile + "_sumstats_50.txt", 'w+')

l_stats = ["thetapi_m","thetaw_m","thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd"]
result_50.write("simID" + '\t' + "f0" + '\t' + "f1" + '\t' + "f2" + '\t' + "f3" + '\t' + "N_anc" + '\t' + "N_cur" + '\t' + "model" + '\t' + "N_fold" + '\t' + "g_factor")

for stat in l_stats:
	result_50.write('\t' + "func_" + stat + '\t' + "link_" + stat + '\t' + "neu_" + stat)
result_50.write('\n')

#get parameters:
d_para = {}
f_para = open("/home/pjohri1/demo_disc_parameters+5_SingExon_osg.txt", 'r')
for line in f_para:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    d_para[line2[0]] = line1
f_para.close()
#f_para = open("/home/pjohri1/demo_disc_parameters+5_SingExon_osg2.txt", 'r')
#for line in f_para:
#    line1 = line.strip('\n')
#    line2 = line1.split('\t')
#    d_para[line2[0]] = line1
#f_para.close()


i = 1
mark = 0
while i <= tot_numsims:
    simID = "sim" + str(i)
    s_para = ""
	#find out if all required files exist first, for this simID:
    s_pi = ""
    s_50 = ""
    #if 1 > 0:
    try:
        #f_para = open("/scratch/pjohri1/" + folder + "/parameters_" + simID + ".list", 'r')
        f_pi = open("/scratch/pjohri1/" + folder + "/" + simID + "_200.recoverypi", 'r')
        f_50 = open("/scratch/pjohri1/" + folder + "/" + simID + "_numbp50.bigwinsummary", 'r')
        #read the parameter file:
        s_para = d_para[simID]
        #for line in f_para:
        #    line1 = line.strip('\n')
        #    line2 = line1.split('\t')
        #    if line2[0] != "sim":
        #        s_para = line1
        #f_para.close()
        #write pi recovery stats:
        s_pi = s_para + '\t'
        for line in f_pi:
            line1 = line.strip('\n')
            line2 = line1.split('\t')
            if line2[0] != "slope":
                s_pi = s_pi + line1
        f_pi.close()

		#write numbp50 stats:
        s_50 = s_para + '\t'
        for line in f_50:
            line1 = line.strip('\n')
            line2 = line1.split('\t')
            if line2[0] != "functional":
                for x in line2[1:]:
                    if x != "NA":
                        s_50 = s_50 + '\t' + str("{:.5f}".format(float(x)))
                    else:
                        s_50 = s_50 + '\t' + x
        f_50.close()
        if s_pi != "" and s_50 != "":
            result_pi.write(s_pi + '\n')
            result_50.write(s_50 + '\n')
    #else:
    except:
        print (" File 50 and pi not found for:" + simID)
    i = i + 1

result_pi.close()
result_50.close()

print ("done")



