#This script is to subsample from an ms file to create another ms file.
#python subsample_ms.py output1.ms 2 
import sys
import random

f_ms = open(sys.argv[1], 'r')
#if len(f_ms.readlines(  )) > 3:
    #result = open(sys.argv[1].replace(".ms", "_" + str(num_samples) + "indv.ms"), 'w+')
result = open(sys.argv[2], 'w+')
num_samples = int(sys.argv[3])#number of samples to output for the psmc analysis

i = 1
l_lines = []
while i <= 100:
    l_lines.append(i)
    i = i + 1
l_samples = random.sample(l_lines, num_samples)

#f_ms.seek(0)
sample_num = 1
for line in f_ms:
    line1 = line.strip('\n')
    line2 = line1.split()
    if "positions" in line or "//" in line or "segsites" in line:
        result.write(line)
    else:
        if sample_num in l_samples:
            result.write(line)
        sample_num = sample_num + 1
result.close()
f_ms.close()

print ("done")
