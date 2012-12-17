# take in the true vcf and any other vcfs and graph false neg/false pos for SNPs and indels

import optparse
import sys
import matplotlib.pyplot as plt
from variant_compare import *


#---------
# GLOBALS
#---------

LEG_LST = []
COLORS = ['b','r','g','m','y','k']

ERROR_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%']
CORRECT_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%']


#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------

all_file_lst = sys.argv[1:]
assert (len(all_file_lst) % 2 == 0)

# separate files out into true and variant files (every other one)
num_files = len(all_file_lst)/2
print('num files ' + str(num_files))
true_file_lst = []
variant_file_lst = []
for i in range(num_files):
    true_file_lst.append(all_file_lst[2*i])
    variant_file_lst.append(all_file_lst[2*i+1])

WIDTH = 1/float(num_files+1)



#---------
# HELPERS
#---------

# graph a file
def graph_vcf(i, snp_pos, snp_neg, c):
    start = WIDTH*i
    p = plt.bar([0+start], [snp_pos], WIDTH, color=c)
    plt.bar([1+start], [snp_neg], WIDTH, color=c)
    LEG_LST.append(p[0])

def graph_common(i, common_snp, c):
    start = WIDTH*i
    plt.bar([0+start], [common_snp], WIDTH, color=c)


#------
# MAIN
#------


common_snp_lst = []

for i in range(num_files):
    true_file_name = true_file_lst[i]
    variant_file_name = variant_file_lst[i]
    print('\ncomparing: ' + true_file_name + ' ' + variant_file_name)
    true_snp_dict = create_variant_dict(true_file_name)
    snp_dict = create_variant_dict(variant_file_name)
    num_common, bad_snps, false_pos, false_neg = compare_snp_dict(true_snp_dict, snp_dict)
    common_snp_lst.append(num_common)
    graph_vcf(i, false_pos, false_neg, COLORS[i])

plt.xticks([x+WIDTH*num_files/2 for x in range(2)],['SNP false pos','SNP false neg'])
plt.legend(LEG_LST,ERROR_LEG[:num_files],loc=1)
plt.axis([-WIDTH,2,0,16000])
plt.savefig('variant_call_results_all.pdf',format='pdf')

# graph correct calls
plt.clf()
for i in range(len(common_snp_lst)):
    graph_common(i, common_snp_lst[i], COLORS[i])

plt.xticks([x+WIDTH*num_files/2 for x in range(1)],['correct SNPs'])
plt.legend(LEG_LST,CORRECT_LEG[:num_files],loc=1)
plt.axis([-WIDTH-0.1,1+0.9,0,67000])
plt.savefig('correct_results_all.pdf',format='pdf')

