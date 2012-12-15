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



#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------

true_vcf_name = sys.argv[1]
variant_file_lst = sys.argv[2:]

num_files = len(variant_file_lst)
WIDTH = 1/float(num_files+1+2)



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
    #plt.bar([0+start], [bad_snp], WIDTH, color='y', bottom=good_snp)


#------
# MAIN
#------

# get truth
true_snp_dict = create_variant_dict(true_vcf_name)
#weird_pos_lst = weird_bases(filter_lst,threshold)

common_snp_lst = []
#bad_snp_lst = []

for i in range(num_files):
    variant_file_name = variant_file_lst[i]
    print('\ncomparing: ' + variant_file_name)
    snp_dict = create_variant_dict(variant_file_name)
    num_common, bad_snps, false_pos, false_neg = compare_snp_dict(true_snp_dict, snp_dict)
    #if i == num_files-1:
    #    false_neg_array = lst_false_negs(true_snp_dict,snp_dict)
    #    weird_bases(filter_lst, threshold, false_neg_array)
    common_snp_lst.append(num_common)#-bad_snps)
    #bad_snp_lst.append(bad_snps)
    graph_vcf(i, false_pos, false_neg, COLORS[i])

# base numbers
graph_vcf(3, 1498, 13809 - 5235, COLORS[3])
common_snp_lst.append(47712)

# region numbers
graph_vcf(4, 1824, 7455, COLORS[4])
common_snp_lst.append(38232)

plt.xticks([x+WIDTH*(num_files+2)/2 for x in range(2)],['SNP false pos','SNP false neg'])
plt.legend(LEG_LST,['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%'],loc=1)
#plt.legend(LEG_LST,['004','005','006','007','008','010'],loc=1)
plt.axis([-WIDTH,2,0,16000])
#plt.show()
plt.savefig('variant_call_results_all.pdf',format='pdf')

# graph correct calls
plt.clf()
for i in range(len(common_snp_lst)):
    graph_common(i, common_snp_lst[i], COLORS[i])

plt.xticks([x+WIDTH*(num_files+2)/2 for x in range(1)],['correct SNPs'])
plt.legend(LEG_LST,['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%'],loc=1)
#plt.legend(LEG_LST,['004','005','006','007','008','010'],loc=1)
plt.axis([-WIDTH-0.1,1+0.9,0,67000])
#plt.show()
plt.savefig('correct_results_all.pdf',format='pdf')

