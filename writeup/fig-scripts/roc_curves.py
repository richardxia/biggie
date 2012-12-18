#!/usr/bin/env python

import matplotlib.pyplot as plt

#---------
# GLOBALS
#---------

LEG_LST = []
COLORS = ['b','r','g','m','y','k']
#COLORS = ['b','r','g','m','c','k']

ERROR_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%']
CORRECT_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=5%']


num_files = 5
WIDTH = 1/float(num_files+1)



#---------
# HELPERS
#---------

# graph a file
def graph_vcf(i, snp_pos, snp_neg, hi_complex, c):
    start = WIDTH*i
    p = plt.bar([0+start], [snp_pos], WIDTH, color=c)
    plt.bar([1+start], [snp_neg-hi_complex], WIDTH, color=c)
    # plot false positives in region we don't call in hatched color
    plt.bar([1+start], [hi_complex], WIDTH, [snp_neg-hi_complex], color=c, hatch='x')
    LEG_LST.append(p[0])

def graph_common(i, common_snp, c):
    start = WIDTH*i
    plt.bar([0+start], [common_snp], WIDTH, color=c)


#------
# MAIN
#------


common_snp_lst = []

# 1.44cm = 10000
series = [
        (60800, 5780, 700, 0),
        (60930, 14670, 605, 0),
        (59380, 1660, 2135, 0),
        #(47710, 1500, 8575, 0), #approximation of old base=0.1
        (47710, 1216, 9082, 6654), #new base=0.1
        (38230, 1830, 7453, 0), #approximation of old region=5%
        ]

for i, data in enumerate(series):
    num_common, false_pos, false_neg, hi_complex = data
    common_snp_lst.append(num_common)
    graph_vcf(i, false_pos, false_neg, hi_complex, COLORS[i])

plt.xticks([x+WIDTH*num_files/2 for x in range(2)],['SNP false pos','SNP false neg'])
plt.legend(LEG_LST,ERROR_LEG[:num_files],loc=1)
plt.axis([-WIDTH,2,0,16000])
plt.ylabel("# of SNPs")
plt.title('Error rate of various SNP callers')
plt.savefig('variant_call_results_all.pdf',format='pdf')

# graph correct calls
plt.clf()
for i in range(len(common_snp_lst)):
    graph_common(i, common_snp_lst[i], COLORS[i])

plt.xticks([x+WIDTH*num_files/2 for x in range(1)],['correct SNPs'])
plt.legend(LEG_LST,CORRECT_LEG[:num_files],loc=1)
plt.axis([-WIDTH-0.1,1+0.9,0,67000])
plt.ylabel("# of SNPs")
plt.title('Accuracy of various SNP callers')
plt.savefig('correct_results_all.pdf',format='pdf')

# graph roc curve for per-base results
plt.clf()

thresholds = [0.01, 0.1, 1, 10, 100]
false_pos = [27, 1216, 3949, 4464, 4480]
false_neg = [24, 2428, 3285, 3778, 3887]
false_neg_hi_complex = [30095, 6654, 1301, 257, 99]
common_snps = [31202, 52439, 56935, 57486, 57535]
plt.xscale('log')
p = plt.plot(thresholds, false_pos, '-o')
p = plt.plot(thresholds, false_neg, '-o')
p = plt.plot(thresholds, false_neg_hi_complex, '-o')
plt.ylim(ymax=7000)
#p = plt.plot(thresholds, common_snps, '-o')
LEG_LST2 = ['false positives', 'false negatives', 'high complexity false negatives']
plt.legend(LEG_LST2)
plt.xlabel('complexity score threshold')
plt.ylabel('# of SNPs')
plt.title('Accuracy vs complexity score threshold')
plt.savefig('base_accuracy_vs_thresh.pdf', format='pdf')
