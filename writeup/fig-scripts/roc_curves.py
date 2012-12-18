#!/usr/bin/env python

import matplotlib.pyplot as plt

#---------
# GLOBALS
#---------

LEG_LST = []
COLORS = ['b','r','g','m','y','k']
#COLORS = ['b','r','g','m','c','k']

ERROR_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=14%']
CORRECT_LEG = ['BWA+GATK','SNAP+GATK','BWA+mpileup','BWA+BIGGIE,base=0.1','BWA+BIGGIE,region=14%']


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

def graph_common(i, common_snp, hi_complex, c):
    start = WIDTH*i
    plt.bar([0+start], [common_snp], WIDTH, color=c)
    plt.bar([0+start], [hi_complex], WIDTH, [common_snp], color=c, hatch='x')


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
        (52439, 1216, 9082, 6654), #new base=0.1
        (39719, 1799, 21802, 19933), #approximation of old region=5%
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
    graph_common(i, common_snp_lst[i], series[i][3], COLORS[i])

plt.xticks([x+WIDTH*num_files/2 for x in range(1)],['correct SNPs'])
plt.legend(LEG_LST,CORRECT_LEG[:num_files],loc=1)
plt.axis([-WIDTH-0.1,1+0.9+0.8,0,67000])
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

# graph roc curve for region results
plt.clf()
#thresholds = [4, 5, 6, 7, 8, 10]
#false_pos = [1430, 1825, 2095, 2265, 2490, 2765]
#false_neg = [6455, 7455, 8030, 8470, 8825, 9240]
#false_neg_hi_complex = [30095, 6654, 1301, 257, 99, 20]

thresholds = [8, 10, 12, 14, 16]
false_pos = [672, 1127, 1480, 1799, 2037]
#false_neg = [41874, 32399, 26169, 21802, 18748]
false_neg = [802, 1253, 1593, 1869, 2073]
false_neg_hi_complex = [41072, 31146, 24576, 19933, 16675]

p = plt.plot(thresholds, false_pos, '-o')
p = plt.plot(thresholds, false_neg, '-o')
#p = plt.plot(thresholds, false_neg_hi_complex, '-o')
#plt.ylim(ymax=10000)
#plt.ylim(ymin=1000)
LEG_LST2 = ['false positives', 'false negatives', 'high complexity false negatives']
plt.legend(LEG_LST2, loc=5)
plt.xlabel('percentage of complex bases in region')
plt.ylabel('# of SNPs')
plt.title('Accuracy vs percentage complex bases')
plt.savefig('region_accuracy_vs_thresh.pdf', format='pdf')

