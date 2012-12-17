# input: a file of snps (in either vcf or our format) and a file of weirdness scores, and a prefix for output files
# output: a set of files (one for each threshold) that filters out those above a certain score

import sys


#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------

all_snp_file_name = sys.argv[1]
weirdness_file_name = sys.argv[2]
prefix = sys.argv[3]

threshold_lst = [0.01,0.1,1,10,100]
file_lst = []
for t in threshold_lst:
    t_str = str(t).replace('.','_')
    file_lst.append(file(prefix + t_str + '.snps','w'))

#------
# MAIN
#------


all_snp_file = file(all_snp_file_name,'r')
weirdness_file = file(weirdness_file_name,'r')

all_pos_dict = {}

# true vcf file
if 'vcf' in all_snp_file_name:
    for line in all_snp_file:
        if line[0] != '#':
            tokens = line.split()
            pos = int(tokens[1])
            all_pos_dict[pos] = line
# our format
else:
    for line in all_snp_file:
        if line[0] == 'C':
            tokens = line.split()
            pos = int(tokens[4])+1
            all_pos_dict[pos] = line

print('finished with snp file')


for line in weirdness_file:
    if line[0] == 'p':
        tokens = line.split()
        pos = int(tokens[1])+1
        score = float(tokens[3])

        for i in range(len(threshold_lst)):
            t = threshold_lst[i]
            if score < t and pos in all_pos_dict:
                file_lst[i].write(all_pos_dict[pos])

    
all_snp_file.close()
weirdness_file.close()

for f in file_lst:
    f.close()
    


