# take in the true vcf, list of weirdness scores, another vcf and it's threshold

import sys
from variant_compare import *


#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------

true_vcf_name = sys.argv[1]
weirdness_file_name = sys.argv[2]
variant_file_name = sys.argv[3]
threshold = float(sys.argv[4])


#------
# MAIN
#------


true_snp_dict = create_variant_dict(true_vcf_name)
snp_dict = create_variant_dict(variant_file_name)
num_common, bad_snps, false_pos, false_neg = compare_snp_dict(true_snp_dict, snp_dict)
false_neg_array = lst_false_negs(true_snp_dict,snp_dict)
weird_bases(weirdness_file_name, threshold, false_neg_array)

