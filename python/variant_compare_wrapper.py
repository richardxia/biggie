# wrapper for comparing true vcf to another

import sys
from variant_compare import *


#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------

true_vcf_name = sys.argv[1]
variant_file_name = sys.argv[2]


#------
# MAIN
#------


print('comparing: ' + true_vcf_name + ' and ' + variant_file_name)
true_snp_dict = create_variant_dict(true_vcf_name)
variant_snp_dict = create_variant_dict(variant_file_name)
num_common, bad_snps, false_pos, false_neg = compare_snp_dict(true_snp_dict, variant_snp_dict)
    


