# compare true vcf and a candidate vcf (or our SNP format)

import numpy as np

# SNP class, holds list of alleles
class SNP:

    def __init__(self, *args):
        if len(args) == 1:
            self.allele_set = self.initialize_from_lst(*args)
        else:
            self.allele_set = self.initialize_from_vars(*args)

    def initialize_from_vars(self, ref, alt, genotype):
        allele_set = set([])
        if genotype[0] == '0' or genotype[2] == '0':
            allele_set.add(ref)
        if genotype[0] == '1' or genotype[2] == '1':
            allele_set.add(alt)
        return allele_set

    def initialize_from_lst(self, allele_lst):
        allele_set = set(allele_lst)
        return allele_set

    def equals_SNP(self, other_SNP):
        if self.allele_set.issubset(other_SNP.allele_set) and other_SNP.allele_set.issubset(self.allele_set):
            return True
        return False
    

# create SNP/indel lists from vcf file name
def create_variant_dict(file_name):
    if 'vcf' in file_name:
        # true vcf format
        return create_variant_dict_vcf(file_name)
    else:
        # our format
        return create_variant_dict_biggie(file_name)


# create a dictionary of positions (keys) and set of alleles at that SNP
def create_variant_dict_vcf(vcf_file_name):
    vcf_file = file(vcf_file_name,'r')
    snp_dict = {}

    for line in vcf_file:
        if line[0] != '#':
            tokens = line.split()
            pos = int(tokens[1])
            ref_seq = tokens[3]
            alt_seq = tokens[4]
            genotype = tokens[9]

            # call SNP 
            if len(ref_seq) == 1 and len(alt_seq) == 1:
                snp_dict[pos] = SNP(ref_seq,alt_seq,genotype)

    vcf_file.close()
    return snp_dict

def create_variant_dict_biggie(biggie_file_name):
    biggie_file = file(biggie_file_name,'r')
    snp_dict = {}

    for line in biggie_file:
        if line[0] == 'C':
            tokens = line.split()
            pos = int(tokens[4])+1
            alt_seq = tokens[2]
            snp_dict[pos] = SNP([alt_seq[0],alt_seq[2]])

    biggie_file.close()
    return snp_dict


# compare two snp dicts
def compare_snp_dict(snp_ref, snp_alt):
    set_ref = set(snp_ref.keys())
    set_alt = set(snp_alt.keys())

    intersection = set_alt.intersection(set_ref)
    inter_lst = list(intersection)
    num_common = len(intersection)

    # double check the SNP caller got these common bases correct
    bad_snps = 0
    for i in range(num_common):
        base = inter_lst[i]
        if not snp_ref[base].equals_SNP(snp_alt[base]):
            bad_snps += 1

    false_pos = len(snp_alt) - num_common
    false_neg = len(set_ref) - num_common

    print('SNPs')
    print('num common ' + str(num_common))
    print('percent bad ' + str(bad_snps/float(num_common)))
    print('false pos ' + str(false_pos))
    print('false neg ' + str(false_neg))
    return num_common, bad_snps, false_pos, false_neg


# generate file of false negs
def lst_false_negs(snp_ref, snp_alt):
    num_common, bad_snps, false_pos, false_neg = compare_snp_dict(snp_ref, snp_alt)
    
    false_neg_array = np.zeros(false_neg)
    snp_ref_lst = list(snp_ref.keys())
    snp_alt_lst = list(snp_alt.keys())
    false_neg_index = 0
    
    for i in range(len(snp_ref_lst)):
        pos = snp_ref_lst[i]
        # false neg 
        if pos not in snp_alt_lst:
            false_neg_array[false_neg_index] = pos
            false_neg_index += 1

    return false_neg_array


# obtain lst of weird bases from file
def weird_bases(file_name, threshold, false_neg_array):
    f = file(file_name,'r')
    num_uncalled_neg = 0

    for line in f:
        if line[0] == 'p':
            tokens = line.split()
            pos = int(tokens[1])+1
            score = float(tokens[3])

            if score > threshold and pos in false_neg_array:
                print('uncalled false neg ' + str(pos))
                num_uncalled_neg += 1
    print('num uncalled false neg ' + str(num_uncalled_neg))
    return num_uncalled_neg
