# takes in snp file, regions file, and returns snps which are not in the regions
# i.e. region file = high complexity regions, then this returns low complexity snps
import sys

snp_file = sys.argv[1]
regions_file = sys.argv[2]

snps = open(snp_file, 'r')
regions = open(regions_file, 'r')
current_line = None
def next_pos():
    global current_line
    current_line = snps.readline()
    return int(current_line.split()[4])+1
current_pos = next_pos()
try:
    for line in regions:
        if line[0] == 'R':
            tokens = line.split()
            region_start = int(tokens[1])
            region_end = int(tokens[3])
            #print(region_start, " ", region_end)
            while current_pos < region_start:
                print current_line,
                current_pos = next_pos()
            while current_pos >= region_start and current_pos <= region_end:
                current_pos = next_pos()
        if line[0] == 'C':
            break
    while True:
        print current_line,
        current_pos = next_pos()
except:
    pass
