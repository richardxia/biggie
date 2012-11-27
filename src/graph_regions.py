# takes in a file specifying high complexity regions and graphs them


import optparse
import matplotlib.pyplot as plt
import sys



#------------------------------
# PARSE COMMAND LINE ARGUMENTS
#------------------------------


parser = optparse.OptionParser(description='creates a plot of the likelihood over the EM iterations')

parser.add_option('-i', '--region_file_name', type='string', help='path to input file of high complexity regions')

(opts, args) = parser.parse_args()

mandatories = ['region_file_name']
for m in mandatories:
    if not opts.__dict__[m]:
        print 'mandatory option ' + m + ' is missing\n'
        parser.print_help()
        sys.exit()



#---------
# HELPERS
#---------


# extract the start and end of the interval (does not include end)
def parse_line(line):
    tokens = line.split()
    start = int(tokens[7][1:-1])
    end = int(tokens[8][:-1]) + 1
    return start, end



#------
# MAIN
#------


region_file = file(opts.region_file_name,'r')
for line in region_file:
    if line[:22] == 'High complexity region':
        start, end = parse_line(line)
        #if start < 20000000:# and start < 40000000:
        plt.plot(range(start,end), [1]*(end-start), lw=4)

region_file.close()
plt.xlabel('genome position')
plt.ylabel('high complexity')
plt.title('High Complexity Regions, chr 20, all')
plt.show()

