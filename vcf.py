from __future__ import print_function

from phylip import Phylip
from popmap import Popmap

import collections
import operator

class VCF():
	'Class for working with VCF files'

	def __init__(self, phy, popmap,outfile):
		lookup = {'Y':'C,T', 'R':'A,G', 'W':'A,T', 'S':'G,C', 'K':'T,G', 'M':'C,A', 'A':'A,A', 'T':'T,T', 'C':'C,C', 'G':'G,G'} #dictionary of ambiguity codesi
		
		f = open(outfile, 'w')

		f.write("##fileformat=VCFv4.1\n")
		f.write("##fileDate=20170603\n")
		f.write("##source=pyRAD.v.3.0.66\n")
		f.write("##reference=common_allele_at_each_locus\n")
		f.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n")
		f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
		f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
		f.write("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n")
		f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
		f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
		f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
		f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

		od = collections.OrderedDict(sorted(phy.alignment.items()))

		for key, value in od.iteritems():
			f.write("\t")
			f.write(key)
		f.write("\n")
		chromcount=1

		for i in xrange(0, phy.alignLength):
			d = collections.defaultdict(int)
			NS=0 #counter for number of samples at which locus is present
			for ind, seq in od.iteritems():
				if seq[i] != 'N' and seq[i] != '-':
					NS+=1
					templocus = lookup[seq[i]].split(',')
					for allele in templocus:
						d[allele]+=1
			
			sorted_d = sorted( d.items(), key=operator.itemgetter(1), reverse=True )
			
			temp = collections.defaultdict(int)
			nuclist = list()
			counter=0
			for key,value in sorted_d:
				nuclist.append(key)
				temp[key] += counter
				counter+=1
			f.write(str(chromcount))
			f.write("\t1\t.\t")
			f.write(nuclist.pop(0))

			remainder = ",".join(nuclist)
			f.write("\t")
			f.write(remainder)
			
			f.write("\t20\tPASS\tNS=")
			f.write(str(NS))
			f.write(";DP=15\tGT")
			
			for ind,seq in od.iteritems():
				f.write("\t")
				if seq[i] == 'N' or seq[i] == '-':
						f.write("./.")
				else:
					alleles = lookup[seq[i]].split(",")
					alleles[0] = str(temp[alleles[0]])
					alleles[1] = str(temp[alleles[1]])
					tempstring = "|".join(alleles)
					f.write(tempstring)
			f.write("\n")
			chromcount+=1
