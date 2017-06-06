from __future__ import print_function

import collections
from collections import defaultdict

from phylip import Phylip
import operator

class VCF():
	'Class for working with VCF files'

	def __init__(self, phy):
		lookup = {'Y':'C,T', 'R':'A,G', 'W':'A,T', 'S':'G,C', 'K':'T,G', 'M':'C,A', 'A':'A,A', 'T':'T,T', 'C':'C,C', 'G':'G,G'} #dictionary of ambiguity codes
		print("##fileformat=VCFv4.1")
		print("##fileDate=20170603")
		print("##source=pyRAD.v.3.0.66")
		print("##reference=common_allele_at_each_locus")
		print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")
		print("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">")
		print("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
		print("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">")
		print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
		print("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">")
		print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")
		print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", end='')

		od = collections.OrderedDict(sorted(phy.alignment.items()))

		for key, value in od.iteritems():
			print("\t", end='')
			print(key, end='')
		print("")

		chromcount=1

		for i in xrange(0, phy.alignLength):
			#print(i)
			d = defaultdict(int)
			NS=0 #counter for number of samples at which locus is present
			for ind, seq in od.iteritems():
				if seq[i] != 'N':
					NS+=1
					if seq[i] == 'A' or seq[i] == 'T' or seq[i] == 'G' or seq[i] == 'C':
						d[seq[i]] +=2
					elif seq[i] == 'Y':
						d['C'] +=1
						d['T'] +=1
					elif seq[i] == 'R':
						d['A'] +=1
						d['G'] +=1
					elif seq[i] == 'W':
						d['A'] +=1
						d['T'] +=1
					elif seq[i] == 'S':
						d['G'] +=1
						d['C'] +=1
					elif seq[i] == 'K':
						d['T'] +=1
						d['G'] +=1
					elif seq[i] == 'M':
						d['C'] +=1
						d['A'] +=1
			
			sorted_d = sorted( d.items(), key=operator.itemgetter(1), reverse=True )
			#print(sorted_d)
			
			temp = defaultdict(int)
			nuclist = list()
			counter=0
			for key,value in sorted_d:
				nuclist.append(key)
				temp[key] += counter
				counter+=1
			print(chromcount, end='')
			print("\t1\t.\t", end='')
			#print(temp)
			print(nuclist.pop(0), end='')

			remainder = ",".join(nuclist)
			print("\t",remainder, end='')
			
			print("\t20\tPASS\tNS=", end='')
			print(NS, end='')
			print(";DP=15\tGT", end='')
			
			for ind,seq in od.iteritems():
				print("\t",end='')
				if seq[i] == 'N' or seq[i] == '-':
						print("./.", end='')
				else:
					alleles = lookup[seq[i]].split(",")
					alleles[0] = str(temp[alleles[0]])
					alleles[1] = str(temp[alleles[1]])
					#print(lookup[seq[i]], end='')
					tempstring = "|".join(alleles)
					print(tempstring, end='')


			chromcount+=1
			print("")

			

