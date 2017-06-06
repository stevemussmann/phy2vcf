from __future__ import print_function

class Phylip():
	'Class for working with phylip files'

	def __init__(self, infile):
		phy = open(infile, 'r')
		content = phy.readlines()
		content = [x.rstrip('\n') for x in content] #remove newline characters
		
		header = content.pop(0).split() #get the header and split on whitespace

		#get the number of individuals and alignment length
		self.numInds = int(header.pop(0))
		self.alignLength = int(header.pop(0))

		#declare a dictionary to hold alignment
		self.alignment = dict()
		
		for line in content:
			temp = line.split()
			self.alignment[temp[0]] = temp[1]
	
