from __future__ import print_function

class Popmap():
	'Class for parsing a popmap'

	def __init__(self, infile):
		data = open(infile, 'r')
		content = data.readlines()
		content = [x.rstrip('\n') for x in content]

		self.popmap = dict()

		for line in content:
			temp = line.split()
			self.popmap[temp[0]] = temp[1]
