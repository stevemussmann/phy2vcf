#!/usr/bin/python

from comline import ComLine
from phylip import Phylip
from popmap import Popmap
from vcf import VCF

import sys



def main():
	input = ComLine(sys.argv[1:])
	phy = Phylip(input.args.phy)
	pops = Popmap(input.args.popmap)
	VCF(phy,pops,input.args.out)
main()

raise SystemExit
