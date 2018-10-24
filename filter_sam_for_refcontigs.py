#!/usr/local/bin/python
# Python 2.7.6
# filter_sam_for_refcontigs.py
# a python script that filters .sam alignments piped from STDIN for the reference contig names also supplied by an input file, outputs .sam format; called by SNP_calling_and_vcf_filtering.sh
# Copyright (C) 2017, ETH Zurich, Mathias Scharmann
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 		
# If you use this code please cite:
#
# "Scharmann M, Grafe TU, Metali F, Widmer A. (2017) Sex-determination 
# and sex chromosomes are shared across the radiation of dioecious 
# Nepenthes pitcher plants. XXX"
# AND/OR
# "Scharmann M, Metali F, Grafe TU, Widmer A. (2017) Divergence with 
# gene flow among multiple sympatric carnivorous Nepenthes pitcher 
# plants is linked to trap morphology. XXX"
# 	
# contact: mathias.scharmann[-at-]env.ethz.ch or msph52[-at-]gmail.com



# usage example
# samtools view sample_TGCAT-RG.bam | python filter_sam_for_refcontigs.py --contiglist xxfu 

"""
reads STDIN from samtools view and returns only lines containing the reference contigs of names species fied in the contiglist

"""


import argparse
import os
import sys


########################## HEAD

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

# breaks script if non-UNIX linebreaks in input file popmap
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--contiglist", required=True, type=extant_file,
		help="name/path of contig list input file", metavar="FILE")
		
	args = parser.parse_args()
	
	# finish
	return args

################################## CORE


def read_contigfile(contiglist):
	
	contigs = set()
	with open(contiglist, "r") as INFILE:
		for line in INFILE:
			contigs.add(line.strip("\n"))
	
	return contigs


				
################################## MAIN


	
args = 	get_commandline_arguments ()

retain_contigs = read_contigfile(args.contiglist)

for line in sys.stdin:
	if line.split("\t")[2] in retain_contigs:
		sys.stdout.write(line)
	elif line.startswith("@"):
		if line.split("\t")[1].split(":")[1] in retain_contigs:
			sys.stdout.write(line)
		elif line.startswith("@RG") or line.startswith("@PG"): # retain read group in header!!
			sys.stdout.write(line)
			
#print args



	
	

