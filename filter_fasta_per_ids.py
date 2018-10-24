# filter_fasta_per_ids.py
# a python script that filters .fasta files for the name of the entries; called by SNP_calling_and_vcf_filtering.sh
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


import sys

#Usage: filter_fasta_per_ids.py input.fasta filter_ids.txt output.fasta

input_file =sys.argv[1]
id_file =sys.argv[2]
output_file =sys.argv[3]
wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

input_fasta = {}
with open(input_file, "r") as INFILE:
	for line in INFILE:
		if line.startswith(">"):
			id = line.strip(">").strip("\n")
		elif line.startswith("N"):
			input_fasta[id] = line.strip("\n")

coreids = {}
for id in wanted:
	coreids[ id.split("_")[0] ] = id

the_order = sorted(coreids.keys(), key = int)	
the_order = [ coreids[x] for x in the_order[:] ]

print the_order
			
outlines = []
for id in the_order:
		outlines.append( ">"+id+"\n"+input_fasta[id] )

with open(output_file, "w") as OUTPUT:
	OUTPUT.write("\n".join(outlines))
