# vcf_filters_stage1_GNU_parallel.sh
# a bash script that performs several quality-filtering steps on raw.vcf files; called by SNP_calling_and_vcf_filtering.sh
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

## the original dDocent comes with an MIT licence which is included here:
# The MIT License (MIT)
# 
# Copyright (c) 2013 jpuritz
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


vcftools --vcf raw.FUFUXX.vcf --minDP 3 --mac 1 --recode --recode-INFO-all --out FUFUXX.tmp1

## remove any sites that were also mapped in the blank sample ( sample_blank-blank-ddH2o ) :
vcftools --vcf FUFUXX.tmp1.recode.vcf --indv sample_blank-blank-ddH2o --site-depth --out FUFUXX
cat FUFUXX.ldepth | awk '{if ($3 > 0) print $1}' | sort | uniq | grep -v "CHROM" > contig_blacklist.FUFUXX.txt

python ./vcf_subset_by_chromosome.py --vcf FUFUXX.tmp1.recode.vcf --exclude_chrom contig_blacklist.FUFUXX.txt
mv FUFUXX.tmp1.recode.vcfsubset.vcf FUFUXX.tmp2.recode.vcf 
vcftools --vcf FUFUXX.tmp2.recode.vcf --remove-indv sample_blank-blank-ddH2o --recode --recode-INFO-all --out FUFUXX.tmp3

#Because RADseq targets specific locations of the genome, we expect that the allele balance in our data (for real loci) should be close to 0.5, or if fixed in some population, it must be close to 0 or close to 1 (depending on what the reference is!!)
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01 | AB > 0.99"  FUFUXX.tmp3.recode.vcf > FUFUXX.tmp4.recode.vcf

# qual/ depth
vcffilter -f "QUAL / DP > 0.25" FUFUXX.tmp4.recode.vcf > FUFUXX.tmp5.recode.vcf

# we must filter for excess heterozygosity, as there appear to be a lot of sites were all or nearly all indivduals (across the species!) are heterozygous == most most likely paralogous seqs!
## NO! For the sex-association, we EXPECT EXCESS HETEROZYGOSITY FOR HEMIZYGOUS LOCI!! So do not drop them in this case!!
# vcftools --vcf FUFUXX.tmp5.recode.vcf --recode --recode-INFO-all --out FUFUXX.tmp6
# vcftools --vcf FUFUXX.tmp6.recode.vcf --hardy --out FUFUXX

#cat FUFUXX.hwe | awk '{ if ($8 <= 0.05 ) print $1 }' | sort | uniq > excessively_heterozygous_refcontigs.FUFUXX.txt
#rm out.hwe

# python $HOME/tools/vcf_subset_by_chromosome.py --vcf FUFUXX.tmp6.recode.vcf --exclude_chrom excessively_heterozygous_refcontigs.FUFUXX.txt

# final mac1 checking:
vcftools --vcf FUFUXX.tmp5.recode.vcf --mac 1 --recode --recode-INFO-all --out FUFUXX.tmp7
mv FUFUXX.tmp7.recode.vcf FUFUXX.stage1_filtered.vcf
# returns: FUFUXX.stage1_filtered.vcf
# cleanup:
rm FUFUXX.tmp* contig_blacklist.FUFUXX.txt FUFUXX.ldepth FUFUXX.hwe excessively_heterozygous_refcontigs.FUFUXX.txt FUFUXX.log

#############

