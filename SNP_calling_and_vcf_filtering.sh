# SNP_calling_and_vcf_filtering.sh
# a bash script that calls SNPs from .bam alignments and a reference.fasta using FREEBAYES then applies several quality filters to the raw.vcf outputs, mostly similar to dDocent
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


#########################################################
# SNP calling
#########################################################


Nintervals=59
NUMProc=20
FB1=20

# dependencies, should be in directory from which executed:
# filter_sam_for_refcontigs.py
# filter_fasta_per_ids.py

#Create list of BAM files
ls *-RG.bam >bamlist.list

#######

echo "preparing intervals"

grep ">" reference.fasta | sed s/">"// > reference.contigs

NumMaps=$(cat reference.contigs | wc -l)
IntervalSize=$(( $NumMaps / $Nintervals )) # original: 50 split makes maximum of 676 suffixes -> no more plitting higher than that 

split -l $IntervalSize reference.contigs splitmap
ls splitmap* > splitlist

COUNTER=0
for interval in $(ls splitmap*);
do
let COUNTER=COUNTER+1 
mv $interval interval.${COUNTER}.txt
done


echo "making interval BAMs"
# extract with python script simply the refrence contigs of an interval; don't do all the BED and mapping stuff:
parallel -j $NUMProc --noswap --no-notice 'samtools view -h {1} | python filter_sam_for_refcontigs.py --contiglist {2} | samtools view -bS - > {1}.{2}.bam' ::: $(ls *-RG.bam) ::: $(ls interval*)

echo "indexing interval BAMs"
# index the BAMs
parallel -j $NUMProc --noswap --no-notice 'samtools index {1}.{2}.bam' ::: $(ls *-RG.bam) ::: $(ls interval*)


### subset also the reference.fasta to the same intervals:

parallel -j $NUMProc --noswap --no-notice python filter_fasta_per_ids.py reference.fasta interval.{}.txt interval.{}.reference.fasta ::: $(seq 1 $(ls interval.*.txt | wc -l))

##### now run freebayes over these really intervalled .bam files;
# each instance needs its own bamlist on top of its own target list!
# this is achived by a process substitution: '<( ls *interval.{}.txt.bam )'
# process substitution has to be in quotes because otherwise it will be executed before the loop over the parallel ::: arguments and the {} will remain empty!

echo "running freebayes without HWE priors"
parallel -j $FB1 --noswap --no-notice --delay 1 freebayes -L '<( ls *interval.{}.txt.bam )' -v raw.{}.vcf -f interval.{}.reference.fasta -m 5 -q 5 -E 3 -G 3 --min-repeat-entropy 1 -V --hwe-priors-off -i -X -u ::: $(seq 1 $(ls interval.*.txt | wc -l))


echo "done with freebayes"

# the original dDocent contains a workaround to the problem that integer variables cannot be passed into the "range" thingy of bash: {1..500} is OK but {1..$somenumber} is not possible; the workaround hardcoded $NumBED to either 50 or 51... this is more flexible now; supply approximately desired $NumBED in line 364 as the divisor of $NumMaps; actual $NumBED will depend on "split" and deviate slightly.



## when done: clean up temp files:
rm interval.*.txt
rm interval.*.reference.fasta
rm interval.*.reference.fasta.fai
find . -maxdepth 1 -name "*bam.interval.*.txt.bam*" -print0 | xargs -0 rm
rm inbamlist.*
rm splitlist
rm reference.contigs
rm bamlist.list


#########################################################
# vcf_filtering stage_1
#########################################################
for i in {1..60} ; do 
echo $i
cat vcf_filters_stage1_GNU_parallel.sh | sed "s/FUFUXX/$i/g" > fu.${i}.sh
done

parallel -j 20 bash fu.{}.sh ::: {1..60}

rm fu.*.sh

#########################################################
# vcf_filtering stage_2
#########################################################
# concatenate and filter

# concatenate all partial vcfs:
mv 1.stage1_filtered.vcf 01.stage1_filtered.vcf
mv 2.stage1_filtered.vcf 02.stage1_filtered.vcf
mv 3.stage1_filtered.vcf 03.stage1_filtered.vcf
mv 4.stage1_filtered.vcf 04.stage1_filtered.vcf
mv 5.stage1_filtered.vcf 05.stage1_filtered.vcf
mv 6.stage1_filtered.vcf 06.stage1_filtered.vcf
mv 7.stage1_filtered.vcf 07.stage1_filtered.vcf
mv 8.stage1_filtered.vcf 08.stage1_filtered.vcf
mv 9.stage1_filtered.vcf 09.stage1_filtered.vcf

cat 01.stage1_filtered.vcf | grep "#" > TotalSNPs.stage1.vcf
cat *.stage1_filtered.vcf | grep -v "#" >> TotalSNPs.stage1.vcf

sleep 3

vcftools --vcf TotalSNPs.stage1.vcf --geno-depth

python ./vcf_drop_genotypes_exceeding_2std_indiv_mean_depth.py --vcf TotalSNPs.stage1.vcf --gdepth out.gdepth
rm out.gdepth

# use "--max-missing-count" argument of vcftools to include only RADtags that are covered in at least 2 individuals (4 chromosomes) of this set; otherwise the vcf file may contain un-informative sites. This is the ONLY site missingness filter applied:
# note that "--max-missing-count" refers to number of chromosomes, NOT number of individuals (found this out by trying & online)
n_chrom_tot=$(( $( head -100 TotalSNPs.stage1.vcf.indiv_cov_filtered.vcf | grep "sample" | sed 's/\t/\n/g' | grep -v "#" | grep "sample" | wc -l | awk '{print $1}')*2 ))
n_chrom_tot_minus_two=$(( $n_chrom_tot - 4 ))

vcftools --vcf TotalSNPs.stage1.vcf.indiv_cov_filtered.vcf --mac 1 --max-missing-count $n_chrom_tot_minus_two --recode 

mv out.recode.vcf TotalSNPs.popgen.final_filtered.not_excess-Het_filtered.2016-11-03.vcf

#############

rm *.stage1_filtered.vcf TotalSNPs.stage1.vcf TotalSNPs.stage1.vcf.indiv_cov_filtered.vcf





