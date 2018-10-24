# assembly_and_mapper_script.sh
# a bash script for de-novo assembly of RAD-tags (dDocent-style) from demultiplexed, single-end fastq files and mapping of samples to the de-novo reference
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


# start dDocent-like assembly:
# concatenate everything
zcat *.fq.gz | mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' > concat.fasta

#Convert fasta format to just sequences
mawk '!/>/' concat.fasta > concat.seq

#Find every unique sequence and count the number of occurrences
perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' concat.seq > uniq.seqs

#Create a data file with the number of unique sequences and the number of occurrences
ONE=$(mawk '$1 >= 2' uniq.seqs | wc -l)
echo -e "2""\t""$ONE" > uniqseq.data

for ((i = 3; i <= 50; i++));
do
J=$(mawk -v x=$i '$1 >= x' uniq.seqs | wc -l)
echo -e "$i""\t""$J" >> uniqseq.data
done



#Function to convert a file of unique sequences to fasta format
uniq2fasta()
{
i=1
cat $1 | while read line
do
echo ">Contig"$i
echo $line
i=$(($i + 1))
done
}


#Main function for assembly
Assemble()
{

# main arguments: $CUTOFF, $vsearch_id, $simC
# secondary arguments: $NUMProc

mawk -v x=$CUTOFF '$1 >= x' uniq.seqs | cut -f 2 > totaluniqseq


#Convert reads to fasta
uniq2fasta totaluniqseq > uniq.fasta

# instead of rainbow which wants PE data, we cluster with vsearch like pyRAD does:
vsearch --fasta_width 0 --threads $NUMProc --id $vsearch_id --cluster_fast uniq.fasta --msaout vsearch_msaout.txt

# instead of rainbow merge and best contig selection (which do not work with single-end data
# use my python script:
python ./vsearch_msaout_global_align_muscle.v2.py --vsearch_msa vsearch_msaout.txt

#cd-hit to cluster reads based on sequence similarity
cd-hit-est -i vsearch_muscle.fasta -o referenceRC.fasta -M 0 -T $NUMProc -c $simC &>cdhit.log

#seqtk seq -r referenceRC.fasta > reference.fasta.original
#rm referenceRC.fasta

mv referenceRC.fasta reference.fasta.original

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta

samtools faidx reference.fasta
bwa index reference.fasta

}


CUTOFF=10
NUMProc=12
simC=0.9
vsearch_id=0.9

Assemble

echo Done!

######

######################
# mapping with bwa and SNP calling in freebayes (ignore invariant sites!)
######################

for i in $( ls *.fq.gz | sed 's/.fq.gz//g' ) ; do

bwa mem -a reference.fasta $i.fq.gz -L 20,5 -t $NUMProc -M -T 10 -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
samtools sort -@$NUMProc $i.bam -o $i.sorted
rm $i.bam
mv $i.sorted $i-RG.bam
samtools index $i-RG.bam

done

rm *.bam.log bwa.*.log
rm concat.fasta concat.seq uniq.seqs totaluniqseq uniq.fasta vsearch_msaout.txt vsearch_muscle.fasta


