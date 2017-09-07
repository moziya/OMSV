#!/bin/bash
if [ $# -ne 4 ]
then
	echo "There are 4 parameters: "
	echo "	1. the file of medium inversions (called by SVDetection in OMBlast),"
	echo "	2. the file of CNVs (called by duplication caller based on Matlab), and "
	echo "	3. the file of large compex SVs (called by align_complex_caller.sh)"
	echo "	4. the output file"
	exit
fi

echo "##Compiled from three complex callers: medium inversion, copy-number variation, and split-alignment caller" > $4

echo "#chr1	point1	chr2	point2	type	size	sample" >> $4

if [ -f $1 ]
then
	sed '/#/d' $1 | awk '$2=$2"\t"$1' OFS="\t" | awk -v var="$sample" '$6=$4-$2"\t"var' OFS="\t" | cut -f1-7 | sed 's/Inversion/Medium-Inversion/g' >> $4
	rm $1
else
	echo "The 1st input file (Medium-size inversion) is not exist."
fi
if [ -f $2 ]
then
	awk '$2=$2"\t"$1' OFS="\t" $2 | sed '/#/d' | sed 's/sizeChange=//g' | sed 's/Duplication/CNV/g' | cut -f1-5,9 | awk -v var="$sample" '$6=$6"\t"var' OFS="\t" >>  $4
	rm $2
else
	echo "The 2nd input file (CNV) is not exist."
fi
if [ -f $3 ]
then
	sed '/#/d' $3 | sed 's/Inversion/Large-Inversion/g' | awk -v var="$sample" '$6=$4-$2"\t"var' OFS="\t" >> $4
	rm $3
else
	echo "The 3rd input file (Large-inversion) is not exist."
fi



sort -n -k1,1 -k3,3 -k2,2 -k4,4 $4 | awk '$5=="Inter-Translocation"||$5=="Intra-Translocation"{$6="nan"}1' OFS="\t" > ${4}_tp && mv ${4}_tp $4

reg_file="density_PAR_region.bed"
len=$(wc -l < $reg_file)
for((i=1;i<=len;i++))
do
	chr="$(awk -v var="$i" -F'\t' 'NR==var {print $1}' $reg_file)"
	start="$(awk -v var="$i" -F'\t' 'NR==var {print $2}' $reg_file)"
	stop="$(awk -v var="$i" -F'\t' 'NR==var {print $3}' $reg_file)"
	awk -v ch=$chr -v st=$start -v sp=$stop '!(($1==ch&&$2>=st&&$2<=sp)||($3==ch&&$4>=st&&$4<=sp))' OFS="\t" $4 > ${4}_tp
	mv ${4}_tp $4
done



