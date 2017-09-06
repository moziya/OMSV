#!/bin/bash
if [ $# -ne 2 ]
then
	echo "(.osv format only) 2 parameters required: 1. input file, and 2. the threshold of SV size"
	exit
fi

./filterSV -inputSVFile $1 -outputFile temp_Filter.bed
cut -f4 temp_Filter.bed > temp.bed
rm -f temp_Filter.bed
paste $1 temp.bed > 1temp.osv
sed '/^#/ d' 1temp.osv > 2temp.osv
sed -i -e 's/sizeChange=-//g' 2temp.osv
sed -i -e 's/sizeChange=//g' 2temp.osv
awk -v size=$2 '( $8>=size && $NF==0)' 2temp.osv >  3temp.osv
awk 'NF{NF--};1' OFS="\t" 3temp.osv > 4temp.osv
awk '$8="sizeChange="$8' OFS="\t" 4temp.osv > 5temp.osv
sed -nr '/#/p' $1 > xtemp.osv
cat 5temp.osv >> xtemp.osv
mv xtemp.osv $1


rm -f temp.bed
rm -f *temp.osv
