#!/bin/bash
if [ $# -ne 3 ]
then
	echo "3 parameters required: input alignment, output folder and output label"
	exit
fi

./breakpoint -optAlignFile $1 -SVoutputFile $2/${3}Complex.bed -highDensityFile density_PAR_region.bed

echo "#chr1	point1	chr2	point2	type	support" > $2/Refine_complex.bed
awk '$2<25&&$4<25' $2/${3}Complex.bed | sed '/#/d' | cut -f2-7 >> $2/Refine_complex.bed && mv $2/Refine_complex.bed $2/${3}Complex.bed
