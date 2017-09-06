#!/bin/bash
if [ $# -eq 0 ]
then
	folder="OMB_C6661_cm2"
	input_cmap="OMBlast_C6661.cmap"
	input_ref="common_hg38_r.cmap"
	output="align_complex.osv"
elif [ $# -eq 4 ]
then
	folder=$1
	input_cmap=$2
	input_ref=$3
	output=$4
elif [ $# -eq 1 ]
then
	folder=$1
	input_cmap="AssemblyAlign/assembly_q.cmap"
	input_ref="common_hg38_r.cmap"
        output="align_complex.osv"
else
	echo "There are 4 parameter requires: 1. input and output folder, 2. input alignment file, 3. input reference file, 4. output file name"
	exit
fi
input="split_alignment.oma"
java -jar OMTools0.25c+.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/$input_cmap --optresout $folder/$input --writeunmap false --filtermode 1 --postjoinmode 2 --clustermode 2 --thread 72

sed -i '/Unmap/d' $folder/$input

sort -n -k1,1 -k12,12 $folder/$input > $folder/Sort_$input
echo "1.Extract unmapped split maps!"
echo "g++ -std=c++11 extract_unmap_split.cpp -o extract_unmap_split"
./extract_unmap_split -optAlignFile $folder/Sort_$input -outputBeforeFile $folder/Split_map_before.cmap -outputAfterFile $folder/Split_map_after.cmap
echo "2.Re-align the unmapped split maps!"
java -jar OMTools0.25c+.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/Split_map_before.cmap --optresout $folder/Before_$input --writeunmap false --filtermode 1 --postjoinmode 2 --clustermode 2 --thread 72
java -jar OMTools0.25c+.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/Split_map_after.cmap --optresout $folder/After_$input --writeunmap false --filtermode 1 --postjoinmode 2 --clustermode 2 --thread 72
sort -n -k1,1 -k12,12 $folder/After_$input > $folder/Sort_After_$input
sort -n -k1,1 -k12,12 $folder/Before_$input > $folder/Sort_Before_$input
echo "3.Adjust the alignment of split maps!"
echo "g++ -std=c++11 adjust_split.cpp -o adjust_split"

./adjust_split -optAlignFile $folder/Sort_$input -optAlignBeforeFile $folder/Sort_Before_$input -optAlignAfterFile $folder/Sort_After_$input -outputAlignBefore $folder/Adjusted_Before_$input -outputAlignAfter $folder/Adjusted_After_$input
cp $folder/$input $folder/Total_$input
sed '/#/d' $folder/Adjusted_Before_$input >> $folder/Total_$input
sed '/#/d' $folder/Adjusted_After_$input >> $folder/Total_$input

echo "4.Call large inverstions and translocations!"
echo "g++ -std=c++11 breakpoint.cpp -o breakpoint"
echo "10	39687138	39935168" > density_PAR_region.bed
echo "12	34820121	37139973" >> density_PAR_region.bed
echo "17	22754856	23194891" >> density_PAR_region.bed
echo "23	10000	2781479" >> density_PAR_region.bed
echo "23	155701382	156030895" >> density_PAR_region.bed
echo "24	56887902	57217415" >> density_PAR_region.bed
echo "24	10000	2781479" >> density_PAR_region.bed
./breakpoint -optAlignFile $folder/Total_$input -SVoutputFile $folder/$output -highDensityFile density_PAR_region.bed

rm -f $folder/Adjusted_Before_$input
rm -f $folder/Adjusted_After_$input
rm -f $folder/Before_$input
rm -f $folder/After_$input
rm -f $folder/Sort_Before_$input
rm -f $folder/Sort_After_$input
rm -f $folder/Sort_$input

sed -i '/Unknown/d' $folder/$output
echo "#chr1	point1	chr2	point2	type	support" > $folder/Refine_$output
awk '$2<25&&$4<25' $folder/$output | sed '/#/d' | cut -f2-7 >> $folder/Refine_$output && mv $folder/Refine_$output $folder/$output
