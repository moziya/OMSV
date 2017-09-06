#!/bin/bash
if [ $# -eq 4 ]
then
	folder=$1
	input_cmap=$2
	input_ref=$3
	output=$4
else
	echo "There are 4 parameter requires: 1. input and output folder, 2. input alignment file, 3. input reference file, 4. output file name"
	exit
fi
input="split_alignment.oma"
java -jar OMTools/OMTools.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/$input_cmap --optresout $folder/$input --writeunmap false --filtermode 1 --alignmentjoinmode 2 --thread 72

sed -i '/Unmap/d' $folder/$input

sort -n -k1,1 -k12,12 $folder/$input > $folder/Sort_$input
echo "1.Extract unmapped split maps!"
./extract_unmap_split -optAlignFile $folder/Sort_$input -outputBeforeFile $folder/Split_map_before.cmap -outputAfterFile $folder/Split_map_after.cmap
echo "2.Re-align the unmapped split maps!"
java -jar OMTools/OMTools.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/Split_map_before.cmap --optresout $folder/Before_$input --writeunmap false --filtermode 1 --alignmentjoinmode 2 --thread 72
java -jar OMTools/OMTools.jar OMBlastMapper -refmapin $input_ref --optmapin $folder/Split_map_after.cmap --optresout $folder/After_$input --writeunmap false --filtermode 1 --alignmentjoinmode 2 --thread 72
sort -n -k1,1 -k12,12 $folder/After_$input > $folder/Sort_After_$input
sort -n -k1,1 -k12,12 $folder/Before_$input > $folder/Sort_Before_$input
echo "3.Adjust the alignment of split maps!"

./adjust_split -optAlignFile $folder/Sort_$input -optAlignBeforeFile $folder/Sort_Before_$input -optAlignAfterFile $folder/Sort_After_$input -outputAlignBefore $folder/Adjusted_Before_$input -outputAlignAfter $folder/Adjusted_After_$input
cp $folder/$input $folder/Total_$input
sed '/#/d' $folder/Adjusted_Before_$input >> $folder/Total_$input
sed '/#/d' $folder/Adjusted_After_$input >> $folder/Total_$input
mv $folder/Total_$input $folder/$output

rm -f $folder/Adjusted_Before_$input
rm -f $folder/Adjusted_After_$input
rm -f $folder/Before_$input
rm -f $folder/After_$input
rm -f $folder/Sort_Before_$input
rm -f $folder/Sort_After_$input
rm -f $folder/Sort_$input
rm -f $folder/$input
