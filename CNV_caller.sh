#!/bin/bash
if [ $# -ne 4 ]
then
	echo "This duplication caller is developed with Matlab, please make sure the matlab is runable."
	echo "	There are three parameters required: 1. alignment file (.oma), 2. reference map (.cmap), 3. output folder and 4. output label."
	exit
fi
align_file=$1
ref_file=$2
outFold=$3
lab=$4

cut -f1,4,6,7,12,13 $1 | sed '/#/d' | sort -n -k1,1 -k5,5 -k6,6 | awk '($3>10&&$4>1)||($3>20&&$4<1)' > $outFold/cov_mat2.txt
cut -f1,6 $2 | sed '/#/d' | sort -n -k1,1 -k2,2 > $outFold/ref_mat2.txt
matlab -nodisplay -nodesktop -r	"try cal_dup('$outFold/ref_mat2.txt','$outFold/cov_mat2.txt','$outFold/${lab}Duplication.bed',2); exit; catch exit; end;" 
rm -f $outFold/cov_mat2.txt
rm -f $outFold/ref_mat2.txt
