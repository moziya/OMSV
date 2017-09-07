step=0
if [ $# -eq 1 ]
then
	step=$1
fi

if [ $step -lt 1 ]
then

echo "1. Start to call large indels and site of NA12878"
mkdir temp

./OMSV -inputLabel 12878 -outputFolder SV_result/ -SVoutputFile Indel -chrMapFile hg38_r.cmap -optAlignFile data/NA12878_700bp_hg38_combRefOMB2.oma -optTempFolder temp/

sed '/.site\|#/!d' SV_result/12878Indel_2000.osv > SV_result/12878Sites_2000.osv
sed '/Homozygous\|Heterozygous/d' SV_result/12878Indel_2000.osv > SV_result/12878Mixed_indel_2000.osv
sed -i '/.site\|HDI\|HI1I2\|HD1D2/d' SV_result/12878Indel_2000.osv

rm -f SV_result/*_List.txt
./postFilter.sh SV_result/12878Indel_2000.osv 2000

rm -r -f temp
fi

echo "2. Start to call medium-size inversions"
if java -jar */OMTools.jar >/dev/null 2>&1
then
        
	java -jar */OMTools.jar SVDetection --refmapin hg38_r.cmap --optresin data/NA12878_700bp_hg38_combRefOMB2.oma --mininvsig 4 --svout 12878Med_inv.osv --flanksig 0 --deg 0 -svmode 2 -minsupport 10 && sed -i '/.site/d' 12878Med_inv.osv
else
        echo "Medium-size inversion caller requires java and OMTools but not installed. Skip this step."
fi

echo "3. Start to call large complex SVs"
./complex_caller.sh data/NA12878_700bp_hg38_OMB_split.oma . 12878

echo "4. Start to call CNVs"
if command -v matlab >/dev/null 2>&1
then
	
	./CNV_caller.sh data/NA12878_700bp_hg38_combRefOMB2.oma hg38_r.cmap . 12878
else
	echo "CNV caller requires 'matlab' but not installed. Skip this step."
fi

echo "5. Compile complex SV list"
./compile_complex.sh 12878Med_inv.osv 12878Duplication_2.bed 12878Complex.bed SV_result/12878Complex_total.bed
