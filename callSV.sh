step=0
if [ $# -eq 1 ]
then
	step=$1
fi

if [ $step -lt 1 ]
then

echo "1. Start to call large indels"
mkdir temp
./OMSV -inputLabel 12878 -outputFolder SV_result/ -SVoutputFile Indel -chrMapFile hg38_r.cmap -optAlignFile data/NA12878_700bp_hg38_combRefOMB2.oma -optTempFolder temp/
rm -f *_List.txt
./postFilter.sh SV_result/12878Indel_2000.osv 2000

echo "2. Start to call mixed indels"
./OMSV_mixedIndel -inputLabel 12878 -outputFolder SV_result/ -SVoutputFile Mixed_indel -chrMapFile hg38_r.cmap -optAlignFile data/NA12878_700bp_hg38_combRefOMB2.oma -optTempFolder temp/
./postFilter.sh SV_result/12878Mixed_indel_2000.osv 2000
sed '/Homozygous/d' SV_result/12878Mixed_indel_2000.osv | sed '/Heterozygous/d' > SV_result/12878Mixed_indel_2000.osv_tp && mv SV_result/12878Mixed_indel_2000.osv_tp SV_result/12878Mixed_indel_2000.osv


rm -r -f temp
fi


echo "2. Start to call medium-size inversions"
java -jar OMTools/OMTools.jar SVDetection --refmapin hg38_r.cmap --optresin data/NA12878_700bp_hg38_combRefOMB2.oma --mininvsig 4 --svout 12878Med_inv.osv --flanksig 0 --deg 0 -svmode 2 -minsupport 10 && sed -i '/.site/d' 12878Med_inv.osv

echo "3. Start to call large complex SVs"
./complex_caller.sh data/NA12878_700bp_hg38_OMB_split.oma . 12878

echo "4. Start to call CNVs"
./CNV_caller.sh data/NA12878_700bp_hg38_combRefOMB2.oma hg38_r.cmap . 12878

echo "5. Compile complex SV list"
./compile_complex.sh 12878Med_inv.osv 12878Duplication_2.bed 12878Complex.bed SV_result/12878Complex_total.bed
