# OMSV TUTORIAL
=================================================================================
## Overall
Here we provide the package of OMSV with all source codes (C++, bash script, and Matlab codes) and binary files. The scrips of OMSV have been tested in Debian GNU/Linux 9.0 (stretch) and CentOS Linux release 7.3.1611 (Core) platform and the matlab code was developed in Matlab R2011b (7.13.0.564) 64-bit (glnxa64). The way to start our pipeline within a few steps and the details of the programs, including input parameters and supported file formats, are shown in "Readme\_OMSV.txt".  
**All data** (raw optical mapping data and alignments) are accessible at [here](https://drive.google.com/drive/folders/0B4T3OTv54s2DbDdLRFdqdERXZkE). The alignment file, as well as the 2-round split alignment, of NA12878 are provided for the quick start. 

=================================================================================
## Quick Start

To quickly try our OMSV Caller, please first download the package and alignment files of OM data (e.g. NA12878\_alignment.tar.bz2 and NA12878\_split\_alignment.tar.bz2) from the data website mentioned above and type following commands (please make sure the [OMTools](https://github.com/aldenleung/OMTools) has been installed):  
>tar -xzvf OMSV.tar.gz  
>tar -jxf NA12878\_alignment.tar.bz2 -C OMSV/data  
>tar -jxf NA12878\_split\_alignment.tar.bz2 -C OMSV/data  
>cd OMSV  
>chmod 777 makefile  
>./makefile  
>./callSV.sh (Run the original version, of which the results are consistent with the published results)

callSV\_old.sh (callSV.sh) provide the examples of calling the OMSV callers (Matlab is required to call CNVs). Then you can find the resulting lists in SV\_result with prefix 12878 (e.g. 12878Indel.osv, 12878Mixed\_indel.osv, and 12878Complex\_total.bed).

=================================================================================
## Parameters and Settings of the tools:  

+ OMSV/OMSV\_mixedIndel (type ./OMSV to check the parameters) to call indels and mixed indels and site variations:  
	- -inputLabel:  Default value: 878. The index/label of genome.  
	- -outputFolder:	Default value: ./. The path of the folder to store the output fils.
	- -SVoutputFile:	Default value: Detected_structual_variants. The prefix of the file name of SVs (.osv).
	- -chrMapFile:	Default value: hg38\_r.cmap. The file name of the reference map(.cmap).
	- -optAlignFile:	The file name of the alignment map file(.oma).
	- -optTempFolder:	Default value: ./. The folder to store the processed alignment maps by chromosomes.
	- -likelihoodRatioCutOff:	Default value: 1e+06. The cutoff of the likelihood ratio for SV all hypothesis (reciprocal of the one in the paper). The default value changes along with the experiment data.
	- -numberOfSupportIndelMolecule:	Default value: 10. The minimum coverage of a segment being called SVs. The default value changes along with the experiment data.
	- -numberOfSupportSignalMolecule:	Default value: 10. The minimum coverage of a segment to call signal variations. The default value changes along with the experiment data.
	- -minIndelSize:	Default value (b): 2000. The minimum length of a segment to call SVs.
	- -minIndelRatio:	Default value: 0.05. The length proportion of a minimum SV could be detected on a segment. E.g. segment = 10000b, then the length of the minimum SV should be larger than 10000\*0.05=500b.
	- -resolutionLimit:	Default value: 1749. The minimum length of a segment to call signal variations.
	- -digestionRate:	Default value: 0.875. The digestion rate of labels (signals) measured in the experiment.
	- -falseCutRate:	Default value: 1e-05. The rate of false cut of a non-label position.
	- -pValueCutOff:	Default value: 1e-09. The cutoff of p-value when call signal variations.
	- -cauchyMean:	Default value: 1.0096. The mean value of cauchy distribution of null hypothesis when calling SVs. Reset a new value only if you have good reason.
	- -cauchyScale:	Default value: 0.0291. The parameter to calculate cauchy distribution. Reset a new value only if you have good reason.
	- -confidenceLimit:	Default value: 9. The lowest alignment confidence for molecules (optical maps) to call SVs or signal variations.
	- -numberOfChromosome:	Default value: 24. The first n chromosomes to detect SVs (chr23 is chrX and chr24 is chrY).

+ complex\_caller.sh has 3 parameters: input alignment (.oma), output folder and output label. Call complex SVs, e.g. translocation, large inversion, and duplications (type ./complex_caller.sh check parameters).

+ CNV\_caller.sh has 4 parameters: 1. alignment file (.oma), 2. reference map (.cmap), 3. output folder and 4. output label. This duplication caller is developed with Matlab, please make sure the matlab is runable. Detect CNV candidates.

+ Medium-size inversions are called by the component of [OMTools](https://github.com/aldenleung/OMTools). Please install the OMTools into the same folder of OMSV, and type "java -jar OMTools.jar SVDetection" to check the parameters.

+ postFilter.sh has 2 parameters: SV list file and cutoff size. This component could remove the false SVs caused by N-gaps, fragile sites, or pseudo-autosomal regions and filter the cases with given SV size.

+ Other utilities and the parameters are posted in our [website](http://yiplab.cse.cuhk.edu.hk/omsv/) (refer to **Section 2. Commands of tools**).

=================================================================================
## File formats:

**Input**: OMSV accepts .oma (OM Alignment) format files as alignment results, which is the output file format of OMBlast, and .cmap (consensus map) file as reference genome maps.  

**Output**: The output file format of OMSV is .osv (OM SV), which generally has 13 columns:  
1. chr: chromosome of SVs.
1. start: start of coordinate of SVs on reference.
1. stop: stop of coordinate of SVs on reference.
1. type: SV types.
1. variant\_id: SV id.
1. obsolete column (preallocated for the next version).
1. sample id.
1. SV size.
1. zygosity of SVs: homozygous or heterozygous.
1. likelihoood: the likelihood of the SV calls.
1. score: the likelihood ratio of the SV calls.
1. coverage: number of molecules covering the SV region.
1. against: number of molecules covering the SV region but support reference type.

=================================================================================
## Citation:
> Le Li, Tsz-Piu Kwok, Alden King-Yung Leung, et.al. **OMSV enables accurate and comprehensive identification of large structural variations from nanochannel-based single-molecule optical maps**. *Manuscript under revision, 2017*.
