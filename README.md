# QUICK START
=================================================================================
The scrips of OMSV have been tested in Debian GNU/Linux 9.0 (stretch) and CentOS Linux release 7.3.1611 (Core) platform. And the matlab code was developed in Matlab R2011b (7.13.0.564) 64-bit (glnxa64).

Quick Start

To quickly try our OMSV Caller, please first download the package(OMSV.tar.gz) and alignment files of OM data (e.g. NA12878\_alignment.tar.bz2 and NA12878_split_alignment.tar.bz2) from our website http://yiplab.cse.cuhk.edu.hk/omsv/ and type following commands:
>tar -xzvf OMSV.tar.gz 
>tar -jxf NA12878_alignment.tar.bz2 -C OMSV/data
>tar -jxf NA12878_split_alignment.tar.bz2 -C OMSV/data
>cd OMSV
>chmod 777 makefile
>./makefile
>./callSV.sh

callSV.sh provide the examples of calling the OMSV callers (Matlab is required to call CNVs). Then you can find the resulting lists in SV\_result with prefix 12878 (e.g. 12878Indel.osv, 12878Mixed\_indel.osv, and 12878Complex\_total.bed).

=================================================================================
Parameters and Settings of the tools:
   OMSV (type ./OMSV to check the parameters) to call indels and mixed indels and site variations:
        -inputLabel:
                 Default value: 878. The index/label of genome.
        -outputFolder:
                 Default value: ./. The path of the folder to store the output fils.
        -SVoutputFile:
                 Default value: Detected_structual_variants. The prefix of the file name of SVs (.osv).
        -chrMapFile:
                 Default value: hg38_r.cmap. The file name of the reference map(.cmap).
        -optAlignFile:
                 Default value: C6661_700bp_hg38_combRefOMB2.oma. The file name of the alignment map file(.oma).
        -optTempFolder:
                 Default value: ./. The folder to store the processed alignment maps by chromosomes.
        -likelihoodRatioCutOff:
                 Default value: 1e+06. The cutoff of the likelihood ratio for SV all hypothesis (reciprocal of the one in the paper). The default value changes along with the experiment data.
        -numberOfSupportIndelMolecule:
                 Default value: 10. The minimum coverage of a segment being called SVs. The default value changes along with the experiment data.
        -numberOfSupportSignalMolecule:
                 Default value: 10. The minimum coverage of a segment to call signal variations. The default value changes along with the experiment data.
        -minIndelSize:
                 Default value (b): 2000. The minimum length of a segment to call SVs.
        -minIndelRatio:
                 Default value: 0.05. The length proportion of a minimum SV could be detected on a segment. E.g. segment = 10000b, then the length of the minimum SV should be larger than 10000*0.05=500b.
        -resolutionLimit:
                 Default value: 1749. The minimum length of a segment to call signal variations.
        -digestionRate:
                 Default value: 0.875. The digestion rate of labels (signals) measured in the experiment.
        -falseCutRate:
                 Default value: 1e-05. The rate of false cut of a non-label position.
        -pValueCutOff:
                 Default value: 1e-09. The cutoff of p-value when call signal variations.
        -cauchyMean:
                 Default value: 1.0096. The mean value of cauchy distribution of null hypothesis when calling SVs. Reset a new value only if you have good reason.
        -cauchyScale:
                 Default value: 0.0291. The parameter to calculate cauchy distribution. Reset a new value only if you have good reason.
        -confidenceLimit:
                 Default value: 9. The lowest alignment confidence for molecules (optical maps) to call SVs or signal variations.
        -numberOfChromosome:
                 Default value: 24. The first n chromosomes to detect SVs.

   complex\_caller.sh has 3 parameters: input alignment (.oma), output folder and output label. Call complex SVs, e.g. translocation, large inversion, and duplications (type ./complex_caller.sh check parameters).

   CNV\_caller.sh has 4 parameters: 1. alignment file (.oma), 2. reference map (.cmap), 3. output folder and 4. output label. This duplication caller is developed with Matlab, please make sure the matlab is runable. Detect CNV candidates.

   Medium-size inversions are called by the component of OMTools.jar. Please type "java -jar OMTools.jar SVDetection" to check the parameters.

   postFilter.sh has 2 parameters: SV list file and cutoff size. This component could remove the false SVs caused by N-gaps, fragile sites, or pseudo-autosomal regions and filter the cases with given SV size.



OMSV accepts .oma (OM Alignment) format files as alignment results, which is the output file format of OMBlast, and .cmap (consensus map) file as reference genome maps.
The resulting file format is .osv (OM SV), which generally has 13 columns:
	1. chr: chromosome of SVs.
	2. start: start of coordinate of SVs on reference.
	3. stop: stop of coordinate of SVs on reference.
	4. type: SV types.
	5. variant\_id: SV id.
	6. obsolete column (preallocated for the next version).
	7. sample id.
	8. SV size.
	9. zygosity of SVs: homozygous or heterozygous.
	10. likelihoood: the likelihood of the SV calls.
	11. score: the likelihood ratio of the SV calls.
	12. coverage: number of molecules covering the SV region.
	13. against: number of molecules covering the SV region but support reference type.
