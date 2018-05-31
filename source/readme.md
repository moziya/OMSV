The c++ source files (c++11, compiled by g++) are put in this folder.  
OMSV\_v3.cpp and OMSV\_mixedIndel\_v3.cpp are the original code (consistent with published work).  

Updates in May31, 2018:
  Added a low-coverage mode in the pipeline to support the cases that the input alignment is low-coverage (e.g. <10x for the alignment of consensus maps to the reference). To trigger this mode, please set "likelihoodRatioCutOff" as 0, "numberOfSupportIndelMolecule" as 1. 
