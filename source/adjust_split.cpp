#include "NoError.h"

FILE *inputAlignmentFile;  
double C[9000][9000]; 
double minIndelRatio; 
double resolutionLimit; 
double digestionRate; 
double falseCutRate; 
double pValueCutOff; 
double cauchyMean; 
double cauchyScale; 
LL numberOfSupportEachSide; 
double confidenceLimit; 
LL numberOfSupportIndelMolecule; 
double likelihoodRatioCutOff; //Start add 

struct	opticalMapType{
        char mapId[100];
	LL numberOfSite;
	char segS[10000];
        LL refStartIndex;
        LL refEndIndex;
        LL refStart;
        LL refEnd;
        LL optStart;
        LL optEnd;
        LL chrId;
        bool orientation;
	double score;
        double confidence;
        char hitEnum[2000];
	bool operator<(const opticalMapType &w) const{
		LL e = strcmp(mapId, w.mapId);
		return ((e<0?true:false) || (e==0 && chrId < w.chrId) || (e == 0 && chrId == w.chrId && optStart < w.optStart));
	}
};


void readSourceFile(char* inputAlignmentFileName, char* inputAlignBefore, char* inputAlignAfter, char* outputBefore=(char*)"Adjusted_before.oma", char* outputAfter=(char*)"Adjusted_after.oma"){
	vector<opticalMapType> opticalMap1;
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
        int numberOfDL = 0;
        char ttt;
        char ttts[100000];
        while(fgetc(inputAlignmentFile) == '#'){
                numberOfDL++;
                fgets(ttts,100000,inputAlignmentFile);
        }
        fclose(inputAlignmentFile);
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
        opticalMap1.clear();
	opticalMap1.resize(5000000);
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 5000, inputAlignmentFile);
        LL numberOfOpticalMap = 0;
        while (fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].mapId) == 1){
                if (numberOfOpticalMap % 100000 == 0) printf("%lld %s\n", numberOfOpticalMap, opticalMap1[numberOfOpticalMap].mapId);
                LL tempNumberOfSites, tempLL;
                fscanf(inputAlignmentFile, "%lld", &opticalMap1[numberOfOpticalMap].numberOfSite);
		tempNumberOfSites = opticalMap1[numberOfOpticalMap].numberOfSite;

                fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].segS);
		fscanf(inputAlignmentFile, "%s\t%s\t%s\t%s\t%s\t%s\t",ttts,ttts,ttts,ttts,ttts,ttts);
		fscanf(inputAlignmentFile, "%lld\t%lld\t",&opticalMap1[numberOfOpticalMap].optStart,&opticalMap1[numberOfOpticalMap].optEnd);
		fgets(ttts,100000,inputAlignmentFile);
                numberOfOpticalMap++;
        }
	opticalMap1.resize(numberOfOpticalMap);
        printf("Number of optical map: %lld\n", (LL)opticalMap1.size());
        for (LL i=0; i<(LL)opticalMap1.size(); i++){
		LL minO = min(opticalMap1[i].optStart,opticalMap1[i].optEnd);
		LL maxO = max(opticalMap1[i].optStart,opticalMap1[i].optEnd);
		opticalMap1[i].optStart = minO;
		opticalMap1[i].optEnd = maxO;
                opticalMap1[i].optEnd--;
        }
        for (LL i=1; i<(LL)opticalMap1.size(); i++){
		if (strcmp(opticalMap1[i].mapId,opticalMap1[i-1].mapId)==0){
			LL minO = min(opticalMap1[i].optStart,opticalMap1[i-1].optStart);
			LL maxO = max(opticalMap1[i].optEnd,opticalMap1[i-1].optEnd);
			opticalMap1[i].optStart = minO;
			opticalMap1[i].optEnd = maxO;
			opticalMap1[i-1].optStart = minO;
			opticalMap1[i-1].optEnd = maxO;
		}
	}
	

        if ((inputAlignmentFile = fopen(inputAlignBefore, "r")) == NULL) puts("ERROR IN READ SOURCE");
        numberOfDL = 0;
        while(fgetc(inputAlignmentFile) == '#'){
                numberOfDL++;
                fgets(ttts,100000,inputAlignmentFile);
        }
        fclose(inputAlignmentFile);
        if ((inputAlignmentFile = fopen(inputAlignBefore, "r")) == NULL) puts("ERROR IN READ SOURCE");
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 5000, inputAlignmentFile);
	char mapId[1000], segS[10000];
	LL numberOfSite;
	FILE* outF;

	outF = fopen(outputBefore,"w");

	fprintf(outF,"#Adjusted alignment for the segments before alignment!\n");
	LL sti = 0;
	while(fscanf(inputAlignmentFile,"%s\t%lld\t%s\t",mapId,&numberOfSite,segS)==3){
		bool tes = false;
		for (;sti<opticalMap1.size();sti++){
			opticalMapType mol = opticalMap1[sti];
			if (strcmp(mol.mapId,mapId)==0){
				tes = true;
				numberOfSite = mol.numberOfSite;
				strcpy(segS,mol.segS);
				break;
			}
		}
		if (!tes)perror("Wrong for the before seg align!");
		memset(ttts,0,sizeof(ttts));
		fgets(ttts,100000,inputAlignmentFile);
		fprintf(outF,"%s\t%lld\t%s\t%s",mapId,numberOfSite,segS,ttts);
	}
	fclose(outF);

        if ((inputAlignmentFile = fopen(inputAlignAfter, "r")) == NULL) puts("ERROR IN READ SOURCE");
        numberOfDL = 0;
        while(fgetc(inputAlignmentFile) == '#'){
                numberOfDL++;
                fgets(ttts,100000,inputAlignmentFile);
        }
        fclose(inputAlignmentFile);
        if ((inputAlignmentFile = fopen(inputAlignAfter, "r")) == NULL) puts("ERROR IN READ SOURCE");
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 5000, inputAlignmentFile);
	char temp1[100],temp2[100],temp3[100],temp4[100],temp5[100],temp6[100];
	LL optS, optE;
	outF = fopen(outputAfter,"w");
	fprintf(outF,"#Adjusted alignment for the segments after alignment!\n");
	sti=0;
	while(fscanf(inputAlignmentFile,"%s\t%lld\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lld\t%lld\t",mapId,&numberOfSite,segS,temp1,temp2,temp3,temp4,temp5,temp6,&optS,&optE)==11){
		bool tes = false;
		for (;sti<opticalMap1.size();sti++){
			opticalMapType mol = opticalMap1[sti];
			if (strcmp(mol.mapId,mapId)==0){
				tes = true;
				numberOfSite = mol.numberOfSite;
				strcpy(segS,mol.segS);
				LL bas = max(mol.optStart,mol.optEnd)+1;
				optS+=bas;
				optE+=bas;
				break;
			}
		}
		if (!tes)perror("Wrong for the before seg align!");
		memset(ttts,0,sizeof(ttts));
		fgets(ttts,100000,inputAlignmentFile);
		fprintf(outF,"%s\t%lld\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lld\t%lld\t%s",mapId,numberOfSite,segS,temp1,temp2,temp3,temp4,temp5,temp6,optS,optE,ttts);
	}
	fclose(outF);	
}


int main(int argv, char *argc[]){
	printf("\n");
	//setDefault();
	LL miniM = 3;
	LL shift = 1000*50;
	int suppCutoff = 5;
	double varPrior = 1e-6;
	char inputAlign[1000];
	char inputBefore[1000];
	char inputAfter[1000];
	char outputBefore[1000];
	char outputAfter[1000];
	
	if (argv == 1){
		printf("\nThe valid parameters are described as follows:\n");
		printf("\t-optAlignFile: \n\t\t  The file name of the alignment map file(.oma).\n");
		printf("\t-optAlignBeforeFile: \n\t\t  The file name of the segments before alignment map file(.oma).\n");
		printf("\t-optAlignAfterFile: \n\t\t The file name of the segments after alignment map file(.oma).\n");
		printf("\t-outputAlignBefore: \n\t\t The output file for before alignment.\n");
		printf("\t-outputAlignAfter: \n\t\t The output file for after alignment.\n");
		return -1;
	}
	for (int i = 1; i < argv; i=i+2){
		string temp(argc[i]);
		//printf("What is the parameter?? %s \n",temp.c_str());
		if (temp.compare("-optAlignFile")==0)
			strcpy(inputAlign, argc[i+1]);
		else if (temp.compare("-optAlignBeforeFile")==0)
			strcpy(inputBefore, argc[i+1]);
		else if (temp.compare("-optAlignAfterFile")==0)
			strcpy(inputAfter, argc[i+1]);
		else if (temp.compare("-outputAlignBefore")==0)
			strcpy(outputBefore,argc[i+1]);
		else if (temp.compare("-outputAlignAfter")==0)
			strcpy(outputAfter,argc[i+1]);
		else
			printf("No such parameter or wrong : %s\n",argc[i]);
	}
	readSourceFile(inputAlign,inputBefore,inputAfter,outputBefore,outputAfter);
	return 0;
}
