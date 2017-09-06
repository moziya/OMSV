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
	double fpr;////
        double fnr;////
        double alignRate;////
        vector<LL> position;
	void calc()////
        {////
                int tot_fp = 0, tot_fn = 0;
                for (int i = 0; hitEnum[i]!='\0'; i++){
                        int curr = 0;
                        int fp, fn;
                        int j;
                        for (j = i; hitEnum[j] >= '0' && hitEnum[j] <= '9'; j++)
                                curr = curr * 10 + (hitEnum[j] - '0');
                        if (hitEnum[j] == 'I'){
                                fp++;
                                tot_fp += curr;
                        }
                        else if (hitEnum[j] == 'D'){
                                fn++;
                                tot_fn += curr;
                        }
                        i = j;
                }
                int mol_len = 0;
                int ali_len = 0;
                int cnt = 0;
                for (auto ele:position){
                        cnt++;
                        mol_len += ele;
                        if (cnt >= min(optStart,optEnd)&&cnt <= max(optStart,optEnd)) ali_len+=ele;
                }
                fpr = tot_fp*1.0/ali_len;
                fnr = tot_fn*1.0/(abs(refEndIndex-refStartIndex)+1);
                alignRate = ali_len*1.0/mol_len;
        }////
	bool operator<(const opticalMapType &w) const{
		LL e = strcmp(mapId, w.mapId);
		return ((e<0?true:false) || (e==0 && chrId < w.chrId) || (e == 0 && chrId == w.chrId && optStart < w.optStart));
	}
};


void readSourceFile(char* inputAlignmentFileName, char* outputBefore, char* outputAfter, LL minM, LL minL){
	vector<opticalMapType> opticalMap1;
	vector<opticalMapType> opticalMap2;
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
        opticalMap2.clear();
	opticalMap1.resize(5000000);
	opticalMap2.resize(5000000);
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 5000, inputAlignmentFile);
        LL numberOfOpticalMap = 0;
        while (fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].mapId) == 1){
                if (numberOfOpticalMap % 100000 == 0) printf("%lld %s\n", numberOfOpticalMap, opticalMap1[numberOfOpticalMap].mapId);
                LL tempNumberOfSites, tempLL;
                fscanf(inputAlignmentFile, "%lld", &tempNumberOfSites);
                LL tempIndex = 0;
                for (LL i=0; i<tempNumberOfSites; i++){
                        if (i != tempNumberOfSites - 1)
                                fscanf(inputAlignmentFile, "%lld;", &tempLL);
                        else fscanf(inputAlignmentFile, "%lld", &tempLL);
                        tempIndex += tempLL;
                        opticalMap1[numberOfOpticalMap].position.push_back(tempIndex);
                }
                fscanf(inputAlignmentFile, "%s", ttts);
                if (ttts[0] == 'c'){
                        if (ttts[3] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (ttts[3] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (ttts[3] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
                        else
                                sscanf(ttts, "chr%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                else{
                        if (ttts[0] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (ttts[0] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (ttts[0] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
			else
                        	sscanf(ttts, "%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                fscanf(inputAlignmentFile, "%s", ttts);
                if (ttts[0] == 'r' || ttts[0] == '-') opticalMap1[numberOfOpticalMap].orientation = false; else opticalMap1[numberOfOpticalMap].orientation = true;
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].score);
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].confidence);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStartIndex, &opticalMap1[numberOfOpticalMap].refEndIndex);
		opticalMap1[numberOfOpticalMap].refStartIndex--;
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].optStart, &opticalMap1[numberOfOpticalMap].optEnd);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStart, &opticalMap1[numberOfOpticalMap].refEnd);
                fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].hitEnum);
		opticalMap1[numberOfOpticalMap].calc();////
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
	FILE* outFB = fopen(outputBefore,"w");
	fprintf(outFB,"#Split maps before alignment\n");
	FILE* outFA = fopen(outputAfter,"w");
	fprintf(outFA,"#Split maps after alignment\n");
	for (LL i = 0; i<opticalMap1.size();i++){
		if (i<opticalMap1.size()-1&&strcmp(opticalMap1[i+1].mapId,opticalMap1[i].mapId)==0)continue;
		if(i==0||strcmp(opticalMap1[i].mapId,opticalMap1[i-1].mapId)!=0){
			LL bef_idx = min(opticalMap1[i].optStart,opticalMap1[i].optEnd);
			if (bef_idx>=minM&&opticalMap1[i].position[bef_idx]>minL){
				LL st_idd = 1;
				for (LL j = 0; j < bef_idx;j++)
					fprintf(outFB,"%s\t%lld\t%lld\t%lld\t1\t%lld\t0.0\t11.0\t11.0\n",opticalMap1[i].mapId,opticalMap1[i].position[bef_idx],bef_idx,st_idd++,max(opticalMap1[i].position[j],(LL)1));
				fprintf(outFB,"%s\t%lld\t%lld\t%lld\t0\t%lld\t0.0\t11.0\t11.0\n",opticalMap1[i].mapId,opticalMap1[i].position[bef_idx],bef_idx,st_idd++,opticalMap1[i].position[bef_idx]);
			}
			LL aft_idx = max(opticalMap1[i].optStart,opticalMap1[i].optEnd);
			LL end_idx = opticalMap1[i].position.size()-1;
			if (end_idx - aft_idx >=minM && opticalMap1[i].position[end_idx]-opticalMap1[i].position[aft_idx]>minL){
				LL st_idd = 1;
				for (LL j = aft_idx+1; j < end_idx;j++)
					fprintf(outFA,"%s\t%lld\t%lld\t%lld\t1\t%lld\t0.0\t11.0\t11.0\n",opticalMap1[i].mapId,opticalMap1[i].position[end_idx]-opticalMap1[i].position[aft_idx],end_idx-aft_idx-1,st_idd++,max((LL)1,opticalMap1[i].position[j]-opticalMap1[i].position[aft_idx]));
				fprintf(outFA,"%s\t%lld\t%lld\t%lld\t0\t%lld\t0.0\t11.0\t11.0\n",opticalMap1[i].mapId,opticalMap1[i].position[end_idx]-opticalMap1[i].position[aft_idx],end_idx-aft_idx-1,st_idd++,opticalMap1[i].position[end_idx]-opticalMap1[i].position[aft_idx]);
				
			}
		}
	}
	
	fclose(outFA);
	fclose(outFB);
}


int main(int argv, char *argc[]){
	printf("\n");
	//setDefault();
	LL minM = 5;
	LL minL = 1000*50;
	int suppCutoff = 5;
	double varPrior = 1e-6;
	char inputMR[1000];
	char outputFile[1000];
	char outputFile2[1000];
	
	if (argv == 1){
		printf("\nThe valid parameters are described as follows:\n");
		printf("\t-outputBeforeFile: \n\t\t The file name of segments before alignment.\n");
		printf("\t-outputAfterFile: \n\t\t The file name of segments after alignment.\n");
		printf("\t-optAlignFile: \n\t\t  The file name of the alignment map file(.oma).\n");
		printf("\t-minSiteNumber: \n\t\t  The minimum number of sites on the split map(.oma). Default: %lld\n",minM);
		printf("\t-minLength: \n\t\t  The minimum length of split map(.oma). Default: %lld\n",minL);
		return -1;
	}
	for (int i = 1; i < argv; i=i+2){
		string temp(argc[i]);
		//printf("What is the parameter?? %s \n",temp.c_str());
		if (temp.compare("-optAlignFile")==0)
			strcpy(inputMR, argc[i+1]);
		else if (temp.compare("-outputBeforeFile")==0)
			strcpy(outputFile,argc[i+1]);
		else if (temp.compare("-outputAfterFile")==0)
			strcpy(outputFile2,argc[i+1]);
		else if (temp.compare("-minSiteNumber")==0)
			minM = atol(argc[i+1]);
		else if (temp.compare("-minLength")==0)
			minL = atol(argc[i+1]);
		else
			printf("No such parameter or wrong : %s\n",argc[i]);
	}
	readSourceFile(inputMR,outputFile,outputFile2,minM,minL);
	return 0;
}
