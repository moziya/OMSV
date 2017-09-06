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

struct bedReg{
	LL chr;
	LL start;
	LL stop;
};
vector<bedReg> highDens;

void readHighDensity(char* inputAlignmentFileName){
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
        highDens.clear();
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 5000, inputAlignmentFile);
	bedReg tpBed;
	while(fscanf(inputAlignmentFile, "%lld\t%lld\t%lld\n", &tpBed.chr,&tpBed.start,&tpBed.stop) == 3){
		highDens.push_back(tpBed);
	}
	fclose(inputAlignmentFile);
}

vector<opticalMapType> readSourceFile(char* inputAlignmentFileName){
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
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].optStart, &opticalMap1[numberOfOpticalMap].optEnd);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStart, &opticalMap1[numberOfOpticalMap].refEnd);
                fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].hitEnum);
		opticalMap1[numberOfOpticalMap].calc();////
                numberOfOpticalMap++;
        }
	opticalMap1.resize(numberOfOpticalMap);
        printf("Number of optical map: %lld\n", (LL)opticalMap1.size());
        for (LL i=0; i<(LL)opticalMap1.size(); i++){
                if (opticalMap1[i].orientation){
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart - 1];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd];
                }else {
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd - 1];
                }
        }
	return opticalMap1;
}

struct breakpoint{//actually is break point pair
	LL chr1; // if chr1 != chr2, let chr1 < chr2
	LL point1; // if chr1 == chr2, let point1 < point2
	LL pointIdx1;
	LL chr2;
	LL point2;
	LL pointIdx2;
	char type[100];
	vector<LL> supportMapId;
	vector<double> supportScores;
	vector<LL> nonSupportMapId;
	vector<LL> nonSupportMapId1;
	vector<LL> nonSupportMapId2;
	vector<double> nonSuppScores;
	vector<double> nonSuppScores1;	
	vector<double> nonSuppScores2;	
	bool operator<(const breakpoint &x) const{
		return (chr1 < x.chr1 || (chr1 == x.chr1 && point1 < x.point1) || (chr1 == x.chr1 && point1 == x.point1 && chr2 < x.chr2) || (chr1 == x.chr1 && point1 == x.point1 && chr2 == x.chr2 && point2 < x.point2));
	}
	double score = 0;
};

bool ss(breakpoint y, breakpoint x){
	return (y.chr1 < x.chr1 || (y.chr1 == x.chr1 && y.chr2 < x.chr2) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && (y.point1+y.point2) < (x.point2+x.point1)));
}
bool st(breakpoint y, breakpoint x){
	return (y.chr1 < x.chr1 || (y.chr1 == x.chr1 && y.chr2 < x.chr2) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) < 0) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) == 0  && (y.point1+y.point2) < (x.point2+x.point1)));
}

vector<breakpoint> candidateSet; 
vector<breakpoint> completeSet; 

double calc_score(char* hitEnum, int miniM){
	int tempCount = 0;
	double N=0, n=0;
	for (LL j = 0; hitEnum[j]; j++)
	{
		if (hitEnum[j] >= '0' && hitEnum[j] <= '9')
	        	tempCount = tempCount * 10 + (hitEnum[j]-'0');
	    	else{
        		if (hitEnum[j] == 'M')
				n+=tempCount;
			N+=tempCount;
			tempCount = 0;
		}
	}
	if (n<miniM) return 0;
	return log(n)*n/N;
}

double calc_score_div(char* hitEnum, int rel_pos, int miniM){
	int tempCount = 0;
	int N1=0, n1=0, N2=0, n2=0, N=0;
	for (LL j = 0; hitEnum[j]; j++){
		if (hitEnum[j] >= '0' && hitEnum[j] <= '9')
        		tempCount = tempCount * 10 + (hitEnum[j]-'0');
		else{
        		if (hitEnum[j] == 'M'){
				N+=tempCount;
				if (N<rel_pos)
				{
					n1+=tempCount;
				}
				else if (N>=rel_pos&&N-tempCount<rel_pos)
				{
					n1+=(rel_pos+tempCount-N);
					n2+=(N-rel_pos);
				}
				else
				{
					n2+=tempCount;
				}
			}
			else if (hitEnum[j] == 'D')
				N+=tempCount;
			if (N<rel_pos)
				N1+=tempCount;
			else if (N>=rel_pos&&N-tempCount<rel_pos){
				N1+=(rel_pos+tempCount-N);
				N2+=(N-rel_pos);
			}
			else
				N2+=tempCount;
			tempCount = 0;
		}
	}
	if (n1<miniM || n2<miniM) return 0;
	return min(log(n1)*double(n1)/N1,log(n2)*double(n2)/N2);
}

void swapNums(LL &x, LL &y)
{
	LL temp = x;
	x = y;
	y = temp;
}

void scanSplitedMap(vector<opticalMapType> opticalMap1, int miniM){
	LL tempCC = (LL)opticalMap1.size();
	sort(opticalMap1.begin(), opticalMap1.end());
	for (LL i=0; i<tempCC-1; i++){
		if (strcmp(opticalMap1[i].mapId, opticalMap1[i+1].mapId) != 0)continue;
		if (opticalMap1[i].chrId == opticalMap1[i+1].chrId && opticalMap1[i].orientation == opticalMap1[i+1].orientation && (opticalMap1[i].orientation == (opticalMap1[i].refStart < opticalMap1[i+1].refStart)) && max(opticalMap1[i+1].refStart - opticalMap1[i].refEnd,opticalMap1[i].refStart - opticalMap1[i+1].refEnd) < 5000000 && max(opticalMap1[i+1].refStart - opticalMap1[i].refEnd,opticalMap1[i].refStart - opticalMap1[i+1].refEnd) > 0){
			//it requires the opticalMaps are all + oriented
			//1. chr same, 2. orientation same, 3. gap is small, 4. gap orientation is same
		}
		else
		{
			breakpoint tempBP;
			tempBP.supportMapId.clear();
			tempBP.supportScores.clear();
			tempBP.nonSupportMapId1.clear();
			tempBP.nonSupportMapId2.clear();
			tempBP.nonSupportMapId.clear();
			tempBP.nonSuppScores1.clear();
			tempBP.nonSuppScores2.clear();
			tempBP.nonSuppScores.clear();
			tempBP.chr1 = opticalMap1[i].chrId;
			tempBP.chr2 = opticalMap1[i+1].chrId;

			if (tempBP.chr1<tempBP.chr2&&max(opticalMap1[i].optStart,opticalMap1[i].optEnd)>max(opticalMap1[i+1].optStart,opticalMap1[i+1].optEnd)){
				if (!opticalMap1[i].orientation){
					tempBP.point1 = opticalMap1[i].refEnd;
					tempBP.pointIdx1 = opticalMap1[i].refEndIndex;
				} else {
					tempBP.point1 = opticalMap1[i].refStart;
					tempBP.pointIdx1 = opticalMap1[i].refStartIndex;
				}
				
				if (!opticalMap1[i+1].orientation){
					tempBP.point2 = opticalMap1[i+1].refStart;
					tempBP.pointIdx2 = opticalMap1[i+1].refStartIndex;
				} else {
					tempBP.point2 = opticalMap1[i+1].refEnd;
					tempBP.pointIdx2 = opticalMap1[i+1].refEndIndex;
				}
			}else{
				if (opticalMap1[i].orientation){
                                        tempBP.point1 = opticalMap1[i].refEnd;
                                        tempBP.pointIdx1 = opticalMap1[i].refEndIndex;
                                } else {
                                        tempBP.point1 = opticalMap1[i].refStart;
                                        tempBP.pointIdx1 = opticalMap1[i].refStartIndex;
                                }

                                if (opticalMap1[i+1].orientation){
                                        tempBP.point2 = opticalMap1[i+1].refStart;
                                        tempBP.pointIdx2 = opticalMap1[i+1].refStartIndex;
                                } else {
                                        tempBP.point2 = opticalMap1[i+1].refEnd;
                                        tempBP.pointIdx2 = opticalMap1[i+1].refEndIndex;
                                }
			}
			if (tempBP.chr1 == tempBP.chr2 && tempBP.point2 < tempBP.point1){
				swapNums(tempBP.point1,tempBP.point2);
				swapNums(tempBP.pointIdx1,tempBP.pointIdx2);
			}


			if (min(calc_score(opticalMap1[i].hitEnum,miniM),calc_score(opticalMap1[i+1].hitEnum,miniM)) > 0)
			{
				tempBP.supportMapId.push_back(atol(opticalMap1[i].mapId));
				tempBP.supportScores.push_back(min(calc_score(opticalMap1[i].hitEnum,miniM),calc_score(opticalMap1[i+1].hitEnum,miniM)));
			}
			double ovlpDist = (double)max(opticalMap1[i+1].refStart - opticalMap1[i].refEnd,opticalMap1[i].refStart - opticalMap1[i+1].refEnd);
			bool cond1 = opticalMap1[i].chrId != opticalMap1[i+1].chrId;
			bool cond2 = max(opticalMap1[i+1].refStart - opticalMap1[i].refEnd,opticalMap1[i].refStart - opticalMap1[i+1].refEnd) >= 5000000;
//			bool cond3 = ovlpDist < 0;
			bool cond3 = -ovlpDist >=  0.8 * min(opticalMap1[i].refEnd - opticalMap1[i].refStart, opticalMap1[i+1].refEnd - opticalMap1[i+1].refStart);
			bool cond31 = -ovlpDist <  0.3 * min(opticalMap1[i].refEnd - opticalMap1[i].refStart, opticalMap1[i+1].refEnd - opticalMap1[i+1].refStart);
			bool cond4 = opticalMap1[i].orientation != opticalMap1[i+1].orientation;
			bool cond5 = opticalMap1[i].orientation != (opticalMap1[i].refStart < opticalMap1[i+1].refStart);
			

			if (cond1)
				strcpy(tempBP.type,"Inter-Translocation");
			else if (cond2){
				strcpy(tempBP.type,"Intra-Translocation");

//				if (opticalMap1[i].orientation != opticalMap1[i+1].orientation)
  //                                      strcat(tempBP.type,"_Inversion");
    //                            else if ((opticalMap1[i].orientation != (opticalMap1[i].refStart < opticalMap1[i+1].refStart)))
      //                                  strcat(tempBP.type,"_Disorder");
	//			else
	//				strcat(tempBP.type,"_")

			}
			else if (cond4)
				strcpy(tempBP.type,"Inversion");
			else if (cond3){
				strcpy(tempBP.type,"Duplication");
			}
//			else if (cond3&&cond4){
//				strcpy(tempBP.type,"Tandem-Inversion");
//			}
			else{
				strcpy(tempBP.type,"Unknown-Breakpoints");				
			}
		/*		printf("\tOptStart\tOptEnd\tChr\tRefStart\tRefEnd\tOrient\t%\t%s\t%s\n",tempBP.type,opticalMap1[i].mapId,opticalMap1[i+1].mapId);
				char ori_i, ori_j;
				if (opticalMap1[i].orientation)ori_i = '+';
				else ori_i = '-';
				if (opticalMap1[i+1].orientation)ori_j = '+';
                                else ori_j = '-';
				printf("\t%lld\t%lld\t%lld\t%lld\t%lld\t%c\n",opticalMap1[i].optStart,opticalMap1[i].optEnd,opticalMap1[i].chrId,opticalMap1[i].refStartIndex,opticalMap1[i].refEndIndex,ori_i);
				printf("\t%lld\t%lld\t%lld\t%lld\t%lld\t%c\n",opticalMap1[i+1].optStart,opticalMap1[i+1].optEnd,opticalMap1[i+1].chrId,opticalMap1[i+1].refStartIndex,opticalMap1[i+1].refEndIndex,ori_j);
*/
			candidateSet.push_back(tempBP);
		}
	}
	printf("So far, there are %lld candidate complicated SV\n",(LL)candidateSet.size());
}
void printCandidate(){
	for(auto tempBP:candidateSet)
		printf("MapId:%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld\t%lld\n",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size());
}

void printCandidateToFile(FILE* outputFile){
	fprintf(outputFile,"#Complex large SVs called with split alignment.\n");
        for(auto tempBP:completeSet)
	{
		bool nUs = false;
		for (auto bed:highDens){
			if ((tempBP.chr1 == bed.chr && tempBP.point1>=bed.start && tempBP.point1<=bed.stop) || (tempBP.chr2 == bed.chr && tempBP.point2>=bed.start && tempBP.point2<=bed.stop))
				nUs = true;
		}
		if (tempBP.supportScores.size()>1&&!nUs)
			fprintf(outputFile,"MapId:%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld\t%lld\n",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size());
	}
	fclose(outputFile);
}


void mergeCandidate(LL shift){
	completeSet.clear();
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	char type[100] = "None";
	LL cnt = 0;
	shift = 0;
	for (LL i = 0; i < candidateSet.size(); i++)
	{
		if (candidateSet[i].supportScores.size()>0 && pID1 == candidateSet[i].chr1 && pID2 == candidateSet[i].chr2 && strcmp(candidateSet[i].type,type)==0 && abs(pP1-candidateSet[i].point1) <= shift && abs(pP2-candidateSet[i].point2) <= shift)
		{
			bool flag = true;
			for (LL k = 0; k < completeSet[cnt-1].supportMapId.size(); k++){
				if (completeSet[cnt-1].supportMapId[k]==candidateSet[i].supportMapId[0]){
					if (completeSet[cnt-1].supportScores[k] < candidateSet[i].supportScores[0]){
						completeSet[cnt-1].supportScores[k] = candidateSet[i].supportScores[0];
					}
					flag = false;
					break;
				}
			}
			if (flag){
				completeSet[cnt-1].supportScores.push_back(candidateSet[i].supportScores[0]);
				completeSet[cnt-1].supportMapId.push_back(candidateSet[i].supportMapId[0]);
			}
		}
		else
		{
			completeSet.push_back(candidateSet[i]);
			cnt++;
			pID1 = candidateSet[i].chr1;
			pID2 = candidateSet[i].chr2;
			pP1 = candidateSet[i].point1;
			pP2 = candidateSet[i].point2;
			strcpy(type,candidateSet[i].type);
		}
	}
	printf("After merge, there are %lld candidate SV\n",(LL)completeSet.size());
}

void mergeCandidateAgain(LL shift){//to merge the candidates with different but close regions
	candidateSet = completeSet;
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	completeSet.clear();
        LL pID1=-1,pID2=-1, pP1=-1, pP2=-1, cov = -1;
	char type[100] = "None";
        LL cnt = 0;
	vector<int> flags(candidateSet.size(),0);
        for (LL ii = 0; ii < candidateSet.size(); ii++){
		if (flags[ii]!=0)continue;
        	completeSet.push_back(candidateSet[ii]);
                cnt++;
		flags[ii]=1;
                pID1 = candidateSet[ii].chr1;
	        pID2 = candidateSet[ii].chr2;
        	pP1 = candidateSet[ii].point1;
                pP2 = candidateSet[ii].point2;
		cov = (LL)candidateSet[ii].supportScores.size();
		strcpy(type,candidateSet[ii].type);
		
		for (LL i = ii+1; i < candidateSet.size(); i++){
			if (flags[i]!=0)continue;
	                if (candidateSet[i].supportScores.size()>0 && pID1 == candidateSet[i].chr1 && pID2 == candidateSet[i].chr2 && strcmp(candidateSet[i].type,type)==0 && abs(pP1-candidateSet[i].point1) <= shift && abs(pP2-candidateSet[i].point2) <= shift)
        	        {
				flags[i] = 1;
				for (LL t = 0; t < candidateSet[i].supportMapId.size(); t++){
                        		bool flag = true;
	                        	for (LL k = 0; k < completeSet[cnt-1].supportMapId.size(); k++){
        	                        	if (completeSet[cnt-1].supportMapId[k]==candidateSet[i].supportMapId[t]){
                	                        	if (completeSet[cnt-1].supportScores[k] < candidateSet[i].supportScores[t]){
                        	                        	completeSet[cnt-1].supportScores[k] = candidateSet[i].supportScores[t];
	                                	        }
        	                                	flag = false;
	        	                                break;
        	        	                }
                	        	}
	                        	if (flag){
        	                        	completeSet[cnt-1].supportScores.push_back(candidateSet[i].supportScores[t]);
	        	                        completeSet[cnt-1].supportMapId.push_back(candidateSet[i].supportMapId[t]);
                		        }
				}
				if ((LL)candidateSet[i].supportScores.size() >= cov){//always keep the regions originally having the most number of supports
					pP1 = candidateSet[i].point1;
					pP2 = candidateSet[i].point2;
					cov = (LL)candidateSet[i].supportScores.size();
				}
        	        }
                	else if (pID1 != candidateSet[i].chr1 || pID2 != candidateSet[i].chr2 || strcmp(candidateSet[i].type,type)!=0)
				break;
		}
        }
        printf("After merge, there are %lld candidate SV\n",(LL)completeSet.size());
}

void mergeCandidateAgainO(LL shift){//to merge the candidates with different but close regions
	candidateSet = completeSet;
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	completeSet.clear();
        LL pID1=-1,pID2=-1, pP1=-1, pP2=-1, cov = -1;
	char type[100] = "None";
        LL cnt = 0;
        for (LL i = 0; i < candidateSet.size(); i++)
        {
                if (candidateSet[i].supportScores.size()>0 && pID1 == candidateSet[i].chr1 && pID2 == candidateSet[i].chr2 && strcmp(candidateSet[i].type,type)==0 && abs(pP1-candidateSet[i].point1) <= shift && abs(pP2-candidateSet[i].point2) <= shift)
                {
			for (LL t = 0; t < candidateSet[i].supportMapId.size(); t++){
                        	bool flag = true;
	                        for (LL k = 0; k < completeSet[cnt-1].supportMapId.size(); k++){
        	                        if (completeSet[cnt-1].supportMapId[k]==candidateSet[i].supportMapId[t]){
                	                        if (completeSet[cnt-1].supportScores[k] < candidateSet[i].supportScores[t]){
                        	                        completeSet[cnt-1].supportScores[k] = candidateSet[i].supportScores[t];
                                	        }
                                        	flag = false;
	                                        break;
        	                        }
                	        }
                        	if (flag){
                                	completeSet[cnt-1].supportScores.push_back(candidateSet[i].supportScores[t]);
	                                completeSet[cnt-1].supportMapId.push_back(candidateSet[i].supportMapId[t]);
                	        }
			}
			if ((LL)candidateSet[i].supportScores.size() >= cov){//always keep the regions originally having the most number of supports
				pP1 = candidateSet[i].point1;
				pP2 = candidateSet[i].point2;
				cov = (LL)candidateSet[i].supportScores.size();
			}
                }
                else
                {
                        completeSet.push_back(candidateSet[i]);
                        cnt++;
                        pID1 = candidateSet[i].chr1;
                        pID2 = candidateSet[i].chr2;
                        pP1 = candidateSet[i].point1;
                        pP2 = candidateSet[i].point2;
//			printf("%lld:%lld\t%lld:%lld\n",pID1,pP1,pID2,pP2);
			cov = (LL)candidateSet[i].supportScores.size();
			strcpy(type,candidateSet[i].type);
                }
        }
        printf("After merge, there are %lld candidate SV\n",(LL)completeSet.size());
}

void addNonSupport(vector<opticalMapType> opticalMap, int miniM){
	LL cnt = 0;
	LL cnt1 = 0;
	printf("The size of OM set: %lld\n",(LL)opticalMap.size());
	for (LL i = 0; i < completeSet.size(); i++)
	{
		for (LL j = 0; j < opticalMap.size(); j++)
		{
			if (completeSet[i].chr2 != opticalMap[j].chrId && completeSet[i].chr1 != opticalMap[j].chrId)continue;
	        	if (completeSet[i].chr1 == opticalMap[j].chrId)
			{
				if (completeSet[i].pointIdx1 > opticalMap[j].refStartIndex && completeSet[i].pointIdx1 < opticalMap[j].refEndIndex)
				{
					//LL bottomMolIdx = 0;
					//LL topMolIdx = (LL)opticalMap[j].position.size()-1;
					double score_div = calc_score_div(opticalMap[j].hitEnum,completeSet[i].pointIdx1-opticalMap[j].refStartIndex+1,miniM);
					if (score_div>0){ 
						bool flag = true;
			                        for (LL k = 0; k < completeSet[i].nonSupportMapId1.size(); k++){
                        			        if (completeSet[i].nonSupportMapId1[k]==atol(opticalMap[j].mapId)){
			                                        if (completeSet[i].nonSuppScores1[k] < score_div){
                        			                        completeSet[i].nonSuppScores1[k] = score_div;
									cnt1++;
			                                        }
                        			                flag = false;
			                                        break;
                        			        }
			                        }
                        			if (flag){
							completeSet[i].nonSuppScores1.push_back(score_div);
							completeSet[i].nonSupportMapId1.push_back(atol(opticalMap[j].mapId));
							cnt++;
							cnt1++;
						}
					}
				}
			}
            		if (completeSet[i].chr2 == opticalMap[j].chrId)
			{
				if (completeSet[i].pointIdx2 > opticalMap[j].refStartIndex && completeSet[i].pointIdx2 < opticalMap[j].refEndIndex)
                		{
					double score_div = calc_score_div(opticalMap[j].hitEnum,completeSet[i].pointIdx2-opticalMap[j].refStartIndex+1,miniM);
					if (score_div>0){ 
						bool flag = true;
                                                for (LL k = 0; k < completeSet[i].nonSupportMapId2.size(); k++){
                                                        if (completeSet[i].nonSupportMapId2[k]==atol(opticalMap[j].mapId)){
                                                                if (completeSet[i].nonSuppScores2[k] < score_div){
                                                                        completeSet[i].nonSuppScores2[k] = score_div;
									cnt1++;
                                                                }
                                                                flag = false;
                                                                break;
                                                        }
                                                }
                                                if (flag){
                                                        completeSet[i].nonSuppScores2.push_back(score_div);
                                                        completeSet[i].nonSupportMapId2.push_back(atol(opticalMap[j].mapId));
                                                        cnt++;
							cnt1++;
                                                }					
					}
                		}
			}
		}
		breakpoint &tempBP = completeSet[i];
		for (LL t1 = 0; t1 < tempBP.nonSupportMapId1.size(); t1++){
			for (LL s = 0; s < tempBP.supportMapId.size(); s++){
				if (tempBP.supportMapId[s] == tempBP.nonSupportMapId1[t1]){
					tempBP.nonSupportMapId1[t1] = -1;
					break;
				}
			}
		}
		for (LL t1 = tempBP.nonSupportMapId1.size()-1; t1 >= 0; t1--){
			if (tempBP.nonSupportMapId1[t1] == -1){
				tempBP.nonSupportMapId1.erase(tempBP.nonSupportMapId1.begin()+t1);
				tempBP.nonSuppScores1.erase(tempBP.nonSuppScores1.begin()+t1);
			}
		}
		for (LL t1 = 0; t1 < tempBP.nonSupportMapId2.size(); t1++){
			for (LL s = 0; s < tempBP.supportMapId.size(); s++){
				if (tempBP.supportMapId[s] == tempBP.nonSupportMapId2[t1]){
					tempBP.nonSupportMapId2[t1] = -1;
					break;
				}
			}
		}
		for (LL t1 = tempBP.nonSupportMapId2.size()-1; t1 >= 0; t1--){
			if (tempBP.nonSupportMapId2[t1] == -1){
				tempBP.nonSupportMapId2.erase(tempBP.nonSupportMapId2.begin()+t1);
				tempBP.nonSuppScores2.erase(tempBP.nonSuppScores2.begin()+t1);
			}
		}
		if (tempBP.nonSupportMapId1.size() > tempBP.nonSupportMapId2.size()){
			tempBP.nonSupportMapId = tempBP.nonSupportMapId2;
			tempBP.nonSuppScores = tempBP.nonSuppScores2;
		} else {
			tempBP.nonSupportMapId = tempBP.nonSupportMapId1;
			tempBP.nonSuppScores = tempBP.nonSuppScores1;
		}			
			
	}
	printf("Add %lld non-breakpoint supporting cases among total %lld\n",cnt,cnt1);
}

double callSVs(breakpoint tempBP, int suppCutoff, double varPrior){
	if (tempBP.supportScores.size() < suppCutoff||tempBP.supportScores.size()+tempBP.nonSuppScores.size()<suppCutoff*2)return 0;
	double sumSupp = accumulate(tempBP.supportScores.begin(),tempBP.supportScores.end(),0.0);
	double sumNonSupp = accumulate(tempBP.nonSuppScores.begin(),tempBP.nonSuppScores.end(),0.0);
	double sumDiff = abs(sumSupp-sumNonSupp);
	double homo = pow(0.01,sumNonSupp)*pow(0.995,sumSupp);//0.995^2~=0.99
//	LL sizeDiff = abs((LL)tempBP.nonSuppScores.size()-(LL)tempBP.supportScores.size());
	double heter = pow(0.1,sumDiff)*pow(0.99,min(sumSupp,sumNonSupp));
	double wildtype = pow(0.01,sumSupp)*pow(0.995,sumNonSupp);
//	printf("Count: (%lld:%g, %lld:%g), Homo: %g, Heter:%g, WT:%g\n",(LL)tempBP.supportScores.size(),sumSupp,(LL)tempBP.nonSuppScores.size(),sumNonSupp,homo,heter,wildtype);
	if (homo > heter && homo*varPrior > wildtype)return -10.0*log(wildtype/homo)/log(10.0);
	else if (heter > homo && heter*varPrior > wildtype)return 10.0*log(wildtype/heter)/log(10.0);
	else return 0.0;
}

void outputFile(char* outputF, int suppCutoff, double varPrior){
	FILE* output = fopen(outputF,"w");
	sort(completeSet.begin(),completeSet.end(),ss);
	fprintf(output,"#chr1\tpos1\tchr2\tpos2\ttype\tSpNum\tNspNum\tZygosity\tTorF\n");
	for (LL i = 0; i < completeSet.size(); i++){
		breakpoint tempBP = completeSet[i];
		char zygo[100];
		tempBP.score = callSVs(tempBP, suppCutoff, varPrior);
		if (tempBP.score==0.0)continue;
		else if (tempBP.score<0){strcpy(zygo,"Heterozygous");tempBP.score = -tempBP.score + 10.0*log(varPrior)/log(10.0);}
		else {strcpy(zygo,"Homozygous");tempBP.score+=10.0*log(varPrior)/log(10.0);}
		fprintf(output,"%lld\t%lld\t%lld\t%lld\t%s\t%lld\t%lld\t%s\t%g\n",tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),zygo,tempBP.score);
	}
	fclose(output);
}



/*void setDefault() {
	digestionRate = 0.875;
	falseCutRate = 0.000010;
	pValueCutOff = 1e-9;
	cauchyMean = 1.0096;
	cauchyScale = 0.0291;
	numberOfSupportEachSide = 1;
	confidenceLimit = 9;
        numberOfSupportIndelMolecule = 10; 
        numberOfSupportSignalMolecule = 10;
        likelihoodRatioCutOff = 1e6;
	minIndelSize = 2000;
	minIndelRatio = 0.05;
	resolutionLimit = 1749;
	inputTrio = 878;
}*/

bool feq(double x, double y){
	//return fabs(x - y) <= eps;
	return fabs(x - y) <= eps*(min(fabs(x),fabs(y)));//should be normalized by the values of x and y
}

bool feq1(double x, double y){
        return fabs(x - y) <= eps;
       // return fabs(x - y) <= eps*fabs(min(x,y));//should be normalized by the values of x and y
}

LL max(LL x, LL y){
	return x > y ? x : y;
}

LL min(LL x, LL y){
	return x > y ? y : x;
}

double average(vector <double> a){
	double ans = 0;
	for (LL i=0; i<(LL)a.size(); i++)
		ans += a[i];
	ans /= a.size();
	return ans;
}

char outputFolder[1000]; bool overlapped(LL x1, LL y1, LL x2, LL y2){
	if (x1 >= x2 && x1 <= y2) return true;
	if (y1 >= x2 && y1 <= y2) return true;
	if (x2 >= x1 && x2 <= y1) return true;
	if (y2 >= x1 && y2 <= y1) return true;
	return false;
}

double cauchyCDF(double x, double mean, double scale){
	return (1/pi)*atan((x-mean)/scale)+1/2;
}

double cauchyPDF(double x, double mean, double scale){
	return (1/pi)*scale/((x-mean)*(x-mean)+scale*scale);
}

void NCRcalculation(){
	for (int i = 0; i < 9000; ++i){
		C[i][0] = C[i][i] = 1;
		for (int j = 1; j < i; ++j)
			C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]);
		for (int j = i + 1; j < 9000; ++j)
			C[i][j] = 0;
	}
}

int main(int argv, char *argc[]){
	printf("\n");
	//setDefault();
	LL miniM = 3;
	LL shift = 1000*50;
	int suppCutoff = 5;
	double varPrior = 1e-6;
	char inputBN_MCR[1000];
	char inputHigh[1000];
	char inputBN_CR[1000];
	char SVoutputFile[1000];
	char SVoutputFile2[1000];
	
	if (argv == 1){
		printf("\nThe valid parameters are described as follows:\n");
		printf("\t-SVoutputFile: \n\t\t The file name of SVs (.osv).\n");
		printf("\t-optAlignFile: \n\t\t  The file name of the total split alignment map file(.oma).\n");
		printf("\t-highDensityFile: \n\t\t  The file of regions with high site density(.oma).\n");
		printf("\t-shift: \n\t\t Default value: %lld. The maximum allowed pos-error.\n",shift);
		printf("\t-miniM: \n\t\t Default value: %lld. The minimum Matches in HitEnum to be tried as useful.\n",miniM);
		return -1;
	}
	for (int i = 1; i < argv; i=i+2){
		string temp(argc[i]);
		//printf("What is the parameter?? %s \n",temp.c_str());
		if (temp.compare("-shift")==0)
			shift = atol(argc[i+1]);
		else if (temp.compare("-optAlignFile")==0)
			strcpy(inputBN_CR, argc[i+1]);
		else if (temp.compare("-highDensityFile")==0)
			strcpy(inputHigh, argc[i+1]);
		else if (temp.compare("-SVoutputFile")==0)
			strcpy(SVoutputFile2,argc[i+1]);
		else if (temp.compare("-miniM")==0)
			miniM = atol(argc[i+1]);
		else
			printf("No such parameter or wrong : %s\n",argc[i]);
	}
	readHighDensity(inputHigh);
	vector<opticalMapType> opticalMap3 = readSourceFile(inputBN_CR);


	completeSet.clear();
	candidateSet.clear();
	scanSplitedMap(opticalMap3, miniM);
	LL cnc = completeSet.size();
        mergeCandidate(shift);
        while(cnc!=completeSet.size()){
                cnc = completeSet.size();
                mergeCandidateAgain(shift);
        }
	FILE* outputCR = fopen(SVoutputFile2,"w");
	printCandidateToFile(outputCR);
	return 0;
}
