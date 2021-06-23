#include "NoError.h"
// Can choose Alden, Ref, CHM1, CHM1_2, Angel1, Angel2, Angel3, Angel4, 
// Ref_Sim, Alden_sim_T0, Sim_Ref_OMDP, Sim_Alden_OMDP 
// Pipeline_Real, Pipeline_Sim, CHM1ToASM, Autism, Schizo 
// Real_Ref, Real_Assembly_Hg38, Real_Molecule_Assembly 
// Sim_Haploid_Pipeline, Sim_Diploid_Pipeline 
// Sim_Hap_Ref, Sim_Dip_Ref 
// YH, Sim2_Haploid_OMBlast 
// Han // Alden1_1, Alden1_1000, Alden2_1, Alden2_1000 
//#define Alden2_1000

FILE* inputChrConsen; 
FILE* inputOptAlign; 
FILE* outputResult1; 
FILE* outputResult2; 
FILE* outputVariantResultFile; 
char inputAlignmentFileName[1000]; 
char outputFileLocation[500]; 
FILE *inputAlignmentFile; 
FILE *outputFileList[500]; 
LL numberOfOpticalMap; 
//#if defined (Real_Assembly_Hg38) || defined(Real_Molecule_Assembly) //FILE* outputMoleculeToAssembly; //FILE* outputAssembly; 
//#endif 
char tempString[10000]; 
LL tempLongLong; 
double tempDouble; 
double tempDistanceDouble[20000]; 
double C[9000][9000]; 
bool canOpenFile = true;


double minIndelSize; 
double minIndelRatio; 
double resolutionLimit; 
LL inputTrio; 
double digestionRate; 
double falseCutRate; 
double pValueCutOff; 
double cauchyMean; 
double cauchyScale; 
LL numberOfSupportEachSide; 
double confidenceLimit; 
LL numberOfSupportIndelMolecule; 
LL numberOfSupportSignalMolecule; 
double likelihoodRatioCutOff; 
//Start add 
struct opticalMapType{
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
//        char hitEnum[1000];
	//char hitEnum[2000];
	string hitEnum;
	double fpr=0;////
        double fnr=0;////
        double alignRate=0;////
        vector<int> position;
/*	void calc()////
        {////
                int tot_fp = 0, tot_fn = 0;
                for (int i = 0; i<hitEnum.length(); i++)
                {
                        int curr = 0;
                        int fp, fn;
                        int j;
                        for (j = i; hitEnum[j] >= '0' && hitEnum[j] <= '9'; j++)
                        {
                                curr = curr * 10 + (hitEnum[j] - '0');
                        }
                        if (hitEnum[j] == 'I')
                        {
                                fp++;
                                tot_fp += curr;
                        }
                        else if (hitEnum[j] == 'D')
                        {
                                fn++;
                                tot_fn += curr;
                        }
                        i = j;
                }
                int mol_len = 0;
                int ali_len = 0;
                int cnt = 0;
                for (auto ele:position)
                {
                        cnt++;
                        mol_len += ele;
                        if (cnt >= min(optStart,optEnd)&&cnt <= max(optStart,optEnd)) ali_len+=ele;
                }
                fpr = tot_fp*1.0/ali_len;
                fnr = tot_fn*1.0/(abs(refEndIndex-refStartIndex)+1);
                alignRate = ali_len*1.0/mol_len;
        }////*/
        void print(){
		printf("%s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld ", mapId, score, confidence, fpr,fnr,alignRate, chrId, optStart, optEnd, refStart, refEnd, hitEnum.c_str(), (LL)position.size());
//                printf("%s %lf %lf %lld %lld %lld %lld %lld %s %lld ", mapId, score, confidence, chrId, optStart, optEnd, refStart, refEnd, hitEnum, (LL)position.size());
                for (LL i=0; i<(LL)position.size(); i++)
                        printf("%d ", position[i]);
                printf("\n");
        }
        void print(FILE* targetFile){
                sprintf(tempString, "%lld", inputTrio);
		fprintf(targetFile, "%s %s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld ", tempString, mapId, score, confidence, fpr, fnr, alignRate, chrId, optStart, optEnd, refStart, refEnd, hitEnum.c_str(), (LL)position.size());
//                fprintf(targetFile, "%s %s %lf %lf %lld %lld %lld %lld %lld %s %lld ", tempString, mapId, score, confidence, chrId, optStart, optEnd, refStart, refEnd, hitEnum, (LL)position.size());
                for (LL i=0; i<(LL)position.size(); i++)
                        fprintf(targetFile, "%d ", position[i]);
                fprintf(targetFile, "\n");
        }
};
bool ss(opticalMapType q, opticalMapType w){
        LL e = strcmp(q.mapId, w.mapId);
        return ((e<0?true:false) || (e==0 && q.chrId < w.chrId) || (e == 0 && q.chrId == w.chrId && q.optStart < w.optEnd));
}

vector<opticalMapType> opticalMap1;

void readSourceFile(){
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
        int numberOfDL = 0;
        char ttt;
        char ttts[100000];
        while(fgetc(inputAlignmentFile) == '#')
        {
                numberOfDL++;
                fgets(ttts,100000,inputAlignmentFile);
        }
	LL totSize = 0;
        while(fscanf(inputAlignmentFile,"%s",ttts)==1){
                totSize++;
                fgets(ttts,100000,inputAlignmentFile);
        }
        fclose(inputAlignmentFile);
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
	opticalMap1.clear();
	opticalMap1.resize(totSize+10000);
//        opticalMap1 = new opticalMapType[20000000];
        for (LL i=0; i<numberOfDL; i++)
                fgets(tempString, 5000, inputAlignmentFile);
        numberOfOpticalMap = 0;
        while (fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].mapId) == 1){
//                if (numberOfOpticalMap<5 || numberOfOpticalMap % 100000 == 0) printf("%lld %s\n", numberOfOpticalMap, opticalMap1[numberOfOpticalMap].mapId);
                LL tempNumberOfSites, tempLL;
                fscanf(inputAlignmentFile, "%lld", &tempNumberOfSites);
                LL tempIndex = 0;
                for (LL i=0; i<tempNumberOfSites; i++){
                        if (i != tempNumberOfSites - 1)
                                fscanf(inputAlignmentFile, "%lld;", &tempLL);
                        else fscanf(inputAlignmentFile, "%lld", &tempLL);
                        tempIndex += tempLL;
                        opticalMap1[numberOfOpticalMap].position.push_back((int)tempIndex);
                }
                fscanf(inputAlignmentFile, "%s", tempString);
                if (tempString[0] == 'c'){
                        if (tempString[3] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (tempString[3] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (tempString[3] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
                        else
                                sscanf(tempString, "chr%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                else{
                        if (tempString[0] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (tempString[0] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (tempString[0] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
			else
                        	sscanf(tempString, "%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                fscanf(inputAlignmentFile, "%s", tempString);
                if (tempString[0] == 'r' || tempString[0] == '-') opticalMap1[numberOfOpticalMap].orientation = false; else opticalMap1[numberOfOpticalMap].orientation = true;
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].score);
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].confidence);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStartIndex, &opticalMap1[numberOfOpticalMap].refEndIndex);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].optStart, &opticalMap1[numberOfOpticalMap].optEnd);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStart, &opticalMap1[numberOfOpticalMap].refEnd);
		char hitE[5000];
                memset(hitE,0,sizeof(hitE));
                fscanf(inputAlignmentFile, "%s", hitE);
		opticalMap1[numberOfOpticalMap].hitEnum = hitE;
		
              //  int siz = sizeof(hitE)/sizeof(char);
                //opticalMap1[numberOfOpticalMap].hitEnum = new char[siz+1];
//                memset(opticalMap1[numberOfOpticalMap].hitEnum,0,sizeof(opticalMap1[numberOfOpticalMap].hitEnum));
  //              strncpy(opticalMap1[numberOfOpticalMap].hitEnum,hitE,siz);
//                fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].hitEnum);
//		opticalMap1[numberOfOpticalMap].calc();////
		
                numberOfOpticalMap++;
		if (numberOfOpticalMap==totSize){
                        totSize+=1000000;
                        opticalMap1.resize(totSize);
                }
        }
	opticalMap1.resize(numberOfOpticalMap);
        printf("Number of optical map: %lld\n", numberOfOpticalMap);
	LL old_i = 0, new_i=0;////
        for (LL i=0; i<numberOfOpticalMap; i++){
                if (opticalMap1[i].orientation){
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart - 1];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd];
                } else {
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd - 1];
                }
                if (!opticalMap1[i].orientation){
                        LL tempMaxDis = opticalMap1[i].position[(int)opticalMap1[i].position.size() - 1];
                        for (LL j=0; j<(LL)opticalMap1[i].position.size(); j++)
                                opticalMap1[i].position[j] = (int)tempMaxDis - opticalMap1[i].position[j];
                        sort(opticalMap1[i].position.begin(), opticalMap1[i].position.end());
                        opticalMap1[i].optStart = tempMaxDis - opticalMap1[i].optStart;
                        opticalMap1[i].optEnd = tempMaxDis - opticalMap1[i].optEnd;
                        opticalMap1[i].orientation = true;
                }
		new_i = i;////
                if (new_i == 0) continue;////
                if (strcmp(opticalMap1[new_i].mapId,opticalMap1[old_i].mapId)==0)////
                {
                        for(LL tm = old_i; tm < new_i; tm++)
                        {
                                opticalMap1[tm].alignRate += opticalMap1[new_i].alignRate;
                        }
                        opticalMap1[new_i].alignRate = opticalMap1[old_i].alignRate;
                }
                else
                {
                        old_i = new_i;
                }
        }
}

void initOutput(int chr, const char *type){
        char buffer[200];
        char nameOfFile[1000];
        strcpy(nameOfFile, outputFileLocation);
        sprintf(buffer, "%lld_%d",inputTrio,chr);
        strcat(nameOfFile, buffer);
        strcat(nameOfFile, ".bmap");
        if ((outputFileList[chr] = fopen(nameOfFile, type)) == NULL)
        {
                perror("error opening file(): nameOfFile");
        }
}

void addSplitedMap(){
        LL tempCC = numberOfOpticalMap;
	sort(opticalMap1.begin(), opticalMap1.end(), ss);
	LL totSize = numberOfOpticalMap+50000;
	opticalMap1.resize(totSize);
        for (LL i=0; i<numberOfOpticalMap-1; i++){
                if (strcmp(opticalMap1[i].mapId, opticalMap1[i+1].mapId) == 0
                                && opticalMap1[i].chrId == opticalMap1[i+1].chrId
                                && opticalMap1[i+1].refStart - opticalMap1[i].refEnd > 0
                                && opticalMap1[i+1].refStart - opticalMap1[i].refEnd < 5000000){
                        opticalMap1[tempCC] = opticalMap1[i];
                        opticalMap1[tempCC].optStart = opticalMap1[i].optEnd;
                        opticalMap1[tempCC].optEnd = opticalMap1[i+1].optStart;
                        opticalMap1[tempCC].refStart = opticalMap1[i].refEnd;
                        opticalMap1[tempCC].refEnd = opticalMap1[i+1].refStart;
			opticalMap1[tempCC].score = 10.0;
                        opticalMap1[tempCC].confidence = 0.1;
                        opticalMap1[tempCC].orientation = true;
//			opticalMap1[tempCC].calc();////
                        opticalMap1[tempCC].fpr = opticalMap1[i].fpr + opticalMap1[i+1].fpr;////
                        opticalMap1[tempCC].fnr = opticalMap1[i].fnr + opticalMap1[i+1].fnr;////
                        //opticalMap1[tempCC].fpr = opticalMap1[i].fpr + opticalMap1[i+1].fpr;////
                        opticalMap1[tempCC].position.clear();
//			opticalMap1[tempCC].hitEnum = new char[4];
                        opticalMap1[tempCC].hitEnum= "FFF";
                        if (opticalMap1[tempCC].optStart > opticalMap1[tempCC].optEnd)
                                tempCC--;
                        tempCC++;
			if (tempCC==totSize){
				totSize+=50000;
				opticalMap1.resize(totSize);
			}
                }
        }
        numberOfOpticalMap = tempCC;
	opticalMap1.resize(numberOfOpticalMap);
}

void outputToBillMapDestinationFile(){
        for (LL i=0; i<numberOfOpticalMap; i++){
//                if (i % 100000 == 0) printf("%lld\n", i);
		//if (opticalMap1[i].chrId > numOfChr) continue;
                if (outputFileList[opticalMap1[i].chrId] == NULL)
                        initOutput(opticalMap1[i].chrId,"w+");
                else
                        initOutput(opticalMap1[i].chrId,"a+");
                opticalMap1[i].print(outputFileList[opticalMap1[i].chrId]);
                fclose(outputFileList[opticalMap1[i].chrId]);
        }
	opticalMap1.clear();
	vector<opticalMapType>().swap(opticalMap1);
}


//stop add

void setDefault() {
	digestionRate = 0.875;
	falseCutRate = 0.000010;
	pValueCutOff = 1e-9;
	cauchyMean = 1.0096;
	cauchyScale = 0.0291;
	numberOfSupportEachSide = 1;
	confidenceLimit = 9;
        numberOfSupportIndelMolecule = 10; //coverage of the indel
        numberOfSupportSignalMolecule = 10;
        likelihoodRatioCutOff = 1e6; // Settings
	minIndelSize = 2000;
	minIndelRatio = 0.05;
	resolutionLimit = 700;
	inputTrio = 0507;
}

LL newChrSite = 0; 
LL numberOfHomoMissingFragment = 0; 
LL numberOfHomoDeletedSNP = 0; 
LL homoMissingSite = 0; 
LL numberOfNoSupportSV = 0; 
LL homoNonMissingSite = 0; 
LL homoInsertSite = 0; 
LL homoNonInsertSite = 0; LL heterMissingSite = 0; LL heterNonMissingSite = 0; LL heterInsertSite = 0; LL heterNonInsertSite = 0; LL numberOfNickingSite = 0; LL numberOfOM878 = 0; LL numberOfOM891 = 0; 
LL numberOfDetailOpticalMap = 0; LL distancePairCount = 0; LL numberOfVariant = 0; LL numberOfSupportedSV = 0;

LL tempPre = -1; LL tempCounter = 0; vector <LL> listOfChromosome;

// shortestDist[10] refers to distance between 9-11 
//double shortestDist[500000]; 
//LL shortestDistCount[500000];

// shortestGapDist[10] refers to distance between 10-11 
//double shortestGapDist[500000]; 
//LL shortestGapDistCount[500000];


// 0 = Homo Ins/Del // 1 = Homo Non-Ins/Del // 2 = Heter Ins/Del // 3 = Heter Non-Ins/Del // 4 = nothing

LL curId; int chrId;


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


struct variantType{
	LL people;
	LL chr;
	LL start;
	LL end;
	LL size;
	LL sizeExt;
	LL support;
//	double score;
//	double confidence;
//	double fpr;////
  //      double fnr;////
//        double alignRate;////
        LL oppoSupp;////
	double ratio;
	double ratioExt;
	double likelihood;
	bool isSignal;
	bool isDel;
	bool isHomo;
	bool isSupported;
	void update(LL inputPeople, LL inputChr, LL inputStart, LL inputEnd, LL inputSize, LL inputSizeExt, LL inputSignal, LL inputDel, LL inputHomo, LL inputSupport, double inputRatio, double inputRatioExt, double likely, LL supp, LL IoppoSupp)
	{//version for SV
//		fpr = Ifpr;////
  //              fnr = Ifnr;////
    //            alignRate = IalignRate;////
                oppoSupp = IoppoSupp;////
//		score = Iscore;
//		confidence = Iconfidence;
		people = inputPeople;
		chr = inputChr;
		start = inputStart;
		end = inputEnd;
		size = inputSize;
		sizeExt = inputSizeExt;
		ratio = inputRatio;
		ratioExt = inputRatioExt;
		likelihood = likely;
		support = supp;
		isSignal = (inputSignal == 1);
		isDel = (inputDel == 1);
		isHomo = (inputHomo == 1);
		isSupported = (inputSupport == 1);
	};
	void update(LL inputPeople, LL inputChr, LL inputStart, LL inputEnd, LL inputSize, LL inputSignal, LL inputDel, LL inputHomo, LL inputSupport, double inputRatio, double likely, LL supp,  LL IoppoSupp)
	{//version for signal changes
//		fpr = Ifpr;////
  //              fnr = Ifnr;////
    //            alignRate = IalignRate;////
                oppoSupp = IoppoSupp;////
//		score = Iscore;
//		confidence = Iconfidence;
                people = inputPeople;
                chr = inputChr;
                start = inputStart;
                end = inputEnd;
                size = inputSize;
                sizeExt = 0;
                ratio = inputRatio;
                ratioExt = 0;
		likelihood = likely;
		support = supp;
                isSignal = (inputSignal == 1);
                isDel = (inputDel == 1);
                isHomo = (inputHomo == 1);
                isSupported = (inputSupport == 1);
        };
	bool operator<(const variantType &x) const{
		return (chr < x.chr || (chr == x.chr && start < x.start) || (chr == x.chr && start == x.start && end < x.end) || (chr == x.chr && start == x.start && end == x.end && isSupported));
	}

	/*void print(){
		printf("%lld %lld %lld %lld %lld %d %d %d\n", people, chr, start, end, size, isSignal, isDel, isHomo);
	};*/
};
struct chrType{
	LL length;
	LL numberOfSites;
//	LL position[50000];
	// distance[0] = position[1] - position[0]
//	LL distance[50000];
//	LL oldDistance[50000];
//	LL gapDistance[50000];
	// coverage[1] refer to nicking site 1
//	LL coverage[50000];
	// occurrence[1] refer to nicking site 1
//	LL occurrence[50000];
	// gapCount[1] refers to gap between site 1 to 2
//	LL gapCount[50000];
//	LL gapCoverage[50000];
	// If it is significant by Kevin's definition
//	bool gapSigni[50000];
//	bool signi[50000];
	// distance Difference
//	double distanceDifference[50000];
//	double gapDistanceDifference[50000];
	LL* position = new LL[10];
        LL* distance = new LL[10];
        int* coverage = new int[10];
        int* occurrence = new int[10];
        int* gapCount = new int[10];
        int* gapCoverage = new int[10];
        bool* gapSigni = new bool[10];
        bool* signi = new bool[10];
};
chrType chromosome;

struct optAlignType{
	LL belongs;
	char mapId[100];
	LL optStart;
	LL optEnd;
	LL refStart;
	LL refEnd;
	LL numberOfSites;
	double score;
	double confidence;
	//char hitEnum[2000];
        string hitEnum;
	// position = distance[i+1] - distance[i]
	int position[2500];
	int oldPosition[2500];
	double fpr;////
        double fnr;////
        double alignRate;////

};

struct distanceType{
	char mapId[20];
	LL start;
	LL end;
//	double score;
//	double confidence;
//	double fpr;////
  //      double fnr;////
    //    double alignRate;////
	double distance;
	LL cnt;
	bool operator<(const distanceType &x) const{
		return (start < x.start || (start == x.start && end < x.end));
	}
	void print(){
		printf("chrId:%d, id:%s start:%lld end:%lld oldDis:%lld, molDis:%lf cnt:%lld\n", chrId, mapId, chromosome.position[start], chromosome.position[end], chromosome.position[end]-chromosome.position[start], distance, cnt);
	}
	void printToFile(FILE* outputFile){
		fprintf(outputFile, "chrId:%d, id:%s start:%lld end:%lld oldDis:%lld, molDis:%lf cnt:%lld\n", chrId, mapId, chromosome.position[start], chromosome.position[end], chromosome.position[end] - chromosome.position[start], distance, cnt);
	}
};
vector<optAlignType> opticalMap;
vector<distanceType> distancePair;
vector<variantType> variant;
void init(){
	// Has reduced by half
//	opticalMap = new optAlignType[3000000];
//	distancePair = new distanceType[100000000];
	variant.clear();
	variant.resize(100000);
}
char outputFolder[1000];

FILE* createInputTrioFile(const char *suffix, LL inputTrio, const char *mode){
	char nameOfFile[1000];
	char buffer[200];
	memset(buffer, 0, sizeof(buffer));
	strcpy(nameOfFile, outputFolder);
	sprintf(buffer, "%lld_%g_", inputTrio,minIndelSize);
	strcat(nameOfFile, buffer);
	strcat(nameOfFile, suffix);
	//printf("The file name is: %s\nAnd mode is %s\n",nameOfFile,mode);
	FILE* retFile = fopen(nameOfFile, mode);
	if (retFile == NULL) perror("error in opening file!");
	return retFile;
}



void initResultFullList(char* SVoutputFile){
	char nameOfFile[1000];
	char buffer[200];
	char extPart[200];
	outputResult1 = createInputTrioFile("Distance_Pairs_List.txt", inputTrio, "w+");
        outputResult2 = createInputTrioFile("Signal_Indels_List.txt", inputTrio, "w+"); /*
#if defined(Real_Molecule_Assembly)
	outputMoleculeToAssembly = createInputTrioFile("/local/shared/bill/BillMapFromREAL_Molecule_ASM_HG38/distanceList/", "moleculeToAssembly_List.txt", inputTrio, "w+");
#endif if defined(Real_Assembly_Hg38)
*/
	sprintf(extPart,"_%g",minIndelSize);
	memset(nameOfFile, 0, sizeof(nameOfFile));
	strcpy(nameOfFile, outputFolder);
	sprintf(buffer, "%lld", inputTrio);
	strcat(nameOfFile, buffer);
	strcat(nameOfFile, SVoutputFile);
	strcat(nameOfFile, extPart);
	strcat(nameOfFile, ".osv");
//	printf("%s\n",nameOfFile);
	if ((outputVariantResultFile = fopen(nameOfFile, "w+")) == NULL){
		perror("error opening file(): nameOfFile");
	}
}

bool overlapped(LL x1, LL y1, LL x2, LL y2){
	if (x1 >= x2 && x1 <= y2) return true;
	if (y1 >= x2 && y1 <= y2) return true;
	if (x2 >= x1 && x2 <= y1) return true;
	if (y2 >= x1 && y2 <= y1) return true;
	return false;
}
void getChromosomeList(char* chrMapFile){
	listOfChromosome.clear();
	char nameOfFile[1000];
	strcpy(nameOfFile, chrMapFile);
	if ((inputChrConsen = fopen(nameOfFile, "r")) == NULL){
		perror("error opening file(): ChromosomeInfo");
	}
	fscanf(inputChrConsen,"%s",tempString);
	while(tempString[0]=='#'){
		fgets(tempString, 10000, inputChrConsen);
		fscanf(inputChrConsen,"%s",tempString);
	}
	LL tempChrId, tempNumberOfSites;
	double tempDouble1, tempDouble2;
	tempChrId = atoll(tempString);
	if (fscanf(inputChrConsen, "%lf %lld %lld %lld %lf %lf %lf %lf", &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) != 8) perror("Wrong in read chromosome maps!");
	listOfChromosome.push_back(tempChrId);
        fgets(tempString, 10000, inputChrConsen);
	while (fscanf(inputChrConsen, "%lld %lf %lld %lld %lld %lf %lf %lf %lf", &tempChrId, &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) == 9){
		listOfChromosome.push_back(tempChrId);
		fgets(tempString, 10000, inputChrConsen);
	}
	fclose(inputChrConsen);
	sort(listOfChromosome.begin(), listOfChromosome.end());
	vector<LL>::iterator tempIt;
	tempIt = unique(listOfChromosome.begin(), listOfChromosome.end());
	listOfChromosome.resize(distance(listOfChromosome.begin(), tempIt));
	printf("Number of Chromosomes: %lld\n", (LL)listOfChromosome.size());
}

void readChromosomeInfo(int chr,char* chrMapFile){
	char nameOfFile[1000];
	strcpy(nameOfFile, chrMapFile);
	if ((inputChrConsen = fopen(nameOfFile, "r")) == NULL){
		perror("error opening file(): ChromosomeInfo");
	}
	fscanf(inputChrConsen,"%s",tempString);
        while(tempString[0]=='#'){
                fgets(tempString, 10000, inputChrConsen);
                fscanf(inputChrConsen,"%s",tempString);
        }
	LL cc = 0;
	LL tempChrId, tempNumberOfSites;
	double tempDouble1, tempDouble2;
	tempChrId = atoll(tempString);
	bool done= false;
	delete[] chromosome.position;
	delete[] chromosome.distance;
	delete[] chromosome.coverage;
        delete[] chromosome.occurrence;
        delete[] chromosome.gapCount;
        delete[] chromosome.gapCoverage;
        delete[] chromosome.gapSigni;
        delete[] chromosome.signi;
	if (fscanf(inputChrConsen, "%lf %lld %lld %lld %lf %lf %lf %lf", &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) != 8) perror("Wrong in read chromosome maps!");
	if (tempChrId == chr){
                chromosome.numberOfSites = tempNumberOfSites;
		chromosome.position = new LL[tempNumberOfSites+10];
                chromosome.distance = new LL[tempNumberOfSites+10];
                chromosome.coverage = new int[tempNumberOfSites+10];
                chromosome.occurrence = new int[tempNumberOfSites+10];
                chromosome.gapCount = new int[tempNumberOfSites+10];
                chromosome.gapCoverage = new int[tempNumberOfSites+10];
                chromosome.gapSigni = new bool[tempNumberOfSites+10];
                chromosome.signi = new bool[tempNumberOfSites+10];
                done = true;
                chromosome.position[cc] = (LL) tempDouble2;
                chromosome.length = (LL) tempDouble1;
                chromosome.gapCount[cc] = 0;
                chromosome.coverage[cc] = 0;
                chromosome.occurrence[cc] = 0;
                cc++;
        }
        fgets(tempString, 10000, inputChrConsen);
	while (fscanf(inputChrConsen, "%lld %lf %lld %lld %lld %lf %lf %lf %lf", &tempChrId, &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) == 9){
		if (tempChrId == chr){
			chromosome.numberOfSites = tempNumberOfSites;
			if (!done){
                                chromosome.position = new LL[tempNumberOfSites+10];
                                chromosome.distance = new LL[tempNumberOfSites+10];
                                chromosome.coverage = new int[tempNumberOfSites+10];
                                chromosome.occurrence = new int[tempNumberOfSites+10];
                                chromosome.gapCount = new int[tempNumberOfSites+10];
                                chromosome.gapCoverage = new int[tempNumberOfSites+10];
                                chromosome.gapSigni = new bool[tempNumberOfSites+10];
                                chromosome.signi = new bool[tempNumberOfSites+10];
                                done = true;
                        }
			chromosome.position[cc] = (LL) tempDouble2;
			chromosome.length = (LL) tempDouble1;
			chromosome.gapCount[cc] = 0;
			chromosome.coverage[cc] = 0;
			chromosome.occurrence[cc] = 0;
			cc++;
		}
		fgets(tempString, 10000, inputChrConsen);
	}
        chromosome.numberOfSites++;
	for (LL i=0; i<chromosome.numberOfSites-1; i++){
		chromosome.distance[i] = chromosome.position[i+1] - chromosome.position[i];
//		chromosome.oldDistance[i] = chromosome.distance[i];
	}
	printf("chromosome's Sites: %lld\n", chromosome.numberOfSites);
	fclose(inputChrConsen);
}
void readOpticalAlign(int chr, char* outputFileLocation){
	char buffer[200];
	char nameOfFile[1000];
	strcpy(nameOfFile, outputFileLocation);
	sprintf(buffer, "%lld_%d",inputTrio,chr);
	strcat(nameOfFile, buffer);
	strcat(nameOfFile, ".bmap");
	if ((inputOptAlign = fopen(nameOfFile, "r")) == NULL){
		canOpenFile = false;
		return;
	}
	char temp[10000],temp2[100000];
	LL totSiz = 0;
	while (fscanf(inputOptAlign,"%s",temp)==1){
		totSiz++;
		fgets(temp2,100000,inputOptAlign);
	}
	fclose(inputOptAlign);
	opticalMap.clear();
	opticalMap.resize(totSiz+1);
	inputOptAlign = fopen(nameOfFile, "r");
	LL cc = 0;
	char ht[200000];
	memset(ht,0,sizeof(ht));
	while (fscanf(inputOptAlign, "%lld %s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld", &opticalMap[cc].belongs, opticalMap[cc].mapId, &opticalMap[cc].score, &opticalMap[cc].confidence, &opticalMap[cc].fpr, &opticalMap[cc].fnr, &opticalMap[cc].alignRate, &tempLongLong, &opticalMap[cc].optStart, &opticalMap[cc].optEnd, &opticalMap[cc].refStart, &opticalMap[cc].refEnd, ht, &opticalMap[cc].numberOfSites) == 14){
		opticalMap[cc].hitEnum = ht;
		memset(ht,0,sizeof(ht));
		int tempPosition[5000];
//		opticalMap[cc].position = new int[opticalMap[cc].numberOfSites];
  //              opticalMap[cc].oldPosition = new int[opticalMap[cc].numberOfSites];
		for (LL i=0; i<opticalMap[cc].numberOfSites; i++){
			fscanf(inputOptAlign, "%d", &opticalMap[cc].position[i]);
			opticalMap[cc].oldPosition[i] = opticalMap[cc].position[i];
			tempPosition[i] = opticalMap[cc].position[i];
		}
		for (LL i=0; i<opticalMap[cc].numberOfSites; i++){
			if (i == 0)
				opticalMap[cc].position[i] = 0;
			else
				opticalMap[cc].position[i] = tempPosition[i] - tempPosition[i-1];
		}
		cc++;
	}
	opticalMap[cc].refStart = -1;
	opticalMap[cc].refEnd = -1;
	printf("Number of optical map(cc) of chr %d: %lld\n", chr, cc);
	fclose(inputOptAlign);
}
void printDistancePair(const char *fileName){
	FILE* once;
	if((once = fopen(fileName,"w+"))==NULL)
	{
		printf("Something wrong when open distance files!\n");
		return;
	}
	fprintf(once,"%lld\n", distancePairCount);
	for (LL i=0; i<distancePairCount; i++){
		distancePair[i].printToFile(once);
	}
}

double cauchyCDF(double x, double mean, double scale){
	return (1/pi)*atan((x-mean)/scale)+1/2;
}

double cauchyPDF(double x, double mean, double scale){
	return (1/pi)*scale/((x-mean)*(x-mean)+scale*scale);
}


double calculateProbDistancePairNormalDistribution(double *x, LL size){
	if (size == 1 || size == 0) return 1.0;
	double standardDeviation;
	double returnProbability = 1;
	double medium = 0;
	if (size % 2 == 0)
		medium = (x[size/2-1] + x[size/2])/2;
	else medium = x[size/2];
	standardDeviation = 0.1*medium;
	for (LL i=0; i<size; i++){
		double tempZScore = (x[i]-medium)/standardDeviation;
		double tempProba = (1/(sqrt(2*pi)))*exp(-(tempZScore*tempZScore)/2);
		tempProba = max(tempProba, 1e-1);
		returnProbability *= tempProba;
	}
	return returnProbability;
};

// If deletion, return 1 // If Insertion, return 2 
const double BOTTOM_LIKELIHOOD = 1e-100; 
void advanceLikelihoodDistanceCalculation(LL start, LL end){
	double referenceDistance = chromosome.position[distancePair[start].end] - chromosome.position[distancePair[start].start];
	if (end-start + 1 < numberOfSupportIndelMolecule || referenceDistance <= resolutionLimit) {//why 2500?
//		variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
//		variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0);
//		variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
//		variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0);
		numberOfNoSupportSV++;
		return;
	}
	double nullHypothesisLikelihood = 1;
	double homoIndelHypothesisLikelihood = 1;
	double heterDelHypothesisLikelihood = BOTTOM_LIKELIHOOD;
	double heterInsHypothesisLikelihood = BOTTOM_LIKELIHOOD;
	double heterDelHypothesisMaxLikelihoodSizeChanged = -1;
	double heterInsHypothesisMaxLikelihoodSizeChanged = -1;
	double heterDelHypothesisMaxLikelihood = -1;
	double heterInsHypothesisMaxLikelihood = -1;
	LL cnt = 0;
	memset(tempDistanceDouble,0,sizeof(tempDistanceDouble));

	// For heterozygous deletion case
	// Find all previous distance pairs that overlapped the current region

//	map<string, LL> tempMap;
//	vector<distanceType> tempVec[10000];
	LL oppoSupp = 0;
/*	bool vectorOverFlow = false;
	for (LL i=start-1; i>=0; i--){
		if (distancePair[i].start >= distancePair[start].start && distancePair[i].end <= distancePair[start].end){
			string mapFrom(distancePair[i].mapId);
			LL mapTo = (LL)tempMap.size() + 1;
			if (tempMap.count(mapFrom) == 0){
				tempMap.insert(make_pair(mapFrom, mapTo));
			}
			LL newIndex = tempMap[mapFrom];
			if (newIndex >= 10000){
				vectorOverFlow = true;
				break;
			}
			tempVec[newIndex].push_back(distancePair[i]);
		} else
			break;
	}
	if (vectorOverFlow) return;
	for (LL i=end+1; i<distancePairCount; i++){
		if (distancePair[i].start >= distancePair[start].start && distancePair[i].end <= distancePair[start].end){
			string mapFrom(distancePair[i].mapId);
			LL mapTo = (LL)tempMap.size() + 1;
			if (tempMap.count(mapFrom) == 0){
				tempMap.insert(make_pair(mapFrom, mapTo));
			}
			LL newIndex = tempMap[mapFrom];
			if (newIndex >= 10000){
				vectorOverFlow = true;
				break;
			}
			tempVec[newIndex].push_back(distancePair[i]);
		} else
			break;
	}
	if (vectorOverFlow) return;*/
	double averageScore = 0;////
        double averageConfidence = 0;////
        double averageFPR = 0;////
        double averageFNR = 0;////
        double averageAlignRate = 0;////
/*	for (LL i=1; i<(LL)tempMap.size()+1; i++){
		double tempMoleCoverSize = 0;
		double tempMoleRealSize = 0;
		double tempScore = 0;////
                double tempConfidence =0;////
                double tempFPR = 0;////
                double tempFNR = 0;////
                double tempAlignRate = 0;////
		for (LL j=0; j<(LL)tempVec[i].size(); j++){
			tempMoleCoverSize += (chromosome.position[tempVec[i][j].end] - chromosome.position[tempVec[i][j].start]);
			tempMoleRealSize += tempVec[i][j].distance;
		}
		if (feq1(tempMoleCoverSize, referenceDistance))
		{
			tempDistanceDouble[cnt++] = tempMoleRealSize/referenceDistance;
			if (abs(tempMoleRealSize - referenceDistance)*2 <= minIndelSize) oppoSupp++;////
		}
	}
	//double averageScore = 0;
	//double averageConfidence = 0;
	for (LL i=start; i<=end; i++)
	{
		tempDistanceDouble[cnt++] = distancePair[i].distance/referenceDistance;
		if (abs(distancePair[i].distance - referenceDistance)*2 <= minIndelSize) oppoSupp++;////
	}
	sort(tempDistanceDouble, tempDistanceDouble + cnt);*/

	vector<LL> maps;
        vector<double> dists;
        maps.clear();
        dists.clear();
//        fprintf(outputResult1, "%d %lld %lld %lld %lld\t", chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], distancePair[start].start+1, distancePair[start].end+1);
        for (LL i=start; i<=end; i++)
        {
                bool fal = true;
                for (LL j = 0; j < (LL)maps.size(); j++)
                {
                        if (maps[j] == atol(distancePair[i].mapId) && abs(dists[j] - distancePair[i].distance) < 1){fal = false;break;}
                }
                if (fal)
                {
                        maps.push_back(atol(distancePair[i].mapId));
                        dists.push_back(distancePair[i].distance);
                        tempDistanceDouble[cnt++] = distancePair[i].distance/referenceDistance;
                       // fprintf(outputResult1, "%lld:%lld ",atol(distancePair[i].mapId),(LL)distancePair[i].distance);
                        //printf("Let me see: %g\t%g\t%lld\n", distancePair[i].distance,referenceDistance,oppoSupp);
                        if (abs(distancePair[i].distance - referenceDistance)*2 <= minIndelSize) oppoSupp++;////
                }
        }
  //      fprintf(outputResult1, "%lld\n",(LL)maps.size());
    //    for (LL i = 0; i < (LL)maps.size(); i++)
      //          fprintf(outputResult1, "%lld:%lld ",maps[i],(LL)dists[i]);
        //fprintf(outputResult1, "\n");
//        fprintf(outputResult1, "%lld\n",cnt);
        sort(tempDistanceDouble, tempDistanceDouble + cnt);
	if (cnt < numberOfSupportIndelMolecule) {//why 2500?
  //              variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    //            variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0);
      //          variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
        //        variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0);
                numberOfNoSupportSV++;
                return;
        }
	//use another vector to remove the outliers in two ends to purify the results
        //int outRange = static_cast<int>(cnt/10)>1?static_cast<int>(cnt/10):1;
	int outlierRange = 0;
        double tempDistanceDoubleNoOutlier[cnt];
        LL cntNoOutlier = 0;
        for (LL i = outlierRange; i < cnt-outlierRange;)
        {
                tempDistanceDoubleNoOutlier[cntNoOutlier++] = tempDistanceDouble[i++];
        }
	double sampleMean = 0;
	for (LL i=0; i<cntNoOutlier; i++)
		sampleMean += tempDistanceDoubleNoOutlier[i];
	sampleMean /= cntNoOutlier;
	// Change sampleMean to sampleMedian
	sampleMean = (cntNoOutlier%2==1)?(tempDistanceDoubleNoOutlier[cntNoOutlier/2]+tempDistanceDoubleNoOutlier[(cntNoOutlier/2)+1])/2:tempDistanceDoubleNoOutlier[cntNoOutlier/2];
	// Null hypothesis
	for (LL i=0; i<cntNoOutlier; i++)
		nullHypothesisLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[i], cauchyMean, cauchyScale);
	// HomoIndel hypothesis
	for (LL i=0; i<cntNoOutlier; i++)
		homoIndelHypothesisLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[i], sampleMean, cauchyScale);
	double tempSizeChange = (sampleMean - cauchyMean) * referenceDistance;
	if (fabs(tempSizeChange) < max(minIndelSize, minIndelRatio * referenceDistance))
	{
		tempSizeChange = 0;
		homoIndelHypothesisLikelihood = nullHypothesisLikelihood/2;
	}
	// HeterDel/HeterIns hypothesis
	numberOfSupportEachSide = max(5, (2*cntNoOutlier/5)); //why this formula?
	for (LL i=numberOfSupportEachSide; i<cntNoOutlier-numberOfSupportEachSide; i++){
		double mean1 = 0, mean2 = 0;

		// Try Median instead of Mean.
		mean1 = (i%2==1)?(tempDistanceDoubleNoOutlier[i/2]+tempDistanceDoubleNoOutlier[(i/2)+1])/2:tempDistanceDoubleNoOutlier[i/2];
		LL tempSecondPartSize = cntNoOutlier - i;
		mean2 = (tempSecondPartSize%2==1)?(tempDistanceDoubleNoOutlier[i+tempSecondPartSize/2]+tempDistanceDoubleNoOutlier[i+tempSecondPartSize/2+1])/2:tempDistanceDoubleNoOutlier[i+tempSecondPartSize/2];


		double tempLikelihood;

		// For HeterDel
		tempLikelihood = 1;
		for (LL j=0; j<i; j++)//remove the most least two
			tempLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[j], mean1, cauchyScale);
		for (LL j=i; j<cntNoOutlier; j++)//remove the most largest two
			tempLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[j], cauchyMean, cauchyScale);

		heterDelHypothesisLikelihood += tempLikelihood;
		if (heterDelHypothesisMaxLikelihood < tempLikelihood)
		{
			//printf("There are really this case? (%f\t%f\n)",heterDelHypothesisLikelihood,tempLikelihood);
			heterDelHypothesisMaxLikelihood = tempLikelihood;
			heterDelHypothesisMaxLikelihoodSizeChanged = (mean1 - cauchyMean)*referenceDistance;
		}

		// For HeterIns
		tempLikelihood = 1;
		for (LL j=0; j<i; j++)
			tempLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[j], cauchyMean, cauchyScale);
		for (LL j=i; j<cntNoOutlier; j++)
			tempLikelihood *= cauchyPDF(tempDistanceDoubleNoOutlier[j], mean2, cauchyScale);

		heterInsHypothesisLikelihood += tempLikelihood;
		if (heterInsHypothesisMaxLikelihood < tempLikelihood){
			heterInsHypothesisMaxLikelihood = tempLikelihood;
			heterInsHypothesisMaxLikelihoodSizeChanged = (mean2 - cauchyMean)*referenceDistance;
		}
	}

	double tempTwoPowerCnt = pow(2, cntNoOutlier-2);
	heterDelHypothesisLikelihood /= tempTwoPowerCnt;
	heterInsHypothesisLikelihood /= tempTwoPowerCnt;
	if (fabs(heterInsHypothesisMaxLikelihoodSizeChanged) < max(minIndelSize, minIndelRatio * referenceDistance))
        {
                heterInsHypothesisMaxLikelihoodSizeChanged = 0;
                heterInsHypothesisLikelihood = nullHypothesisLikelihood/2;
        }
	if (fabs(heterDelHypothesisMaxLikelihoodSizeChanged) < max(minIndelSize, minIndelRatio * referenceDistance))
        {
                heterDelHypothesisMaxLikelihoodSizeChanged = 0;
                heterDelHypothesisLikelihood = nullHypothesisLikelihood/2;
        }





//	fprintf(outputResult1, "%d %lld %lld %.2E %.2E %.2E %.2E Total: %lld, Ref: %.0lf, OutlierRange: %d, AveScore: %g, AveConf: %g, AveFPR: %g, AveFNR: %g, AveAlignRate: %g, OppoSupp: %lld\n", chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], nullHypothesisLikelihood, homoIndelHypothesisLikelihood, heterDelHypothesisLikelihood, heterInsHypothesisLikelihood, cnt, referenceDistance, outlierRange, averageScore, averageConfidence, averageFPR, averageFNR, averageAlignRate, oppoSupp);
//	for (LL i=0; i<cnt; i++)
//		tempDistanceDouble[i] = tempDistanceDouble[i] * referenceDistance;
	sort(tempDistanceDouble, tempDistanceDouble + cnt);
//	for (LL i=0; i<cnt; i++){
//		fprintf(outputResult1, "%.0lf ", tempDistanceDouble[i]);
//	}
//	fprintf(outputResult1, "\n");


	double hetIHL = heterInsHypothesisLikelihood;
	double hetDHL = heterDelHypothesisLikelihood;
	double nullHL = nullHypothesisLikelihood;
	double tempMaxLikelihood = max(homoIndelHypothesisLikelihood, max(heterDelHypothesisLikelihood, heterInsHypothesisLikelihood));
	int posFlag = 0;
	if (tempMaxLikelihood/likelihoodRatioCutOff > nullHypothesisLikelihood)
	{
		// Update Homo Indel Variant: also need to consider HI1I2, HD1D2, HDI
		if (feq(tempMaxLikelihood, homoIndelHypothesisLikelihood))
		{
			if (homoIndelHypothesisLikelihood/nullHL > hetIHL/nullHL*hetDHL/nullHL)
			{ // if homoHL_ratio > hetIHL_ratio*hetDHL_ratio, we consider homo Indel
				//double tempSizeChange = (sampleMean - cauchyMean) * referenceDistance;
				if (tempSizeChange < 0)
				{
					posFlag = 1;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 1, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag, tempMaxLikelihood, cnt, oppoSupp);
				}
				else //if (tempSizeChange >= max(minIndelSize, minIndelRatio * referenceDistance))
				{
					posFlag = 2;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 0, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag, tempMaxLikelihood, cnt, oppoSupp);
				}
			}
			else if ( ( hetIHL/likelihoodRatioCutOff>nullHL )&&( hetDHL/likelihoodRatioCutOff>nullHL ) )//(hetIHL/nullHL>likelihoodRatioCutOff&&hetDHL/nullHL>likelihoodRatioCutOff)
			{ // if homoHL_ratio > hetIHL_ratio*hetDHL_ratio, we consider the potential HI1I2, HD1D2, HDI
				double disRight = heterInsHypothesisMaxLikelihoodSizeChanged;
				double disLeft = heterDelHypothesisMaxLikelihoodSizeChanged;
				double disBias = disRight - disLeft;//&&(max(hetIHL/hetDHL,hetDHL/hetIHL)<1e10)
				if ( ( cnt>=numberOfSupportIndelMolecule*1.5)&&\
					( fabs(disBias)>=max(minIndelSize,minIndelRatio*max(fabs(disRight),fabs(disLeft))) ) )
					//( hetIHL/likelihoodRatioCutOff/likelihoodRatioCutOff>nullHL )&&( hetDHL/likelihoodRatioCutOff/likelihoodRatioCutOff>nullHL ) )
				{// if all three differences are significant, then call HI1I2, HD1D2, HDI resp., and store them into separate two variants with the same ref region first
					if (disRight > 0)
					{
						posFlag = 3;
					//	printf("Find one case with position flag %d\n",posFlag);
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 0, 0, 1, hetIHL/nullHL, posFlag, hetIHL, cnt, oppoSupp);
					}
					else
					{
						posFlag = 4;
					//	printf("Find one case with position flag %d\n",posFlag);
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 1, 0, 1, hetIHL/nullHL, posFlag, hetIHL, cnt, oppoSupp);
					}
					if (disLeft > 0)
					{
						posFlag = 5;
					//	printf("Find one case with position flag %d\n",posFlag);
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 0, 0, 1, hetDHL/nullHL, posFlag, hetDHL, cnt, oppoSupp);
					}
					else
					{
						posFlag = 6;
					//	printf("Find one case with position flag %d\n",posFlag);
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 1, 0, 1, hetDHL/nullHL, posFlag, hetDHL, cnt, oppoSupp);
					}
				}
				else
				{// else also output homo Indel
					//double tempSizeChange = (sampleMean - cauchyMean) * referenceDistance;
					if (tempSizeChange < 0)
					{
						posFlag = 7;
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 1, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag, homoIndelHypothesisLikelihood, cnt, oppoSupp);
					}
					else //if (tempSizeChange >= max(minIndelSize, minIndelRatio * referenceDistance))
					{
						posFlag = 8;
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 0, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag,homoIndelHypothesisLikelihood, cnt, oppoSupp);
					}
				}
			}
			else
			{
				if (tempSizeChange < 0)
				{
					posFlag = 9;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 1, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag,homoIndelHypothesisLikelihood, cnt, oppoSupp);
				}
				else //if (tempSizeChange >= max(minIndelSize, minIndelRatio * referenceDistance))
				{
					posFlag = 10;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 0, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag,homoIndelHypothesisLikelihood, cnt, oppoSupp);
				}
			}
		}
		// Update Heter InDel Variant, need to modified tomorrow
		else if (feq(tempMaxLikelihood, hetIHL)||feq(tempMaxLikelihood, hetDHL))
		{
			double disRight = heterInsHypothesisMaxLikelihoodSizeChanged;
                        double disLeft = heterDelHypothesisMaxLikelihoodSizeChanged;
                        double disBias = disRight - disLeft;

			if ( (hetDHL/likelihoodRatioCutOff > nullHL) && (hetIHL/likelihoodRatioCutOff > nullHL) )
			{//if two sides are both confident, then we need to check the possible HI1I2, HD1D2, and HDI
				//printf("Actually no such cases?\n");
                                if ( (cnt>=numberOfSupportIndelMolecule*1.5 )&&\
					( fabs(disBias)>=max(minIndelSize,minIndelRatio*max(fabs(disRight),fabs(disLeft))) ) )
					//( min(hetDHL,hetIHL)/likelihoodRatioCutOff/likelihoodRatioCutOff > nullHL ) )
                                {// if all three differences are significant, then call HI1I2, HD1D2, HDI resp., and store them into separate two variants with the sam$
                                        if (disRight > 0)
					{
						posFlag = -1;
				//		printf("Find one case with position flag %d\n",posFlag);
                                                variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 0, 0, 1, hetIHL/nullHL, posFlag,hetIHL, cnt, oppoSupp);
                                        }
					else
					{
						posFlag = -2;
				//		printf("Find one case with position flag %d\n",posFlag);
                                                variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 1, 0, 1, hetIHL/nullHL, posFlag, hetIHL,cnt, oppoSupp);
                                        }
					if (disLeft > 0)
                                        {
						posFlag = -3;
				//		printf("Find one case with position flag %d\n",posFlag);
					        variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 0, 0, 1, hetDHL/nullHL, posFlag, hetDHL, cnt, oppoSupp);
                                        }
					else
                                        {
						posFlag = -4;
				//		printf("Find one case with position flag %d\n",posFlag);
					        variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 1, 0, 1, hetDHL/nullHL, posFlag,hetDHL,cnt, oppoSupp);
                                	}
				}
                                else if ((disRight*disLeft>0)&&(homoIndelHypothesisLikelihood/nullHL > max(max(hetIHL,hetDHL)/min(hetIHL,hetDHL),likelihoodRatioCutOff)))//is it appropriate? No, I need to
                                {// else if homo is more significant when adjust the heter case (reduce significance if both side support the same case), then output homo Indel
                                       //double tempSizeChange = (sampleMean - cauchyMean) * referenceDistance;
                                        if (tempSizeChange < 0)
					{
						posFlag = -5;
                                                variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 1, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag,homoIndelHypothesisLikelihood, cnt, oppoSupp);
                                        }
					else //if (tempSizeChange >= max(minIndelSize, minIndelRatio * referenceDistance))
                                        {
						posFlag = -6;
					        variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], tempSizeChange, 0, 0, 0, 1, 1, homoIndelHypothesisLikelihood/nullHL, posFlag,homoIndelHypothesisLikelihood, cnt, oppoSupp);
					}
                                }
				else
					if (hetIHL>hetDHL)
					{
						posFlag = -7;
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 0, 0, 1, hetIHL/nullHL, posFlag,hetIHL,cnt, oppoSupp);
					}
					else //if ((hetIHL<hetDHL)&&(fabs(disLeft) >= max(minIndelSize, minIndelRatio * referenceDistance)))
					{
						posFlag = -8;
						variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 1, 0, 1, hetDHL/nullHL, posFlag,hetDHL,cnt, oppoSupp);
					}
			}
			else
			{
				if (hetIHL>hetDHL)
				{
					posFlag = -9;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disRight, 0, 0, 0, 0, 1, hetIHL/nullHL, posFlag,hetIHL,cnt, oppoSupp);
                                }
				else //if (fabs(disLeft) >= max(minIndelSize, minIndelRatio * referenceDistance))
				{
					posFlag = -10;
					variant[numberOfVariant++].update(curId, chrId, chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], disLeft, 0, 0, 1, 0, 1, hetDHL/nullHL, posFlag,hetDHL,cnt,oppoSupp);
				}
			}
		}
	}
// if (chromosome.position[distancePair[start].start]==9810098) // { // printf("Let me see, the segment if (%lld, %lld), and the position flag is %d.\n",chromosome.position[distancePair[start].start], chromosome.position[distancePair[start].end], posFlag); // }
}

void completePairs()
{
//	sort(distancePair, distancePair + distancePairCount);
	sort(distancePair.begin(), distancePair.end());
	vector<distanceType> extraDP;
	extraDP.clear();
	LL ost = -1, oed = -1;
	for (LL i = 0; i < distancePairCount; i++)
	{
		if (distancePair[i].start!= ost || distancePair[i].end!=oed)//this locate the non-redundant reference regions (NRRR)
		{
			distanceType tempPair;
			ost = distancePair[i].start;
			oed = distancePair[i].end;
			if (oed-ost == 1)continue;
			for (LL j = 0; j <= i-1; j++)//this search the possible sub-regions
			{
				//LL ted = distancePair[j].end;
				if (distancePair[j].start < ost) continue;
				else if (distancePair[j].start > ost) break;
				else//for each sub-region holding the same start point with NRRR, we explore the consecutive regions on the same molecule
				{
					tempPair = distancePair[j];
					bool rep = true;
					for (LL k = j+1; distancePair[k].start<oed&&k<distancePairCount; k++)
					{
						if (distancePair[k].start == tempPair.end&&strcmp(distancePair[k].mapId, distancePair[j].mapId)==0)
						{
							tempPair.end = distancePair[k].end;
							tempPair.distance += distancePair[k].distance;
						}
						if (tempPair.end >= oed) break;
					}
					if (tempPair.end == oed) //no need
					{
						rep = false;
						for (LL l = i; l < distancePairCount; l++)
						{
							if (distancePair[i] < distancePair[l])break;
							else if (strcmp(distancePair[l].mapId, tempPair.mapId)==0) {rep = true; break;}
						}
					}
					if (!rep) extraDP.push_back(tempPair);
				}
			}
		}
		else
			continue;
	}
	printf("There are %lld extra distance pairs are added.\n",(LL)extraDP.size());
	// also need to include the new distance pairs into the original set
	distancePair.resize(distancePairCount+extraDP.size());
	for (LL i = 0; i < (LL)extraDP.size(); i++){
		distancePair[distancePairCount++]=extraDP[i];
	}
}
void detectByDistancePair(){
	sort(distancePair.begin(), distancePair.end());//	sort(distancePair, distancePair + distancePairCount);
	distancePair.resize(distancePairCount+1);// printDistancePair();
	distancePair[distancePairCount].start = -1;
	distancePair[distancePairCount].end = -1;
	distancePairCount++;
	LL tempHead = 0, tempTail = 0, tempPreNumVar, tempPreNumNoSuppVar;
	for (LL i=0; i<distancePairCount-1; i++){
		if (distancePair[i].start == distancePair[i+1].start && distancePair[i].end == distancePair[i+1].end){
		} else {
			tempTail = i;
			tempPreNumVar = numberOfVariant;
			tempPreNumNoSuppVar = numberOfNoSupportSV;
			advanceLikelihoodDistanceCalculation(tempHead, tempTail);
			/*
							Print distance distribution for 878 only
							if (inputTrio == 878){
							if (numberOfVariant != tempPreNumVar && numberOfNoSupportSV == tempPreNumNoSuppVar)
							for (LL j=tempHead; j<=tempTail; j++)
							distancePair[j].printToFile(outputResult1);
							}
			 */
			tempHead = tempTail + 1;
		}
	}
	distancePairCount--;
	distancePair.resize(distancePairCount);
}
void processDistancePair(){
	sort(distancePair.begin(), distancePair.end());//	sort(distancePair, distancePair + distancePairCount);
	distancePair.resize(distancePairCount+1);//	printDistancePair();
	distancePair[distancePairCount].start = -1;
	distancePair[distancePairCount].end = -1;
	distancePairCount++;
	LL tempCnt = 1, cc = 0;
	double tempDis = distancePair[0].distance;
	distancePair[0].cnt = 1;
	for (LL i=0; i<distancePairCount-1; i++)
		if (distancePair[i].start == distancePair[i+1].start && distancePair[i].end == distancePair[i+1].end){
			tempDis += distancePair[i+1].distance;
			tempCnt++;
		} else {
			distancePair[cc].start = distancePair[i].start;
			distancePair[cc].end = distancePair[i].end;
			distancePair[cc].distance = tempDis/tempCnt;
			distancePair[cc].cnt = tempCnt;
			cc++;
			tempCnt = 1;
			tempDis = (i+1)>=distancePairCount?0:distancePair[i+1].distance;
		}
	distancePairCount = cc;
	distancePair.resize(cc);
//	memset(shortestDist, 0, sizeof(shortestDist));
//	memset(shortestDistCount, 0, sizeof(shortestDistCount));
//	memset(shortestGapDist, 0, sizeof(shortestGapDist));
//	memset(shortestGapDistCount, 0, sizeof(shortestGapDistCount));
//	for (LL i=0; i<distancePairCount; i++){
//		if (distancePair[i].end - distancePair[i].start == 2){
//			shortestDist[distancePair[i].start + 1] = distancePair[i].distance;
//			shortestDistCount[distancePair[i].start + 1] = distancePair[i].cnt;
//		}
//		if (distancePair[i].end - distancePair[i].start == 1){
//			shortestGapDist[distancePair[i].start] = distancePair[i].distance;
//			shortestGapDistCount[distancePair[i].start] = distancePair[i].cnt;
//		}
//	}	//	printDistancePair();
}
void parseHitEnum(){
	LL tempCount = 0;
	distancePair.clear();
	distancePair.resize(5000000);
	for (LL i=0; opticalMap[i].refStart != -1 || opticalMap[i].refEnd != -1; i++){
		if (opticalMap[i].score < confidenceLimit) continue; 
/*#if defined(CHM1) || defined(CHM1_2) || defined(CHM1ToASM) here score is just the old confidence
		if (opticalMap[i].belongs != 11) continue;
		//#elif defined CHM1_2
		// if (opticalMap[i].belongs != 11) continue;
#elif defined(Alden) || defined(Ref) || defined(Angel1) || defined(Angel2) || \
		defined(Angel3) || defined(Angel4) || defined(Pipeline_Real) || defined(Real_Ref) || \
		defined(Real_Assembly_Hg38) || defined(Real_Molecule_Assembly)
		if (inputTrio == 878)
			if (opticalMap[i].belongs != 78) continue;
		if (inputTrio == 891)
			if (opticalMap[i].belongs != 91) continue;
		if (inputTrio == 892)
			if (opticalMap[i].belongs != 92) continue;
#elif defined(Ref_Sim) || \
		defined(Sim_Ref_OMDP) || \
		defined(Sim_Alden_OMDP) || \
		defined(Alden_sim_T0) || \
		defined(Sim_Haploid_Pipeline) || \
		defined(Sim_Diploid_Pipeline) || \
		defined(Sim_Hap_Ref) || \
		defined(Sim_Dip_Ref) || \
		defined(Sim2_Haploid_OMBlast)
		if (opticalMap[i].belongs != 99) continue;
#elif defined(Schizo) || defined(Autism)
		if (opticalMap[i].belongs != 100) continue;
#elif defined(YH)
		if (opticalMap[i].belongs != 36) continue;
#endif	*/
		//numberOfOM878++;
		LL RefIndex = upper_bound(chromosome.position, chromosome.position + chromosome.numberOfSites, opticalMap[i].refStart) - chromosome.position - 1;
		LL OpticalIndex = upper_bound(opticalMap[i].oldPosition, opticalMap[i].oldPosition + opticalMap[i].numberOfSites, opticalMap[i].optStart) - opticalMap[i].oldPosition - 1;
		// OpticalIndex = 0;
		LL BeginRefIndex = RefIndex;
		LL BeginOpticalIndex = OpticalIndex;
		LL EndRefIndex = lower_bound(chromosome.position, chromosome.position + chromosome.numberOfSites, opticalMap[i].refEnd - 1) - chromosome.position;
		LL countM = 0;
		LL countD = 0;
		LL countI = 0;
		char previousChar = 0;
		double tempDistance = 0;
		LL PreviousIndex = 0;
		if (opticalMap[i].hitEnum[0] == 'F'){
			distancePair[distancePairCount].start = BeginRefIndex;
			distancePair[distancePairCount].end = EndRefIndex;
			distancePair[distancePairCount].distance = opticalMap[i].optEnd - opticalMap[i].optStart;
			strcpy(distancePair[distancePairCount].mapId, opticalMap[i].mapId);;
			distancePairCount++;
			continue;
		}
		for (LL j=0; j<opticalMap[i].hitEnum.length(); j++){
			if (opticalMap[i].hitEnum[j] >= '0' && opticalMap[i].hitEnum[j] <= '9'){
				tempCount = tempCount * 10 + (opticalMap[i].hitEnum[j]-'0');
			}
			else {
				if (opticalMap[i].hitEnum[j] == 'M'){
					for (LL k=0; k<tempCount; k++){ //#if defined(Real_Assembly_Hg38) // fprintf(outputAssembly, "%s %lld %d %lld\n", opticalMap[i].mapId, opticalMap[i].oldPosition[OpticalIndex], chrId, chromosome.position[RefIndex]); //#endif
						if (OpticalIndex != BeginOpticalIndex){
							tempDistance += opticalMap[i].position[OpticalIndex];
							distancePair[distancePairCount].start = PreviousIndex;
							distancePair[distancePairCount].end = RefIndex;
							distancePair[distancePairCount].distance = tempDistance;
							strcpy(distancePair[distancePairCount].mapId, opticalMap[i].mapId);
							distancePairCount++;
							tempDistance = 0;
						}
						chromosome.coverage[RefIndex]++;
						chromosome.occurrence[RefIndex]++;
						chromosome.gapCoverage[RefIndex]++;
						previousChar = 'M';
						OpticalIndex++;

						PreviousIndex = RefIndex;
						countM += tempCount;
						RefIndex++;
					}
				} else if (opticalMap[i].hitEnum[j] == 'D'){
					for (LL k=0; k<tempCount; k++){
						chromosome.coverage[RefIndex+k]++;
						chromosome.gapCoverage[RefIndex+k]++;
						previousChar = 'D';
					}
					RefIndex += tempCount;
					countD += tempCount;
				} else if (opticalMap[i].hitEnum[j] == 'I'){
					chromosome.gapCount[RefIndex-1]++;
					double tempDouble = 0;
					for (LL k=0; k<tempCount; k++){
						if (previousChar != 'D'){
							tempDouble += opticalMap[i].position[OpticalIndex];
						}			//is there any problem? The length of this Inserted segment will be padded if followed by 'M'
						tempDistance += opticalMap[i].position[OpticalIndex];
						previousChar = 'I';
						OpticalIndex++;
					}
					countI += tempCount;
				}
				tempCount = 0;
			}
		}
		if (i < 10){			// printf("%lld %s %lld %lld\n", i, opticalMap[i].mapId, RefIndex - BeginRefIndex - countD + countI, opticalMap[i].numberOfSites);
		}
		chromosome.gapCoverage[RefIndex-1]--;
		if (distancePairCount>distancePair.size()-10000)distancePair.resize(distancePair.size()+100000);
	}
	distancePair.resize(distancePairCount);
	opticalMap.clear();
	vector<optAlignType>().swap(opticalMap);
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

void deletionTest(int chrId){
//	if (chrId == 1){
//		fprintf(outputResult2, "ID\tChr\tSite_Index\tPosition\tCoverage\tOccurence\tP\tU\tP-value\tH0\tHa\tHb\tRef_Dis\tMolecule_Dis\tCount\tType\n");
//	}
	for (LL i=0; i<chromosome.numberOfSites; i++){
		// For not considering close site <= resolutionLimit
		if (i != 0 && i != chromosome.numberOfSites - 1 && (chromosome.position[i] - chromosome.position[i-1] <= resolutionLimit || chromosome.position[i+1] - chromosome.position[i] <= resolutionLimit)) continue;
		if (chromosome.coverage[i] < numberOfSupportSignalMolecule) {
//			variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i], 0, 1, 1, 0, 0, 0, 0, 0, 0);
//			variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i], 0, 1, 1, 1, 0, 0, 0, 0, 0);
			numberOfNoSupportSV++;
			continue;
		} 
		// Filter the extremely high depth region
		if (chromosome.coverage[i] > 5000) continue;

		// printf("%lld %lld %lld %lld\n", i, chromosome.position[i], chromosome.coverage[i], chromosome.occurrence[i]);
		numberOfNickingSite++;
		double pValue = 0;
		double likelihoodh0 = 0;
		double likelihoodha = 0;
		double likelihoodhb = 0;
		double likelihoodRatio = 0;
		for (LL j=chromosome.coverage[i] - chromosome.occurrence[i]; j<=chromosome.coverage[i]; j++){
			pValue += C[chromosome.coverage[i]][j]*pow(1-digestionRate, j)*pow(digestionRate, chromosome.coverage[i] - j);
		}
		likelihoodh0 = C[chromosome.coverage[i]][chromosome.occurrence[i]]*pow(1-digestionRate, chromosome.coverage[i] - chromosome.occurrence[i])*pow(digestionRate, chromosome.occurrence[i]);

		double sizeDifference = (i == chromosome.numberOfSites - 1 ? 10000 : chromosome.position[i+1] - chromosome.position[i]);
		double atLeastOneFalseCutRate = 1 - pow(2.718281828, -falseCutRate*sizeDifference);

		atLeastOneFalseCutRate = 0.1;
		likelihoodha = C[chromosome.coverage[i]][chromosome.occurrence[i]]*pow(atLeastOneFalseCutRate, chromosome.occurrence[i])*pow(1 - atLeastOneFalseCutRate, chromosome.coverage[i] - chromosome.occurrence[i]);
		// double likelihoodRatio2 = likelihoodh0/likelihoodha;
		likelihoodRatio = pow((1-digestionRate)/(1-atLeastOneFalseCutRate), chromosome.coverage[i] - chromosome.occurrence[i]) * pow(digestionRate/atLeastOneFalseCutRate, chromosome.occurrence[i]);

		LL x = chromosome.coverage[i] - chromosome.occurrence[i];
		LL k = chromosome.coverage[i];
		double p = digestionRate;
		double u = atLeastOneFalseCutRate;
		for (LL k1 = 0; k1 <= k; k1++){
			double tempLhb = 0;
			LL k2 = k - k1;
			for (LL y=max(0, x - k2); y <= min(k1, x); y++){
				tempLhb += C[k1][y]*C[k2][x-y]*pow(1-p, y)*pow(p, k1-y)*pow(u, y+k2-x)*pow(1-u, x-y);
			}

			likelihoodhb += C[k][k1]*tempLhb;
		}		
		likelihoodhb /= pow(2, k);


//		fprintf(outputResult2, "%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%lld\t%lld\t%lld\t%s\n", curId, chrId, i+1, chromosome.position[i], chromosome.coverage[i], chromosome.occurrence[i], p, u, pValue, likelihoodh0, likelihoodha, likelihoodhb, (i+1<chromosome.numberOfSites?chromosome.position[i+1]:chromosome.position[i]) - (i-1>=0?chromosome.position[i-1]:0), (LL)shortestDist[i], shortestDistCount[i], "Del");


//		chromosome.distanceDifference[i] = (i+1<chromosome.numberOfSites?chromosome.position[i+1]:chromosome.position[i]) - (i-1>=0?chromosome.position[i-1]:0) - shortestDist[i];


		if (pValue <= pValueCutOff && ((likelihoodh0 * likelihoodRatioCutOff <= likelihoodha) || (likelihoodh0 * likelihoodRatioCutOff <= likelihoodhb)))
			chromosome.signi[i] = true;
		else chromosome.signi[i] = false;

		if (chromosome.signi[i]){
			if (likelihoodha > likelihoodhb){
				variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i], 1, 1, 1, 1, 1, pValue, likelihoodha, x, 0);
			} else {
				variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i], 1, 1, 1, 0, 1, pValue, likelihoodhb, x, 0);
			}
			if (numberOfVariant>variant.size()-1000)variant.resize(variant.size()+10000);
		}
		// fprintf(outputResult1, "%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%s\n", curId, chrId, i+1, chromosome.position[i], chromosome.coverage[i], chromosome.occurrence[i], p, u, pValue, likelihoodh0, likelihoodha, likelihoodhb, "Homo_Miss");
	}
}

void insertionTest(int chrId){
	for (LL i=0; i<chromosome.numberOfSites; i++){
		// For not considering close site <= 1750
		if (i != 0 && i != chromosome.numberOfSites - 1 && (chromosome.position[i] - chromosome.position[i-1] <= resolutionLimit || chromosome.position[i+1] - chromosome.position[i] <= resolutionLimit)) continue;
		if (chromosome.gapCoverage[i] < numberOfSupportSignalMolecule) {
//			variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i+1], 0, 1, 0, 0, 0, 0, 0, 0, 0);
//			variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i+1], 0, 1, 0, 1, 0, 0, 0, 0, 0);
			numberOfNoSupportSV++;
			continue;
		} 
		// Filter the extremely high depth region
		if (chromosome.gapCoverage[i] > 5000) continue;

		double pValue = 0;
		double likelihoodh0 = 0;
		double likelihoodha = 0;
		double likelihoodhb = 0;
		double likelihoodRatio = 0;
		double sizeDifference = (i == chromosome.numberOfSites - 1 ? 10000 : chromosome.position[i+1] - chromosome.position[i]);
		double atLeastOneFalseCutRate = 1 - pow(2.718281828, -falseCutRate*sizeDifference);


		/*
			 TEST
			 atLeastOneFalseCutRate = 0.1;
			 TEST
		 */
		for (LL j=chromosome.gapCount[i]; j<=chromosome.gapCoverage[i]; j++){
			pValue += C[chromosome.gapCoverage[i]][j]*pow(atLeastOneFalseCutRate, j)*pow(1-atLeastOneFalseCutRate, chromosome.gapCoverage[i] - j);
		}
		likelihoodh0 = C[chromosome.gapCoverage[i]][chromosome.gapCount[i]]*pow(atLeastOneFalseCutRate, chromosome.gapCount[i])*pow(1-atLeastOneFalseCutRate, chromosome.gapCoverage[i] - chromosome.gapCount[i]);
		likelihoodha = C[chromosome.gapCoverage[i]][chromosome.gapCount[i]]*pow(digestionRate, chromosome.gapCount[i])*pow(1-digestionRate, chromosome.gapCoverage[i] - chromosome.gapCount[i]);
		likelihoodRatio = pow((1-atLeastOneFalseCutRate)/(1-digestionRate), chromosome.gapCoverage[i] - chromosome.gapCount[i]) * pow(atLeastOneFalseCutRate/digestionRate, chromosome.gapCount[i]);

		// printf("%10lE %10lE\n", likelihoodh0 / likelihoodha, likelihoodRatio);

		LL x = chromosome.gapCount[i];
		LL k = chromosome.gapCoverage[i];
		double p = digestionRate;
		double u = atLeastOneFalseCutRate;
		for (LL k1 = 0; k1 <= k; k1++){
			double tempLhb = 0;
			LL k2 = k - k1;
			for (LL y=max(0, x - k2); y <= min(k1, x); y++){
				tempLhb += C[k1][y]*C[k2][x-y]*pow(u, y)*pow(1-u, k1-y)*pow(1-p, y+k2-x)*pow(p, x-y);
			}
			likelihoodhb += C[k][k1]*tempLhb;
		}
		likelihoodhb /= pow(2, k);


//		fprintf(outputResult2, "%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%lld\t%lld\t%lld\t%s\n", curId, chrId, i+1, chromosome.position[i], chromosome.gapCoverage[i], chromosome.gapCount[i], p, u, pValue, likelihoodh0, likelihoodha, likelihoodhb, ((i+1<chromosome.numberOfSites)?chromosome.position[i+1]:chromosome.position[i]) - chromosome.position[i], (LL)shortestGapDist[i], shortestGapDistCount[i], "Ins");

//		chromosome.gapDistanceDifference[i] = ((i+1<chromosome.numberOfSites)?chromosome.position[i+1]:chromosome.position[i]) - chromosome.position[i] - shortestGapDist[i];

		if (chromosome.gapCoverage[i] >= numberOfSupportSignalMolecule && pValue <= pValueCutOff && ((likelihoodh0 * likelihoodRatioCutOff <= likelihoodha) || (likelihoodh0 * likelihoodRatioCutOff <= likelihoodhb)))
			chromosome.gapSigni[i] = true;
		else chromosome.gapSigni[i] = false;

		if (chromosome.gapSigni[i]){
			if (likelihoodha > likelihoodhb){
				variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i+1], 1, 1, 0, 1, 1, pValue, likelihoodha, x, 0);
			} else {
				variant[numberOfVariant++].update(curId, chrId, chromosome.position[i], chromosome.position[i+1], 1, 1, 0, 0, 1, pValue, likelihoodha, x, 0);
			}
			if (numberOfVariant>variant.size()-1000)variant.resize(variant.size()+10000);
		}
		// fprintf(outputResult1, "%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%10lE\t%s\n", curId,chrId, i+1, chromosome.position[i], chromosome.gapCoverage[i], chromosome.gapCount[i], p, u, pValue, likelihoodh0, likelihoodha, likelihoodhb, "Homo_Ins");
	}
}

void initStatData(){
	homoMissingSite = 0;
	homoNonMissingSite = 0;
	homoInsertSite = 0;
	homoNonInsertSite = 0;
	heterMissingSite = 0;
	heterNonMissingSite = 0;
	heterInsertSite = 0;
	heterNonInsertSite = 0;
	numberOfNickingSite = 0;
	numberOfHomoDeletedSNP = 0;
	numberOfHomoMissingFragment = 0;
	numberOfVariant = 0;
	numberOfNoSupportSV = 0;
}
void initData(){
	memset(chromosome.position, 0, sizeof(chromosome.position));
	memset(chromosome.coverage, 0, sizeof(chromosome.coverage));
	memset(chromosome.occurrence, 0, sizeof(chromosome.occurrence));
	memset(chromosome.gapCoverage, 0, sizeof(chromosome.gapCoverage));
	memset(chromosome.gapCount, 0, sizeof(chromosome.gapCount));
	chromosome.length = 0;
	chromosome.numberOfSites = 0;
	//memset(opticalMap, 0, sizeof(opticalMap));
	//memset(distancePair, 0, sizeof(distancePair));
	distancePairCount = 0;
}
void correctOverlapVariant(){
//	delete opticalMap;
//	delete distancePair;
//	variantType *tempVariant = new variantType[500000];
//	variantType *tempRealVariant = new variantType[500000];
	vector<variantType> tempVariant;
	tempVariant.resize(variant.size()+10000);
	vector<variantType> tempRealVariant;
	tempRealVariant.resize(variant.size()+10000);
	LL tempCnt = 0, tempRealCnt = 0;
	// Consider Segmental Only
	for (LL i=0; i<numberOfVariant; i++)
		if (variant[i].isSupported && !variant[i].isSignal)
			tempRealVariant[tempRealCnt++] = variant[i];

	sort(tempRealVariant.begin(), tempRealVariant.begin() + tempRealCnt);//	sort(tempRealVariant, tempRealVariant + tempRealCnt);
//	printf("numberOfRealVariant before filter overlapping: %lld\n", tempRealCnt);


	//here to combine the HI1I2, HD1D2, HDI
	LL tempLen = 0;
	for (LL sRep = 0; sRep < tempRealCnt - 1;)
	{
		if (tempRealVariant[sRep+1].start == tempRealVariant[sRep].start && tempRealVariant[sRep+1].end == tempRealVariant[sRep].end)
		{
			LL chooseX = tempRealVariant[sRep].ratio>=tempRealVariant[sRep+1].ratio?sRep:sRep+1;
			LL chooseY = tempRealVariant[sRep].ratio<tempRealVariant[sRep+1].ratio?sRep:sRep+1;
			variantType tempVS = tempRealVariant[chooseX];// always store the max ratio of the "Both" SVs in the first slot (e.g. "size", "ratio"), and left the smaller one in the second
							              // slot (e.g. "sizeExt", "ratioExt")
			tempVS.ratioExt = tempRealVariant[chooseY].ratio;
			tempVS.sizeExt = tempRealVariant[chooseY].size;
			tempRealVariant[tempLen++] = tempVS;
			sRep++;
		}
		else
			tempRealVariant[tempLen++] = tempRealVariant[sRep];
		sRep++;
	}
	tempRealCnt = tempLen;
	
	//this three rows are used to include the overlapping variants, then we need to delete the following long sentence
/*	numberOfVariant = tempRealCnt;
        for (LL i=0; i<numberOfVariant; i++)
        {
                variant[i] = tempRealVariant[i];
                //tempVariant[tempCnt++] = tempRealVariant[i];
        }
*/



	LL tempStart = 0, tempEnd = 0;
        while (tempStart < tempRealCnt && tempEnd < tempRealCnt){
                LL tempIndex = tempStart;
                while (tempIndex + 1 < tempRealCnt && \
                                ((tempRealVariant[tempIndex+1].start <= tempRealVariant[tempIndex].start && tempRealVariant[tempIndex+1].end >= tempRealVariant[tempIndex].end) || \
                                 (tempRealVariant[tempIndex].start <= tempRealVariant[tempIndex+1].start && tempRealVariant[tempIndex].end >= tempRealVariant[tempIndex+1].end)
                                ))
                        tempIndex++;
                tempEnd = tempIndex;
                LL tempMaxLength = -1, tempMaxLengthIndex = -1;
                LL tempDelCnt = 0, tempInsCnt = 0, tempHomoCnt = 0, tempHeterCnt = 0;

////////////////////////S1:find the maximum length////////////////////////////////
/*
                // Find non-signal variants first
                for (LL i=tempStart; i<=tempEnd; i++){
                        if (!tempRealVariant[i].isSignal && tempRealVariant[i].end - tempRealVariant[i].start > tempMaxLength){
                                tempMaxLength = tempRealVariant[i].end - tempRealVariant[i].start;
                                tempMaxLengthIndex = i;
                        }
                }
                // If can't find, find others
                if (tempMaxLengthIndex!=-1) tempVariant[tempCnt++] = tempRealVariant[tempMaxLengthIndex];
*/
////////////////////////S2:find the cases without completely covering any others///////////////////////////////////////////////
                for (LL i=tempStart; i<=tempEnd; i++)
                {
			if (tempRealVariant[i].isSignal)continue;
                        bool flag = true;
                        for (LL j = tempStart; j <=tempEnd; j++)
                        {
				if (i==j)continue;
                                if (tempRealVariant[i].end >= tempRealVariant[j].end && tempRealVariant[i].start <= tempRealVariant[j].start)
                                {
                                        flag = false;
                                        break;
                                }
                        }
                        if (flag)  tempVariant[tempCnt++] = tempRealVariant[i];
                }	
                //tempCnt++;
                tempStart = tempEnd+1;
        }



	printf("numberOfRealVariant after filter overlapping: %lld\n", tempCnt);
	
	for (LL i=0; i<numberOfVariant; i++)
		if (variant[i].isSupported && variant[i].isSignal)
			tempVariant[tempCnt++] = variant[i];

	numberOfVariant = tempCnt;
	sort(tempVariant.begin(), tempVariant.begin() + numberOfVariant);
	//numberOfVariant = tempRealCnt;
	for (LL i=0; i<numberOfVariant; i++)
		variant[i] = tempVariant[i];

	tempVariant.clear();//	free(tempVariant);
	tempRealVariant.clear();//	free(tempRealVariant);
}

void outputVariantResult(int argv, char* argc[]){
	fprintf(outputVariantResultFile, "#");
	for (int i = 0; i < argv; i++)
		fprintf(outputVariantResultFile,"%s ",argc[i]);
//	fprintf(outputVariantResultFile, "\n#chr\tstart\tstop\tvariant_call_type\tvariant_call_id\tAverage_align_score\tsample_id\tSV_attributes\tzygosity\tlikelihood\tscores\tCoverage\tAverage_align_confidence\n");
	fprintf(outputVariantResultFile, "\n#chr\tstart\tstop\tvariant_call_type\tsize_extra\tAverage_align_score\tsample_id\tSV_attributes\tzygosity\tlikelihood\tscore\tCoverage\tOppoSupp\n");
	LL variantCallId = 1;
	for (LL i=0; i<numberOfVariant; i++){
		char tempH[20];
		char tempIndel[20];
		memset(tempH,0,sizeof(tempH));
		memset(tempH,0,sizeof(tempIndel));
		if (variant[i].size>0&&variant[i].sizeExt>0)
		{
			strcpy(tempH,"HI1I2");
			strcpy(tempIndel,"Insertions");
		}
		else if (variant[i].size<0&&variant[i].sizeExt<0)
		{
			strcpy(tempH,"HD1D2");
			strcpy(tempIndel,"Deletions");
		}
		else if (variant[i].size*variant[i].sizeExt<0)
		{
			strcpy(tempH,"HDI");
			strcpy(tempIndel,"Both");
		}
		else
		{
			strcpy(tempH,(variant[i].isHomo?"Homozygous":"Heterozygous"));
			strcpy(tempIndel,(variant[i].isDel?"Deletion":"Insertion"));
		}
//		if (strcmp(tempH,"Homozygous")!=0&&strcmp(tempH,"Heterozygous")!=0){
			fprintf(outputVariantResultFile, "%lld\t%lld\t%lld\t%s%s\t%lld\t%lf\t%lld\tsizeChange=%lld\t%s\t%g\t%g\t%lld\t%lld\n", variant[i].chr, variant[i].start, variant[i].end, tempIndel, variant[i].isSignal?".site":"", variant[i].sizeExt, 0.0, inputTrio, variant[i].size, tempH, variant[i].likelihood, variant[i].ratio,variant[i].support, variant[i].oppoSupp);
			numberOfSupportedSV++;
//		}
	}
}

int main(int argv, char *argc[]){
	printf("\n");
	setDefault();
	char SVoutputFile[1000];
	memset(SVoutputFile,0,sizeof(SVoutputFile));
	strcpy(SVoutputFile,"Detected_structual_variants");
	memset(outputFolder,0,sizeof(outputFolder));
	char chrMapFile[1000];
	int numOfChr = 24;
	memset(chrMapFile,0,sizeof(chrMapFile));
	strcpy(inputAlignmentFileName,"Alden_1000_c666_1/combined_1000_c1.oma");
	strcpy(outputFileLocation,"Alden_1000_c666_1/");
	strcpy(chrMapFile,"/home/lil/OM/RefMaps_c666_1/GRCh38_25chr.condense.1000.cmap");
		if (argv == 1)
		{
			printf("\nThe valid parameters are described as follows:\n");
                        printf("\t-inputLabel: \n\t\t Default value: %lld. The index of genome in the trio data. No effect if only investigate one genome.\n",inputTrio);
                        printf("\t-outputFolder: \n\t\t The path of the folder to store the output fils. This folder must be exist!\n");
			printf("\t-SVoutputFile: \n\t\t Default value: %s. The file name of SVs (.osv).\n",SVoutputFile);
			printf("\t-chrMapFile: \n\t\t Default value: %s. The file name of the reference map(.cmap).\n",chrMapFile);
			printf("\t-optAlignFile: \n\t\t Default value: %s. The file name of the alignment map file(.oma).\n",inputAlignmentFileName);
			printf("\t-optTempFolder: \n\t\t Default value: %s. The folder to store the processed alignment maps by chromosomes.\n",outputFileLocation);
                        printf("\t-likelihoodRatioCutOff: \n\t\t Default value: %g. The cutoff of the likelihood ratio for SV all hypothesis. The default value changes along with the experiment data.\n",likelihoodRatioCutOff);
                        printf("\t-numberOfSupportIndelMolecule: \n\t\t Default value: %lld. The minimum coverage of a segment being called SVs. The default value changes along with the experiment data.\n",numberOfSupportIndelMolecule);
                        printf("\t-numberOfSupportSignalMolecule: \n\t\t Default value: %lld. The minimum coverage of a segment to call signal variations. The default value changes along with the experiment data.\n",numberOfSupportSignalMolecule);
                        printf("\t-minIndelSize: \n\t\t Default value (b): %g. The minimum length of a segment to call SVs.\n",minIndelSize);
                        printf("\t-minIndelRatio: \n\t\t Default value: %g. The length proportion of a minimum SV could be detected on a segment. E.g. segment = 10000b, then the length of the minimum SV should be larger than 10000*0.05=500b.\n",minIndelRatio);
                        printf("\t-resolutionLimit: \n\t\t Default value: %g. The minimum length of a segment to call signal variations.\n",resolutionLimit);
                        printf("\t-digestionRate: \n\t\t Default value: %g. The digestion rate of labels (signals) measured in the experiment.\n",digestionRate);
                        printf("\t-falseCutRate: \n\t\t Default value: %g. The rate of false cut of a non-label position.\n",falseCutRate);
                        printf("\t-pValueCutOff: \n\t\t Default value: %g. The cutoff of p-value when call signal variations.\n",pValueCutOff);
                        printf("\t-cauchyMean: \n\t\t Default value: %g. The mean value of cauchy distribution of null hypothesis when calling SVs. Reset a new value only if you have good reason.\n",cauchyMean);
                        printf("\t-cauchyScale: \n\t\t Default value: %g. The parameter to calculate cauchy distribution. Reset a new value only if you have good reason.\n",cauchyScale);
                        printf("\t-confidenceLimit: \n\t\t Default value: %g. The lowest alignment confidence for molecules (optical maps) to call SVs or signal variations.\n",confidenceLimit);
			printf("\t-numberOfChromosome: \n\t\t Default value: %d. The first n chromosomes to detect SVs.\n",numOfChr);
                        return -1;
		}
	bool paraWrongFlag = false;
	for (int i = 1; i < argv; i=i+2)
	{
		string temp(argc[i]);
		//printf("What is the parameter?? %s \n",temp.c_str());
		if (temp.compare("-inputLabel")==0)
			inputTrio = atol(argc[i+1]);
		else if (temp.compare("-numberOfChromosome")==0)
                        numOfChr = atoi(argc[i+1]);
		else if (temp.compare("-optAlignFile")==0)
		{
			memset(inputAlignmentFileName,0,sizeof(inputAlignmentFileName));
			strcpy(inputAlignmentFileName, argc[i+1]);
		}
		else if (temp.compare("-optTempFolder")==0)
		{
			memset(outputFileLocation,0,sizeof(outputFileLocation));
			strcpy(outputFileLocation, argc[i+1]);
		}
		else if (temp.compare("-SVoutputFile")==0)
		{
			memset(SVoutputFile,0,sizeof(SVoutputFile));
			strcpy(SVoutputFile,argc[i+1]);
		}
		else if (temp.compare("-chrMapFile")==0)
		{
			memset(chrMapFile,0,sizeof(chrMapFile));
			strcpy(chrMapFile,argc[i+1]);
		}
		else if (temp.compare("-outputFolder")==0)
			strcpy(outputFolder,argc[i+1]);
		else if (temp.compare("-likelihoodRatioCutOff")==0)
		{	likelihoodRatioCutOff = atof(argc[i+1]);
			if (likelihoodRatioCutOff < 0)
			{
				printf("Error! An negative likelihoodRatioCutOff is inputted!\n");
				paraWrongFlag = true;
			}
		}
		else if (temp.compare("-numberOfSupportIndelMolecule")==0)
		{
			numberOfSupportIndelMolecule = atoi(argc[i+1]);
			if (numberOfSupportIndelMolecule < 0)
                        {
                                printf("Error! An negative numberOfSupportIndelMolecule is inputted!\n");
				paraWrongFlag = true;
                        }
			if (numberOfSupportIndelMolecule >= 1000)
                        {
                                printf("Warning! An improper large (>=1000) numberOfSupportIndelMolecule is inputted!\n");
                        }
		}
		else if (temp.compare("-minIndelSize")==0)
		{
			minIndelSize = atof(argc[i+1]);
			if (minIndelSize < 1)
                        {
                                printf("Error! Please make sure the minIndelSize is a positive integer!\n");
				paraWrongFlag = true;
                        }
			if (minIndelSize >= 100000)
                        {
                                printf("Warning! An improper large number (>=100000b) Of minIndelSize is inputted!\n");
                        }
		}
		else if (temp.compare("-minIndelRatio")==0)
		{
			minIndelRatio = atof(argc[i+1]);
			if (minIndelRatio <= 0||minIndelRatio >= 1)
                        {
                                printf("Error! An improper number Of minIndelRatio is inputted, should be in between 0 and 1!\n");
				paraWrongFlag = true;
                        }
		}
		else if (temp.compare("-resolutionLimit")==0)
		{
			resolutionLimit = atof(argc[i+1]);
			if (resolutionLimit < 0)
                        {
                                printf("Error! An negative number Of resolutionLimit is inputted!\n");
				paraWrongFlag = true;
                        }
			if (resolutionLimit >= 10000)
                        {
                                printf("Warning! An improper large number (>=10000) of resolutionLimit is inputted!\n");
                        }
		}
		else if (temp.compare("-digestionRate")==0)
		{
			digestionRate = atof(argc[i+1]);
			if (digestionRate <= 0||digestionRate>1)
                        {
                                printf("Error! An improper digestionRate is inputted, should be in between 0 and 1!\n");
				paraWrongFlag = true;
                        }
			else if (digestionRate < 0.5)
                        {
                                printf("Warning! An very low level of digestionRate (<0.5) is adapted!\n");
                        }
		}
		else if (temp.compare("-falseCutRate")==0)
		{
			falseCutRate = atof(argc[i+1]);
			if (falseCutRate < 0||falseCutRate>=1)
                        {
                                printf("Error! An improper falseCutRate is inputted, should be in between 0 and 1!\n");
				paraWrongFlag = true;
                        }
			else if (falseCutRate > 0.8)
				printf("Warning! An very large (>0.8) falseCutRate is inputted!\n");
		}
		else if (temp.compare("-pValueCutOff")==0)
		{
			pValueCutOff = atof(argc[i+1]);
			if (pValueCutOff <= 0||pValueCutOff>=1)
                        {
                                printf("Error! An improper pValueCutOff is inputted, should be in between 0 and 1!\n");
				paraWrongFlag = true;
                        }
			else if (pValueCutOff > 0.1)
                        {
                                printf("Warning! An non-significant (>0.1) pValueCutOff is inputted!\n");
                        }
		}
		else if (temp.compare("-cauchyMean")==0)
		{
			cauchyMean = atof(argc[i+1]);
			if (cauchyMean < 0.8||cauchyMean > 1.2)
                        {
                                printf("Error! cauchyMean should be around 1!\n");
				paraWrongFlag = true;
                        }
			else
				printf("You have changed the value of cauchyMean to be %g, please make sure it is your intention!\n",cauchyMean);
		}
		else if (temp.compare("-cauchyScale")==0)
		{
			cauchyScale = atof(argc[i+1]);
			if (cauchyScale < 0)
                        {
                                printf("Warning! An negative cauchyScale is inputted!\n");
				paraWrongFlag = true;
                        }
			else
				printf("You have changed the value of cauchyScale to be %g, please make sure it is your intention!\n",cauchyScale);
		}
		else if (temp.compare("-confidenceLimit")==0)
		{
			confidenceLimit = atof(argc[i+1]);
			if (cauchyScale < 0)
                        {
                                printf("Error! An negative confidenceLimit is inputted!\n");
				paraWrongFlag = true;
                        }
		}
		else if (temp.compare("-numberOfSupportSignalMolecule")==0)
		{
			numberOfSupportSignalMolecule = atol(argc[i+1]);
			if (numberOfSupportSignalMolecule < 0)
                        {
                                printf("Error! An negative numberOfSupportSignalMolecule is inputted!\n");
				paraWrongFlag = true;
                        }
			if (numberOfSupportSignalMolecule > 1000)
                        {
                                printf("Warning! An very large (>1000) numberOfSupportSignalMolecule is inputted!\n");
                        }
		}
		else
		{
			printf("No such parameter or wrong : %s\n",argc[i]);
			paraWrongFlag = true;
		}
	}
	if (paraWrongFlag)
	{
		printf("\nThe format to arrange parameters is:\n\t./makeRefine_en -param_1 val_1 -param_2 val_2 ... -param_n val_n\n");
		printf("\nPlease check your input parameters!\n\n");
                return -1;
	}
	
	init();
	NCRcalculation();
	//	findResolutionDistribution(); return 0;
	initResultFullList(SVoutputFile);
	getChromosomeList(chrMapFile);
	initStatData();
	//new content
	bool canOF = true;
	numOfChr = 24;
        for (LL tempChr = 0; tempChr < (LL)listOfChromosome.size(); tempChr++){
                LL chrID = listOfChromosome[tempChr];
                if (chrID > numOfChr) continue;
                char buff[200];
                char nameOF[1000];
                strcpy(nameOF, outputFileLocation);
                sprintf(buff, "%lld_%d",inputTrio,chrID);
                strcat(nameOF, buff);
                strcat(nameOF, ".bmap");
                if ((inputOptAlign = fopen(nameOF, "r")) == NULL)
                {
                        canOF = false;
                        break;
                }
        }
        if (!canOF)
        {
                readSourceFile();
                puts("DONE readSourceFile");
                addSplitedMap();
                puts("DONE addSplitedMap");
                outputToBillMapDestinationFile();
        }
	//readSourceFile();
        //puts("DONE readSourceFile");
        //addSplitedMap();
        //puts("DONE addSplitedMap");
        //outputToBillMapDestinationFile(numOfChr);
	//content end
	for (LL tempChrId = 0; tempChrId < (LL)listOfChromosome.size(); tempChrId++){
		chrId = listOfChromosome[tempChrId];
		//if (chrId > numOfChr) continue;
		// for (chrId = 12; chrId <= 12; chrId++){
		canOpenFile = true;
		initData();
		// printf("Doing %d\n", chrId);
		readChromosomeInfo(chrId,chrMapFile);
		// printf("Done Read Chromosome\n");
		readOpticalAlign(chrId,outputFileLocation);
		if (!canOpenFile)
		{
			printf("Open File fails?\n");
			continue;
		}
		// printf("Done Read Optical\n");
		//printDistancePair("test_output/Distance_before_parsing");
		printf("Distance Pair Count-1: %lld\n", distancePairCount);
		parseHitEnum();
//								printf("Done parsehitEnum\n");
		printf("Distance Pair Count-2: %lld\n", distancePairCount);
		completePairs();
		// printf("pair count: %lld\n", distancePairCount);
		//printDistancePair("test_output/Distance_after_parsing");
		detectByDistancePair();
		//printDistancePair("test_output/Distance_after_detection");
		printf("Distance Pair Count-3: %lld\n", distancePairCount);
//				printf("Done detectByDistancePair\n");

//		processDistancePair();
//				printf("Done processDistancePair\n");
		//printDistancePair("test_output/Distance_after_processing");
//		printf("Distance Pair Count-4: %lld\n", distancePairCount);
//		printf("Number of the sites: %lld\n",chromosome.numberOfSites);
		deletionTest(chrId);
//				printf("Done Deletion\n");
		insertionTest(chrId);
		// printf("Done Insertion\n");
	}
	correctOverlapVariant();
	outputVariantResult(argv,argc);
//	printf("Number of No support SV: %lld\n", numberOfNoSupportSV);
	printf("Number of support SV: %lld\n", numberOfSupportedSV);
	printf("\n");
	return 0;
}
