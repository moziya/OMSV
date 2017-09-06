#include "NoError.h"

FILE* inputBed;
FILE* inputFragileSite;
FILE* inputPseudo;
FILE* inputNRegion;
FILE* outputF;


struct stdSV
{
	char chr[10];
	LL start;
	LL end;
};
vector<stdSV> SVbed;
vector<stdSV> fragileSite;
vector<stdSV> Nregion;
vector<stdSV> pseudoAuto;

struct outRow
{
	char chr[10];
	LL start;
	LL end;
	LL score = 0;
};
vector<outRow> resList;

void swapIfRev(LL &num1, LL &num2)
{
	LL temp;
	if (num1 > num2)
	{
		printf("Swap two values: (%lld, %lld).\n",num1,num2);
		temp = num2;
		num2 = num1;
		num1 = temp;
	}
}

int relation(LL start1, LL stop1, LL start2, LL stop2, int shift)
{//1 is the SV bed, 2 is the FS bed
	swapIfRev(start1,stop1);
	swapIfRev(start2,stop2);
	int typSame = 0;
//	if (strcmp(type1,type2)==0) typSame += 20;
//	if (strcmp(zygosity1,zygosity2)==0) typSame += 10;
	if (start2 > stop1+shift || start1 > stop2+shift) return 0; // no relation
	shift = 0;
	if (abs(start1-start2)<=shift && abs(stop1-stop2)<=shift) return 4+typSame; // complete
	if ((start2 <= start1-shift && stop2 > stop1+shift)||(start2 < start1-shift && stop2 >= stop1+shift)) return 3+typSame; // covered
	if ((start1 <= start2-shift && stop1 > stop2+shift)||(start1 < start2-shift && stop1 >= stop2+shift)) return 2+typSame; // covering
	return 1+typSame; // overlap
}


void readFS(char* inputFS)
{
	fragileSite.clear();
	int numberOfDL = 0;
	char temp[10000];
	while(fgetc(inputFragileSite)=='#')
	{
		numberOfDL++;
		fgets(temp,10000,inputFragileSite);
	}
	fclose(inputFragileSite);
	if ((inputFragileSite = fopen(inputFS,"r"))==NULL) perror("Wrong in reading standard SVs!");
	for (int i = 0; i < numberOfDL; i++)
	{
		fgets(temp,10000,inputFragileSite);
	}
	stdSV tempSV;
	while(fscanf(inputFragileSite,"%s\t%lld\t%lld\n",tempSV.chr, &tempSV.start, &tempSV.end)==3)
	{
		if (strcmp(tempSV.chr,"X")==0)strcpy(tempSV.chr,"23");
                else if (strcmp(tempSV.chr,"Y")==0)strcpy(tempSV.chr,"24");
		fragileSite.push_back(tempSV);
	}
	fclose(inputFragileSite);
	printf("There are %lld fragile site regions.\n", (LL)fragileSite.size());
}

void readNR(char* inputNR)
{
        Nregion.clear();
        int numberOfDL = 0;
        char temp[10000];
        while(fgetc(inputNRegion)=='#')
        {
                numberOfDL++;
                fgets(temp,10000,inputNRegion);
        }
        fclose(inputNRegion);
        if ((inputNRegion = fopen(inputNR,"r"))==NULL) perror("Wrong in reading standard SVs!");
        for (int i = 0; i < numberOfDL; i++)
        {
                fgets(temp,10000,inputNRegion);
        }
        stdSV tempSV;
        while(fscanf(inputNRegion,"%s\t%lld\t%lld\n",tempSV.chr, &tempSV.start, &tempSV.end)==3)
        {
		if (strcmp(tempSV.chr,"X")==0)strcpy(tempSV.chr,"23");
                else if (strcmp(tempSV.chr,"Y")==0)strcpy(tempSV.chr,"24");
                Nregion.push_back(tempSV);
        }
        fclose(inputNRegion);
	printf("There are %lld N-regions.\n", (LL)Nregion.size());
}

void readPAR(char* inputPAR)
{
        pseudoAuto.clear();
        int numberOfDL = 0;
        char temp[10000];
        while(fgetc(inputPseudo)=='#')
        {
                numberOfDL++;
                fgets(temp,10000,inputPseudo);
        }
        fclose(inputPseudo);
        if ((inputPseudo = fopen(inputPAR,"r"))==NULL) perror("Wrong in reading standard SVs!");
        for (int i = 0; i < numberOfDL; i++)
        {
                fgets(temp,10000,inputPseudo);
        }
        stdSV tempSV;
        while(fscanf(inputPseudo,"%s\t%lld\t%lld\n",tempSV.chr, &tempSV.start, &tempSV.end)==3)
        {
                pseudoAuto.push_back(tempSV);
        }
        fclose(inputPseudo);
	printf("There are %lld PARs.\n", (LL)pseudoAuto.size());
}


void readSVBed(char* inputSV)
{
        SVbed.clear();
        int numberOfDL = 0;
        char temp[100000];
        while(fgetc(inputBed)=='#')
        {
                numberOfDL++;
                fgets(temp,100000,inputBed);
        }
        fclose(inputBed);
        if ((inputBed = fopen(inputSV,"r"))==NULL) perror("Wrong in reading standard SVs!");
        for (int i = 0; i < numberOfDL; i++)
        {
                fgets(temp,100000,inputBed);
        }
        stdSV tempSV;
        while(fscanf(inputBed,"%s\t%lld\t%lld\n",tempSV.chr, &tempSV.start, &tempSV.end)==3)
        {
                SVbed.push_back(tempSV);
		fgets(temp,100000,inputBed);
        }
        fclose(inputBed);
}

void crossValidate(int shift)
{
	resList.clear();
	for (LL i = 0; i < (LL)SVbed.size(); i++)
	{
		outRow tempOR;
		strcpy(tempOR.chr,SVbed[i].chr);
		tempOR.start = SVbed[i].start;
		tempOR.end = SVbed[i].end;
		for (LL j = 0; j < (LL)fragileSite.size(); j++)
		{
			if (strcmp(SVbed[i].chr,fragileSite[j].chr)!=0)continue;
			if (relation(SVbed[i].start, SVbed[i].end, fragileSite[j].start, fragileSite[j].end, shift)>0)tempOR.score = 1;
		}
		for (LL j = 0; j < (LL)Nregion.size(); j++)
                {
                        if (strcmp(SVbed[i].chr,Nregion[j].chr)!=0)continue;
                        if (relation(SVbed[i].start, SVbed[i].end, Nregion[j].start, Nregion[j].end, shift)>0)tempOR.score = 2;
                }
		for (LL j = 0; j < (LL)pseudoAuto.size(); j++)
                {
                        if (strcmp(SVbed[i].chr,pseudoAuto[j].chr)!=0)continue;
                        if (relation(SVbed[i].start, SVbed[i].end, pseudoAuto[j].start, pseudoAuto[j].end, shift)>0)tempOR.score = 2;
                }
		resList.push_back(tempOR);
	}
}

void printFile()
{
	fprintf(outputF,"#chr\tstart\tend\tValidation_Result\n");
	for (LL i = 0; i < (LL)resList.size(); i++)
	{
		fprintf(outputF,"%s\t%lld\t%lld\t%lld\n",resList[i].chr, resList[i].start, resList[i].end,resList[i].score);
	}
}



int main(int argc, char* argv[])
{
        char inputSV[500];
	char inputFS[500];
	char inputNR[500];
	char inputPAR[500];
        char outFile[500];
	int shift = -1;
        memset(inputSV,0,sizeof(inputSV));
        memset(outFile,0,sizeof(outFile));
	strcpy(inputFS,"fragile_site_list.bed");
        strcpy(inputNR,"hg38_Nregion.bed");
        strcpy(inputPAR,"pseudo_autosomal.bed");
                if (argc==1)//||temp.compare("--help")==0||temp.compare("--Help")==0||temp.compare("--HELP")==0)
                {
                        printf("\nThe valid parameters are described as follows:\n");
                        printf("\t-inputSVFile: \n\t\t The input file name of the detected SVs.\n");
			printf("\t-inputFSFile: \n\t\t The input fragial site list on hg38.\n");
			printf("\t-inputNRFile: \n\t\t The input N-regions on hg38.\n");
			printf("\t-inputPARFile: \n\t\t The input pseudo-autosomal regions on hg38.\n");
                        printf("\t-outputFile: \n\t\t The output file name.\n");
                        printf("\t-shift: \n\t\t The tolerated position error of SVs. The default value is: %lf.\n",shift);
                        return -1;
                }
        //}
        for (int i = 1; i < argc; i=i+2)
        {
                string temp(argv[i]);
                if (temp.compare("-inputSVFile")==0)
                        strcpy(inputSV,argv[i+1]);
		else if (temp.compare("-inputFSFile")==0)
                        strcpy(inputFS,argv[i+1]);
		else if (temp.compare("-inputNRFile")==0)
                        strcpy(inputNR,argv[i+1]);
		else if (temp.compare("-inputPARFile")==0)
                        strcpy(inputPAR,argv[i+1]);
                else if (temp.compare("-outputFile")==0)
                        strcpy(outFile,argv[i+1]);
                else if (temp.compare("-shift")==0)
                        shift = atoi(argv[i+1]);
                else
                {
                        printf("No such parameter or wrong : %s\n",argv[i]);
                        printf("\nThe format to arrange parameters is:\n\t./convertFromCaoToStandard -param_1 val_1 -param_2 val_2 ... -param_n val_n\n");
                        return -1;
                }
        }
        int ccnt = 0;
        int dotPos = -1;
        while(outFile[ccnt]!='\0')
        {
                if (outFile[ccnt]=='.') dotPos = ccnt;
                ccnt++;
        }
	char outFiles[500];
        memset(outFiles,0,sizeof(outFiles));
        if (dotPos == -1) dotPos = ccnt;
        for (int i = 0; i < dotPos; i++)
        {
                outFiles[i] = outFile[i];
        }
        strcat(outFiles,".bed");
        if (((inputBed = fopen(inputSV,"r"))==NULL)||((inputFragileSite = fopen(inputFS,"r"))==NULL)||((inputPseudo = fopen(inputPAR,"r"))==NULL)||((outputF = fopen(outFiles,"w+"))==NULL)||((inputNRegion = fopen(inputNR,"r"))==NULL))
                perror("Failed to open files!");
	fprintf(outputF,"#");
	for (int par = 0; par < argc; par++)
		fprintf(outputF,"%s  ",argv[par]);
	fprintf(outputF,"\n");
	readSVBed(inputSV);
	printf("Done 1.\n");
	readFS(inputFS);
	readNR(inputNR);
	readPAR(inputPAR);
	printf("Done 2.\n");
	crossValidate(shift);
	printf("Done 3.\n");
	printFile();
        return 0;
}

