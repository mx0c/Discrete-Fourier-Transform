#include "MMSP.h"

t_mmsignal* generateSinusSignal(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude){
	double *pdValues = (double*)malloc(sizeof(double)*iNumSamples);
	pdValues = generateSinus(iNumSamples, iNumSamplesPerPeriod, dAmplitude);

	t_mmsignal* signal = (t_mmsignal*)malloc(sizeof(t_mmsignal));
	signal->iNumSamples = iNumSamples;
	signal->pdValues = pdValues;

	return signal;
};

double* generateSinus(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude){
	double *pdValues = (double*)malloc(sizeof(double)*iNumSamples);

	double dSampleStep = (2*M_PI)/iNumSamplesPerPeriod;
	double dStep=0;

	//printf("samplestep: %f\n",dSampleStep);

	int i = 0;
	for(i;i<iNumSamples;i++){
		dStep+=dSampleStep;
		double dVal = dAmplitude*sin(dStep);
		//printf("%f\n",dVal);
		pdValues[i] = dVal;
	}

	return pdValues;
};

t_mmsignal* initSignal(void){
	t_mmsignal* signal = (t_mmsignal*)malloc(sizeof(t_mmsignal));
	return signal;
};

void deleteSignal(t_mmsignal* pmmsvSignal){
	if(pmmsvSignal != NULL){
		if(pmmsvSignal->pdValues != NULL){
			free(pmmsvSignal->pdValues);
		}
		free(pmmsvSignal);
	}
};

//saves Signalstruct (s_MMSignal) in file
void writeSignal(t_mmsignal* pmmsvSignal,char* pcFilename){
	FILE *fp = fopen(pcFilename, "w");

	int i = 0;
	for(i=0;i<pmmsvSignal->iNumSamples;i++){
		fprintf(fp,"%f\n",pmmsvSignal->pdValues[i]);
	}
	fflush(fp);
	fclose(fp);
};

t_mmsignal* readSignal(char* pcFilename){
	FILE *fp = fopen(pcFilename, "r");
	double* tmpValues = (double*)malloc(sizeof(double));
	int tmpNumSamples = 0;

	if(fp != NULL){
		while (fscanf(fp, "%lf",&tmpValues[tmpNumSamples]) != EOF) {
			tmpNumSamples++;
			tmpValues = (double*)realloc(tmpValues,sizeof(double)*(tmpNumSamples+1)); //+1 da tmpNumSamples arrayIndex ist und dieser bei 0 startet
		}
		fclose(fp);
		t_mmsignal* signal = (t_mmsignal*)malloc(sizeof(t_mmsignal));
		signal->iNumSamples = tmpNumSamples;
		signal->pdValues = tmpValues;

		return signal;
	}else{
		printf("\nFehler: Datei nicht gefunden!");
	}
};

double getMeanvalue(t_mmsignal* pmmsvSignal){
	double sum=0;
	int i=0;
	for(i=0;i<pmmsvSignal->iNumSamples;i++){
		sum += pmmsvSignal->pdValues[i];
	}

	sum = sum/pmmsvSignal->iNumSamples;
	return sum;
};

double getArea(t_mmsignal *pmmsvSignal){
	double sum=0;
	int i = 0;

	for(i=0;i<pmmsvSignal->iNumSamples;i++){
		sum+=fabs(pmmsvSignal->pdValues[i]); //fabs() for absolute value of a double
	}

	return sum;
};

double getVariance(t_mmsignal *pmmsvSignal){
	int mv = getMeanvalue(pmmsvSignal);
	double erg=0;

	int i=0;
	for(i=0;i<pmmsvSignal->iNumSamples;i++){
		erg += pow((pmmsvSignal->pdValues[i]-mv),2);
	}

	erg = erg/pmmsvSignal->iNumSamples;
	return erg;
};

//compare function for qsort (in getMedian())
int cmpfx(const void *x, const void *y){
	double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
};

//TODO even odd difference
double getMedian(t_mmsignal *pmmsvSignal){
	double* sortedSignal = pmmsvSignal->pdValues;
	qsort(sortedSignal,pmmsvSignal->iNumSamples,sizeof(double),cmpfx);
	return (sortedSignal[pmmsvSignal->iNumSamples/2]);
};

t_Extrema *getExtrema(t_mmsignal *pmmsvSignal){
	//init
  t_LocalExtrema* minExtrema = (t_LocalExtrema*)malloc(sizeof(t_LocalExtrema));
	t_LocalExtrema* maxExtrema = (t_LocalExtrema*)malloc(sizeof(t_LocalExtrema));
	maxExtrema->piLocalExtremaPos = (int*)malloc(sizeof(int));
	minExtrema->piLocalExtremaPos = (int*)malloc(sizeof(int));

  double tmp=0;
	int i=0,maxAmount=0,minAmount=0;

	for(i=1;i<(pmmsvSignal->iNumSamples-1);i++){
		tmp = pmmsvSignal->pdValues[i];
		if(pmmsvSignal->pdValues[i-1]<=tmp && pmmsvSignal->pdValues[i+1]<=tmp){
			//local max
			maxAmount++;
		  maxExtrema->iNumLocalExtrema = maxAmount;
			maxExtrema->piLocalExtremaPos = (int*)realloc(maxExtrema->piLocalExtremaPos,sizeof(int)*maxAmount);
			maxExtrema->piLocalExtremaPos[maxAmount-1] = i;
		}else if(pmmsvSignal->pdValues[i-1]>=tmp && pmmsvSignal->pdValues[i+1]>=tmp){
			//local min
			minAmount++;
			minExtrema->iNumLocalExtrema = minAmount;
			minExtrema->piLocalExtremaPos = (int*)realloc(minExtrema->piLocalExtremaPos,sizeof(int)*minAmount);
			minExtrema->piLocalExtremaPos[minAmount-1] = i;
		}
	}

	//find globalMax
	double globalMax=0;
	for(i=0;i<(maxAmount-1);i++){
		if(pmmsvSignal->pdValues[maxExtrema->piLocalExtremaPos[i]] > globalMax){
			globalMax = maxExtrema->piLocalExtremaPos[i];
		}
	}

	//find globalMin
	double globalMin=0;
	for(i=0;i<(minAmount-1);i++){
		if(fabs(pmmsvSignal->pdValues[minExtrema->piLocalExtremaPos[i]]) > globalMin){
			globalMin = minExtrema->piLocalExtremaPos[i];
		}
	}

	//create t_Extrema pointer and return it
	t_Extrema* extrema = (t_Extrema*)malloc(sizeof(t_Extrema));
	extrema->ptleMinima = minExtrema;
	extrema->ptleMaxima = maxExtrema;
	extrema->iGlobalMinPos = globalMin;
	extrema->iGlobalMaxPos = globalMax;
	return(extrema);
};

t_Histogram *getHistogram(t_mmsignal *pmmsvSignal, int iNumBins){
	t_Extrema* extrema = getExtrema(pmmsvSignal);
	double binSize = fabs(pmmsvSignal->pdValues[extrema->iGlobalMaxPos]) + fabs(pmmsvSignal->pdValues[extrema->iGlobalMinPos]);
  binSize = binSize/iNumBins;
	double minValue = pmmsvSignal->pdValues[extrema->iGlobalMinPos];
	printf("\n\nbinSize = %lf minValue: %lf\n",binSize,minValue);

	//Init of the Bins Array
	int* bins = (int*)malloc(iNumBins*sizeof(int));
	int i = 0;
	for(i=0;i<iNumBins;i++){
		bins[i] = 0;
	}

	//filling the Bins
	int sampleNum = 0;
	int binNum = 1;
	for(binNum=1;binNum<=iNumBins;binNum++){
		for(sampleNum=0;sampleNum<pmmsvSignal->iNumSamples;sampleNum++){
			if(binNum!=iNumBins){
				if((pmmsvSignal->pdValues[sampleNum] >= minValue+binSize*(binNum-1)) && (pmmsvSignal->pdValues[sampleNum] < minValue+binSize*binNum)){
					bins[binNum-1] += 1;
				}
			}else{
				if((pmmsvSignal->pdValues[sampleNum] >= minValue+binSize*(binNum-1)) && (pmmsvSignal->pdValues[sampleNum] <= minValue+binSize*binNum)){
					bins[binNum-1] += 1;
				}
			}
		}
	}

	//Debug Ausgabe
	/*for(i=1;i<=iNumBins;i++){
		printf("Bin: %i Start: %lf Ende: %lf Value: %i\n",i,minValue+binSize*(i-1),minValue+binSize*i,bins[i-1]);
	}*/

	//create and construct t_Histogram pointer and return it
	t_Histogram* hist = (t_Histogram*)malloc(sizeof(t_Histogram));
	hist->iNumBins = iNumBins;
	hist->piBinValues = bins;
	hist->dMinValue = minValue;
	hist->dMaxValue = pmmsvSignal->pdValues[extrema->iGlobalMaxPos];
	hist->dEntropy = getEntropy(hist);
	return(hist);
};

double getEntropy(t_Histogram *ptHistogram){
	double probabilitys[ptHistogram->iNumBins];
	int amount=0;
	int i=0;

	//get total amount of samples
	for(i=0;i<ptHistogram->iNumBins;i++){
		amount += ptHistogram->piBinValues[i];
	}

	//calculating the probabilitys
	for(i=0;i<ptHistogram->iNumBins;i++){
		probabilitys[i] = ptHistogram->piBinValues[i]/(double)amount;
		printf("\nprobability:%lf  TotalAmount:%i   BinValue:%i",probabilitys[i],amount,ptHistogram->piBinValues[i]);
	}

	//calculating entropy
	double entropy=0;
	for(i=0;i<ptHistogram->iNumBins;i++){
		entropy += probabilitys[i]*log2(probabilitys[i]);
	}
	//printf("\n%lf",entropy);
	entropy = -entropy;
	return entropy;
};

void deleteHistogram(t_Histogram *ptHistogram){
	if(ptHistogram != NULL){
		if(ptHistogram->piBinValues != NULL){
			free(ptHistogram->piBinValues);
		}
		free(ptHistogram);
	}
};

t_mmsignal *convoluteSignals(t_mmsignal *pmmsvSignalA, t_mmsignal *pmmsvSignalB){
	double* resArray = (double*)malloc(sizeof(double));
	t_mmsignal* greaterSignal = pmmsvSignalA->iNumSamples<pmmsvSignalB->iNumSamples?pmmsvSignalB:pmmsvSignalA;
	t_mmsignal* lesserSignal = pmmsvSignalA->iNumSamples<=pmmsvSignalB->iNumSamples?pmmsvSignalA:pmmsvSignalB;

	int n;
  for (n = 0; n < greaterSignal->iNumSamples + lesserSignal->iNumSamples - 1; n++){
    int min,max,j;

		//dynamically reallocating and initializing the result array
		resArray = (double*)realloc(resArray,sizeof(double)*(n+1));
    resArray[n] = 0;

		//for smooth "fade-in" and "fade-out" purposes
    min = (n >= lesserSignal->iNumSamples - 1) ? n - (lesserSignal->iNumSamples - 1) : 0;
    max = (n < greaterSignal->iNumSamples - 1) ? n : greaterSignal->iNumSamples - 1;

    for (j = min; j <= max; j++){
			printf("xIndex%i yIndex%i\n",j,n-j);
      resArray[n] += greaterSignal->pdValues[j] * lesserSignal->pdValues[n - j];
    }

		//scaling the result
		resArray[n] /= lesserSignal->iNumSamples;
  }

	t_mmsignal* re = initSignal();
	re->pdValues = resArray;
	re->iNumSamples = n+1;

	return(re);
};

//generate Gaussglocke
t_mmsignal *generateBinomial(int iPascalLine){

};

void dft(int iNumSamples, double *pdRealIn, double *pdImgIn, double *pdRealOut, double *pdImgOut, int iDir){
	//http://www.dspguide.com/ch8/6.htm
	double reSum,imSum = 0;
	int i = 0;
	int n = 0;

	//Forward (Converting to Frequency Domain)
	if(iDir == 1){
		for(i=0;i<iNumSamples;i++){
			reSum = 0;
			imSum = 0;
			for(n=0;n<iNumSamples;n++){
				double angle = 2.0 * M_PI * n * i / iNumSamples; //exponent von e
				reSum +=  pdRealIn[n] * cos(angle) + pdImgIn[n] * sin(angle);
				imSum += -pdRealIn[n] * sin(angle) + pdImgIn[n] * cos(angle);
			}
			pdRealOut[i] = reSum;
			pdImgOut[i] = imSum;
		}
	//Backward (Converting to Time Domain)
	}else if(iDir == -1){
		for(i=0;i<iNumSamples;i++){
			reSum = 0;
			imSum = 0;
			for(n=0;n<iNumSamples;n++){
				double angle = -(2.0 * M_PI * n * i / iNumSamples);
				reSum +=  pdRealIn[n] * cos(angle) + pdImgIn[n] * sin(angle);
				imSum += -pdRealIn[n] * sin(angle) + pdImgIn[n] * cos(angle);
			}
			pdRealOut[i] = reSum/iNumSamples;
			pdImgOut[i] = imSum/iNumSamples;
		}
	}
};

//http://www.dspguide.com/ch8/8.htm
void cart2polar(double *pdRealIn, double *pdImgIn,double *pdAmpOut, double *pdFreqOut){
	//iterate through pdRealIn[]
	int counter = 0;
	while(*pdRealIn){
		pdAmpOut[counter] = sqrt(pow(pdRealIn[counter],2)+pow(pdImgIn[counter],2));
		pdFreqOut[counter] = atan(pdImgIn[counter]/pdRealIn[counter]);
		counter++;
		++pdRealIn;
	}
};

void polar2cart( double *pdAmpIn, double *pdFreqIn,double *pdRealOut, double *pdImgOut){
	//iterate through pdRealIn[]
	int counter = 0;
	while(*pdAmpIn){
		pdRealOut[counter] = pdAmpIn[counter]*cos(pdFreqIn[counter]);
		pdImgOut[counter] = pdAmpIn[counter]*sin(pdFreqIn[counter]);
		counter++;
		++pdAmpIn;
	}
};
