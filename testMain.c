#include "MMSP.h"

int main (int iArgc, char **ppcArgs){

	int iNumSamples = 100;
	int iNumSamplesPerPeriod = 25;
	double dAmplitude = 10;

	t_mmsignal* signal = generateSinusSignal(iNumSamples, iNumSamplesPerPeriod, dAmplitude);
	writeSignal(signal,(char*)".\\data\\signal.sig");

	t_mmsignal* test = readSignal((char*)".\\data\\signal.sig");
	//writeSignal("signal.dat",test);

	printf("\n\nAverage: %lf",getMeanvalue(test));
	printf("\nVarianz: %lf",getVariance(test));
	printf("\nArea: %lf",getArea(test));

	//Test for Extrema
	printf("\n\nGlobalMaxPos = %i",getExtrema(test)->iGlobalMaxPos);
	printf("\nGlobalMinPos = %i",getExtrema(test)->iGlobalMinPos);

	int i = 0;
	//prints all local Minimas
	/*for(i=0;i<getExtrema(test)->ptleMinima->iNumLocalExtrema;i++){
		printf("\n%lf",test->pdValues[getExtrema(test)->ptleMinima->piLocalExtremaPos[i]]);
	}*/

	//Test for histogram
	t_Histogram* h = getHistogram(test,10); //Test with 5 Bins
	for(i=0;i<h->iNumBins;i++){
		printf("\nBinNum:%i  Data:%i",i+1,h->piBinValues[i]);
	}

	//prints median
	printf("\nMedian: %lf",getMedian(test));

	//prints entropy
	printf("\nEntropy: %lf\n\n",getEntropy(h));

	//test for Convolution
	t_mmsignal* kernel = initSignal();
	double tmpA[] = {0,1,2,1,0};
	kernel->pdValues = tmpA;
	kernel->iNumSamples = 5;

	t_mmsignal* kernel2 = initSignal();
	double tmpC[] = {4,3,1,3,4};
	kernel2->pdValues = tmpC;
	kernel2->iNumSamples = 5;

	t_mmsignal* orig = initSignal();
	double tmpB[] = {1,2,3,4,5,3,2,1,3,6};
	orig->pdValues = tmpB;
	orig->iNumSamples = 10;

  	t_mmsignal* faltung1 = convoluteSignals(signal,kernel);
	t_mmsignal* inv_faltung1 = convoluteSignals(kernel,signal);
	t_mmsignal* faltung2 = convoluteSignals(signal,kernel2);

	/*for(i=0;i<faltung->iNumSamples;i++){
		printf("\n%lf",faltung->pdValues[i]);
	}*/

	writeSignal(faltung1,(char*)".\\data\\conv.sig");
  	//writeSignal(orig,"orig.sig");
	writeSignal(inv_faltung1,(char*)".\\data\\inv_conv.sig");
	writeSignal(faltung2,(char*)".\\data\\conv2.sig");

	//Display signal in gnuplot
	//system("gnuplot -persist -e \"plot 'signal.sig' with lines, 'faltung1.sig' with lines, 'faltung2.sig' with lines\"");
	system("gnuplot -persist -e \"plot '.\\data\\signal.sig' with lines, '.\\data\\conv.sig' with lines, '.\\data\\inv_conv.sig' with lines, '.\\data\\conv2.sig' with lines\"");
	return 0;
}
