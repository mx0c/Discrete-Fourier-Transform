#ifndef _MMSP_H_
#define _MMSP_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>

typedef struct s_Histogram{
	int iNumBins; //Number of intervalls values are distributed to
	int *piBinValues;
	double dMinValue;
	double dMaxValue;
	double dEntropy;
} t_Histogram;

typedef struct s_LocalExtrema{
	int iNumLocalExtrema;
	int *piLocalExtremaPos;
} t_LocalExtrema;

typedef struct s_Extrema{
	int iGlobalMaxPos;
	int iGlobalMinPos;
	t_LocalExtrema *ptleMaxima;
	t_LocalExtrema *ptleMinima;
} t_Extrema;

typedef struct s_MMSignal{
	int iNumSamples;
	double *pdValues;
	double dMean;
	double dVariance;
	double dArea;
	double dMedian;
	t_Extrema *pteExtrema;
	t_Histogram *ptHistogram;
} t_mmsignal;

double *generateSinus(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude);
t_mmsignal *generateSinusSignal(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude);
t_mmsignal *initSignal(void);
void deleteSignal (t_mmsignal *pmmsvSignal);
void deleteExtrema (t_Extrema *pteExtrema);
void writeSignal(t_mmsignal *pmmsvSignal, char *pcFilename);
t_mmsignal *readSignal(char *pcFilename);
double getArea(t_mmsignal *pmmsvSignal);
double getMeanvalue(t_mmsignal *pmmsvSignal);
double getVariance(t_mmsignal *pmmsvSignal);
double getMedian(t_mmsignal *pmmsvSignal);
t_Extrema *getExtrema(t_mmsignal *pmmsvSignal);
t_Histogram *getHistogram(t_mmsignal *pmmsvSignal, int iNumBins);
double getEntropy(t_Histogram *ptHistogram);
void deleteHistogram(t_Histogram *ptHistogram);
t_mmsignal *convoluteSignals(t_mmsignal *pmmsvSignalA, t_mmsignal *pmmsvSignalB);
void dft(int iNumSamples, double *pdRealIn, double *pdImgIn, double *pdRealOut, double *pdImgOut, int iDir);
void cart2polar( double *pdRealIn, double *pdImgIn,double *pdAmpOut, double *pdFreqOut);
void polar2cart( double *pdAmpIn, double *pdFreqIn,double *pdRealOut, double *pdImgOut);
#endif
