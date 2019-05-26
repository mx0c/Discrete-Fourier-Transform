# Signal Processing Library
A Multimodal Signal Processing Library in C i wrote in my 3rd Semester of studying CS. The Testmain uses GNUPlot to display the Signals.

| Method|returns|Description|
|-|-|-|  
|*generateSinus(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude)|double||
|*generateSinusSignal(int iNumSamples, int iNumSamplesPerPeriod, double dAmplitude)|t_mmsignal||
|*initSignal(void)|t_mmsignal||
|deleteSignal (t_mmsignal *pmmsvSignal)|void||
|deleteExtrema (t_Extrema *pteExtrema)|void||
|writeSignal(t_mmsignal *pmmsvSignal, char *pcFilename)|void||
|*readSignal(char *pcFilename)|t_mmsignal||
|getArea(t_mmsignal *pmmsvSignal)|double||
|getMeanvalue(t_mmsignal *pmmsvSignal)|double||
|getVariance(t_mmsignal *pmmsvSignal)|double||
|getMedian(t_mmsignal *pmmsvSignal)|double||
|*getExtrema(t_mmsignal *pmmsvSignal)|t_Extrema||
|*getHistogram(t_mmsignal *pmmsvSignal, int iNumBins)|t_Histogram||
|getEntropy(t_Histogram *ptHistogram)|double||
|deleteHistogram(t_Histogram *ptHistogram)|void||
|*convoluteSignals(t_mmsignal *pmmsvSignalA, t_mmsignal *pmmsvSignalB)|t_mmsignal||
|dft(int iNumSamples, double *pdRealIn, double *pdImgIn, double *pdRealOut, double *pdImgOut, int iDir)|void||
|cart2polar( double *pdRealIn, double *pdImgIn,double *pdAmpOut, double *pdFreqOut)|void||
|polar2cart( double *pdAmpIn, double *pdFreqIn,double *pdRealOut, double *pdImgOut)|void||

# Compiling the Testmain
- $ gcc *.c -o out
