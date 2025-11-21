#ifndef BSCAN_H
#define BSCAN_H


#include <cmath>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <thread>
#include <sstream>
#include <vector>
#include <optional>

#include <kiss_fft.h>


#include "Settings.h"
#include "UtilityMathFunctions.h"
#include "IO.h"

using namespace std;

class BScan {
public:
    BScan();
    BScan(const string filePath);
    Settings BScanSettings;
    void preprocessSpectrumInPlace();
    vector<vector<float>> fftBScan();
    virtual ~BScan();
    uint64_t getTime();
    static void stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower); //old function

    void stretchSpectraInPlace(float minimumReferencePower);
    static void FIAALoop(vector<vector<float>>& spectra,BScan* scan,int fromIndex, int toIndex,
                              size_t N, int K, int numberOfPartitions,int numberOfIterations, double vt, vector<float>& startingColumn,vector<vector<float>>& processedImage );

    vector<vector<float>> processBScan(size_t M,const size_t N, int K,int numberOfIterationsFirstColumn,int numberOfIterations, double vt,int NThreads) ;
    vector<vector<float>> processBScan();
    void FIAAPartitioned(const vector<float>& x, vector<float>& powerSpectrum,int numberOfIterations = settings->numberOfIterations );
    void FIAAPartitioned(const vector<float>& x,
                              size_t N, int K, int numberOfPartitions,int numberOfIterations, double vt, vector<float>& powerSpectrum ) ;
    pair<vector<float>,vector<float>> FIAA(const vector<float>& x, int K, int numberOfIterations, double vt, vector<float>& powerSpectrum, int powerSpectrumIndex);
    pair<vector<float>,vector<float>> FIAA(const vector<float>& x, int K, int numberOfIterations, double vt);
    tuple<vector<vector<float>>,int,int> getProcessedBScan();
    void calculateLowResBitmap();
    //float** imageRIAA = nullptr;
    //float** imageFFT = nullptr;

    vector<vector<float>> imageRIAA;
    vector<vector<float>> imageFFT;

    unsigned char* lowResBitmap = nullptr; //no alpha channel



protected:

private:
    vector<vector<float>> spectra;
    vector<float> offset;
    vector<float> chirp;
    vector<float> referenceSpectrum;
    vector<float> intensity;
    vector<float> window;

    static void fftPartBScan(float** spectra, float** image, int Nz,const int startXIndex, const int stopXIndex);
    static void fftPartBScan(vector<vector<float>>& spectra, vector<vector<float>>& image, const int startXIndex, const int stopXIndex);


};

#endif // BSCAN_H
