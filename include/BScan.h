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

#include <kiss_fft.h>


#include "Settings.h"
#include "UtilityMathFunctions.h"
#include "IO.h"

using namespace std;

class BScan
{
    public:
        BScan(const string filePath);
        Settings BScanSettings;
        void preprocessSpectrumInPlace();
        float** fftBScan();
        virtual ~BScan();
        uint64_t getTime();
        static void stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower);
        static void fiaa_oct_loop(float** spectra,BScan* scan,int fromIndex, int toIndex,
        size_t N, int K, int numberOfPartitions,int q_i, double vt, float* startingColumn,float** processedImage );

        float** processBScan(size_t M,const size_t N, int K,int q_init,int q_i, double vt,int NThreads) ;
        float** processBScan();
        void fiaa_oct_partitioned(const float* x, float* diaaf_floatingPoint,int numberOfIterations = settings->numberOfIterations );
        void fiaa_oct_partitioned(const float* x,
            size_t N, int K, int numberOfPartitions,int q_i, double vt, float* diaaf_floatingPoint ) ;
        pair<float*, float*> fiaa_oct(const float* x, size_t N, int K, int q_i, double vt, float* diaaf_floatingPoint = nullptr);
        tuple<float**,int,int> getProcessedBScan();
        void calculateLowResBitmap();


    protected:

    private:
        float** spectra = nullptr;
        float** processedBScan = nullptr;
        float** imageFFT = nullptr;
        unsigned char* lowResBitmap = nullptr; //no alpha channel
        float* offset = nullptr;
        float* chirp = nullptr;
        float* referenceSpectrum = nullptr;
        float* intensity = nullptr;
        float* window = nullptr;

};

#endif // BSCAN_H
