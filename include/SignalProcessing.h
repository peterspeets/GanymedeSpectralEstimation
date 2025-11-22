#ifndef SIGNALPROCESSING_H
#define SIGNALPROCESSING_H

#include <complex>
#include <cmath>
#include <string>
#include <vector>

#include "Settings.h"
#include "globals.h"
#include "UtilityMathFunctions.h"

/*

This class includes some leftover signal processing functions that are no longer in use. All
functions that are actually used for signal processing is included in the BScan class.

*/

class SignalProcessing {
public:
    SignalProcessing();



    static vector<double> hannWindow(const int N) ;
    static vector<double> tukeyWindow(const int N, const double alpha) ;
    static complex<float>* hilbert(const float* signal, const int N) ;
    static complex<float>* hilbert(const complex<float>* signal, const int N) ;

    static float* calculateEnvelope(float *spectrum,size_t length) ;

    static void equalizeEnvelopeInPlace(float** spectra,size_t N,size_t M) ;

    static void stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower) ;







    static void preprocessSpectrumInPlace(float** spectra, float* offsetSpectrum, float* chirp, float* referenceSpectrum,float* window   ) ;


    virtual ~SignalProcessing();

protected:

private:
};

#endif // SIGNALPROCESSING_H
