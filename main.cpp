#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <list>
#include <map>
#include <complex>
#include <memory>
#include <utility>
#include <cmath>
#include <vector>
#include <kiss_fft.h>
#include <kiss_fftr.h>
#include <png.h>


#include <cstdio>
#include <cstdlib>

#include <miniz.h>
#include "Tag.h"
#include "XML_Interpreter.h"
#include "Settings.h"
#include "UtilityMathFunctions.h"
#include "IO.h"
#include "BScan.h"
#include "Window.h"
#include "globals.h"

using namespace std;

shared_ptr<Settings> settings = nullptr;
shared_ptr<BScan> scan = nullptr;


class SignalProcessing {

public:

    static vector<double> hannWindow(const int N) {
        /*
        N: window length

        This function calculates a Hann window and returns the window as a vector of size N.
        */

        vector<double> window(N);
        for(int i = 0; i< N; i++) {
            window[i] = sin(M_PI*i/N)*sin(M_PI*i/N);
        }
        return window;
    }

    static vector<double> tukeyWindow(const int N, const double alpha) {
        /*
        N: length of the vector to return
        alpha: the tapering parameter of the Tukey window, a lower value means a
        steeper edge, and a higher value means a very smooth window.
        */

        vector<double> window(N);
        for(int i = 0; i< N/2+1; i++) {
            if(i < alpha*N/2 ) {
                window[i] = 0.5*(1.0 - cos(2*M_PI*i/(alpha * (N-1))));
            } else {
                window[i] = 1;
            }

            window[N-i-1] = window[i];
        }
        return window;
    }

    template<typename T>
    static complex<T>* hilbert(const T* signal, const int N ) {
        /*
        signal: real valued signal to calculate the Hilbert transform for
        N: length of the signal.

        Hilbert function for a real function. This function calculates the Hilbert transform
        by passing it to the complex Hilbert transform.
        */

        complex<T> temp[N] ;
        complex<T>* output;
        for(int i=0; i<N; i++) {

            temp[i] = complex<T>(signal[i],0.0);
        }


        output = hilbert(temp,N);

        return output;
    }

    template<typename T>
    static complex<T>* hilbert(const complex<T>* signal, const int N) {
        /*
        signal: pointer to an array of complex numbers
        N: length of the array.

        This function calculates the Hilbert transform of a complex signal.
        The Hilbert transforms is calculated by Fourier transforming the
        signal, rotating it 90 degrees in the complex plane, and inverting
        the sign of the negative frequencies.
        */

        kiss_fft_cfg cfg = kiss_fft_alloc(N, 0, NULL, NULL);

        kiss_fft_cpx ft1[N];
        kiss_fft_cpx ft2[N];
        complex<T>* output = new complex<T>[N];

        for(int i = 0; i < N; i++) {
            ft1[i].r = signal[i].real();
            ft1[i].i = signal[i].imag();
        }

        kiss_fft(cfg,ft1,ft2);
        kiss_fft_free(cfg);


        T temp;


        for(int i = 0; i < N/2; i++) {
            temp = ft2[i].r;
            ft2[i].r = ft2[i].i;
            ft2[i].i = -temp;
        }


        for(int i = N/2; i < N; i++) {
            temp = ft2[i].r;
            ft2[i].r = -ft2[i].i;
            ft2[i].i = temp;
        }



        cfg = kiss_fft_alloc(N, 1, NULL, NULL);
        kiss_fft(cfg,ft2,ft1);

        for(int i = 0; i < N; i++) {
            output[i] = complex<T>(signal[i].real() - ft1[i].i/N,signal[i].imag() + ft1[i].r /N) ;

        }
        kiss_fft_free(cfg);
        return output;
    }


    static float* calculateEnvelope(float *spectrum,size_t length) {
        /*
        spectrum: pointer to array of the interference spectrum to calculate
        the envelope for.
        length: the size of the array spectrum.

        This function calculates the envelope of a function based on a spline
        interpolation of the peaks in the spectrum. This is good for finding
        a visually appealing envelope, but for many signal processing applications
        use the absolute value of the Hilbert transform instead.
        */
        float* env = new float[length];
        int NChunks = settings->NChunksEnvelopeSubtraction;
        int chunkLength = length/NChunks;
        float maximumValuesPerChunk[NChunks+2];
        float indexOfMaximumPerChunk[NChunks+2];
        int index;
        int i, j;

        for(i = 0; i< length; i++) {
            env[i] =  abs(spectrum[i]);
        }

        for(i = 0; i< NChunks; i++) {
            int maxIndex = (i+1)*chunkLength;
            if(maxIndex>= length) {
                maxIndex = length-1;
            }

            index = distance(env, max_element(env + i*chunkLength, env+maxIndex));
            maximumValuesPerChunk[i+1] = env[index];
            indexOfMaximumPerChunk[i+1] = static_cast<float>(index);
        }


        for(i = 1; i< NChunks+2; i++) {
            if( indexOfMaximumPerChunk[i] - indexOfMaximumPerChunk[i-1] < 5 ) {
                indexOfMaximumPerChunk[i] = (indexOfMaximumPerChunk[i] + indexOfMaximumPerChunk[i+1])/2;
                maximumValuesPerChunk[i] = (maximumValuesPerChunk[i] + maximumValuesPerChunk[i+1])/2;
            }
        }

        indexOfMaximumPerChunk[0] = -1;
        maximumValuesPerChunk[0] = maximumValuesPerChunk[1];
        //create two virtual chucks at indices -1 and length. When set to 0 or length-1, a double entry can occur if the
        //first or the last element of the list is the maximum.
        indexOfMaximumPerChunk[NChunks + 1] = static_cast<float>(length);
        maximumValuesPerChunk[NChunks + 1] = maximumValuesPerChunk[NChunks];

        float maxValue = maximumValuesPerChunk[0];
        for(int i = 1; i< NChunks + 1; i++) {
            if(maxValue < maximumValuesPerChunk[i]) {
                maxValue = maximumValuesPerChunk[i];
            }
        }

        UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(indexOfMaximumPerChunk,
                maximumValuesPerChunk,NChunks+2);


        for(i = 0; i < length; i++) {
            env[i] = spline.evaluate(  static_cast<float>( i ) );
            if(env[i] < 0.01*maxValue) {
                env[i] = 0.01*maxValue;
            }
        }
        //delete spline;
        return env;
    }






    static void equalizeEnvelopeInPlace(float** spectra,size_t N,size_t M) {
        /*
        spectra: 2D array to remove the envelope for
        N: array length
        M: array length

        This function overwrites the array it receives with a signal with the envelope
        removed such that it resembles a block wave. This does work well visually, but
        do not use for data extraction. For example, a signal with beating will be
        distorted.
        */
        int i, j;
        float* env;


        for(i = 0 ; i<N; i++) {




            env = calculateEnvelope(spectra[i],M);

            for(j = 0; j< M; j++) {
                spectra[i][j] /= (0.99*env[j] + 0.01);
            }



            float mean = 0.0;
            for(j = 0; j< M; j++) {
                mean += spectra[i][j];
            }
            mean /= M;
            for(j = 0; j< M; j++) {
                spectra[i][j] -= mean;
            }

            complex<float>* hb = hilbert(spectra[i],M);
            for(j = 0; j< M; j++) {
                env[j] = abs( hb[j]);
                spectra[i][j] /= (0.9999*env[j] + 0.0001);
            }


            delete[] hb;
            delete[] env;

        }

    }

    static void stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower) {
        /*
        spectra: values that will be replaced with a stretched version.
        referenceSpectrum: measured reference spectrum
        minimumReferencePower: minimum fraction of the reference spectrum that need to be present in order to keep the data


        The FFT works well over any OCT data, however the RIAA algorithm needs an apodized interference spectum. However,
        the measured OCT interference is zero or close to zero near the edges. This function takes the reference spectrum,
        and keeps only all values that are between the edges as given by the minimum reference power.

        In order to keep most information in the data, but also to keep the data in a power of 2 array length, the
        center data is stretched and interpolated to fit again on the same array shape.

        This function processes in place, and overwrites float** spectra.
        */
        int minIndex, maxIndex;
        float maximumReference =0.0;
        int indexOfMaximum = 0;


        for(int i = 0; i < settings->sizeZSpectrum; i++) {
            if(referenceSpectrum[i] > maximumReference) {
                maximumReference = referenceSpectrum[i];
                indexOfMaximum = i;
            }
        }

        for(minIndex = indexOfMaximum; referenceSpectrum[minIndex] > minimumReferencePower*maximumReference && minIndex > 0; minIndex--) {}
        for(maxIndex = indexOfMaximum; referenceSpectrum[maxIndex] > minimumReferencePower*maximumReference
                && maxIndex < settings->sizeZSpectrum; maxIndex++) {}

        float xrange[settings->sizeZSpectrum];
        for(int i = 0; i < settings->sizeZSpectrum; i++) {
            xrange[i] = 1.0*i/settings->sizeZSpectrum;
        }

        for(int i = 0; i < settings->sizeXSpectrum; i++) {

            UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(xrange,spectra[i],
                    settings->sizeZSpectrum);
            for(int j = 0; j < settings->sizeZSpectrum; j++) {
                float x = 1.0*j*( maxIndex - minIndex ) / (settings->sizeZSpectrum*settings->sizeZSpectrum) + (1.0*minIndex)/settings->sizeZSpectrum;
                spectra[i][j] = spline.evaluate(x);
            }


            //delete spline;
        }

        return;
    }






    static void preprocessSpectrumInPlace(float** spectra, float* offsetSpectrum, float* chirp, float* referenceSpectrum,float* window   ) {
        /*
        spectra: data to process
        offsetSpectrum: offset spectrum that is subtracted from the measured spectrum before processing
        chirp: nonlinearity of the Ganymede spectrometer.
        referenceSpectrum: spectrum of the reference arm
        window: window function

        This function preprocesses the spectra before the FFT or RIAA. The noisy low signal data is
        removed from the edges based on the reference power. If a reference spectrum is not a nullpointer,
        this is removed from the envelope. The signal is compensated for dispersion based on the
        polynomial coefficients stored in the global settings.
        */
        int i, j;


        if(offsetSpectrum) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                for(i = 0; i < settings->sizeXSpectrum ; i++) {

                    spectra[i][j] -= offsetSpectrum[j]  ;
                }
            }
            if(referenceSpectrum) {
                for(j = 0; j < settings->sizeZSpectrum; j++) {
                    referenceSpectrum[j] -= offsetSpectrum[j];
                }
            }
        }
        if(!referenceSpectrum) {
            referenceSpectrum = new float[settings->sizeZSpectrum];
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                referenceSpectrum[j] = spectra[0][j];
                for(i = 1; i < settings->sizeXSpectrum ; i++) {
                    referenceSpectrum[j] += abs(spectra[i][j]);
                }
                referenceSpectrum[j] /= settings->sizeXSpectrum;
            }
        }


        for(j = 0; j < settings->sizeZSpectrum; j++) {
            if(referenceSpectrum[j] <= 1e-2 ) {
                if(j == 0) {
                    referenceSpectrum[j] = referenceSpectrum[j+1] ;
                } else if(j == settings->sizeZSpectrum-1) {
                    referenceSpectrum[j] = referenceSpectrum[j-1] ;
                } else {
                    referenceSpectrum[j] = 0.5*(referenceSpectrum[j+1] + referenceSpectrum[j-1]) ;
                }
            }
        }


        for(j = 0; j < settings->sizeZSpectrum; j++) {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {
                spectra[i][j] -= referenceSpectrum[j]  ;
            }
        }




        kiss_fft_cfg cfg = kiss_fft_alloc(settings->sizeZSpectrum, 0, NULL, NULL);
        kiss_fft_cfg icfg = kiss_fft_alloc(settings->sizeZSpectrum, 1, NULL, NULL);


        kiss_fft_cpx positiveFrequencies[settings->sizeZSpectrum];
        kiss_fft_cpx negativeFrequencies[settings->sizeZSpectrum];

        kiss_fft_cpx positiveSignal[settings->sizeZSpectrum];
        kiss_fft_cpx negativeSignal[settings->sizeZSpectrum];



        double temp[settings->sizeZSpectrum];

        for(i = 0; i < settings->sizeXSpectrum ; i++) {
            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                positiveSignal[j].r = spectra[i][j];
                positiveSignal[j].i = 0.0;
            }



            kiss_fft(cfg, positiveSignal, positiveFrequencies);

            for(j = 0; j < settings->sizeZSpectrum /2; j++) {
                negativeFrequencies[j].r = 0.0;
                negativeFrequencies[j].i = 0.0;
            }


            for(j = settings->sizeZSpectrum/2; j < settings->sizeZSpectrum ; j++) {
                negativeFrequencies[j].r = positiveFrequencies[j].r;
                negativeFrequencies[j].i = positiveFrequencies[j].i;
                positiveFrequencies[j].r = 0.0;
                positiveFrequencies[j].i = 0.0;
            }
            kiss_fft(icfg, positiveFrequencies,positiveSignal);
            kiss_fft(icfg, negativeFrequencies,negativeSignal);



            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                double dispersionPhase = j*settings->dispersionCoefficients[0];
                for(int k = 1; k < settings->numberOfDispersionCoefficients-1; k++) {
                    dispersionPhase = j*(settings->dispersionCoefficients[k] + dispersionPhase);
                }
                if(settings->numberOfDispersionCoefficients > 1) {
                    dispersionPhase += settings->dispersionCoefficients[settings->numberOfDispersionCoefficients-1];
                }
                dispersionPhase  *= 2*M_PI;
                positiveSignal[j].r = positiveSignal[j].r  * cos(dispersionPhase ) + positiveSignal[j].i*sin(dispersionPhase) ;
                negativeSignal[j].r = negativeSignal[j].r  * cos(dispersionPhase ) - negativeSignal[j].i*sin(dispersionPhase) ;
                temp[j] = negativeSignal[j].r;
            }

            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                spectra[i][j] = positiveSignal[j].r + negativeSignal[j].r;
            }
        }




        kiss_fft_free(cfg);
        kiss_fft_free(icfg);



        for(j = 0; j < settings->sizeZSpectrum; j++) {
            double meanAmplitude= abs(spectra[0][j]);
            for(i = 1; i < settings->sizeXSpectrum ; i++) {
                meanAmplitude += abs(spectra[i][j]);
            }

            meanAmplitude /= settings->sizeXSpectrum;
            if(window) {
                for(i = 0; i < settings->sizeXSpectrum ; i++) {
                    spectra[i][j] *= window[j]/meanAmplitude;
                }
            } else {
                for(i = 0; i < settings->sizeXSpectrum ; i++) {
                    spectra[i][j] /= meanAmplitude;
                }
            }
        }




        if(chirp) {
            /*
            This chirp is the chirp caused by the spectrometer of the Ganymede system.
            The array consists of floating point index numbers of the CCD of the spectrometer, these are the actual
            index numbers of of the measured data. For example, a number 11.1 on index 11 of the chirp array, would
            mean that the actual index of the measured spectra is 11.1. You cannot select index 11.1, so an interpolation
            is done here.
            */
            for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
                UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(chirp,spectra[i],
                        settings->sizeZSpectrum);
                for(j = 0; j < settings->sizeZSpectrum; j++) {
                    spectra[i][j] = spline.evaluate(  static_cast<float>( j) );
                }
                //delete spline;
            }
        }



        stretchSpectraInPlace(spectra,referenceSpectrum,0.01);

        for(i = 0; i < settings->sizeXSpectrum  ; i++) {
            float meanValue = 0.0;
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                meanValue += spectra[i][j];
            }
            meanValue /= settings->sizeZSpectrum;
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectra[i][j] -= meanValue;
            }
        }

    }






    static void preprocessSpectrumInPlaceOud(float** spectra, float* offsetSpectrum, float* chirp, float* referenceSpectrum   ) {
        /*
        TODO: remove in next commit.
        */
        int i, j;


        for(j = 0; j < settings->sizeZSpectrum; j++) {
            if(referenceSpectrum[j] <= 1e-2 ) {
                if(j == 0) {
                    referenceSpectrum[j] = referenceSpectrum[j+1] ;
                } else if(j == settings->sizeZSpectrum-1) {
                    referenceSpectrum[j] = referenceSpectrum[j-1] ;
                } else {
                    referenceSpectrum[j] = 0.5*(referenceSpectrum[j+1] + referenceSpectrum[j-1]) ;
                }

            }
        }

        if(offsetSpectrum) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                referenceSpectrum[j] -= offsetSpectrum[j];
            }
        }

        for(j = 0; j < settings->sizeZSpectrum; j++) {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {
                if(offsetSpectrum) {
                    spectra[i][j] -= offsetSpectrum[j]  ;
                }
                spectra[i][j] -= referenceSpectrum[j]  ;
                //spectra[i][j] *= (window[j]) /referenceSpectrum[j] ;
                spectra[i][j] *= 1.0 /(referenceSpectrum[j]) ;
            }
        }


        kiss_fft_cfg cfg = kiss_fft_alloc(settings->sizeZSpectrum, 0, NULL, NULL);
        kiss_fft_cfg icfg = kiss_fft_alloc(settings->sizeZSpectrum, 1, NULL, NULL);


        kiss_fft_cpx positiveFrequencies[settings->sizeZSpectrum];
        kiss_fft_cpx negativeFrequencies[settings->sizeZSpectrum];

        kiss_fft_cpx positiveSignal[settings->sizeZSpectrum];
        kiss_fft_cpx negativeSignal[settings->sizeZSpectrum];



        double temp[settings->sizeZSpectrum];

        for(i = 0; i < settings->sizeXSpectrum ; i++) {
            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                positiveSignal[j].r = spectra[i][j];
                positiveSignal[j].i = 0.0;
            }



            kiss_fft(cfg, positiveSignal, positiveFrequencies);

            for(j = 0; j < settings->sizeZSpectrum /2; j++) {
                negativeFrequencies[j].r = 0.0;
                negativeFrequencies[j].i = 0.0;
            }


            for(j = settings->sizeZSpectrum/2; j < settings->sizeZSpectrum ; j++) {
                negativeFrequencies[j].r = positiveFrequencies[j].r;
                negativeFrequencies[j].i = positiveFrequencies[j].i;
                positiveFrequencies[j].r = 0.0;
                positiveFrequencies[j].i = 0.0;
            }
            kiss_fft(icfg, positiveFrequencies,positiveSignal);
            kiss_fft(icfg, negativeFrequencies,negativeSignal);



            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                double dispersionPhase = j*settings->dispersionCoefficients[0];
                for(int k = 1; k < settings->numberOfDispersionCoefficients-1; k++) {
                    dispersionPhase = j*(settings->dispersionCoefficients[k] + dispersionPhase);
                }
                if(settings->numberOfDispersionCoefficients > 1) {
                    dispersionPhase += settings->dispersionCoefficients[settings->numberOfDispersionCoefficients-1];
                }
                dispersionPhase  *= 2*M_PI;
                positiveSignal[j].r = positiveSignal[j].r  * cos(dispersionPhase ) + positiveSignal[j].i*sin(dispersionPhase) ;
                negativeSignal[j].r = negativeSignal[j].r  * cos(dispersionPhase ) - negativeSignal[j].i*sin(dispersionPhase) ;
                temp[j] = negativeSignal[j].r;
            }

            for(j = 0; j < settings->sizeZSpectrum ; j++) {
                spectra[i][j] = positiveSignal[j].r + negativeSignal[j].r;
            }
        }

        kiss_fft_free(cfg);
        kiss_fft_free(icfg);


        float averagePower[settings->sizeXSpectrum];
        for(i = 0; i < settings->sizeXSpectrum; i++) {
            averagePower[i] = 0.0;
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                averagePower[i] += spectra[i][j]*spectra[i][j];
            }
            averagePower[i] /= settings->sizeZSpectrum;
        }



        equalizeEnvelopeInPlace(spectra, settings->sizeXSpectrum, settings->sizeZSpectrum);
        if(chirp) {
            for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
                UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(chirp,spectra[i],
                        settings->sizeZSpectrum);
                for(j = 0; j < settings->sizeZSpectrum; j++) {
                    spectra[i][j] = spline.evaluate(  static_cast<float>( j) );
                }
                //delete spline;
            }
        }
        for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
            float meanValue = 0.0;
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                meanValue += spectra[i][j];
            }
            meanValue /= settings->sizeZSpectrum;
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectra[i][j] -= meanValue;
            }
        }
        stretchSpectraInPlace(spectra,referenceSpectrum,0.01);

        for(i = 0; i < settings->sizeXSpectrum; i++) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectra[i][j] *= sqrt(averagePower[i]);
            }
        }
    }
};





int main(int argc, char *argv[]) {

    QApplication app(argc, argv);
    //TODO: header is utf-8 in the Ganymede software, here ASCI. If somebody puts an emoticon into the filename, the code might fail.
    settings = make_shared<Settings>();
    //The path is required for the import of the objective settings file.
    settings->pathToExecutable = QCoreApplication::applicationDirPath().toStdString();

    map<string, vector<double>> mapFromData;
    mapFromData = IO<double>::loadObjectiveDispersionData(
                      settings->pathToExecutable + "\\..\\..\\settings\\objectives.yaml");
    settings->objectiveDispersionData.insert(mapFromData.begin(), mapFromData.end());

    Window window;
    window.show();
    //the scan object contains all data and some non static processing functions.


    string filePath = "D:\\data\\ThorlabsCppTestData\\wedgeBscan\\sk.txt";
    /*


    tuple<float**, int, int> ctuple = IO<float>::load2DArrayFromFile("D:\\data\\ThorlabsCppTestData\\cLevinson.txt");
    int N = get<2>(ctuple);
    complex<float>* carray = new complex<float>[N];


    cout << N <<" " <<get<1>(ctuple ) << " " << get<2>(ctuple) << endl;

    for(int i = 0; i < N;i++){
        carray[i] = complex<float>(get<0>(ctuple)[0][i] , get<0>(ctuple)[1][i]);
    }


    complex<float>* Aarray = new complex<float>[N];

    vector<complex<float>> cvector = vector<complex<float>>(N);
    vector<complex<float>> Avector = vector<complex<float>>(N);
    for(int i = 0; i< N; i++){
        cvector[i] = carray[i];
    }


    UtilityMathFunctions<float>::levinson(carray,N,Aarray);
    UtilityMathFunctions<float>::levinson(cvector,Avector);

    for(int i = 0; i < N;i++){
        cout << i << " " << Aarray[i] << " " << Avector[i] << endl;
    }

    */



    scan = make_shared<BScan>();

    //set the image. When not debugging, settings->sizeXSpectrum and settings->sizeZSpectrum will be zero.
    window.setImage(scan->fftBScan(),settings->sizeXSpectrum,settings->sizeZSpectrum);
    cout << "end" << endl;

    //return 0;
    return app.exec();
}


