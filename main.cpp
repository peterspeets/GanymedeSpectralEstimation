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


using namespace std;

shared_ptr<Settings> settings;


class SignalProcessing {

public:

    static double* hannWindow(const int N) {
        int i;
        double* window = new double[N];
        for(i = 0; i< N; i++) {
            window[i] = sin(M_PI*i/N)*sin(M_PI*i/N);
        }
        return window;
    }

    static double* tukeyWindow(const int N, const double alpha) {
        double* window = new double[N];
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

        //cout << "\ntest: " << signal[500] << endl;
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
            output[i] = complex<T>(signal[i].real() - ft1[i].i/N ,signal[i].imag() + ft1[i].r /N) ;

        }

        //cout << "\ntest: " << output[500] << endl;

        kiss_fft_free(cfg);


        return output;


    }


    static float* calculateEnvelope(float *spectrum,size_t length) {
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
            if( indexOfMaximumPerChunk[i] - indexOfMaximumPerChunk[i-1] < 5 ){
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
        UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(indexOfMaximumPerChunk,maximumValuesPerChunk,NChunks+2);
        for(i = 0; i < length; i++) {
            env[i] = spline->evaluate(  static_cast<float>( i ) );

        }
        delete spline;
        return env;
    }





    static void equalizeEnvelopeInPlace(float** spectra,size_t N,size_t M) {
        int i, j;
        float* env;


        for(i = 0 ; i<N; i++) {




            env = calculateEnvelope(spectra[i],M);

            if(i == 121) {

                IO<float>::saveArrayToFile( spectra[i], M, "D:\\data\\spe0.txt");
                IO<float>::saveArrayToFile( env, M, "D:\\data\\env0.txt");
            }

            for(j = 0; j< M; j++) {
                spectra[i][j] /= (0.99*env[j] + 0.01);
            }


            complex<float>* hb = hilbert(spectra[i],M);
            for(j = 0; j< M; j++) {
                env[j] = abs( hb[j]);
                spectra[i][j] /= (0.9999*env[j] + 0.0001);
            }





            if(i == 121) {

                IO<float>::saveArrayToFile( spectra[i], M, "D:\\data\\spe2.txt");
                IO<float>::saveArrayToFile( env, M, "D:\\data\\env2.txt");
            }
            delete[] hb;
            delete[] env;

        }

    }


    static void preprocessSpectrumInPlace(float** spectra, float* offsetSpectrum, float* chirp, float* referenceSpectrum   ) {
        int i, j;
        //double* window = hannWindow(settings->sizeZSpectrum);
        double* window = tukeyWindow(settings->sizeZSpectrum,0.1);
        //for(i = 0; i < settings->sizeZSpectrum; i++){
        //   window[i] = 1.0;
        //}


        IO<float>::saveArrayToFile( spectra[settings->sizeXSpectrum/2], settings->sizeZSpectrum, "D:\\data\\rawSpectrum.txt");

        IO<float>::saveArrayToFile( referenceSpectrum, settings->sizeZApodization, "D:\\data\\reference.txt");

        IO<float>::saveArrayToFile( window, settings->sizeZSpectrum, "D:\\data\\window.txt");
        if(chirp)
            IO<float>::saveArrayToFile( chirp, settings->sizeZChirp, "D:\\data\\chirp.txt");
        if(offsetSpectrum)
            IO<float>::saveArrayToFile( offsetSpectrum, settings->sizeZOffset, "D:\\data\\offset.txt");

        for(j = 0; j < settings->sizeZSpectrum; j++) {
            if(referenceSpectrum[j] <= 1e-2) {
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


        for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                if(offsetSpectrum) {
                    spectra[i][j] -= offsetSpectrum[j]  ;
                }
                spectra[i][j] -= referenceSpectrum[j]  ;
                spectra[i][j] *= (window[j]) /referenceSpectrum[j] ;
                //spectra[i][j] *= 1.0 /referenceSpectrum[j] ;
            }
        }

        IO<float>::saveArrayToFile( spectra[0], settings->sizeZSpectrum, "D:\\data\\tst.txt");

        equalizeEnvelopeInPlace(spectra, settings->sizeXSpectrum, settings->sizeZSpectrum);

        IO<float>::saveArrayToFile(spectra[250], settings->sizeZSpectrum, "D:\\data\\env.txt");

        if(chirp) {
            for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
                UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(chirp,spectra[i],settings->sizeZSpectrum);
                for(j = 0; j < settings->sizeZSpectrum; j++) {
                    spectra[i][j] = spline->evaluate(  static_cast<float>( j) );
                }
                delete spline;
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

        IO<float>::saveArrayToFile( spectra[settings->x_px/2], settings->sizeZSpectrum, "D:\\data\\processedSpectrum.txt");

        delete[] window;
    }

};





int main() {

    //string filePath = "C:\\data\\ThorlabsCppTestData\\MicSlideTest\\MicSlideTest_0004_Mode2D.oct";
    //string filePath = "C:\\data\\ThorlabsCppTestData\\MilkTest\\Milk flow measurement_0028_Mode2D.oct";
    //string filePath = "C:\\data\\ThorlabsCppTestData\\MilkTest\\Milk flow measurement_0001_ModeDoppler.oct";
    //string filePath = "C:\\cpp\\onionBscan\\";
    string filePath = "C:\\cpp\\wedgeBscan\\";


    //TODO: header is utf-8 in the Ganymede software, here ASCI. If somebody puts an emoticon into the filename, the code might fail.
    settings = make_shared<Settings>();


    cout << "load test data."<< endl;
    pair<float*,int> xpair = IO<float>::loadArrayFromFile("D:\\data\\wedgeSpectrum.txt");
    float* x = xpair.first;
    int N = xpair.second;
    int q_i = 2;
    int K = 16*N;
    double vt = 1.0;




    pair<float*, float*> riaa_res;
    //riaa_res = UtilityMathFunctions<float>::fiaa_oct(x, N, K, q_i, vt);



    IO<float>::saveArrayToFile(x, N, "D:\\data\\x.txt");


    float* offset = nullptr;
    float* chirp = nullptr;
    float* referenceSpectrum = nullptr;
    float* intensity = nullptr;
    float** spectra = nullptr;
    float** processedBscan = nullptr;

    if(filePath.substr(filePath.length() - 4) == ".oct") {
        IO<float>::GanymedeFileLoader fileLoader =  IO<float>::GanymedeFileLoader(filePath);
        settings = fileLoader. loadSettings() ;
        offset = fileLoader.loadCalibrationSpectrum(settings->pathOffset);
        chirp = fileLoader.loadCalibrationSpectrum(settings->pathChirp);
        referenceSpectrum = fileLoader.loadCalibrationSpectrum(settings->pathApodization);
        intensity = fileLoader.loadCalibrationSpectrum(settings->pathIntensity);
        spectra = fileLoader.loadSpectrum(0);

    } else {
        settings = make_shared<Settings>(filePath);
        tuple<float**,int,int> loadingSpectraTuple =  IO<float>::load2DArrayFromFile(filePath + "rawData.txt");
        pair<float*,int> loadingReferencePair =  IO<float>::loadArrayFromFile(filePath + "sk.txt");
        spectra = get<0>(loadingSpectraTuple);
        referenceSpectrum = loadingReferencePair.first;
    }


    cout << "Spectrum size: " << settings->sizeXSpectrum << "x" << settings->sizeZSpectrum  <<endl;

    //spectra[0] = spectra[settings->sizeXSpectrum-10];

    //settings->sizeXSpectrum = 1;
    //settings->sizeZSpectrum /= 2;
    //for(int i = 0; i < settings->sizeXSpectrum; i++){
    //        spectra[i] += settings->sizeZSpectrum /2 ;
    //}
    //referenceSpectrum += settings->sizeZSpectrum /2 ;

    cout << "Preprocessing:";
    SignalProcessing::preprocessSpectrumInPlace(spectra,offset, chirp,referenceSpectrum);
    cout << " done" << endl;
    //processedBscan = UtilityMathFunctions<float>::processBScan(spectra,  settings->sizeXSpectrum,settings->sizeZSpectrum,  K, q_i, 1.0);




    float** image_fft = new float*[settings->sizeXSpectrum];
    for (int i = 0; i < settings->sizeXSpectrum; i++) {
        image_fft[i] = new float[settings->sizeZSpectrum];
    }

    kiss_fft_cfg icfg = kiss_fft_alloc(settings->sizeZSpectrum, 1, NULL, NULL);

    kiss_fft_cpx in[settings->sizeZSpectrum];
    kiss_fft_cpx out[settings->sizeZSpectrum];

    for (int j = 0; j < settings->sizeZSpectrum; j++) {
        in[j].i = 0.0;
    }


    for (int i = 0; i < settings->sizeXSpectrum; i++) {
        for (int j = 0; j < settings->sizeZSpectrum; j++) {
            in[j].r = spectra[i][j];
        }
        kiss_fft(icfg, in, out);
        for (int j = 1; j < settings->sizeZSpectrum+1; j++) {
            image_fft[i][j-1] = out[j].r*out[j].r + out[j].i*out[j].i;
        }
    }


    kiss_fft_free(icfg);

    IO<float>::saveArrayToFile(image_fft[settings->sizeXSpectrum-10], settings->sizeZSpectrum, "D:\\data\\fft.txt");

    //return 0;


    //IO<float>::savePng("D:\\data\\spectra.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum,  spectra);
    //IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2,  image_fft,true );


    processedBscan = UtilityMathFunctions<float>::processBScan(spectra,  settings->sizeXSpectrum,settings->sizeZSpectrum,  K, q_i, 1.0);




    IO<float>::saveArrayToFile(processedBscan[settings->sizeXSpectrum-10], K, "D:\\data\\riaa.txt");
    IO<float>::saveArrayToFile(processedBscan[settings->sizeXSpectrum-1], K, "D:\\data\\riaalast.txt");

    float** image = processedBscan ;

    //IO<float>::savePng("D:\\data\\testImage3.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, settings->sizeXSpectrum-0,  789,  image,true );
    //IO<float>::savePng("D:\\data\\testImage4.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, 789,  settings->sizeZSpectrum/2,  image,true );


    double x_res = 0.1*settings->x_mm/(settings->sizeXSpectrum - 0);
    double z_res = settings->z_mm/(settings->sizeZSpectrum/2);


    if(x_res > z_res) {
        //x_res poorer,reshape to keep z_res
        double newXShape_dbl = (settings->sizeXSpectrum - 0)*(x_res/z_res);
        int newXShape = static_cast<int>(round(newXShape_dbl));
        cout << "saving image (keeping z shape)" << newXShape << "x" <<  settings->sizeZSpectrum/2 <<  endl;
        IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, 2*newXShape,  2*settings->sizeZSpectrum/2,  image_fft,true );
        IO<float>::savePng("D:\\data\\testImageRFIAA.png", settings->sizeXSpectrum-1,  K/2, 2*newXShape,  2*settings->sizeZSpectrum/2,  image,true );
    } else {
        //z_res poorer,reshape to keep x_res
        double newZShape_dbl = (settings->sizeZSpectrum/2)*(z_res/x_res);
        int newZShape = static_cast<int>(round(newZShape_dbl));
        cout << "saving image (keeping x shape) " << (settings->sizeXSpectrum - 0) << "x" <<  newZShape<<  endl;
        IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0, settings->sizeZSpectrum/2,  2*settings->sizeXSpectrum-0, 2*newZShape,  image_fft,true );
        IO<float>::savePng("D:\\data\\testImageRFIAA.png", settings->sizeXSpectrum-1, K/2,  2*settings->sizeXSpectrum-0, 2*newZShape,  image,true );
    }


    double newXShape_dbl = (settings->sizeXSpectrum - 0)*(x_res/z_res);
    int newXShape = static_cast<int>(round(newXShape_dbl));
    //cout << "saving image (keeping z shape)" << newXShape << "x" <<  settings->sizeZSpectrum/2 <<  endl;




    for (int i = 0; i < settings->sizeXSpectrum; i++) {
        delete[] image[i];
    }
    delete[] image;






    //system("pause");
    cout << "end" << endl;
    return 0;


}
