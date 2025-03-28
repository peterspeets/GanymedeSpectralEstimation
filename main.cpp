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


    static float* calculateEnvelope(float *spectrum,size_t length){
        float* env = new float[length];
        int NChunks = settings->NChunksEnvelopeSubtraction;
        int chunkLength = length/NChunks;
        float maximumValuesPerChunk[NChunks+2];
        float indexOfMaximumPerChunk[NChunks+2];
        int index;
        int i, j;

        for(i = 0; i< length; i++){
            env[i] =  abs(spectrum[i]);
        }

        for(i = 0; i< NChunks; i++){
            int maxIndex = (i+1)*chunkLength;
            if(maxIndex>= length){
               maxIndex = length-1;
            }

            index = distance(env, max_element(env + i*chunkLength, env+maxIndex));
            maximumValuesPerChunk[i+1] = env[index];
            indexOfMaximumPerChunk[i+1] = static_cast<float>(index);
        }

        indexOfMaximumPerChunk[0] = 0;
        maximumValuesPerChunk[0] = maximumValuesPerChunk[1];
        indexOfMaximumPerChunk[NChunks + 1] = static_cast<float>(length-1);
        maximumValuesPerChunk[NChunks + 1] = maximumValuesPerChunk[NChunks];

        IO<float>::saveArrayToFile( indexOfMaximumPerChunk, NChunks+2, "D:\\data\\inds.txt");
        IO<float>::saveArrayToFile( maximumValuesPerChunk, NChunks+2, "D:\\data\\maxs.txt");


        UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(indexOfMaximumPerChunk,maximumValuesPerChunk,NChunks+2);
        for(i = 0; i < length; i++) {
            env[i] = spline->evaluate(  static_cast<float>( i ) );
        }
        delete spline;
        return env;
    }


    static void equalizeEnvelopeInPlace(float** spectra,size_t N,size_t M){
        int i, j;
        float* env;


        for(i = 0 ; i<N;i++){

            env = calculateEnvelope(spectra[i],M);
            for(j = 0; j< M; j++){
                spectra[i][j] /= (0.5*env[j] + 0.5);
            }
            delete[] env;
        }

    }


    static void preprocessSpectrumInPlace(float** spectra, float* offsetSpectrum, float* chirp, float* referenceSpectrum   ) {
        int i, j;
        double* window = hannWindow(settings->sizeZSpectrum);
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
                spectra[i][j] *= window[j] /referenceSpectrum[j] ;
            }
        }

        equalizeEnvelopeInPlace(spectra, settings->sizeXSpectrum, settings->sizeZSpectrum);

        IO<float>::saveArrayToFile(spectra[250], settings->sizeZSpectrum, "D:\\data\\env.txt");

        if(chirp){
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
    pair<float*,int> xpair = IO<float>::loadArrayFromFile("D:\\data\\spectrum.txt");
    float* x = xpair.first;
    int N = xpair.second;
    int q_i = 10;
    int K = 16*N;
    double vt = 1.0;


    pair<float*, float*> riaa_res;
    //riaa_res = UtilityMathFunctions<float>::fiaa_oct(x, N, K, q_i, vt);
    //IO<float>::saveArrayToFile(riaa_res.first, K, "D:\\data\\riaa.txt");



    float* offset = nullptr;
    float* chirp = nullptr;
    float* referenceSpectrum = nullptr;
    float* intensity = nullptr;
    float** spectra = nullptr;
    float** processedBscan = nullptr;

    if(filePath.substr(filePath.length() - 4) == ".oct"){
        IO<float>::GanymedeFileLoader fileLoader =  IO<float>::GanymedeFileLoader(filePath);
        settings = fileLoader. loadSettings() ;
        offset = fileLoader.loadCalibrationSpectrum(settings->pathOffset);
        chirp = fileLoader.loadCalibrationSpectrum(settings->pathChirp);
        referenceSpectrum = fileLoader.loadCalibrationSpectrum(settings->pathApodization);
        intensity = fileLoader.loadCalibrationSpectrum(settings->pathIntensity);
        spectra = fileLoader.loadSpectrum(0);

    }else{
        settings = make_shared<Settings>(filePath);
        tuple<float**,int,int> loadingSpectraTuple =  IO<float>::load2DArrayFromFile(filePath + "rawData.txt");
        pair<float*,int> loadingReferencePair =  IO<float>::loadArrayFromFile(filePath + "sk.txt");
        spectra = get<0>(loadingSpectraTuple);
        referenceSpectrum = loadingReferencePair.first;
    }


    cout << "Spectrum size: " << settings->sizeXSpectrum << "x" << settings->sizeZSpectrum  <<endl;



    cout << "Preprocessing:";
    SignalProcessing::preprocessSpectrumInPlace(spectra,offset, chirp,referenceSpectrum);
    cout << " done" << endl;

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


    //IO<float>::savePng("D:\\data\\spectra.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum,  spectra);
    //IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2,  image_fft,true );


    processedBscan = UtilityMathFunctions<float>::processBScan(spectra,  settings->sizeXSpectrum,settings->sizeZSpectrum,  K, q_i, 1.0);
    float** image = processedBscan ;

    //IO<float>::savePng("D:\\data\\testImage3.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, settings->sizeXSpectrum-0,  789,  image,true );
    //IO<float>::savePng("D:\\data\\testImage4.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, 789,  settings->sizeZSpectrum/2,  image,true );


    double x_res = settings->x_mm/(settings->sizeXSpectrum - 0);
    double z_res = settings->z_mm/(settings->sizeZSpectrum/2);


    if(x_res > z_res) {
        //x_res poorer,reshape to keep z_res
        double newXShape_dbl = (settings->sizeXSpectrum - 0)*(x_res/z_res);
        int newXShape = static_cast<int>(round(newXShape_dbl));
        cout << "saving image (keeping z shape)" << newXShape << "x" <<  settings->sizeZSpectrum/2 <<  endl;
        IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0,  settings->sizeZSpectrum/2, 2*newXShape,  2*settings->sizeZSpectrum/2,  image_fft,true );
        IO<float>::savePng("D:\\data\\testImageRFIAA.png", settings->sizeXSpectrum-0,  K/2, 2*newXShape,  2*settings->sizeZSpectrum/2,  image,true );
    } else {
        //z_res poorer,reshape to keep x_res
        double newZShape_dbl = (settings->sizeZSpectrum/2)*(z_res/x_res);
        int newZShape = static_cast<int>(round(newZShape_dbl));
        cout << "saving image (keeping x shape) " << (settings->sizeXSpectrum - 0) << "x" <<  newZShape<<  endl;
        IO<float>::savePng("D:\\data\\testImageFFT.png", settings->sizeXSpectrum-0, settings->sizeZSpectrum/2,  2*settings->sizeXSpectrum-0, 2*newZShape,  image_fft,true );
        IO<float>::savePng("D:\\data\\testImageRFIAA.png", settings->sizeXSpectrum-0, K/2,  2*settings->sizeXSpectrum-0, 2*newZShape,  image,true );
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
