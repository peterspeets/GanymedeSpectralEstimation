#include "BScan.h"

BScan::BScan() {
    /*
    Constructor for an empty B-scan. The settings are set to default values.
    */

    shared_ptr<Settings> oldSettings = settings;
    settings = make_shared<Settings>();

    settings->pathToExecutable = oldSettings->pathToExecutable;
    map<string, vector<double>> mapFromData = IO<double>::loadObjectiveDispersionData(
                                 settings->pathToExecutable + "\\..\\..\\settings\\objectives.yaml");
    settings->objectiveDispersionData.insert(mapFromData.begin(), mapFromData.end());

    settings->bytesPerPixelIntensity = 0;
    settings->bytesPerPixelSpectrum = 0;
    settings->bytesPerPixelChirp = 0;
    settings->bytesPerPixelApodization = 0;
    settings->bytesPerPixelOffset = 0;
    settings->sizeXIntensity = 0;
    settings->sizeZIntensity = 0;
    settings->sizeXSpectrumRaw = 0;
    settings->sizeXSpectrum = 0;
    settings->sizeZSpectrum = 0;
    settings->sizeXChirp = 0;
    settings->sizeZChirp = 0;
    settings->sizeXApodization = 0;
    settings->sizeZApodization = 0;
    settings->sizeXOffset = 0;
    settings->sizeZOffset = 0;
    settings->spectrumAveraging = 0;
    settings->numberOfAScans = 0;
    settings->numberOfBScans = 0;
    settings->x_px = 0;
    settings->z_px = 0;

}


BScan::BScan(const string filePath) {
    /*
    filePath: path to the spectral data.

    Constructor for a B-scan with a given filename. If the file name ends in .oct, it is assumed to
    be Thorlabs .oct data, otherwise it is considered test data stored in .txt files.
    */
    shared_ptr<Settings> oldSettings = settings;

    if(filePath.substr(filePath.length() - 4) == ".oct") {
        //If the file ends in .oct, the file loading is slightly more involved. Use the GanymedeFileLoader.

        cout << "file loader" << endl;

        IO<float>::GanymedeFileLoader fileLoader =  IO<float>::GanymedeFileLoader(filePath);
        cout << "load settigns" << endl;
        settings = fileLoader.loadSettings() ;
        cout << "load calibratio spectra:" << endl;
        cout << "offset" << endl;
        offset = fileLoader.loadCalibrationSpectrum(settings->pathOffset);
        cout << "chirp" << endl;
        chirp = fileLoader.loadCalibrationSpectrum(settings->pathChirp);
        cout << "reference" << endl;
        referenceSpectrum = fileLoader.loadCalibrationSpectrum(settings->pathApodization);
        cout << "intensity" << endl;
        intensity = fileLoader.loadCalibrationSpectrum(settings->pathIntensity,2048,4);

        settings->numberOfStoredSpectra = settings->pathsSpectra.size();

        vector<vector<float>> singleAScan;
        cout << "Reserve mem for spectra." << endl;
        spectra = vector<vector<float>>(settings->sizeXSpectrum, vector<float>(settings->sizeZSpectrum, 0));

        cout << "paths to spectra" << endl;
        for (pair<int,string> const& pathSpectraPair : settings->pathsSpectra) {
            if(pathSpectraPair.first  == 0) {
                singleAScan = fileLoader.loadSpectrum(pathSpectraPair.first);
                for(int i = 0; i< settings->sizeXSpectrum; i++) {
                    for(int j = 0; j< settings->sizeZSpectrum; j++) {
                        spectra[i][j] = singleAScan[i][j]/(1.0*settings->numberOfStoredSpectra);
                    }
                }
            } else {
                singleAScan = fileLoader.loadSpectrum(pathSpectraPair.first);

                for(int i = 0; i< settings->sizeXSpectrum; i++) {
                    for(int j = 0; j< settings->sizeZSpectrum; j++) {
                        spectra[i][j] += singleAScan[i][j]/(1.0*settings->numberOfStoredSpectra);
                    }
                }
            }

            cout <<pathSpectraPair.first << " " << pathSpectraPair.second << endl;
        }
        cout << "..." << endl;
        if (settings->objectiveDispersionData.count("Native")) {
            settings->objectiveDispersionData.erase("Native");
        }

        settings->objectiveLabel = oldSettings->objectiveLabel;


    } else {
        string directoryPath = filePath;
        if(filePath.substr(filePath.length() - 4) == ".txt") {
            for(int i = filePath.length()-1; i > 0; i--) {
                if(filePath[i] == '\\' || filePath[i] == '/') {
                    directoryPath = filePath.substr(0,i+1);
                    break;
                }
            }
        }

        settings = make_shared<Settings>(directoryPath );

        spectra = IO<float>::load2DVectorFromFile(directoryPath + "rawData.txt");
        vector<float> loadingReference =  IO<float>::loadVectorFromFile(directoryPath  + "sk.txt");

        vector<double> loadingDispersion =  IO<double>::loadVectorFromFile(directoryPath  + "phase.txt");



        bool hasNativeDispersionCorrection = false;
        for(pair<string,vector<double>> const objectiveDispersionDataPair : oldSettings->objectiveDispersionData) {
            if(objectiveDispersionDataPair.first == "Native") {
                cout << "Found native label." << endl;
                hasNativeDispersionCorrection = true;
                break;
            }
        }

        settings->objectiveLabel = oldSettings->objectiveLabel;
        if(!hasNativeDispersionCorrection) {
            cout << "Set to native, only once." << endl;
            settings->objectiveLabel = "Native";
        }
        settings->objectiveDispersionData["Native"] = loadingDispersion;



        /*
        Settings is now overwritten by all variables that are loaded by the BScan constructor.However,
        the the problem is that some settings, such as the executable path and some other settings
        should be kept. Therefore, this needs to change to a BScan settings settings, such as the SizeXSpectrum,
        and the global settings, such as the dispersion correction. For now, do this with this temporary solution
        of copying some old settings here.
        */

    }


    settings->pathToExecutable = oldSettings->pathToExecutable;
    map<string, vector<double>> mapFromData = IO<double>::loadObjectiveDispersionData(
                                 settings->pathToExecutable + "\\..\\..\\settings\\objectives.yaml");
    settings->objectiveDispersionData.insert(mapFromData.begin(), mapFromData.end());
    settings->numberOfDispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].size();

    settings->dispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].data();
    settings->alwaysRedoPreprocessing = oldSettings->alwaysRedoPreprocessing;
    preprocessSpectrumInPlace();
    BScanSettings = *settings;

    imageFFT = vector<vector<float>>(settings->sizeXSpectrum, vector<float>(settings->sizeZSpectrum, 0));
    imageRIAA = vector<vector<float>>(settings->sizeXSpectrum, vector<float>(settings->upscalingFactor*settings->sizeZSpectrum, 0));

    fftBScan();
    cout << "Process B scan" << endl;
    processBScan();
    cout << "objective label "<< settings->objectiveLabel << endl;
}

uint64_t BScan::getTime() {
    /*
    Timing function for code optimization.
    */
    using namespace std::chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


void BScan::fftPartBScan(float** spectra, float** image, int Nz,const int startXIndex,const int stopXIndex) {
    /*
    spectra: spectral data
    image: the absolute value of the FFT data. This means the imaginary part is removed.
    Nz: the length of each spectrum.
    startIndex: startIndex of the A-scan. The FFT will be done between startIndex and stopIndex
    stopIndex: last index of the A-scan. The FFT will be done between startIndex and stopIndex

    Do an FFT of only part of the spectra. This function started by several threads to divide the spectrum up
    into different chunks to do the FFT of.
    */
    kiss_fft_cfg icfg = kiss_fft_alloc(Nz, 1, NULL, NULL);
    kiss_fft_cpx in[Nz];
    kiss_fft_cpx out[Nz];
    for (int j = 0; j < Nz; j++) {
        in[j].i = 0.0;
    }


    for (int i = startXIndex; i < stopXIndex; i++) {
        for (int j = 0; j < Nz; j++) {
            in[j].r = spectra[i][j];
        }
        kiss_fft(icfg, in, out);
        for (int j = 1; j < Nz+1; j++) {
            image[i][j-1] = (out[j].r*out[j].r + out[j].i*out[j].i)/(Nz*Nz) ;
        }
    }


    kiss_fft_free(icfg);
}



void BScan::fftPartBScan(vector<vector<float>>& spectra, vector<vector<float>>& image, const int startXIndex, const int stopXIndex) {
    /*
    spectra: spectral data
    image: the absolute value of the FFT data. This means the imaginary part is removed.
    startIndex: startIndex of the A-scan. The FFT will be done between startIndex and stopIndex
    stopIndex: last index of the A-scan. The FFT will be done between startIndex and stopIndex

    Do an FFT of only part of the spectra. This function started by several threads to divide the spectrum up
    into different chunks to do the FFT of.
    */

    size_t Nz = spectra[0].size();

    kiss_fft_cfg icfg = kiss_fft_alloc(Nz, 1, NULL, NULL);
    kiss_fft_cpx in[Nz];
    kiss_fft_cpx out[Nz];
    for (int j = 0; j < Nz; j++) {
        in[j].i = 0.0;
    }

    for (int i = startXIndex; i < stopXIndex; i++) {
        for (int j = 0; j < Nz; j++) {
            in[j].r = spectra[i][j];
        }
        kiss_fft(icfg, in, out);
        for (int j = 1; j < Nz+1; j++) {
            image[i][j-1] = (out[j].r*out[j].r + out[j].i*out[j].i)/(Nz*Nz) ;
        }
    }

    kiss_fft_free(icfg);
}



vector<vector<float>> BScan::fftBScan() {
    /*
    This function Fourier transforms the BScan spectra into an OCT image. This
    function uses threads and fftPartBScan to speed up the process.

    */
    if(settings->sizeXSpectrum == 0) {
        return imageFFT;
    }

    static uint64_t stopTime = 0;
    static uint64_t startingTime = getTime();


    int startIndex = 0;
    int stopIndex = 0;
    int NThreads = settings->NThreads;

    vector<thread> threads;

    for(int threadIndex = 0; threadIndex < NThreads; threadIndex++) {

        startIndex = threadIndex  *settings->sizeXSpectrum/NThreads;
        stopIndex = (threadIndex + 1) *settings->sizeXSpectrum/NThreads;
        threads.emplace_back([=] {BScan::fftPartBScan(spectra, imageFFT,startIndex, stopIndex);});
    }

    for(thread& t : threads) {
        if(t.joinable()) {
            t.join();
        }
    }

    return imageFFT;

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
            imageFFT[i][j-1] = (out[j].r*out[j].r + out[j].i*out[j].i)/(settings->sizeZSpectrum*settings->sizeZSpectrum) ;
        }
    }


    kiss_fft_free(icfg);


    return imageFFT;
}

tuple<vector<vector<float>>,int,int> BScan::getProcessedBScan() {
    /*
    This function gets the image as calculated with the RIAA algorithm, together with
    its dimensions and returns it as a tuple.
    */
    tuple<vector<vector<float>>,int,int> output = make_tuple(imageRIAA,settings->sizeXSpectrum,
                             settings->upscalingFactor*settings->sizeZSpectrum);
    return output;
}


void BScan::preprocessSpectrumInPlace() {
    /*
    This function preprocesses the B-scan spectra in place. This means that for changes to the dispersion,
    the spectra need to be reloaded. First the offset is subtracted, then the reference spectrum is compensated,
    then the dispersion correction is made. After the dispersion correction, the signal is apodized, and the
    chirp of the spectrometer is removed. Since the RIAA algorithm does not work well with spectra that go to 0, the edges where the signal is not present or
    very low are removed. Then, in order to continue to have a 2^N bin size, the spectra are interpolated and
    stretched to fill the old array again.
    */
    int i, j;
    if(offset.size() > 0) {
        for(j = 0; j < settings->sizeZSpectrum; j++) {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {

                spectra[i][j] -= offset[j]  ;
            }
        }
        if(referenceSpectrum.size() > 0) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                referenceSpectrum[j] -= offset[j];
            }
        }
    }
    if(referenceSpectrum.size() == 0) {
        referenceSpectrum.reserve(settings->sizeZSpectrum);

        for(j = 0; j < settings->sizeZSpectrum; j++) {
            referenceSpectrum.emplace_back(spectra[0][j]);
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

            double dispersionPhase = j*settings->objectiveDispersionData[settings->objectiveLabel][0];


            for(int k = 1; k < settings->objectiveDispersionData[settings->objectiveLabel].size()-1; k++) {
                dispersionPhase = j*(settings->objectiveDispersionData[settings->objectiveLabel][k] + dispersionPhase);
            }

            if(settings->objectiveDispersionData[settings->objectiveLabel].size() > 1) {
                dispersionPhase +=
                    settings->objectiveDispersionData[settings->objectiveLabel][settings->objectiveDispersionData[settings->objectiveLabel].size()-1];
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
        if(window.size() > 0) {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {
                spectra[i][j] *= window[j]/meanAmplitude;
            }
        } else {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {
                spectra[i][j] /= meanAmplitude;
            }
        }
    }


    if(chirp.size() > 0) {
        for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
            UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(chirp,spectra[i]);
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectra[i][j] = spline.evaluate(  static_cast<float>( j) );
            }
        }
    }



    stretchSpectraInPlace(0.01);

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

void BScan::stretchSpectraInPlace(float minimumReferencePower) {
    /*
    minimumReferencePower: cutoff value for the reference. Lower means more signal is retained, but RIAA might fail, higher
    means less signal is retained, but RIAA works better.

    Since the signal needs to have a box shape for the RIAA to work, the edges with very low signal need to be removed.
    Since this removal causes the array length to be unequal to a power of 2, the spectrum is interpolated, and placed
    into the original array.
    */

    int minIndex, maxIndex;
    float maximumReference = 0.0;
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

    vector<float> xrange;
    xrange.reserve(settings->sizeZSpectrum);
    for (int i = 0; i < settings->sizeZSpectrum; ++i) {
        xrange.push_back(1.0*i / settings->sizeZSpectrum);
    }

    for(int i = 0; i < settings->sizeXSpectrum; i++) {

        UtilityMathFunctions<float>::SplineInterpolation spline = UtilityMathFunctions<float>::SplineInterpolation(xrange,spectra[i]);


        for(int j = 0; j < settings->sizeZSpectrum; j++) {
            float x = 1.0*j*( maxIndex - minIndex ) / (settings->sizeZSpectrum*settings->sizeZSpectrum) + (1.0*minIndex)/settings->sizeZSpectrum;
            spectra[i][j] = spline.evaluate(x);
        }

    }

    return;
}



void BScan::stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower) {
    /*
    spectra: 2D array that contains the spectra which are to be processed.
    referenceSpectrum: array that contains the reference spectrum. This is used to identify low and high signal regions
    minimumReferencePower: cutoff value for the reference. Lower means more signal is retained, but RIAA might fail, higher
    means less signal is retained, but RIAA works better.

    Since the signal needs to have a box shape for the RIAA to work, the edges with very low signal need to be removed.
    Since this removal causes the array length to be unequal to a power of 2, the spectrum is interpolated, and placed
    into the original array.
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

void BScan::FIAALoop(vector<vector<float>>& spectra,BScan* scan,int fromIndex, int toIndex,
                     size_t N, int K, int numberOfPartitions,int numberOfIterations, double vt,
                     vector<float>& startingColumn, vector<vector<float>>& processedImage ) {
    /*
    spectra: the spectra that are to be processed
    scan: the stance that contains all data
    fromIndex: the starting A-scan index of the chunk
    toIndex: the stopping A-scan index of the chunk
    N: the number of A-scans
    K: The length of the A-scans (the size in the z direction)
    numberOfPartitions: the chunks over which each A-scan is chopped
    numberOfIterations: Length of the upscaled Ascan
    vt: noise parameter. Usually set to 1.0.
    startingColumn: Preprocessed starting column. This is an A-scan that is already processed.
    processedImage: store here the RIAA processed image.


    This function is meant to be used with multiple threads processing the image.
    The spectra are first divided of chunks such that multiple threads can work simultaneously. This
    function processes the spectra from fromIndex to toIndex. Since the RIAA algorithm scales quadratically
    with the length of the A-scan, but longer A-scans do not significantly improve the image, quite a speed-up
    can be achieved by processing each A-scan in chunks as well.

    */

    int sign = 1;
    if(fromIndex > toIndex) {
        sign = -1;
    }

    int i = fromIndex;
    for(i = fromIndex; i != toIndex ; i+=sign) {
        for(int j = 0; j < K; j++) {
            if(i == fromIndex) {
                processedImage[i][j] = startingColumn[j];
            } else {
                processedImage[i][j] = processedImage[i-sign][j];
            }
        }
        scan->FIAAPartitioned(spectra[i],N,K,4,numberOfIterations,vt,processedImage[i]);
    }
    cout << "ended at " << i-sign << endl;
}




void BScan::FIAAPartitioned(const vector<float>& x, vector<float>& powerSpectrum,int numberOfIterations) {
    /*

    Wrapper function that loads some relevant data from the global settings.
    */
    BScan::FIAAPartitioned(x,
                           settings->sizeZSpectrum, settings->sizeZSpectrum*settings->upscalingFactor,settings->NChunksRIAA,
                           numberOfIterations, settings->RIAA_NoiseParameter, powerSpectrum );
}

void BScan::FIAAPartitioned(const vector<float>& x,
                            size_t N, int K, int numberOfPartitions,int numberOfIterations, double vt, vector<float>& powerSpectrum ) {
    /*
    x: Initial value
    N: Number of A-scans
    K: length of the A-scan
    numberOfPartitions: The number of partitions of the B-scan
    numberOfIterations: number of iterations
    vt: noise parameter (set to 1.0)
    powerSpectrum: initial value for the power spectrum (square of the FFT). This is the image.


    This function processes the spectra according to the RIAA algorithm
    */

    kiss_fft_cfg cfg = kiss_fft_alloc(N, 0, NULL, NULL);
    kiss_fft_cpx FT[N];
    kiss_fft_cpx signal[N];

    for(int i = 0; i < N; i++) {
        signal[i].r = x[i];
        signal[i].i = 0.0;
    }

    kiss_fft( cfg, signal, FT);
    kiss_fft_free(cfg);
    kiss_fft_cpx partialSignal[N/numberOfPartitions];
    vector<float> partialSignalReal = vector<float>(N/numberOfPartitions);

    kiss_fft_cpx partialFT[N/numberOfPartitions];

    if(powerSpectrum.size() == 0 ) {
        cout << "Populating powerSpectrum vector." << endl;
        cfg = kiss_fft_alloc(K, 0, NULL, NULL);
        kiss_fft_cpx temp[K];
        kiss_fft_cpx powerSpectrumLocal[K];
        powerSpectrum.reserve(K);

        for(int i = 0; i<N; i++) {
            temp[i].r = x[i];
            temp[i].i = 0.0;
        }
        for(int i = N; i<K; i++) {
            temp[i].r = 0.0;
            temp[i].i = 0.0;
        }
        kiss_fft( cfg, temp, powerSpectrumLocal);
        for(int i = 0; i<K; i++) {
            powerSpectrumLocal[i].r = (powerSpectrumLocal[i].r *powerSpectrumLocal[i].r  + powerSpectrumLocal[i].i *powerSpectrumLocal[i].i)/(N*N);
            powerSpectrumLocal[i].i=0;
            powerSpectrum.emplace_back(powerSpectrumLocal[i].r);
        }
        kiss_fft_free(cfg);
    }

    kiss_fft_cfg icfg = kiss_fft_alloc(N/numberOfPartitions, 1, NULL, NULL);

    for(int chunkIndex = 0; chunkIndex < numberOfPartitions; chunkIndex++) {
        //cout << "test " << chunkIndex << "  "<< numberOfPartitions << endl;
        for(int i = 0; i < N/(numberOfPartitions*2); i++) {
            partialFT[i].r = FT[ i + N*chunkIndex/(2*numberOfPartitions)  ].r;
            partialFT[i].i = FT[ i + N*chunkIndex/(2*numberOfPartitions) ].i;
            partialFT[N/(numberOfPartitions) - i - 1].r = FT[N- (i + N*chunkIndex/(2*numberOfPartitions) ) -1].r;
            partialFT[N/(numberOfPartitions) - i - 1].i = FT[N- (i + N*chunkIndex/(2*numberOfPartitions) ) -1].i;
        }

        kiss_fft( icfg, partialFT, partialSignal);

        for(int i = 0; i < N/numberOfPartitions; i++) {
            partialSignalReal[i] = partialSignal[i].r /N ;
        }

        pair<vector<float>, vector<float>> fiaa_output;
        if(chunkIndex != -1) {
            fiaa_output = FIAA(partialSignalReal,  K/numberOfPartitions,numberOfIterations,
                               vt,powerSpectrum, chunkIndex*K/(2*numberOfPartitions) );
        } else {
            for(int i = chunkIndex*K/(numberOfPartitions*2); i < K; i++) {
                powerSpectrum[i] = 0.0;
            }
        }
    }

    kiss_fft_free(icfg);
}

pair<vector<float>,vector<float>> BScan::FIAA(const vector<float>& x, int K, int numberOfIterations, double vt) {

    kiss_fft_cpx powerSpectrumLocal[K];
    kiss_fft_cpx temp[K];
    kiss_fft_cfg cfg = kiss_fft_alloc(K, 0, NULL, NULL);

    vector<float> powerSpectrum;
    powerSpectrum.reserve(K);
    size_t N = x.size();



    for(int i = 0; i<N; i++) {
        temp[i].r = x[i];
        temp[i].i = 0.0;
    }
    for(int i = N; i<K; i++) {
        temp[i].r = 0.0;
        temp[i].i = 0.0;
    }
    kiss_fft( cfg, temp, powerSpectrumLocal);
    for(int i = 0; i<K; i++) {
        powerSpectrum.emplace_back( (powerSpectrumLocal[i].r *powerSpectrumLocal[i].r  + powerSpectrumLocal[i].i *powerSpectrumLocal[i].i)/(N*N) );
    }

    kiss_fft_free(cfg);
    //cout << "Do FIAA with initial FFT power spectrum" << endl;
    return FIAA(x, K, numberOfIterations, vt, powerSpectrum,0);
}

pair<vector<float>,vector<float>> BScan::FIAA(const vector<float>& x, int K, int numberOfIterations, double vt,
vector<float>& powerSpectrum,int powerSpectrumIndex) {
    /*
    x: initial value
    K: length of the A-scan
    numberOfIterations: number of iterations
    vt: noise parameter (set to 1.0)
    powerSpectrum: place the result in this vector.

    */
    int i, j;
    static int evaluated = 0;
    static uint64_t fiaaTime = 0;
    static uint64_t startingTime = getTime();
    size_t N = x.size();


    vector<float> Eta(numberOfIterations + 1);
    float eta = 0.0;
    float af;

    for(i = 0; i < N; i++) {
        eta += abs(x[i]*x[i]);
    }
    eta /= N;


    Eta[0] = eta;
    kiss_fft_cpx powerSpectrumLocal[K];
    for(i = 0; i<K; i++) {
        powerSpectrumLocal[i].r = powerSpectrum[i+powerSpectrumIndex];
        powerSpectrumLocal[i].i=0;
    }


    kiss_fft_cpx temp[K];
    kiss_fft_cpx temp2[K];
    kiss_fft_cpx q[K];
    kiss_fft_cpx Fa1[K];



    vector<complex<float>> c(N);
    kiss_fft_cpx diaa_num[K];
    complex<float> diaa_den[K];
    float diag_a[N];
    vector<complex<float>> A(N);
    vector<complex<float>> y(N);
    vector<complex<float>> fa1(2*N-1);



    uint64_t time0;
    uint64_t time1;
    uint64_t time3;
    uint64_t time4;


    kiss_fft_cfg cfg = kiss_fft_alloc(K, 0, NULL, NULL);
    kiss_fft_cfg icfg = kiss_fft_alloc(K, 1, NULL, NULL);

    for(int k = 0; k < numberOfIterations; k++) {
        evaluated++;
        startingTime = getTime();

        kiss_fft(icfg,powerSpectrumLocal,q);

        for(i = 0; i<N; i++) {
            c[i] = complex<float>(q[i].r, 0.0*q[i].i);
        }
        c[0] += vt*eta;

        tuple<vector<complex<float>>, float> levinsonOut = UtilityMathFunctions<float>::levinson(c,A);

        af = sqrt(get<1>(levinsonOut));
        for(i = 0; i < N; i++) {
            A[i] /= af;
        }

        UtilityMathFunctions<float>::gohberg(A,x,y);


        for(i = 0; i<N; i++) {
            temp[i].r = y[i].real();
            temp[i].i = y[i].imag();
        }
        for(i =  N; i<K; i++) {
            temp[i].r = 0;
            temp[i].i = 0;
        }
        kiss_fft(cfg,temp,diaa_num);

        UtilityMathFunctions<float>::polynomialEstimation(A,fa1);//fa1 has size 2*N-1

        for(i = 0; i < 2*N-1; i++) {
            temp[i].r = fa1[i].real();
            temp[i].i = fa1[i].imag();
        }


        kiss_fft(cfg,temp,Fa1);


        diaa_den[0] = complex<float>(Fa1[0].r, Fa1[0].i);

        for(i = 1; i<K; i++) {
            diaa_den[i] = complex<float>(Fa1[ K - i ].r, Fa1[K - i  ].i);
        }

        for(i = 0; i<K; i++) {
            powerSpectrumLocal[i].r  = abs( (diaa_num[i].r*diaa_num[i].r + diaa_num[i].i*diaa_num[i].i  )  /(diaa_den[i]*conj(diaa_den[i])) );
            powerSpectrumLocal[i].i = 0.0;
        }


        for(i = 1; i < N; i++) {
            temp[i].r = A[N-i].real();
            temp[i].i = -A[N-i].imag();
        }
        temp[0].r = 0.0;
        temp[0].i = 0.0;

        diag_a[0] = A[0].real()*A[0].real() + A[0].imag()*A[0].imag()  -(temp[0].r * temp[0].r +  temp[0].i*temp[0].i);

        for(i = 1; i < N; i++) {
            diag_a[i] = diag_a[i-1] +  A[i].real()*A[i].real() + A[i].imag()*A[i].imag()  -(temp[i].r * temp[i].r +  temp[i].i*temp[i].i);
        }

        eta = 0.0;
        for(i = 0; i < N; i++) {
            eta +=  abs(y[i]*y[i] / (diag_a[i]*diag_a[i])) ;
        }
        eta /= N;
        Eta[k+1] = eta;
        fiaaTime += getTime() - startingTime;

    }


    kiss_fft_free(cfg);
    kiss_fft_free(icfg);

    pair<vector<float>, vector<float>> result;

    for(i = 0; i < K; i++) {
        powerSpectrum[i+powerSpectrumIndex] = powerSpectrumLocal[i].r;
    }

    result.first = powerSpectrum;
    result.second = Eta;
    return result;
}



vector<vector<float>> BScan::processBScan() {
    /*
    Start the RIAA processing of the BScan with the settings saved in the global settings.
    */
    return BScan::processBScan(settings->sizeXSpectrum,settings->sizeZSpectrum,
                               settings->upscalingFactor *settings->sizeZSpectrum,
                               settings->initialNumberOfIterations,settings->numberOfIterations,
                               settings->RIAA_NoiseParameter,settings->NThreads);
}



vector<vector<float>> BScan::processBScan(size_t M,const size_t N, int K,int numberOfIterationsFirstColumn,int numberOfIterations,
double vt,int NThreads) {

    /*
    M: number of A-scans
    N: length of A-scan
    K: length of A-scan including upscaling
    numberOfIterationsFirstColumn: number of iterations for the first A-scan (usually about 15)
    numberOfIterations: number of iterations (can be lower than 5)
    vt: noise parameter, set to about 1.0
    NThreads: number of threads.

    This function creates a number of threads to process parts of the B-scan. For each chunk, one
    a-scan is analyzed first with numberOfIterationsFirstColumn number of iterations. Then the next A-scan in the chunk
    is analyzed with a much lower number of iterations, but the initial value is the previous A-scan,
    or, if it is the second A-scan, the A-scan that is iterated over numberOfIterationsFirstColumn times.
    */

    int i, j;
    uint64_t time0;
    uint64_t time1;
    pair<vector<float>,vector<float>> fiaa_output;

    time0 = getTime();

    for(int i = 0; i < M; i++) {
        if( (i) % (M/NThreads) == 0 && i != 0) {
            //cout << "do fiaa " << endl;
            fiaa_output = FIAA(spectra[i],K,numberOfIterationsFirstColumn,vt);
            //cout << "swap " << endl;
            imageRIAA[i].swap(fiaa_output.first);
            //cout << "swapped." << endl;
        }
    }

    int startIndex;
    int stopIndex;

    cout << "startIndex " << startIndex << " stopIndex " << stopIndex << endl;

    vector<thread> threads;

    for(int threadIndex = 0; threadIndex < NThreads; threadIndex++) {

        if(threadIndex % 2 == 0) {


            if(threadIndex >= NThreads -1 ) {
                startIndex = M-1;
            } else {
                startIndex = (threadIndex+1) * (M/NThreads);
            }

            stopIndex = threadIndex * (M/NThreads);
            cout << "backwards " << startIndex << "  " << stopIndex << endl;
            threads.emplace_back(FIAALoop,ref(spectra),this,startIndex-1,stopIndex-1,N,K,4,numberOfIterations,vt,ref(imageRIAA[startIndex]),
                                 ref(imageRIAA));
        } else {
            startIndex = (threadIndex) * (M/NThreads);

            if(threadIndex >= NThreads -1 ) {
                stopIndex =  M-1;
            } else {
                stopIndex = (threadIndex+1) * (M/NThreads)-1;
            }
            cout << "forwards " << startIndex << "  " << stopIndex << endl;
            threads.emplace_back(FIAALoop, ref(spectra),this,startIndex+1,stopIndex+1,N,K,4,numberOfIterations,vt,ref(imageRIAA[startIndex]),
                                 ref(imageRIAA));
        }
    }



    for(thread& t : threads) {
        if(t.joinable()) {
            t.join();
        }
    }



    time1 = getTime();
    cout << "done in " << 1e-3*(time1 - time0) << " s." << endl;

    return imageRIAA;
}



BScan::~BScan() {

}
