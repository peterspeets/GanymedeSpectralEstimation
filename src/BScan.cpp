#include "BScan.h"

BScan::BScan(){

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


BScan::BScan(const string filePath)
{
    shared_ptr<Settings> oldSettings = settings;

    if(filePath.substr(filePath.length() - 4) == ".oct") {
        IO<float>::GanymedeFileLoader fileLoader =  IO<float>::GanymedeFileLoader(filePath);
        settings = fileLoader.loadSettings() ;
        offset = fileLoader.loadCalibrationSpectrum(settings->pathOffset);
        chirp = fileLoader.loadCalibrationSpectrum(settings->pathChirp);
        referenceSpectrum = fileLoader.loadCalibrationSpectrum(settings->pathApodization);
        intensity = fileLoader.loadCalibrationSpectrum(settings->pathIntensity);

        settings->numberOfStoredSpectra = settings->pathsSpectra.size();

        float** singleAScan;
        spectra = new float*[settings->sizeXSpectrum];

        for(int i = 0; i< settings->sizeXSpectrum; i++) {
            spectra[i] = new float[settings->sizeZSpectrum];
        }



        cout << "paths to spectra" << endl;
        for (pair<int,string> const& pathSpectraPair : settings->pathsSpectra){
            if(pathSpectraPair.first  == 0){
                singleAScan = fileLoader.loadSpectrum(pathSpectraPair.first);
                for(int i = 0; i< settings->sizeXSpectrum; i++) {
                    for(int j = 0; j< settings->sizeZSpectrum; j++) {
                        spectra[i][j] = singleAScan[i][j]/(1.0*settings->numberOfStoredSpectra);
                    }
                }
            }else{
                singleAScan = fileLoader.loadSpectrum(pathSpectraPair.first);

                for(int i = 0; i< settings->sizeXSpectrum; i++) {
                    for(int j = 0; j< settings->sizeZSpectrum; j++) {
                        spectra[i][j] += singleAScan[i][j]/(1.0*settings->numberOfStoredSpectra);
                    }
                }
            }

            cout <<pathSpectraPair.first << " " << pathSpectraPair.second << endl;
            for(int i = 0; i< settings->sizeXSpectrum; i++) {
                delete[] singleAScan[i];
            }
            delete[] singleAScan;
        }
        cout << "..." << endl;





        if (settings->objectiveDispersionData.count("Native")){
            settings->objectiveDispersionData.erase("Native");
        }

        settings->objectiveLabel = oldSettings->objectiveLabel;


    } else {
        string directoryPath = filePath;
        if(filePath.substr(filePath.length() - 4) == ".txt"){
            for(int i = filePath.length()-1; i > 0; i--){
                if(filePath[i] == '\\' || filePath[i] == '/'){
                    directoryPath = filePath.substr(0,i+1);
                    break;
                }
            }
        }

        settings = make_shared<Settings>(directoryPath );
        tuple<float**,int,int> loadingSpectraTuple =  IO<float>::load2DArrayFromFile(directoryPath + "rawData.txt");
        pair<float*,int> loadingReferencePair =  IO<float>::loadArrayFromFile(directoryPath  + "sk.txt");
        spectra = get<0>(loadingSpectraTuple);
        //referenceSpectrum = loadingReferencePair.first;
        pair<double*,int> loadingDispersionPair =  IO<double>::loadArrayFromFile(directoryPath  + "phase.txt");



        bool hasNativeDispersionCorrection = false;
        for(pair<string,vector<double>> const objectiveDispersionDataPair : oldSettings->objectiveDispersionData){
            if(objectiveDispersionDataPair.first == "Native"){
                cout << "Found native label." << endl;
                hasNativeDispersionCorrection = true;
                break;
            }
        }

        settings->objectiveLabel = oldSettings->objectiveLabel;
        if(!hasNativeDispersionCorrection){
            cout << "Set to native, only once." << endl;
            settings->objectiveLabel = "Native";
        }
        settings->objectiveDispersionData["Native"] = vector<double>(loadingDispersionPair.first,
                                                                                     loadingDispersionPair.first + loadingDispersionPair.second);

        if(settings->objectiveLabel == "Native"){
            cout << "Settings already where on native.";
            settings->dispersionCoefficients = loadingDispersionPair.first;
            settings->numberOfDispersionCoefficients = loadingDispersionPair.second;
        }

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
    imageRIAA = new float*[settings->sizeXSpectrum];
    fftBScan();
    cout << "Process B scan" << endl;
    processBScan();
    cout << "objective label "<< settings->objectiveLabel << endl;
}

uint64_t BScan::getTime() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}


void BScan::fftPartBSscan(float** spectra, float** image, int Nz,int startXIndex, int stopXIndex) {
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


float** BScan::fftBScan(){
    if(settings->sizeXSpectrum == 0){
        return imageFFT;
    }

    if(!imageFFT ){
        imageFFT = new float*[settings->sizeXSpectrum];
        for (int i = 0; i < settings->sizeXSpectrum; i++) {
            imageFFT[i] = new float[settings->sizeZSpectrum];
        }
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
        threads.emplace_back(fftPartBSscan,spectra, imageFFT, settings->sizeZSpectrum,startIndex,stopIndex);
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

tuple<float**,int,int> BScan::getProcessedBScan(){
    tuple<float**,int,int> output = make_tuple(imageRIAA,settings->sizeXSpectrum,settings->upscalingFactor*settings->sizeZSpectrum);
    return output;
}


void BScan::preprocessSpectrumInPlace() {
    int i, j;
    if(offset) {
        for(j = 0; j < settings->sizeZSpectrum; j++) {
            for(i = 0; i < settings->sizeXSpectrum ; i++) {

                spectra[i][j] -= offset[j]  ;
            }
        }
        if(referenceSpectrum) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                referenceSpectrum[j] -= offset[j];
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
            //cout << "dispersionPhase 0 "<<dispersionPhase  <<endl;

            for(int k = 1; k < settings->numberOfDispersionCoefficients-1; k++) {
                dispersionPhase = j*(settings->dispersionCoefficients[k] + dispersionPhase);
            }
            //cout << "dispersionPhase 1 "<<dispersionPhase  <<endl;
            if(settings->numberOfDispersionCoefficients > 1) {
                dispersionPhase += settings->dispersionCoefficients[settings->numberOfDispersionCoefficients-1];
            }
            //cout << "dispersionPhase 2 "<<dispersionPhase  <<endl;
            //cout << "dispersion phase: " << dispersionPhase <<endl;

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
        for(i = 0; i < settings->sizeXSpectrum -0 ; i++) {
            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(chirp,spectra[i],settings->sizeZSpectrum);
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectra[i][j] = spline->evaluate(  static_cast<float>( j) );
            }
            delete spline;
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


void BScan::stretchSpectraInPlace(float** spectra, float* referenceSpectrum, float minimumReferencePower) {
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
    for(maxIndex = indexOfMaximum; referenceSpectrum[maxIndex] > minimumReferencePower*maximumReference && maxIndex < settings->sizeZSpectrum; maxIndex++) {}

    float xrange[settings->sizeZSpectrum];
    for(int i = 0; i < settings->sizeZSpectrum; i++) {
        xrange[i] = 1.0*i/settings->sizeZSpectrum;
    }

    for(int i = 0; i < settings->sizeXSpectrum; i++) {

        UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(xrange,spectra[i],settings->sizeZSpectrum);
        for(int j = 0; j < settings->sizeZSpectrum; j++) {
            float x = 1.0*j*( maxIndex - minIndex ) / (settings->sizeZSpectrum*settings->sizeZSpectrum) + (1.0*minIndex)/settings->sizeZSpectrum;
            spectra[i][j] = spline->evaluate(x);
        }


        delete spline;
    }

    return;
}

void BScan::fiaa_oct_loop(float** spectra,BScan* scan,int fromIndex, int toIndex,
        size_t N, int K, int numberOfPartitions,int q_i, double vt, float* startingColumn,float** processedImage ) {
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
        scan->fiaa_oct_partitioned(spectra[i],N,K,4,q_i,vt,processedImage[i]);
    }
    cout << "ended at " << i-sign << endl;
}

void BScan::fiaa_oct_partitioned(const float* x, float* diaaf_floatingPoint,int numberOfIterations){
    BScan::fiaa_oct_partitioned(x,
        settings->sizeZSpectrum, settings->sizeZSpectrum*settings->upscalingFactor,settings->NChunksRIAA ,
        numberOfIterations, settings->RIAA_NoiseParameter, diaaf_floatingPoint );
}

void BScan::fiaa_oct_partitioned(const float* x,
        size_t N, int K, int numberOfPartitions,int q_i, double vt, float* diaaf_floatingPoint ) {

    //TODO: limit scope of stack intensive arrays, before calling FIAA function.
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
    float partialSignalReal[N/numberOfPartitions];
    kiss_fft_cpx partialFT[N/numberOfPartitions];

    if(!diaaf_floatingPoint ) {
        cout << "Creating new OPL array." << endl;
        cfg = kiss_fft_alloc(K, 0, NULL, NULL);
        kiss_fft_cpx temp[K];
        kiss_fft_cpx diaaf[K];
        diaaf_floatingPoint = new float[K];

        for(int i = 0; i<N; i++) {
            temp[i].r = x[i];
            temp[i].i = 0.0;
        }
        for(int i = N; i<K; i++) {
            temp[i].r = 0.0;
            temp[i].i = 0.0;
        }
        kiss_fft( cfg, temp, diaaf);
        for(int i = 0; i<K; i++) {
            diaaf[i].r = diaaf_floatingPoint[i] = (diaaf[i].r *diaaf[i].r  + diaaf[i].i *diaaf[i].i)/(N*N);
            diaaf[i].i=0;
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

        pair<float*, float*> fiaa_output;


        if(chunkIndex != -1) {
            fiaa_output = fiaa_oct(partialSignalReal, N/numberOfPartitions,  K/numberOfPartitions,q_i,
                                   vt,&diaaf_floatingPoint[chunkIndex*K/(2*numberOfPartitions)] );
            delete[] fiaa_output.second;
        } else {
            for(int i = chunkIndex*K/(numberOfPartitions*2); i < K; i++) {
                diaaf_floatingPoint[i] = 0.0;
            }
        }
    }

    kiss_fft_free(icfg);
}




pair<float*, float*> BScan::fiaa_oct(const float* x,
        size_t N, int K, int q_i, double vt, float* diaaf_floatingPoint ) {
    int i, j;

    ostringstream  filename;
    static int evaluated = 0;
    static uint64_t fiaaTime = 0;
    static uint64_t startingTime = getTime();



    float* Eta = new float[q_i+1];
    float eta = 0.0;
    float af;
    for(i = 0; i < N; i++) {
        eta += abs(x[i]*x[i]);
    }
    eta /= N;


    Eta[0] = eta;
    kiss_fft_cpx diaaf[K];



    kiss_fft_cpx temp[K];
    kiss_fft_cpx temp2[K];
    kiss_fft_cpx q[K];
    kiss_fft_cpx Fa1[K];



    complex<float> c[N];
    kiss_fft_cpx diaa_num[K];
    complex<float> diaa_den[K];
    float diag_a[N];
    complex<float>* A = new complex<float>[N];
    complex<float>* y = new complex<float>[N];
    complex<float>* fa1 = new complex<float>[2*N-1] ;



    uint64_t time0;
    uint64_t time1;
    uint64_t time3;
    uint64_t time4;


    kiss_fft_cfg cfg = kiss_fft_alloc(K, 0, NULL, NULL);
    kiss_fft_cfg icfg = kiss_fft_alloc(K, 1, NULL, NULL);


    if(!diaaf_floatingPoint ) {

        diaaf_floatingPoint = new float[K];

        for(i = 0; i<N; i++) {
            temp[i].r = x[i];
            temp[i].i = 0.0;
        }
        for(i = N; i<K; i++) {
            temp[i].r = 0.0;
            temp[i].i = 0.0;
        }
        kiss_fft( cfg, temp, diaaf);
        for(i = 0; i<K; i++) {
            diaaf[i].r = (diaaf[i].r *diaaf[i].r  + diaaf[i].i *diaaf[i].i)/(N*N);
            diaaf[i].i=0;
        }

    } else {

        for(i = 0; i<K; i++) {
            diaaf[i].r = diaaf_floatingPoint[i];
            diaaf[i].i=0;
        }
    }




    for(int k = 0; k < q_i; k++) {
        evaluated++;
        startingTime = getTime();

        kiss_fft(icfg,diaaf,q);


        for(i = 0; i<N; i++) {
            c[i] = complex<float>(q[i].r, 0.0*q[i].i);
        }
        c[0] += vt*eta;

        tuple<complex<float>*, float> levinsonOut = UtilityMathFunctions<float>::levinson(c,N,A);

        af = sqrt(get<1>(levinsonOut));

        for(i = 0; i < N; i++) {
            A[i] /= af;
        }

        UtilityMathFunctions<float>::tvec_gs_i(A,x,N,y);



        for(i = 0; i<N; i++) {
            temp[i].r = y[i].real();
            temp[i].i = y[i].imag();
        }
        for(i =  N; i<K; i++) {
            temp[i].r = 0;
            temp[i].i = 0;
        }
        kiss_fft(cfg,temp,diaa_num);



        UtilityMathFunctions<float>::polynomialEstimation(A,N,fa1);//fa1 has size 2*N-1


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
            diaaf[i].r  = abs( (diaa_num[i].r*diaa_num[i].r + diaa_num[i].i*diaa_num[i].i  )  /(diaa_den[i]*conj(diaa_den[i])) );
            diaaf[i].i = 0.0;
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
        Eta[k] = eta;

        fiaaTime += getTime() - startingTime;
    }


    kiss_fft_free(cfg);
    kiss_fft_free(icfg);
    delete[] A;
    delete[] y;
    delete[] fa1;
    pair<float*, float*> result;

    for(i = 0; i < K; i++) {
        diaaf_floatingPoint[i] = diaaf[i].r;
    }

    result.first = diaaf_floatingPoint;
    result.second = Eta;
    return result;
}


float** BScan::processBScan(){
    return BScan::processBScan(settings->sizeXSpectrum,settings->sizeZSpectrum,
                               settings->upscalingFactor *settings->sizeZSpectrum,
                               settings->initialNumberOfIterations,settings->numberOfIterations,
                               settings->RIAA_NoiseParameter,settings->NThreads);
}


calculateLowResBitmap(){
}





float** BScan::processBScan(size_t M,const size_t N, int K,int q_init,int q_i, double vt,int NThreads) {
    int i, j;


    for(i = 0 ; i < M; i++){
        if(!imageRIAA[i]){
            delete[] imageRIAA[i];
        }

    }


    uint64_t time0;
    uint64_t time1;
    pair<float*, float*> fiaa_output;

    time0 = getTime();

    for(int i = 0; i < M; i++) {
        if( (i) % (M/NThreads) == 0 && i != 0) {
            fiaa_output = fiaa_oct(spectra[i],N,K,q_init,vt);
            imageRIAA[i] = fiaa_output.first;
            delete[] fiaa_output.second;
        } else {

            imageRIAA[i] = new float[K];
        }
    }

    int startIndex;
    int stopIndex;


    vector<thread> threads;

    for(int threadIndex = 0; threadIndex < NThreads; threadIndex++) {

        if(threadIndex % 2 == 0) {


            if(threadIndex >= NThreads -1 ){
                startIndex = M-1;
            }else{
                startIndex = (threadIndex+1) * (M/NThreads);
            }

            stopIndex = threadIndex * (M/NThreads);
            cout << "backwards " << startIndex << "  " << stopIndex << endl;
            threads.emplace_back(fiaa_oct_loop,spectra,this,startIndex-1,stopIndex-1,N,K,4,q_i,vt,imageRIAA[startIndex],imageRIAA);
        } else {
            startIndex = (threadIndex) * (M/NThreads);

            if(threadIndex >= NThreads -1 ){
                stopIndex =  M-1;
            }else{
                stopIndex = (threadIndex+1) * (M/NThreads)-1;
            }
            cout << "forwards " << startIndex << "  " << stopIndex << endl;
            threads.emplace_back(fiaa_oct_loop,spectra,this,startIndex+1,stopIndex+1,N,K,4,q_i,vt,imageRIAA[startIndex],imageRIAA);
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



BScan::~BScan()
{  /*Todo, change class logic such that spectra are deleted here.*/
}
