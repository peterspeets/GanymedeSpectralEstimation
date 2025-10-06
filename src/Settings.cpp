#include "Settings.h"
#include "IO.h"

Settings::Settings() {
    return;
}




Settings::Settings(string pathToDirectory) {


    bytesPerPixelIntensity = 1;
    bytesPerPixelSpectrum = 1;
    bytesPerPixelChirp = 1;
    bytesPerPixelApodization = 1;
    bytesPerPixelOffset = 1;

    pair<double*, int> sizePair = IO<double>::loadArrayFromFile(pathToDirectory + "Size.txt");
    pair<float*, int> referenceSpectrumPair = IO<float>::loadArrayFromFile(pathToDirectory + "sk.txt");
    tuple<float**, int,int> spectraTuple = IO<float>::load2DArrayFromFile(pathToDirectory + "rawData.txt");

    sizeXIntensity = get<1>(spectraTuple);
    sizeZIntensity= get<2>(spectraTuple);
    sizeXSpectrumRaw= get<1>(spectraTuple);
    sizeXSpectrum= get<1>(spectraTuple);
    sizeZSpectrum= get<2>(spectraTuple);
    sizeXChirp= get<1>(spectraTuple);
    sizeZChirp= get<2>(spectraTuple);;
    sizeXApodization= get<1>(spectraTuple);
    sizeZApodization= get<2>(spectraTuple);
    sizeXOffset= get<1>(spectraTuple);
    sizeZOffset= get<2>(spectraTuple);


    numberOfAScans = get<1>(spectraTuple);
    x_px = get<1>(spectraTuple);
    z_px = get<2>(spectraTuple);

    spectrumAveraging  = 1;

    numberOfBScans = 1;

    x_mm = sizePair.first[0];
    z_mm = sizePair.first[1];

    wavelength_nm = 930.0;
    bandwidth_nm = 100.0;
    refractiveIndex = 1.33;
    referenceIntensity = 1.0;
    electronCountScaling = 1.0;


    pathsSpectra[0] = pathToDirectory + "rawData.txt";

    scanStartIndices[0] = 0;

    pathIntensity = pathToDirectory + "sk.txt";
    pathChirp = pathToDirectory + "";
    pathApodization= pathToDirectory + "sk.txt";
    pathOffset = pathToDirectory + "";

    if(NThreads < 0) {
        NThreads = thread::hardware_concurrency()/2;
    }


}

Settings::Settings(list<Tag*> tags)  {



    for(Tag* tag : tags) {
        if(tag->label == "DataFile" && tag->content.substr(0,13) == "data/Spectral" && tag->isOpeningTag) {
            bytesPerPixelSpectrum = stod(tag->keys["BytesPerPixel"]) ;

            sizeXSpectrumRaw = stoi(tag->keys["SizeX"]);
            sizeZSpectrum = stoi(tag->keys["SizeZ"]);

            bool isInList = false;
            for (const pair<int, string>& pair : pathsSpectra) {

                if(pair.second == tag->content ) {
                    isInList  = true;
                    break;
                }
            }
            if(!isInList) {
                int spectrumNumber = stoi(tag->content.substr(13,tag->content.size()-18));
                pathsSpectra[spectrumNumber]=tag->content;
                scanStartIndices[spectrumNumber]=stoi(tag->keys["ScanRegionStart0"]);
            }
            continue;
        }
        if(tag->label == "DataFile" && tag->content == "data/Intensity.data" && tag->isOpeningTag) {
            bytesPerPixelIntensity = stod(tag->keys["BytesPerPixel"]) ;
            sizeXIntensity = stoi(tag->keys["SizeX"]);
            sizeZIntensity = stoi(tag->keys["SizeZ"]);
            pathIntensity = tag->content;
            continue;
        } else if(tag->label == "DataFile" && tag->content == "data/Chirp.data" && tag->isOpeningTag) {
            bytesPerPixelChirp= stod(tag->keys["BytesPerPixel"]) ;
            sizeZChirp = stoi(tag->keys["SizeZ"]);
            pathChirp = tag->content;
            continue;
        } else if(tag->label == "DataFile" && tag->content == "data/ApodizationSpectrum.data" && tag->isOpeningTag) {
            bytesPerPixelApodization= stod(tag->keys["BytesPerPixel"]) ;
            sizeZApodization= stoi(tag->keys["SizeZ"]);
            pathApodization = tag->content;
            continue;
        } else if(tag->label == "DataFile" && tag->content == "data/OffsetErrors.data" && tag->isOpeningTag) {
            bytesPerPixelOffset= stod(tag->keys["BytesPerPixel"]) ;
            sizeZOffset = stoi(tag->keys["SizeZ"]);
            pathOffset= tag->content;
            continue;
        } else if(tag->label == "IntensityAveraging" && tag->isOpeningTag) {
            for(Tag* innerTag : tag->innerTags) {
                if(innerTag->label == "Spectra" && innerTag->isOpeningTag ) {
                    spectrumAveraging = stoi(innerTag->content);
                    continue;
                } else if(innerTag->label == "AScans" && innerTag->isOpeningTag ) {
                    numberOfAScans = stoi(innerTag->content);
                    continue;
                } else if(innerTag->label == "BScans" && innerTag->isOpeningTag ) {
                    numberOfBScans = stoi(innerTag->content);
                    continue;
                }

            }
            continue;
        } else if(tag->label == "Image" && tag->isOpeningTag) {

            for(Tag* innerTagImage : tag->innerTags) {
                if(innerTagImage->label == "SizeReal" && innerTagImage->isOpeningTag ) {
                    for(Tag* innerTagSizeReal: innerTagImage->innerTags) {
                        if(innerTagSizeReal->label == "SizeX" && innerTagSizeReal->isOpeningTag ) {
                            x_mm = stod(innerTagSizeReal->content);
                            continue;
                        } else if(innerTagSizeReal->label == "SizeZ" && innerTagSizeReal->isOpeningTag ) {
                            z_mm = stod(innerTagSizeReal->content);
                            continue;
                        }
                    }

                    continue;
                }


                if(innerTagImage->label == "SizePixel" && innerTagImage->isOpeningTag ) {
                    for(Tag* innerTagSizePixel: innerTagImage->innerTags) {
                        if(innerTagSizePixel->label == "SizeX" && innerTagSizePixel->isOpeningTag ) {
                            x_px = stoi(innerTagSizePixel->content);
                            continue;
                        } else if(innerTagSizePixel->label == "SizeZ" && innerTagSizePixel->isOpeningTag ) {
                            z_px = stoi(innerTagSizePixel->content);
                            continue;
                        }
                    }

                    continue;
                }

            }
            continue;
        } else if(tag->label == "Acquisition" && tag->isOpeningTag) {
            for(Tag* innerTagAcquisition : tag->innerTags) {
                if(innerTagAcquisition->label == "RefractiveIndex" && innerTagAcquisition->isOpeningTag ) {
                    refractiveIndex = stod(innerTagAcquisition->content);
                    continue;
                } else if(innerTagAcquisition->label == "ReferenceIntensity" && innerTagAcquisition->isOpeningTag ) {
                    referenceIntensity = stod(innerTagAcquisition->content);
                    continue;
                }

            }
        } else if(tag->label == "Instrument" && tag->isOpeningTag) {
            for(Tag* innerTagInstrument: tag->innerTags) {
                if(innerTagInstrument->label == "SourceBandwidth" && innerTagInstrument->isOpeningTag ) {
                    bandwidth_nm = stod(innerTagInstrument->content);
                    continue;
                } else if(innerTagInstrument->label == "CentralWavelength" && innerTagInstrument->isOpeningTag ) {
                    wavelength_nm = stod(innerTagInstrument->content);
                    continue;
                } else if(innerTagInstrument->label == "BinaryToElectronCountScaling" && innerTagInstrument->isOpeningTag ) {
                    electronCountScaling = stod(innerTagInstrument->content);
                    continue;
                }

            }
        }
    }

    sizeXSpectrum = sizeXSpectrumRaw-scanStartIndices[0];

    if(NThreads < 0) {
        NThreads = thread::hardware_concurrency()/2;
    }
}

void Settings::copyMembers(const Settings& other){



    numberOfDispersionCoefficients = other.numberOfDispersionCoefficients;
    if(dispersionCoefficients){
        dispersionCoefficients = new double[numberOfDispersionCoefficients];
        for(int i = 0; i < numberOfDispersionCoefficients ; i++){
            dispersionCoefficients[i] = other.dispersionCoefficients[i];
        }
    }else{
        dispersionCoefficients = nullptr;
    }


    NChunksEnvelopeSubtraction = other.NChunksEnvelopeSubtraction;
    NThreads = other.NThreads;
    bytesPerPixelIntensity = other.bytesPerPixelIntensity;
    bytesPerPixelSpectrum = other.bytesPerPixelSpectrum;
    bytesPerPixelChirp= other.bytesPerPixelChirp;
    bytesPerPixelApodization= other.bytesPerPixelApodization;
    bytesPerPixelOffset= other.bytesPerPixelOffset;
    sizeXIntensity= other.sizeXIntensity;
    sizeZIntensity= other.sizeZIntensity;
    sizeXSpectrumRaw= other.sizeXSpectrumRaw;
    sizeXSpectrum= other.sizeXSpectrum;
    sizeZSpectrum= other.sizeZSpectrum;
    sizeXChirp= other.sizeXChirp;
    sizeZChirp= other.sizeZChirp;
    sizeXApodization= other.sizeXApodization;
    sizeZApodization= other.sizeZApodization;
    sizeXOffset= other.sizeXOffset;
    sizeZOffset= other.sizeZOffset;

    initialNumberOfIterations = other.initialNumberOfIterations;
    numberOfIterations = other.numberOfIterations;
    upscalingFactor = other.upscalingFactor;
    NChunksRIAA = other.NChunksRIAA;

    spectrumAveraging = other.spectrumAveraging;
    numberOfAScans= other.numberOfAScans;
    numberOfBScans= other.numberOfBScans;
    x_px= other.x_px;
    z_px= other.z_px;
    x_mm= other.x_mm;
    z_mm= other.z_mm;



    wavelength_nm= other.wavelength_nm;
    bandwidth_nm= other.bandwidth_nm;
    refractiveIndex= other.refractiveIndex;
    referenceIntensity= other.referenceIntensity;
    electronCountScaling= other.electronCountScaling;
    RIAA_NoiseParameter = other.RIAA_NoiseParameter;

    numberOfDispersionCoefficients = other.numberOfDispersionCoefficients;

    pathsSpectra= other.pathsSpectra;
    scanStartIndices= other.scanStartIndices;
    pathIntensity= other.pathIntensity;
    pathChirp= other.pathChirp;
    pathApodization= other.pathApodization;
    pathOffset= other.pathOffset;


    return;

}


Settings::Settings(Settings& other) {

    copyMembers(other);

}


Settings& Settings::operator=(const Settings& other) {

    copyMembers(other);

    return *this;
}


Settings::~Settings() {
    if(dispersionCoefficients){
        delete[] dispersionCoefficients;
    }
}
