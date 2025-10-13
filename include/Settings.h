#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <thread>
#include "Tag.h"
//#include "IO.h"



class Settings {
public:
    Settings();
    Settings(list<Tag*>);
    Settings(string);
    Settings(Settings& other);
    Settings& operator=(const Settings& other);
    void writeToYAML(string);

    virtual ~Settings();
    //TODO: create helper functions to be able to make these variables consts.

    int NChunksEnvelopeSubtraction = 64;
    int NThreads = -1;
    int bytesPerPixelIntensity;
    int bytesPerPixelSpectrum;
    int bytesPerPixelChirp;
    int bytesPerPixelApodization;
    int bytesPerPixelOffset;
    int sizeXIntensity;
    int sizeZIntensity;
    int sizeXSpectrumRaw;
    int sizeXSpectrum;
    int sizeZSpectrum;
    int sizeXChirp;
    int sizeZChirp;
    int sizeXApodization;
    int sizeZApodization;
    int sizeXOffset;
    int sizeZOffset;

    int initialNumberOfIterations = 15;
    int numberOfIterations = 3;
    int upscalingFactor = 8;
    int NChunksRIAA = 16;

    int spectrumAveraging;
    int numberOfAScans;
    int numberOfBScans;
    int x_px;
    int z_px;
    double x_mm;
    double z_mm;

    double wavelength_nm;
    double bandwidth_nm;
    double refractiveIndex;
    double referenceIntensity;
    double electronCountScaling;
    double RIAA_NoiseParameter = 1.0;

    size_t numberOfDispersionCoefficients = 1;
    double* dispersionCoefficients = nullptr;
    map<string, vector<double>> objectiveDispersionData;
    string objectiveLabel = "No correction";
    string pathToExecutable;
    string pathToData;

    map<int, string> pathsSpectra;
    map<int, int> scanStartIndices;
    string pathIntensity;
    string pathChirp;
    string pathApodization;
    string pathOffset;
    int colorMap = 1;
    bool useDecibelColorScale = true;
    double decibelFloor = -120;
    double percentageFloor = 20.0;
    double decibelCeil = 0;
    double percentageCeil= 100;
    bool useRIAA = false;
    bool alwaysRedoPreprocessing = false;




protected:
private:
    void copyMembers(const Settings&);

};

#endif // SETTINGS_H
