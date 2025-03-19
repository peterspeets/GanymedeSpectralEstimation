#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "Tag.h"
//#include "IO.h"

class Settings
{
    public:
        Settings();
        Settings(list<Tag*>);
        Settings(string);
        virtual ~Settings();
        //TODO: create helper functions to be able to make these variables consts.

        int NChunksEnvelopeSubtraction = 64;
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

        map<int, string> pathsSpectra;
        map<int, int> scanStartIndices;
        string pathIntensity;
        string pathChirp;
        string pathApodization;
        string pathOffset;




    protected:
    private:
};

#endif // SETTINGS_H
