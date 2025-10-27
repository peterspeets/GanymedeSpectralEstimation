#ifndef IO_H
#define IO_H
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <map>


#include <miniz.h>
#include <kiss_fft.h>
#include "Tag.h"
#include <png.h>

#include "UtilityMathFunctions.h"
#include "XML_Interpreter.h"
#include "Settings.h"



using namespace std;

class Settings;

extern shared_ptr<Settings> settings;




template <typename floatingPointType>
class IO {
public:
    IO();
    virtual ~IO();
    template <typename T>
    static void saveArrayToFile(const T* array, const int N, const string& filename);
    template <typename T>
    static void saveArrayToFile(const complex<T>* cpx, const int N, const string& filename);
    static void saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename);

    template <typename T>
    static void saveVectorToFile(const vector<T>& vectorToSave,const string& filePath);
    template <typename T>
    static void saveVectorToFile(const vector<complex<T>>& vectorToSave,const string& filePath);
    template <typename T>
    static void save2DVectorToFile(const vector<vector<T>>& vectorToSave,const string& filePath, char separator = ',');

    template <typename T>
    static void save2DArrayToFile(T** array, const int M, const int N, const std::string& filename, char separator = ',');
    template <typename T>
    static void save2DArrayToFile(T** array, const int M,const int N, const int newM, const int newN,const string& filename,char separator);

    template <typename T>
    static void test(T** array, int);

    static void savePng(const string filename, const int width, const int height, const unsigned char* image);
    static void savePng(const string filename, const int width, const int height, const floatingPointType* const* image,
                        bool decibelColorScale = false,double decibelFloor = -120 );
    static void savePng(const string filename, const int width, const int height, const int imageWidth,
                        const int imageHeight, const floatingPointType* const* image, bool decibelColorScale = false,double decibelFloor = -120 );

    static pair<floatingPointType*, int> loadArrayFromFile(const string& filename);
    static tuple<floatingPointType**, int, int> load2DArrayFromFile(const string& filename);
    static map<string, vector<double>> loadObjectiveDispersionData(const string& filename);

    static vector<floatingPointType> loadVectorFromFile(const string& filePath);
    static vector<vector<floatingPointType>> load2DVectorFromFile(const string& filePath);


    class FileLoader {
    public:
        FileLoader(const string filePath);
        ~FileLoader();
    };

    class GanymedeFileLoader {
    public:

        shared_ptr<Settings> loadSettings() ;
        floatingPointType** loadSpectrum(int spectrumIndex, floatingPointType* referenceSpectrum) ;
        floatingPointType** loadSpectrum(int spectrumIndex) ;
        floatingPointType* loadCalibrationSpectrum( string spectrumFileName,int arrayLength = 2048,int bytesPerPixel = 4 );
        GanymedeFileLoader(const string filePath);
        ~GanymedeFileLoader();

    private:
        mz_zip_archive OCTArchive;
        string loadHeader();
        int getFileIndex( string fileName);

    };


protected:

private:

    static string trimString(string);
};


#include "IO.tpp"

#endif // IO_H
