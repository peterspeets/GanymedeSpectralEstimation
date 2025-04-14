#ifndef IO_H
#define IO_H
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <kiss_fft.h>
#include <miniz.h>
#include "Tag.h"
#include <png.h>

#include "UtilityMathFunctions.h"
#include "XML_Interpreter.h"
#include "Settings.h"



using namespace std;

class Settings;

extern shared_ptr<Settings> settings;




template <typename floatingPointType>
class IO
{
    public:
        IO();
        virtual ~IO();
        template <typename T>
        static void saveArrayToFile(const T* array, const int N, const string& filename);

        template <typename T>
        static void saveArrayToFile(const complex<T>* cpx, const int N, const string& filename);

        template <typename T>
        static void test(T);

        static void saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename);


        template <typename T>
        static void save2DArrayToFile(T** array, const int M, const int N, const std::string& filename, char separator = ',');
        template <typename T>
        static void save2DArrayToFile(T** array, const int M,const int N, const int newM, const int newN ,const string& filename,char separator);

        template <typename T>
        static void test(T** array, int);

        static void savePng(const string filename, const int width, const int height, const unsigned char* image);
        static void savePng(const string filename, const int width, const int height, const floatingPointType* const* image,
                            bool decibelColorScale = false,double decibelFloor = -120 );
        static void savePng(const string filename, const int width, const int height, const int imageWidth,
                            const int imageHeight, const floatingPointType* const* image, bool decibelColorScale = false,double decibelFloor = -120 );

        static pair<floatingPointType*, int> loadArrayFromFile(const string& filename);
        static tuple<floatingPointType**, int, int> load2DArrayFromFile(const string& filename);

        class GanymedeFileLoader{
        public:

        shared_ptr<Settings> loadSettings() ;
        floatingPointType** loadSpectrum(int spectrumIndex, floatingPointType* referenceSpectrum) ;
        floatingPointType** loadSpectrum(int spectrumIndex) ;
        floatingPointType* loadCalibrationSpectrum( string spectrumFileName,int arrayLength = 2048,int bytesPerPixel = 4 );
        GanymedeFileLoader(const string filePath);
        ~GanymedeFileLoader() ;

        private:
        mz_zip_archive OCTArchive;
        char* loadHeader();
        int getFileIndex( string fileName);

        };


    protected:

    private:
};


#include "IO.tpp"

#endif // IO_H
