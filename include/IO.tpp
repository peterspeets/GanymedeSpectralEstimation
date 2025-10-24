#include "IO.h"
#include "Settings.h"



template <typename floatingPointType>
IO<floatingPointType>::IO() {
    //This class contains static functions
}


template<typename T>
void IO<T>::savePng(const string filename, const int width, const int height, const unsigned char* image) {

    /*
    filename: path of the file to store the png to.
    width: width of the image in pixels
    height: height of the image in pixels
    image: flattened array of the bitmap, row for row.

    This function saves a bytestring that corresponds to a bitmap image into the PNG format.
    It uses the LibPNG library, which requires raw pointers.
    */

    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        perror("Failed to open file for writing");
        exit(EXIT_FAILURE);
    }
    //The nullpointers are what would usually be the error handling. Setting this as a nullpointer means no error handling.
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png) {
        perror("Failed to create PNG write struct");
        exit(EXIT_FAILURE);
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        perror("Failed to create PNG info struct");
        png_destroy_write_struct(&png, nullptr);
        exit(EXIT_FAILURE);
    }

    if (setjmp(png_jmpbuf(png))) {
        perror("Error during PNG creation");
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    png_init_io(png, fp);

    // Write header
    png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    // Write image data
    for (int i = 0; i < height; i++) {
        png_write_row(png, &image[i * width * 3]);
    }

    // End write
    png_write_end(png, nullptr);

    fclose(fp);
    png_destroy_write_struct(&png, &info);
}

template<typename T>
void IO<T>::savePng(const string filename, const int width, const int height, const T* const* image,
                    bool decibelColorScale,double decibelFloor ) {

    /*
    filename: path of the png file
    width: width in pixels
    height: height in pixels
    image: pointer to 2D bitmap
    decibelColorScale: if true, then use logarithmic scale, otherwise linear
    decibelFloor: lowest pixel value in decibel, usually -120. Scale is normalized to maximum pixel value, so this value is negative.

    This is a wrapper function of the IO<T>::savePng(const string filename, const int width, const int height, const unsigned char* image)
    function. It takes a floating point (templated type T) image, performs scaling, creates a bitmap and passes it to the other savePng
    overloaded function.

    */

    T maxValue = image[0][0];
    T minValue = image[0][0];
    int i, j;

    for(i = 0; i < width; i++) {
        for(j = 0; j < height; j++) {
            if(image[i][j] > maxValue) {
                maxValue = image[i][j];
            } else if(image[i][j] < minValue) {
                minValue = image[i][j];
            }
        }
    }

    unsigned char* image_png = new unsigned char[(width * height * 3)];


    double absoluteFloor = pow(10,0.1*decibelFloor);
    unsigned char pixelValue = 0;
    for(i = 0; i < width; i++) {
        for(j = 0; j < height; j++) {
            if(decibelColorScale) {
                if(  ((image[i][j] - minValue)/(maxValue - minValue)) <= absoluteFloor ) {
                    pixelValue = 0.0;

                } else {
                    pixelValue = static_cast<unsigned char>( 255*(10*log10(((image[i][j] - minValue)/(maxValue - minValue))) -decibelFloor )/abs(
                                     decibelFloor) );
                }

            } else {
                pixelValue = static_cast<unsigned char>( 255*(((image[i][j] - minValue)/(maxValue - minValue)) ));
            }
            if(image[i][j] >= maxValue) {
                pixelValue = 255;
            }
            image_png[(j * width + i) * 3] = pixelValue;
            image_png[(j * width + i) * 3+1] = pixelValue;
            image_png[(j * width + i) * 3+2] = pixelValue;

        }

    }

    savePng(filename, width, height, image_png);
    delete[] image_png;
}



template<typename T>
void IO<T>::testa(float a) {
    //remove next commit
    cout << "Test function. " << a << endl;
    return;
}


template<typename T>
void IO<T>::savePng(const string filename, const int width, const int height, const int imageWidth, const int imageHeight,
                    const T* const* image, bool decibelColorScale,double decibelFloor ) {
    /*
    filename: path of file to save to
    width: width in pixels of the original data
    height: height in pixels of the original data
    imageWidth: width in pixels of the image that will be saved.
    imageHeight: height in pixels of the image that will be saved.
    image: 2D array oft the data
    decibelColorScale: if true, then use decibel scale
    decibelFloor: lowest pixel value in decibel, if decibelColorScale is true. Usually this is -120. The decibel is scaled to the maximum pixel value

    This function is a wrapper function for the other overloaded savePng functions. This function interpolates the image
    with a spline interpolation to match the new image size, and passes the result to lower level savePng() overloads.

    */
    int i, j;
    double minValue = image[0][0];
    double maxValue = image[0][0];

    if(decibelColorScale) {
        T** logImage = new T*[width ];
        for(i = 0; i < width; i++) {
            logImage[i] = new T[height];
        }

        for(i = 0; i < width; i++) {
            for(j = 0; j < height; j++) {
                if(image[i][j] > maxValue) {
                    maxValue = image[i][j];
                } else if(image[i][j] < minValue) {
                    minValue = image[i][j];
                }
            }
        }


        double absoluteFloor = pow(10,0.1*decibelFloor);
        T pixelValue = 0;
        for(i = 0; i < width; i++) {
            for(j = 0; j < height; j++) {

                if(  ((image[i][j] - minValue)/(maxValue - minValue)) <= absoluteFloor ) {
                    pixelValue = 0.0;

                } else {
                    pixelValue =  255*(10*log10(((image[i][j] - minValue)/(maxValue - minValue))) -decibelFloor )/abs(decibelFloor) ;
                }


                if(image[i][j] >= maxValue) {
                    pixelValue = 255;
                }
                logImage[i][j] = pixelValue;
            }

        }
        savePng(filename, width, height,imageWidth,imageHeight, logImage,false);

        for(i = 0; i < width; i++) {
            delete[] logImage[i];
        }
        delete[] logImage;

        return;
    }

    //If the desired width and height already as desired, no rescaling is needed.
    if(width == imageWidth && height == imageHeight) {
        savePng(filename, width, height,  image, decibelColorScale,decibelFloor);

    } else if(width == imageWidth && height != imageHeight) {//otherwise, rescale:

        float** resizedImage = new float*[imageWidth];
        float* zLinspace = new float[height];
        for (int i = 0; i < height; i++) {
            zLinspace [i] = 1.0*i/height;
        }

        for (int i = 0; i < imageWidth; i++) {
            resizedImage[i] = new float[imageHeight];
        }
        for(i = 0; i < width ; i++) {

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(zLinspace,image[i],height);
            for(j = 0; j < imageHeight; j++) {
                resizedImage[i][j] = spline->evaluate(  static_cast<float>( j)/imageHeight );
            }
            delete spline;
        }
        savePng(filename, imageWidth, imageHeight,  resizedImage, decibelColorScale,decibelFloor);

        for (int i = 0; i < imageWidth; i++) {
            delete[] resizedImage[i];
        }
        delete[] resizedImage;
        delete[] zLinspace ;
    } else if(width != imageWidth && height == imageHeight) {

        float** resizedImage = new float*[imageWidth];
        float* xLinspace = new float[width];
        for (int i = 0; i < width; i++) {
            xLinspace[i] = 1.0*i/width;
        }

        for (int i = 0; i < imageWidth; i++) {
            resizedImage[i] = new float[imageHeight];
        }
        float* xProfile = new float[height];
        for(i = 0; i < height ; i++) {
            for(j = 0; j < width; j++) {
                xProfile[j] = image[j][i];
            }

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(xLinspace,xProfile,width);
            for(j = 0; j < imageWidth; j++) {
                resizedImage[j][i] = spline->evaluate(  static_cast<float>( j)/imageWidth);
            }
            delete spline;
        }
        savePng(filename, imageWidth, imageHeight,  resizedImage, decibelColorScale,decibelFloor);
        for (int i = 0; i < imageWidth; i++) {
            delete[] resizedImage[i];
        }
        delete[] resizedImage;
        delete[] xLinspace;
        delete[] xProfile;
    } else if(width != imageWidth && height != imageHeight) {
        float** shorterImage = new float*[width];
        float** resizedImage = new float*[imageWidth];
        float* zLinspace = new float[height];
        float* xLinspace = new float[width];
        for (int i = 0; i < width; i++) {
            xLinspace[i] = 1.0*i/width;
        }
        for (int i = 0; i < height; i++) {
            zLinspace [i] = 1.0*i/height;
        }

        for (int i = 0; i < imageWidth; i++) {
            resizedImage[i] = new float[imageHeight];
        }
        for (int i = 0; i < width; i++) {
            shorterImage[i] = new float[imageHeight];
        }


        for(i = 0; i < width ; i++) {

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(zLinspace,image[i],height);
            for(j = 0; j < imageHeight; j++) {
                shorterImage[i][j] = spline->evaluate(  static_cast<float>( j)/imageHeight);
            }
            delete spline;
        }

        float* xProfile = new float[height];
        for(i = 0; i < imageHeight ; i++) {
            for(j = 0; j < width; j++) {
                xProfile[j] = shorterImage[j][i];
            }

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<float>::splineInterpolation(xLinspace,xProfile,width);
            for(j = 0; j < imageWidth; j++) {
                resizedImage[j][i] = spline->evaluate(  static_cast<float>( j)/imageWidth);
            }
            delete spline;
        }



        savePng(filename, imageWidth, imageHeight,  resizedImage, decibelColorScale,decibelFloor);

        for (int i = 0; i < imageWidth; i++) {
            delete[] resizedImage[i];

        }
        for (int i = 0; i < width; i++) {
            delete[] shorterImage[i];
        }

        delete[] shorterImage;
        delete[] resizedImage;
        delete[] xLinspace;
        delete[] zLinspace;
        delete[] xProfile;

    }
}

template <typename T>
map<string, vector<double>> IO<T>::loadObjectiveDispersionData(const string& filename) {
    /*
    filename: path to the dispersion data file.

    This function loads the dispersion data YAML file. This file should contain the polynomial
    coefficients of the dispersion calibration. The x-axis of the polynomial is not the wavenumber,
    but the index, after chirp (the chirp of the spectrometer.) correction. This is objective
    dependent.
    */
    ifstream file(filename);
    if (!file) {
        cout << "Cannot open: " << filename << endl;
    } else {

        cout << "opening " << filename << endl;
    }
    string line;

    string keyword;
    string stringOfDoubles;
    string tempString;
    map<string, vector<double>> objectiveLensSettings;
    int previousIndex = 0;

    while(getline(file,line)) {
        for(int i = 0; i < line.size(); i++) {
            if(line[i] == ':') {
                keyword = line.substr(0,i);
                stringOfDoubles = line.substr(i+1,line.size());
                break;
            }

        }

        previousIndex = 0;
        for(int i = 0; i < stringOfDoubles.size(); i++) {


            if(stringOfDoubles[i] == ',') {
                tempString = stringOfDoubles.substr(previousIndex,i-previousIndex);

                previousIndex = i+1;
                objectiveLensSettings[keyword].push_back(stod(tempString));

            }

        }

        tempString = stringOfDoubles.substr(previousIndex,stringOfDoubles.size()-previousIndex);

        if(tempString.size() > 0 &&  ((tempString[0] >= 48 && tempString[0] < 58) || tempString[0] == ' ' || tempString[0] == '-') ) {
            objectiveLensSettings[keyword].push_back(stod(tempString));
        }




    }

    return objectiveLensSettings;


}


template <typename T>
string IO<T>::trimString(string inputString) {
    /*
    inputString: string to trim

    This utility function removes leading and trailing spaces from a string.
    It does not do so in place, but creates a new string.
    */
    string outputString;
    for(int i = 0; i < inputString.size(); i++) {
        if( (inputString[i] >= 48 && inputString[i] < 58) || inputString[i] == '.'  ) {
            outputString = inputString.substr(i,inputString.size());
            break;
        }
    }
    for(int i = outputString.size()-1; i > 0; i--) {
        if( (outputString[i] >= 48 && outputString[i] < 58) || outputString[i] == '.'  ) {
            outputString = outputString.substr(0,i+1);
            break;
        }
    }

    return outputString;
}

template <typename T>
tuple<T**, int, int> IO<T>::load2DArrayFromFile(const string& filename) {
    /*
    filename: path to the text file that contains a 2D array.

    This function takes a path to a 2D text file and stores it in a 2D array.
    The rows are the first dimension, the columns the second: array[row][column].
    */
    ifstream file(filename);
    if (!file) {
        cout << "Cannot open: " << filename << endl;
    }
    string line;
    int N = 1;
    int M = 1;
    char separationCharacter = ',';

    vector<vector<T>> tempVector;

    getline(file,line);
    for(int i = 0; i < line.size(); i++) {
        if(line[i] == ',') {
            separationCharacter = ',';
            break;
        } else if(line[i] == '\t') {
            separationCharacter = '\t';
            break;
            // Spaces may be added for readability. This 30 size offset is  a (temporary ?) heuristic to see if no other separation character has been found.
        } else if(line[i] == ' ' && i > 30) {
            separationCharacter = ' ';
            break;
        }
    }

    vector<T> row;
    int previousIndex = 0;


    for (int i = 0; i < line.size(); i++) {
        if(line[i] == separationCharacter) {
            row.push_back(stof(line.substr(previousIndex,i-previousIndex)));
            N++;
            previousIndex = i;
        }
    }

    tempVector.push_back(row);

    while(getline(file,line)) {
        row.clear();
        previousIndex = 0;
        for (int i = 0; i < line.size(); i++) {
            if(line[i] == separationCharacter) {
                row.push_back(stof(line.substr(previousIndex,i-previousIndex)));
                previousIndex = i;
            }
        }
        tempVector.push_back(row);
        M++;
    }



    T** array = new T*[M];
    for(int i = 0; i<M; i++) {
        array[i] = new T[N];
    }

    for(int i = 0; i<M; i++) {
        row = tempVector[i];
        copy(row.begin(),row.end(),array[i]);
    }


    return make_tuple(array, M,N);
}

template <typename T>
pair<T*, int> IO<T>::loadArrayFromFile(const string& filename) {
    /*
    filename: path to the text file

    This function loads an array from a text file.
    */
    ifstream file(filename);
    if (!file) {
        cout << "Cannot open: " << filename << endl;
    }

    int N = 1; // number of entries. The array in the file should have at least 1 element.
    string line;
    char separationCharacter = '\n';


    getline(file,line);
    for(int i = 0; i < line.size(); i++) {
        if(line[i] == ',') {
            separationCharacter = ',';
            break;
        } else if(line[i] == '\t') {
            separationCharacter = '\t';
            break;
            // Spaces may be added for readability. This 25 size offset is  a (temporary ?) heuristic to see if no other separation character has been found.
        } else if(line[i] == ' ' && i > 25) {
            separationCharacter = ' ';
            break;
        }
    }

    T* loadedArray;
    if(separationCharacter == '\n') {
        vector<T> tempVector;
        tempVector.push_back( stof(line));
        while(getline(file,line)) {
            tempVector.push_back( stof(line));
            N++;
        }
        loadedArray = new T[N];
        copy(tempVector.begin(),tempVector.end(),loadedArray);
    } else {
        for(int i = 0; i < line.size(); i++) {
            if(line[i] == separationCharacter) {
                N++;
            }
        }
        loadedArray = new T[N];
        int Nindex = 0;
        int previousIndex = 0;
        for (int i = 0; i < line.size(); i++) {
            if(line[i] == separationCharacter) {
                loadedArray[Nindex] = stof(line.substr(previousIndex,i-previousIndex));
                Nindex++;
                previousIndex = i;
            }
        }
    }
    return make_pair(loadedArray, N);
}





template<typename floatingPointType>
template <typename T>
void IO<floatingPointType>::saveArrayToFile(const T* array, const int N, const string& filename) {
    /*
    array: array to save
    N: number of elements
    filename: path to the text file

    This function saves an array to a text file. It is saved with a line break separation (\n).
    */
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << array[i] << endl;
    }
    outFile.close();
}

template <typename T>
void IO<T>::saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename) {
    /*
    array: array to save
    N: number of elements
    filename: path to the text file

    This function saves an array to a text file. It is saved with a line break separation (\n).
    Overloaded function to work with complex numbers. It uses the kiss_fft struct complex numbers.
    */
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].r ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile << endl;
    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].i ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile.close();
}

template <typename floatingPointType>
template <typename T>
void IO<floatingPointType>::saveArrayToFile(const complex<T>* cpx, const int N, const string& filename) {
    /*
    array: array to save
    N: number of elements
    filename: path to the text file

    This function saves an array to a text file. It is saved with a line break separation (\n).
    Overloaded function to work with complex numbers. It uses the standard library complex numbers.
    */
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].real() ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile << endl;
    for (int i = 0; i < N; ++i) {
        outFile << cpx[i].imag() ;
        if(i < N-1) {
            outFile << ",";
        }
    }
    outFile.close();

}


template<typename floatingPointType>
template <typename T>
void IO<floatingPointType>::save2DArrayToFile(T** array, const int M,const int N, const string& filename,char separator) {

    /*
    array: pointer to a 2D array
    M: number of rows
    N: number of columns
    filename: path to the file
    separator: column separator.

    This function saves a 2D array to a text file.
    */

    ofstream outFile(filename);
    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }
    for(int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            outFile << array[i][j];
            if(j != N-1) {
                outFile << separator;
            }
        }
        if(i != M-1) {
            outFile  << endl;
        }
    }

    outFile.close();
}



template<typename floatingPointType>
template <typename T>
void IO<floatingPointType>::save2DArrayToFile(T** array, const int M,const int N, const int newM, const int newN,const string& filename,
        char separator) {
    /*
    array: 2D array
    M: number of rows
    N: number of columns
    newM: the new number of colunms
    newN: the new number of rows
    filename: the path of the file
    separator: the column separation character

    This function saves an interpolated array. The original uninterpolated array is stored with size M,N,
    and a new array shape is interpolated such that it has a new shape (newM,newN). This function is written
    for a quick debug comparison between RIAA optimized images and smaller FFT images.
    */
    T** newArray = new T*[M];
    for(int i = 0; i < newM; i++) {
        newArray[i] = new T[N];
    }

    if(newM == M && newN == N) {
        save2DArrayToFile(array, M, N,  filename, separator);
    } else if(newM == M && newN != N) {
        T** newArray = new T*[M];
        T* zLinspace = new T[N];
        for (int i = 0; i < N; i++) {
            zLinspace[i] = 1.0*i/N;
        }

        for (int i = 0; i < newM; i++) {
            newArray[i] = new T[newN];
        }


        for(int i = 0; i < M ; i++) {

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<T>::splineInterpolation(zLinspace,array[i],N);


            for(int j = 0; j < newN; j++) {
                newArray[i][j] = spline->evaluate(  static_cast<T>( j)/newN );
            }
            delete spline;
        }

        save2DArrayToFile(newArray, newM, newN, filename, separator);

        for (int i = 0; i < newM; i++) {
            delete[] newArray[i];
        }
        delete[] newArray;
        delete[] zLinspace ;
        return;
    } else if(newM != M && newN==N) {

        T** newArray = new T*[newM];
        T* xLinspace = new T[M];
        for (int i = 0; i < M; i++) {
            xLinspace[i] = 1.0*i/M;
        }

        for (int i = 0; i < newM; i++) {
            newArray[i] = new T[newN];
        }
        T* column = new T[N];
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                column[j] = array[j][i];
            }

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<T>::splineInterpolation(xLinspace,column,M);
            for(int j = 0; j < newM; j++) {
                newArray[j][i] = spline->evaluate(  static_cast<T>( j)/newM);
            }
            delete spline;
        }
        save2DArrayToFile(newArray, newM, newN, filename, separator);
        for (int i = 0; i < newM; i++) {
            delete[] newArray[i];
        }
        delete[] newArray;
        delete[] xLinspace;
        delete[] column;
        return;
    } else if(M != newM && N != newN) {
        T** shorterArray = new T*[M];
        T** newArray = new T*[newM];
        T* zLinspace = new T[N];
        T* xLinspace = new T[M];
        for (int i = 0; i < M; i++) {
            xLinspace[i] = 1.0*i/M;
        }
        for (int i = 0; i < N; i++) {
            zLinspace [i] = 1.0*i/N;
        }

        for (int i = 0; i < newM; i++) {
            newArray[i] = new T[newN];
        }
        for (int i = 0; i < M; i++) {
            shorterArray[i] = new T[newN];
        }


        for(int i = 0; i < M ; i++) {

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<T>::splineInterpolation(zLinspace,array[i],N);
            for(int j = 0; j < newN; j++) {
                shorterArray[i][j] = spline->evaluate(  static_cast<T>( j)/newN);
            }
            delete spline;
        }

        T* xProfile = new T[N];
        for(int i = 0; i < newN; i++) {
            for(int j = 0; j < M; j++) {
                xProfile[j] = shorterArray[j][i];
            }

            UtilityMathFunctions<float>::SplineInterpolation* spline = UtilityMathFunctions<T>::splineInterpolation(xLinspace,xProfile,M);
            for(int j = 0; j < newM; j++) {
                newArray[j][i] = spline->evaluate(  static_cast<T>( j)/newM);
            }
            delete spline;
        }



        save2DArrayToFile(newArray,  newM, newM,  filename, separator);

        for (int i = 0; i < newM; i++) {
            delete[] newArray[i];

        }
        for (int i = 0; i < M; i++) {
            delete[] shorterArray[i];
        }

        delete[] shorterArray;
        delete[] newArray;
        delete[] xLinspace;
        delete[] zLinspace;
        delete[] xProfile;
        return;
    }
}





template <typename T>
void save2DArrayToFile(const T* array, const int M,const int N, const string& filename, char separator = ',') {

    /*
    array: 2D array to save
    M: number of rows
    N: number of columns
    filename: path of the file to save.
    separator: column separation character.

    This function saves a 2D array.
    */

    ofstream outFile(filename);
    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }
    for(int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            outFile << array[i][j];
            if(j != N-1) {
                outFile << separator;
            }
        }
        if(i != M-1) {
            outFile  << endl;
        }
    }

    outFile.close();
}









template <typename T>
char* IO<T>::GanymedeFileLoader::loadHeader() {
    /*
    This function loads the header from the header from the headerFileName as stored in the GanymedeFileLoader object.
    It returns the file contents as a byte string.
    */
    string headerFileName= "header.xml";
    size_t uncomp_size;
    int file_index = mz_zip_reader_locate_file(&OCTArchive, headerFileName.c_str(), NULL, MZ_ZIP_FLAG_IGNORE_PATH);
    void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
    char* file_contents = static_cast<char*>(p);
    return file_contents;
}



template <typename T>
int IO<T>::GanymedeFileLoader::getFileIndex( string fileName) {
    /*
    fileName: path to the file to load.

    The minizip library indices every file from a zip file. This function returns that index
    of the file name given. The zip file itself is already loaded in the GanymedeFileLoader object.
    */
    int i;
    int file_index = mz_zip_reader_locate_file(&OCTArchive, fileName.c_str(), NULL, MZ_ZIP_FLAG_IGNORE_PATH);

    /*mz_zip_reader_locate_file fails to locate directories, so here a workaround.
    If no file has been found, it loops over all file indices, and tries to load it.
    If that also fails, it replaces the slashes with backslashes.
    */

    if(file_index < 0) {
        int numFiles = mz_zip_reader_get_num_files(&OCTArchive);
        mz_zip_archive_file_stat fileStat;
        for(i = 0; i < numFiles; i++ ) {
            mz_zip_reader_file_stat(&OCTArchive, i, &fileStat);
            if(fileName ==  fileStat.m_filename) {
                file_index = i;
            }
        }

        string newFileName = "";
        for(i = 0; i < fileName.size(); i++) {
            if(fileName[i] == '/' ) {
                newFileName += '\\';
            } else {
                newFileName += fileName[i];
            }

        }
        for(i = 0; i < numFiles; i++ ) {
            mz_zip_reader_file_stat(&OCTArchive, i, &fileStat);
            if(newFileName ==  fileStat.m_filename) {
                file_index = i;
            }
        }
    }
    return file_index;
}


template <typename T>
T** IO<T>::GanymedeFileLoader::loadSpectrum(int spectrumIndex, T* referenceSpectrum) {
    /*
    spectrumIndex: index of the file that contains the spectrum
    referenceSpectrum: pointer to the reference spectrum

    This function loads the B scan data from the file given by spectrumIndex. Since each
    B scan (unless disabled in the settings file of the Ganymede software) has a reference spectrum,
    this spectrum is loaded and stored into referenceSpectrum.
    */
    size_t uncomp_size;
    int i, j, k;
    string spectrumFileName = settings->pathsSpectra[spectrumIndex];

    int file_index = getFileIndex(spectrumFileName);
    void* p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
    unsigned char* file_contents = static_cast<unsigned char*>(p);


    T** spectrum = new T*[settings->sizeXSpectrum];
    for(i = 0; i< settings->sizeXSpectrum; i++) {
        spectrum[i] = new T[settings->sizeZSpectrum];
    }

    for(i = settings->scanStartIndices[spectrumIndex]; i < settings->sizeXSpectrumRaw; i++) {
        for(j = 0; j < settings->sizeZSpectrum; j++) {
            spectrum[i-settings->scanStartIndices[spectrumIndex]][j] = 0.0;
            for(k = 0; k < settings->bytesPerPixelSpectrum; k++) {
                spectrum[i-settings->scanStartIndices[spectrumIndex]][j] += settings->electronCountScaling*(file_contents[settings->bytesPerPixelSpectrum*
                        (i*settings->sizeXSpectrum+j)+k]<<(8*k));

            }
        }
    }

    for(i = 0; i < settings->scanStartIndices[spectrumIndex]; i++) {
        for(j = 0; j < settings->sizeZSpectrum; j++) {
            referenceSpectrum[0] = 0.0;
            for(k = 0; k < settings->bytesPerPixelSpectrum; k++) {
                referenceSpectrum[j] += settings->electronCountScaling*(file_contents[settings->bytesPerPixelSpectrum*(i*settings->sizeXSpectrum+j)+k]<<
                                        (8*k));
            }
        }
    }

    for(j = 0; j < settings->sizeZSpectrum; j++) {
        referenceSpectrum[j] /= settings->scanStartIndices[spectrumIndex];
    }

    /*
    miniz uses malloc to allocate memory, not new.
    */

    free(file_contents);


    return spectrum;
}

template <typename T>
T** IO<T>::GanymedeFileLoader::loadSpectrum( int spectrumIndex) {
    /*
    spectrumIndex: index of the B-scan file to load.

    This function loads a 2D array of B-scan spectral data.
    */
    size_t uncomp_size;
    int i, j, k;


    string spectrumFileName = settings->pathsSpectra[spectrumIndex];
    int file_index = getFileIndex(spectrumFileName);
    void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
    unsigned  char* file_contents = static_cast<unsigned char*>(p);
    T** spectrum = new T*[settings->sizeXSpectrum];

    for(i = 0; i< settings->sizeXSpectrum; i++) {
        spectrum[i] = new T[settings->sizeZSpectrum];
    }
    for(i = settings->scanStartIndices[spectrumIndex]; i < settings->sizeXSpectrumRaw; i++) {
        for(j = 0; j < settings->sizeZSpectrum; j++) {
            spectrum[i-settings->scanStartIndices[spectrumIndex]][j] = 0.0;
            for(k = settings->bytesPerPixelSpectrum-1; k >= 0; k--) { //most significant bit in back
                spectrum[i-settings->scanStartIndices[spectrumIndex]][j] += (settings->electronCountScaling*
                        (file_contents[settings->bytesPerPixelSpectrum*(i*settings->sizeZSpectrum+j)+k  ]<<(8*k))
                                                                            );
            }
        }

    }



    /*
    miniz uses malloc to allocate memory, not new.
    */
    free(file_contents);
    return spectrum;

}


template <typename T>
T* IO<T>::GanymedeFileLoader::loadCalibrationSpectrum( string spectrumFileName,int arrayLength, int bytesPerPixel  ) {
    /*
    spectrumFileName: file name of the B-scan file to load.
    arrayLength: length of the calibration spectrum.
    bytesPerPixel: the number of bytes of the spectrum. usually 2 for bulk B-scan data, and 4 for
    all other (reference,offset) spectra
    */
    int i;
    size_t uncomp_size;
    int file_index = getFileIndex(spectrumFileName);
    void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
    T* spectrum = reinterpret_cast<T*>(p);
    return spectrum;
}

template <typename T>
shared_ptr<Settings> IO<T>::GanymedeFileLoader::loadSettings() {
    /*
    This function loads the settings from the Thorlabs Ganymede .oct file and creates
    a new settings object.
    */
    char* file_contents = loadHeader();
    XML_Interpreter interpreter = XML_Interpreter(file_contents);
    //settings = new Settings(interpreter.tags);
    shared_ptr<Settings> settings = make_shared<Settings>(interpreter.tags);
    free(file_contents);
    return settings;
}



template <typename T>
IO<T>::GanymedeFileLoader::GanymedeFileLoader(const string filePath) {
    /*
    filePath: path of the .oct file to load.

    This function loads the Thorlabs .oct file. A Thorlabs .oct file is a zip file that contains files with
    the spectra for the B-scan, the offset spectrum, the reference spectrum, the spectrometer nonlinearity and
    the metadata.
    */
    mz_bool status;
    int i;
    memset(&OCTArchive, 0, sizeof(OCTArchive)); //when OCTArchive is not zeroed out, mz_zip_reader_init_mem() fails.
    status = mz_zip_reader_init_file(&OCTArchive, filePath.c_str(), 0);
    if(status == 0) {
        cout << "File not found." << endl;
    }
    int numFiles = mz_zip_reader_get_num_files(&OCTArchive);
    mz_zip_archive_file_stat fileStat;
    for(i = 0; i < numFiles; i++ ) {
        mz_zip_reader_file_stat(&OCTArchive, i, &fileStat);
        //cout << "File " << i << ": " << "\"" << fileStat.m_filename << "\"" << endl;
    }
    settings = loadSettings();


    for (const pair<int, string>& pair : settings->pathsSpectra) {
        //cout << pair.first << " " << pair.second << endl;
    }
    cout << "Loading spectrum: " ;
    cout << settings->pathsSpectra[0] << endl;





    //preprocessSpectrumInPlace(spectra, offset, chirp, referenceSpectrum);





    // Call this here instead of the destructor, since it is no longer needed.
    /*
    miniz uses malloc to allocate memory, not new.
    */
    //delete[] file_contents;

};


template <typename T>
IO<T>::GanymedeFileLoader::~GanymedeFileLoader() {

}


