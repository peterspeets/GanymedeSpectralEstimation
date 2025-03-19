
#include "IO.h"

IO::IO()
{
    //ctor
}


template <typename T>
void IO::saveArrayToFile(const T* array, const int N, const string& filename) {
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
void IO::test(T** array, int a ){

}

//template void IO::test<float>(float**, int);
//template void IO::test<double>(double**, int);



template <typename T>
void IO::save2DArrayToFile(T** array, const int M,const int N, const string& filename,char separator) {
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }
    for(int i = 0; i < M; i++){
        for (int j = 0; j < N; j++) {
            outFile << array[i][j];
            if(j != N-1){
                outFile << separator;
            }
        }
        if(i != M-1){
            outFile  << endl;
        }
    }

    outFile.close();
}


//template void IO::save2DArrayToFile<float>(float** array, const int M, const int N, const std::string& filename, char separator);
//template void IO::save2DArrayToFile<double>(double** array, const int M, const int N, const std::string& filename, char separator);




template <typename T>
void saveArrayToFile(const T* array, const int N, const string& filename) {
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
void save2DArrayToFile(const T* array, const int M,const int N, const string& filename, char separator = ',') {
    ofstream outFile(filename);

    if (!outFile) {
        cout << "Cannot open: " << filename << endl;
        return;
    }
    for(int i = 0; i < M; i++){
        for (int j = 0; j < N; j++) {
            outFile << array[i][j];
            if(j != N-1){
                outFile << separator;
            }
        }
        if(i != M-1){
            outFile  << endl;
        }
    }

    outFile.close();
}

void saveArrayToFile(const kiss_fft_cpx* cpx, const int N, const string& filename) {
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



pair<float*, int> loadArrayFromFile(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cout << "Cannot open: " << filename << endl;
    }

    int N = 1; // number of entries. The array in the file should have at least 1 element.
    string line;
    char separationCharacter = '\n';


    getline(file,line);
    for(int i = 0; i < line.size();i++){
        if(line[i] == ','){
            separationCharacter = ',';
            break;
        }else if(line[i] == '\t'){
            separationCharacter = '\t';
            break;
        // Spaces may be added for readability. This 25 size offset is  a (temporary ?) heuristic to see if no other separation character has been found.
        }else if(line[i] == ' ' && i > 25){
            separationCharacter = ' ';
            break;
        }
    }

    float* loadedArray;
    if(separationCharacter == '\n'){
        vector<float> tempVector;
        tempVector.push_back( stof(line));
        while(getline(file,line)){
            tempVector.push_back( stof(line));
            N++;
        }
        loadedArray = new float[N];
        copy(tempVector.begin(),tempVector.end(),loadedArray);
    }else{
        for(int i = 0; i < line.size();i++){
            if(line[i] == separationCharacter){
                N++;
            }
        }
        loadedArray = new float[N];
        int Nindex = 0;
        int previousIndex = 0;
        for (int i = 0; i < line.size(); i++) {
            if(line[i] == separationCharacter){
                loadedArray[Nindex] = stof(line.substr(previousIndex,i-previousIndex));
                Nindex++;
                previousIndex = i;
            }
        }
    }
    return make_pair(loadedArray, N);
}


tuple<float**, int, int> load2DArrayFromFile(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cout << "Cannot open: " << filename << endl;
    }
    string line;
    int N = 1;
    int M = 1;
    char separationCharacter = ',';

    vector<vector<float>> tempVector;

    getline(file,line);
    for(int i = 0; i < line.size();i++){
        if(line[i] == ','){
            separationCharacter = ',';
            break;
        }else if(line[i] == '\t'){
            separationCharacter = '\t';
            break;
        // Spaces may be added for readability. This 25 size offset is  a (temporary ?) heuristic to see if no other separation character has been found.
        }else if(line[i] == ' ' && i > 30){
            separationCharacter = ' ';
            break;
        }
    }



    vector<float> row;
    int previousIndex = 0;


    for (int i = 0; i < line.size(); i++) {
        if(line[i] == separationCharacter){
            row.push_back(stof(line.substr(previousIndex,i-previousIndex)));
            N++;
            previousIndex = i;
        }
    }

    tempVector.push_back(row);

    while(getline(file,line)){
        row.clear();
        previousIndex = 0;
        for (int i = 0; i < line.size(); i++) {
            if(line[i] == separationCharacter){
                row.push_back(stof(line.substr(previousIndex,i-previousIndex)));
                previousIndex = i;
            }
        }
        tempVector.push_back(row);
        M++;
    }



    float** array = new float*[M];
    for(int i = 0; i<M;i++){
        array[i] = new float[N];
    }

    for(int i = 0; i<M;i++){
        row = tempVector[i];
        copy(row.begin(),row.end(),array[i]);
    }


    return make_tuple(array, M,N);
}



class GanymedeFileLoader {

private:
    mz_zip_archive OCTArchive;

    char* loadHeader() {
        string headerFileName= "header.xml";
        size_t uncomp_size;
        int file_index = mz_zip_reader_locate_file(&OCTArchive, headerFileName.c_str(), NULL, MZ_ZIP_FLAG_IGNORE_PATH);
        void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
        char* file_contents = static_cast<char*>(p);
        return file_contents;
    }

    int getFileIndex( string fileName) {
        int i;
        int file_index = mz_zip_reader_locate_file(&OCTArchive, fileName.c_str(), NULL, MZ_ZIP_FLAG_IGNORE_PATH);

        /*mz_zip_reader_locate_file fails to locate directories, so here a workaround.
        If no file has been found, it looks for .
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

public:

    float** loadSpectrum(int spectrumIndex, float* referenceSpectrum) {
        size_t uncomp_size;
        int i, j, k;
        string spectrumFileName = settings->pathsSpectra[spectrumIndex];
        int file_index = getFileIndex(spectrumFileName);
        void* p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
        unsigned char* file_contents = static_cast<unsigned char*>(p);


        float** spectrum = new float*[settings->sizeXSpectrum];
        for(i = 0; i< settings->sizeXSpectrum; i++) {
            spectrum[i] = new float[settings->sizeZSpectrum];
        }

        for(i = settings->scanStartIndices[spectrumIndex]; i < settings->sizeXSpectrumRaw; i++) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                spectrum[i-settings->scanStartIndices[spectrumIndex]][j] = 0.0;
                for(k = 0; k < settings->bytesPerPixelSpectrum; k++) {
                    spectrum[i-settings->scanStartIndices[spectrumIndex]][j] += settings->electronCountScaling*(file_contents[settings->bytesPerPixelSpectrum*(i*settings->sizeXSpectrum+j)+k]<<(8*k));

                }
            }
        }



        for(i = 0; i < settings->scanStartIndices[spectrumIndex]; i++) {
            for(j = 0; j < settings->sizeZSpectrum; j++) {
                referenceSpectrum[0] = 0.0;
                for(k = 0; k < settings->bytesPerPixelSpectrum; k++) {
                    referenceSpectrum[j] += settings->electronCountScaling*(file_contents[settings->bytesPerPixelSpectrum*(i*settings->sizeXSpectrum+j)+k]<<(8*k));
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


    float** loadSpectrum( int spectrumIndex) {
        size_t uncomp_size;
        int i, j, k;
        string spectrumFileName = settings->pathsSpectra[spectrumIndex];
        int file_index = getFileIndex(spectrumFileName);
        void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
        unsigned  char* file_contents = static_cast<unsigned char*>(p);
        float** spectrum = new float*[settings->sizeXSpectrum];



        for(i = 0; i< settings->sizeXSpectrum; i++) {
            spectrum[i] = new float[settings->sizeZSpectrum];
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



        saveArrayToFile( spectrum[0], settings->sizeZSpectrum, "D:\\data\\rawSpectruma.txt");




        /*
        miniz uses malloc to allocate memory, not new.
        */
        //delete[] file_contents;
        free(file_contents);
        return spectrum;

    }

    float* loadCalibrationSpectrum( string spectrumFileName,int arrayLength = 2048,int bytesPerPixel = 4 ) {
        int i;
        size_t uncomp_size;
        int file_index = getFileIndex(spectrumFileName);
        void * p = mz_zip_reader_extract_to_heap(&OCTArchive, file_index, &uncomp_size, NULL);
        float* spectrum = reinterpret_cast<float*>(p);
        return spectrum;
    }

    Settings* loadSettings() {
        char* file_contents = loadHeader();
        XML_Interpreter interpreter = XML_Interpreter(file_contents);
        settings = new Settings(interpreter.tags);
        free(file_contents);
        return settings;
    }




    GanymedeFileLoader(const string filePath) {
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
        Settings* settings = loadSettings();


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


    ~GanymedeFileLoader() {
        delete settings;
    }

};






IO::~IO()
{
    //dtor
}
