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
#include "BScan.h"
#include "Window.h"
#include "globals.h"

using namespace std;

shared_ptr<Settings> settings = nullptr;
shared_ptr<BScan> scan = nullptr;




int main(int argc, char *argv[]) {

    QApplication app(argc, argv);
    //TODO: header is utf-8 in the Ganymede software, here ASCI. If somebody puts an emoticon into the filename, the code might fail.
    settings = make_shared<Settings>();
    //The path is required for the import of the objective settings file.
    settings->pathToExecutable = QCoreApplication::applicationDirPath().toStdString();

    map<string, vector<double>> mapFromData;
    mapFromData = IO<double>::loadObjectiveDispersionData(
                      settings->pathToExecutable + "\\..\\..\\settings\\objectives.yaml");
    settings->objectiveDispersionData.insert(mapFromData.begin(), mapFromData.end());

    Window window;
    window.show();
    //the scan object contains all data and some non static processing functions.
    scan = make_shared<BScan>();

    //set the image. When not debugging, settings->sizeXSpectrum and settings->sizeZSpectrum will be zero.
    window.setImage(scan->fftBScan(),settings->sizeXSpectrum,settings->sizeZSpectrum);
    cout << "end" << endl;

    //return 0;
    return app.exec();
}


