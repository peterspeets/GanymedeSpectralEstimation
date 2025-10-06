#include "window.h"




ObjectiveSettingsWindow::ObjectiveSettingsWindow(){
    this->setVisible(true);



}


Window::Window() {
    resize(960, 500);
    setWindowTitle("Ganymede spectral estimation");



    QIcon icon("../../icons/icontrans.png");

    setWindowIcon(icon);


    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    mainLayout = new QVBoxLayout(centralWidget);
    centralWidget->setLayout(mainLayout);




    buildMenuBar();
    buildSideBar();

    sideBar = new SideBar(this);






    topToolBarLayout = new QHBoxLayout();
    sideBarAndCanvasLayout = new QHBoxLayout();
    mainCanvasLayout = new QGridLayout();

    mainLayout->addLayout(topToolBarLayout);
    mainLayout->addLayout(sideBarAndCanvasLayout);

    sideBarAndCanvasLayout->addLayout(sideBar);
    sideBarAndCanvasLayout->addLayout(mainCanvasLayout);

    if(imageLabel ){
        delete imageLabel ;
    }
    imageLabel = new QLabel();
    imageLabel->setMinimumSize(600, 300);
    imageLabel->setAlignment(Qt::AlignCenter);
    mainCanvasLayout->addWidget(imageLabel);

    /*
    Build settings window:
    */

    objectiveSettingsWindow = new ObjectiveSettingsWindow();
    objectiveSettingsWindow->setTitle("Objective lens settings");


    show();
}

void Window::displayImage(unsigned char* bitMapImage,size_t width, size_t height){

    QImage image(bitMapImage,  width ,  height ,3*width, QImage::Format_RGB888);


    QPixmap pixmap = QPixmap::fromImage(image);

    imageLabel->setPixmap(pixmap);



    return;

}

void Window::buildSideBar() {

}

void Window::buildMenuBar() {
    fileMenu = menuBar()->addMenu(tr("&File"));
    QAction* loadFileAction = fileMenu->addAction(tr("&Load File"));
    connect(loadFileAction, &QAction::triggered, this, &Window::loadFileWithDialog);
    QAction* saveFileAction = fileMenu->addAction(tr("&Save File"));
    connect(saveFileAction, &QAction::triggered, this, &Window::saveFileWithDialog);
    QAction* exitAction = fileMenu->addAction(tr("&Quit"));
    connect(exitAction, &QAction::triggered, this, &Window::confirmExit);

    settingsMenu = menuBar()->addMenu(tr("&Settings"));
    QAction* objectiveLensSettingsAction = settingsMenu->addAction(tr("Objective lens settings"));
    connect(objectiveLensSettingsAction , &QAction::triggered, this, &Window::openchangeObjectiveSettingsWindow);
}

void Window::openchangeObjectiveSettingsWindow(){
    //changeObjectiveSettingsWindow.visible = true;
    return;
}




tuple<unsigned char,unsigned char,unsigned char> Window::applyColorMap(double pixelValue){
    switch (settings->colorMap)
    {
        case 0: return ColorMaps::greyScale(pixelValue);
        case 1: return ColorMaps::viridis(pixelValue);
        case 2: return ColorMaps::plasma(pixelValue);
        case 3: return ColorMaps::inferno(pixelValue);
        case 4: return ColorMaps::magma(pixelValue);
        case 5: return ColorMaps::cividis(pixelValue);
        case 6: return ColorMaps::jet(pixelValue);
    }
}





void Window::confirmExit() {
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "Exit Confirmation", "Are you sure you want to exit?",
                                  QMessageBox::Yes | QMessageBox::No);
    if (reply == QMessageBox::Yes) {
        QApplication::quit();
    }
}

void Window::setImage(float** floatImage,const size_t floatImageWidth, const size_t floatImageHeight_){

    if(floatImageWidth == 0 || floatImageHeight_ == 0){
        return;
    }

    size_t floatImageHeight = floatImageHeight_/2;

    cout << imageWidth << "x" << imageHeight << endl;
    double absoluteFloor = pow(10,0.1*settings->decibelFloor);
    double aspectRatio = settings->x_mm/settings->z_mm;
    if(aspectRatio >= 1.0*maxImageWidth/maxImageHeight ){
        //Requested width is too high: using full width

        imageWidth = maxImageWidth;
        imageHeight= static_cast<int>(1.0*imageWidth/aspectRatio);
    }else{
        //requested width is too high. use full height
        imageHeight= maxImageHeight;
        imageWidth = static_cast<int>(aspectRatio*imageHeight);
    }

    //float floatResizedImage[imageWidth][imageHeight];
    float** floatResizedImage;
    floatResizedImage = new float*[imageWidth];
    for(int i = 0; i < imageWidth; i++){
        floatResizedImage[i] = new float[imageHeight];
    }

    int x1, x2, y1,y2;
    float x, y;
    float temp1,temp2, newValue;
    tuple<unsigned char, unsigned char, unsigned char> RGBValue;



    float minValue = floatImage[0][0];
    float maxValue = floatImage[0][0];
    for(int i = 0; i < floatImageWidth; i++){
        for(int j = 0; j < floatImageHeight; j++){
            if(floatImage[i][j] < minValue ){
                minValue = floatImage[i][j];
            }
            if(floatImage[i][j] > maxValue ){
                maxValue = floatImage[i][j];
            }

        }

    }


    if(bitMapImage){
        delete[] bitMapImage;
    }
    bitMapImage = new unsigned char[3*imageWidth*imageHeight];



    floatResizedImage[0][0] = floatImage[0][0];
    floatResizedImage[0][imageHeight-1] = floatImage[0][floatImageHeight-1];
    floatResizedImage[imageWidth-1][0] = floatImage[floatImageWidth-1][0];
    floatResizedImage[imageWidth-1][imageHeight-1] = floatImage[floatImageWidth-1][floatImageHeight-1];

    for(int i = 0; i < 3; i++){


        if(settings->useDecibelColorScale){
                bitMapImage[i] =  static_cast<unsigned char>(255*(10*log10(floatImage[0][0]) -settings->decibelFloor  )/abs(settings->decibelFloor));
                bitMapImage[imageHeight-1 + i] =  static_cast<unsigned char>(255*(10*log10(floatImage[0][floatImageHeight-1] ) -settings->decibelFloor  )/abs(settings->decibelFloor));
                bitMapImage[ imageHeight*imageWidth - imageHeight + i] =  static_cast<unsigned char>(255*(10*log10(floatImage[floatImageWidth-1][0]) -settings->decibelFloor  )/abs(settings->decibelFloor));
                bitMapImage[ imageHeight*imageWidth - 1 + i] =  static_cast<unsigned char>(255*(10*log10(floatImage[floatImageWidth-1][floatImageHeight-1]) -settings->decibelFloor  )/abs(settings->decibelFloor));


        }else{
            bitMapImage[i] =  static_cast<unsigned char>(255*(floatImage[0][0] - minValue)/(maxValue - minValue));
            bitMapImage[imageHeight-1 + i] =  static_cast<unsigned char>(255*(floatImage[0][floatImageHeight-1] - minValue)/(maxValue - minValue));
            bitMapImage[ imageHeight*imageWidth - imageHeight + i] =  static_cast<unsigned char>(255*(floatImage[floatImageWidth-1][0] - minValue)/(maxValue - minValue));
            bitMapImage[ imageHeight*imageWidth - 1 + i] =  static_cast<unsigned char>(255*(floatImage[floatImageWidth-1][floatImageHeight-1] - minValue)/(maxValue - minValue));

        }
    }



    for(int i =0; i < imageWidth; i++){

        for(int j = 0; j < imageHeight; j++){
            if( (i == 0 && j == 0) || (i == imageWidth-1 && j ==0)
               || (i == 0 && j == imageHeight-1) || (i == imageWidth-1 && j == imageHeight-1) ){
                continue;
            }

            x = i * floatImageWidth/imageWidth;
            y = j * floatImageHeight/imageHeight;
            x1 = floor(x);
            y1 = floor(y);
            x2 =  x1 + 1;
            y2 =  y1 + 1;

            if(x2 >= imageWidth){
                temp1 =  (x2-x)*(floatImage[x1][y1]*(y2-y) + floatImage[x1][y2]*(y-y1) );
                temp2 =  (x-x1)*(floatImage[x1][y1]*(y2-y) + floatImage[x1][y2]*(y-y1) );
                newValue = (((temp1 + temp2)/((x2-x1)*(y2-y1))) - minValue)/(maxValue-minValue) ;
                if(settings->useDecibelColorScale){
                    if(  newValue <= absoluteFloor ){
                        newValue = 0.0;
                    }else{
                        newValue =  (10*log10(newValue) -settings->decibelFloor )/abs(settings->decibelFloor) ;
                    }
                }

            }else if(y2 >= imageHeight){
                temp1 =  (x2-x)*(floatImage[x1][y1]*(y2-y) + floatImage[x1][y1]*(y-y1) );
                temp2 =  (x-x1)*(floatImage[x2][y1]*(y2-y) + floatImage[x2][y1]*(y-y1) );
                newValue = (((temp1 + temp2)/((x2-x1)*(y2-y1))) - minValue)/(maxValue-minValue) ;

                if(settings->useDecibelColorScale){
                    if(  newValue <= absoluteFloor ){
                        newValue = 0.0;
                    }else{
                        newValue =  (10*log10(newValue) -settings->decibelFloor  )/abs(settings->decibelFloor) ;
                    }
                }

            }else{
                temp1 =  (x2-x)*(floatImage[x1][y1]*(y2-y) + floatImage[x1][y2]*(y-y1) );
                temp2 =  (x-x1)*(floatImage[x2][y1]*(y2-y) + floatImage[x2][y2]*(y-y1) );
                newValue = (((temp1 + temp2)/((x2-x1)*(y2-y1))) - minValue)/(maxValue-minValue) ;


                if(settings->useDecibelColorScale){
                    if(  newValue <= absoluteFloor ){
                        newValue = 0.0;
                    }else{
                        newValue =  (10*log10(newValue) -settings->decibelFloor  )/abs(settings->decibelFloor) ;
                    }
                }

            }

            //cout << newValue << endl;

            floatResizedImage[i][j] = newValue;
            //newValue = 1.0*j/imageHeight;
            RGBValue = applyColorMap(newValue);

            //bitMapImage[ 3*(j*imageWidth + i )  ] = static_cast<unsigned char>( 255*newValue );
            //bitMapImage[ 3*(j*imageWidth + i ) +1 ] = static_cast<unsigned char>( 255*newValue );
            //bitMapImage[ 3*(j*imageWidth + i )  +2] = static_cast<unsigned char>( 255*newValue );


            bitMapImage[ 3*(j*imageWidth + i )  ] = get<0>(RGBValue);
            bitMapImage[ 3*(j*imageWidth + i ) +1 ] = get<1>(RGBValue);
            bitMapImage[ 3*(j*imageWidth + i )  +2] = get<2>(RGBValue);


            if(isnan(floatResizedImage[i][j])){
                cout << i  <<"," << j << ": " << x << " " << y << " " << x1 << "  " << y1 << " " << temp1 << " " <<temp1 << endl;

            }

        }
    }


    displayImage(bitMapImage,imageWidth, imageHeight);
}





void Window::loadFileWithDialog() {
    QString filePath = QFileDialog::getOpenFileName(nullptr,"Open File","C:\\","OCT Files (*.oct);;Text Files (*.txt)");

    if (filePath.isEmpty()) {
        return;
    }else {
        cout << "Loaded " << filePath.toStdString() << endl;

        scan =  std::make_shared<BScan>(filePath.toStdString());

        cout << scan->imageFFT << endl;
        cout << scan->imageFFT[0][0] << endl;

        cout << scan->BScanSettings.sizeXSpectrum << "x" << scan->BScanSettings.sizeZSpectrum << endl << endl;

        setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
        cout << "Loaded scan" << endl;
    }

    return;

}



void Window::saveFileWithDialog() {
    QString filePath = QFileDialog::getOpenFileName(nullptr,"Save PNG","C:\\",";;PNG (*.png)");
    int width, height;
    float** image;
    width = settings->sizeXSpectrum;
    if(settings->useRIAA){
        width = settings->sizeZSpectrum * settings->upscalingFactor;
        image = scan->imageRIAA;
    }else{
        width = settings->sizeZSpectrum ;
        image = scan->imageFFT;
    }

    if (filePath.isEmpty()) {
        return;
    }else {
        cout << "Saving " << filePath.toStdString() << endl;
        IO<float>::savePng(filePath.toStdString(), width ,height, width,height,image, settings->useDecibelColorScale,settings->decibelFloor);


    }

    return;

}



Window::~Window() {
    //dtor
}


