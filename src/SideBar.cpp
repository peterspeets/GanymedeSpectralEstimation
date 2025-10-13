#include <iostream>
#include "SideBar.h"
#include "Window.h"

SideBar::SideBar(Window* parent){
    parentWindow = parent;



    selectColorMapComboBoxLabel = new QLabel("Color map");
    selectColorMapComboBox = new QComboBox();
    selectColorMapComboBox->addItem("Greyscale");
    selectColorMapComboBox->addItem("Viridis");
    selectColorMapComboBox->addItem("Plasma");
    selectColorMapComboBox->addItem("Inferno");
    selectColorMapComboBox->addItem("Magma");
    selectColorMapComboBox->addItem("Cividis");
    selectColorMapComboBox->addItem("Jet");
    selectColorMapComboBox->setCurrentIndex(settings->colorMap);
    selectColorMapComboBox->setMinimumWidth(100);
    selectColorMapComboBox->setMinimumWidth(100);
    connect(selectColorMapComboBox, &QComboBox::currentIndexChanged, this, &selectColorMapComboBoxToggled);



    floorPixelValueSpinBoxLabel = new QLabel();
    floorPixelValueSpinBox = new QDoubleSpinBox();
    floorPixelValueSpinBoxLayout = new QVBoxLayout;
    floorPixelValueSpinBox->setMinimumWidth(100);
    floorPixelValueSpinBox->setMaximumWidth(200);
    floorPixelValueSpinBoxLayout->addWidget(floorPixelValueSpinBoxLabel);
    floorPixelValueSpinBoxLayout->addWidget(floorPixelValueSpinBox);
    floorPixelValueSpinBoxLayout->setSpacing(1);
    floorPixelValueSpinBoxLayout->setContentsMargins(0, 0, 0, 0);


    ceilPixelValueSpinBoxLabel = new QLabel();
    ceilPixelValueSpinBox = new QDoubleSpinBox();
    ceilPixelValueSpinBoxLayout = new QVBoxLayout;
    ceilPixelValueSpinBox ->setMinimumWidth(100);
    ceilPixelValueSpinBox ->setMaximumWidth(200);
    ceilPixelValueSpinBoxLayout->addWidget(ceilPixelValueSpinBoxLabel);
    ceilPixelValueSpinBoxLayout->addWidget(ceilPixelValueSpinBox );
    ceilPixelValueSpinBoxLayout->setSpacing(1);
    ceilPixelValueSpinBoxLayout->setContentsMargins(0, 0, 0, 0);

    connect(floorPixelValueSpinBox, &QDoubleSpinBox::valueChanged, this, &floorPixelValueSpinBoxChanged);
    connect(ceilPixelValueSpinBox, &QDoubleSpinBox::valueChanged, this, &ceilPixelValueSpinBoxChanged);



    setColorScalingSpinboxes();




    startButton = new QPushButton("Process image");

    connect(startButton, &QPushButton::released, this, &startButtonClicked);

    selectFFTRadioButton= new QRadioButton("FFT");
    selectRIAARadioButton= new QRadioButton("RIAA");
    selectFFTRadioButton->setChecked(!settings->useRIAA);
    selectRIAARadioButton->setChecked(settings->useRIAA);

    upscalingFactorComboBoxLayout = new QVBoxLayout;
    upscalingFactorComboBoxLabel = new QLabel("Select upscaling factor:");
    upscalingFactorComboBox = new QComboBox();

    upscalingFactorComboBox->addItem("1");
    upscalingFactorComboBox->addItem("2");
    upscalingFactorComboBox->addItem("4");
    upscalingFactorComboBox->addItem("8");
    upscalingFactorComboBox->addItem("16");
    upscalingFactorComboBox->setMinimumWidth(100);
    upscalingFactorComboBox->setMaximumWidth(200);
    upscalingFactorComboBoxLayout->addWidget(upscalingFactorComboBoxLabel);
    upscalingFactorComboBoxLayout->addWidget(upscalingFactorComboBox);
    upscalingFactorComboBoxLayout->setSpacing(1);

    upscalingFactorComboBoxLayout->setContentsMargins(0, 0, 0, 0);

    bool comboBoxElementNotFound = true;
    for(int i = 0; i < upscalingFactorComboBox->count(); i++ )
    {
        if(settings->upscalingFactor == upscalingFactorComboBox->itemText( i ).toInt()){
            comboBoxElementNotFound  = false;

            upscalingFactorComboBox->setCurrentIndex(i);
            break;
        }
    }
    if(comboBoxElementNotFound){
        cout << "Nonviable upscaling factor. Setting upscaling to 1." << endl;
        settings->upscalingFactor = 1;
        upscalingFactorComboBox->setCurrentIndex(0);

    }
    connect(upscalingFactorComboBox, &QComboBox::currentIndexChanged, this, &upscalingFactorComboBoxChanged);

    QButtonGroup *selectAlgorithmRadioButtonsGroup = new QButtonGroup(this);
    selectAlgorithmRadioButtonsGroup->addButton(selectFFTRadioButton);
    selectAlgorithmRadioButtonsGroup->addButton(selectRIAARadioButton);


    QHBoxLayout *selectAlgorithmRadioButtonsLayout = new QHBoxLayout;
    selectAlgorithmRadioButtonsLayout->addWidget(selectFFTRadioButton);
    selectAlgorithmRadioButtonsLayout->addWidget(selectRIAARadioButton);


    connect(selectAlgorithmRadioButtonsGroup, QOverload<QAbstractButton *>::of(&QButtonGroup::buttonClicked),
        this, selectAlgorithmRadioButtonClicked);

    useDecibelColorScaleCheckBox = new QCheckBox("Use log scale");
    useDecibelColorScaleCheckBox->setChecked(settings->useDecibelColorScale);
    connect(useDecibelColorScaleCheckBox, &QCheckBox::checkStateChanged, this, &useDecibelColorScaleCheckBoxToggled);

    redoPreprocessingCheckBox = new QCheckBox("Force reload\n(turn on to change dispersion correction))");
    redoPreprocessingCheckBox->setChecked(settings->alwaysRedoPreprocessing);
    connect(redoPreprocessingCheckBox, &QCheckBox::checkStateChanged, this, &redoPreprocessingCheckBoxToggled);

    objectiveSelectionComboBoxLayout = new QVBoxLayout;
    objectiveSelectionComboBoxLabel = new QLabel("Select objective lens for dispersion correction");

    objectiveSelectionComboBox = new QComboBox();

    populateObjectiveSelectionComboBox();

    objectiveSelectionComboBox->setMinimumWidth(100);
    objectiveSelectionComboBox->setMaximumWidth(200);

    objectiveSelectionComboBoxLayout->addWidget(objectiveSelectionComboBoxLabel);
    objectiveSelectionComboBoxLayout->addWidget(objectiveSelectionComboBox);
    objectiveSelectionComboBoxLayout->setSpacing(1);
    connect(objectiveSelectionComboBox, &QComboBox::currentIndexChanged, this, &objectiveSelectionComboBoxChanged);

    addWidget(selectColorMapComboBox);
    addWidget(useDecibelColorScaleCheckBox);
    addLayout(ceilPixelValueSpinBoxLayout);
    addLayout(floorPixelValueSpinBoxLayout);

    addStretch();
    addWidget(redoPreprocessingCheckBox);
    addLayout(selectAlgorithmRadioButtonsLayout);

    addLayout(upscalingFactorComboBoxLayout);
    addStretch();

    addLayout(objectiveSelectionComboBoxLayout);

    addStretch();
    addStretch();
    addWidget(startButton);

}

void SideBar::populateObjectiveSelectionComboBox(){
    objectiveSelectionComboBox->clear();
    for (pair<const string, vector<double>> &dispersionDataPair: settings->objectiveDispersionData) {
        objectiveSelectionComboBox->addItem(QString::fromStdString(dispersionDataPair.first));
    }

    int index = objectiveSelectionComboBox->findText("Native");
    if (index != -1) {
        objectiveSelectionComboBox->setCurrentIndex(index);
    }
    return;
}

void SideBar::redoPreprocessingCheckBoxToggled(){

    settings->alwaysRedoPreprocessing = redoPreprocessingCheckBox->isChecked();

    return;
}

void SideBar::selectAlgorithmRadioButtonClicked(QAbstractButton* button){
    settings->useRIAA = selectRIAARadioButton->isChecked();

    if(settings->useRIAA){
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    }else{
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }


    return;
}


void SideBar::objectiveSelectionComboBoxChanged(int index){
    settings->objectiveLabel =  objectiveSelectionComboBox->itemText(index).toStdString();
    //double* dispersionCoefficients = nullptr;
    //map<string, vector<double>> objectiveDispersionData;

    settings->dispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].data();
    settings->numberOfDispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].size();
    return;
}

void SideBar::floorPixelValueSpinBoxChanged(double newValue){

    if(settings->useDecibelColorScale){
        settings->decibelFloor = newValue;

        ceilPixelValueSpinBox->setRange(settings->decibelFloor, 0.0);
    }else{
        settings->percentageFloor = newValue;
        ceilPixelValueSpinBox->setRange(settings->percentageFloor, 100);
    }
    if(scan && scan->imageFFT){
        if(settings->useRIAA){
            parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
        }else{
            parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
        }
    }
    return;
}

void SideBar::ceilPixelValueSpinBoxChanged(double newValue){

    if(settings->useDecibelColorScale){
        settings->decibelCeil = newValue;
        floorPixelValueSpinBox->setRange(-120, settings->decibelCeil);
    }else{
        settings->percentageCeil= newValue;
        floorPixelValueSpinBox->setRange(0, settings->percentageCeil);
    }
    if(scan){
        if(settings->useRIAA && scan->imageRIAA){
            parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
        }else if(scan->imageFFT){
            parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
        }
    }
    return;
}


void SideBar::setColorScalingSpinboxes(){
    /*
    An update of the spinbox due to the change of the range also counts
    as a change in the value. Therefore, before changing to the new range,
    the spinBox range needs to be widened to encompass both the decibel
    and the percentage values;
    */
    floorPixelValueSpinBox->setRange(-9999, 9999);
    ceilPixelValueSpinBox->setRange(-9999, 9999);


    if(settings->useDecibelColorScale){
        /*
        The value needs to be set first to not to conflict with the new range settings.
        Since the new range is set after each new range is set, a dummy new range
        is necessary. Setting the range also triggers the connected slot, and therefore
        also changes the settings-> values. This is why three calls to setValue are needed,
        and a temp newCeil and new Floor.
        */
        double newFloor = settings->decibelFloor;
        double newCeil = settings->decibelCeil;
        floorPixelValueSpinBox->setValue(newFloor);
        ceilPixelValueSpinBox->setValue(newCeil);
        floorPixelValueSpinBox->setValue(newFloor);

        floorPixelValueSpinBoxLabel->setText("Floor (decibel)");
        floorPixelValueSpinBox->setSingleStep(1.0);
        floorPixelValueSpinBox->setRange(-120.0, settings->decibelCeil);
        ceilPixelValueSpinBoxLabel->setText("Ceiling (decibel)");
        ceilPixelValueSpinBox->setSingleStep(1.0);
        ceilPixelValueSpinBox->setRange(settings->decibelFloor, 0);
    }else{
        double newFloor = settings->percentageFloor;
        double newCeil = settings->percentageCeil;
        floorPixelValueSpinBox->setValue(newFloor);
        ceilPixelValueSpinBox->setValue(newCeil);
        floorPixelValueSpinBox->setValue(newFloor);

        floorPixelValueSpinBoxLabel->setText("Floor (percent)");
        floorPixelValueSpinBox->setSingleStep(0.1);
        floorPixelValueSpinBox->setRange(0.0, settings->percentageCeil);
        ceilPixelValueSpinBoxLabel->setText("Ceiling (percent)");
        ceilPixelValueSpinBox->setSingleStep(0.1);



        ceilPixelValueSpinBox->setRange(settings->percentageFloor, 100.0);

    }

    return;
}



void SideBar::upscalingFactorComboBoxChanged(int index){



    settings->upscalingFactor =  upscalingFactorComboBox->itemText(index).toInt();


    return;
}



void SideBar::selectColorMapComboBoxToggled(int index){

    settings->colorMap = index;



    if(settings->useRIAA){
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    }else{
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }
    return;
}


void SideBar::useDecibelColorScaleCheckBoxToggled(){
    settings->useDecibelColorScale = useDecibelColorScaleCheckBox->isChecked();
    setColorScalingSpinboxes();



    if(settings->useRIAA){
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    }else{
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }



    return;
}


void SideBar::startButtonClicked(){
    if(settings->alwaysRedoPreprocessing){
        parentWindow->loadFile(settings->pathToData);
    }
    if(settings->useRIAA){
        scan->processBScan();
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    }else{
        scan->fftBScan();
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }


    return;
}


SideBar::~SideBar()
{
    //dtor
}


