#include <iostream>
#include "SideBar.h"
#include "Window.h"

SideBar::SideBar(Window* parent) {
    /*
    parent: parent window

    This constructor loads all buttons and boxes in the left sidebar.
    */

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
    for(int i = 0; i < upscalingFactorComboBox->count(); i++ ) {
        if(settings->upscalingFactor == upscalingFactorComboBox->itemText( i ).toInt()) {
            comboBoxElementNotFound  = false;
            upscalingFactorComboBox->setCurrentIndex(i);
            break;
        }
    }
    if(comboBoxElementNotFound) {
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

    redoPreprocessingCheckBox = new QCheckBox("Force reload\n(turn on to change dispersion correction)");
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

void SideBar::populateObjectiveSelectionComboBox() {
    /*
    This function repopulates the objective lens selection box based on the data in the settings. This is needed,
    since the data from the text files have their own dispersion data called "Native", whereas the Thorlabs .oct
    data doe not.
    */
    string objectiveLabel = settings->objectiveLabel;//local copy since the global settings one gets overwritten.
    objectiveSelectionComboBox->clear();

    for (pair<const string, vector<double>> &dispersionDataPair: settings->objectiveDispersionData) {
        objectiveSelectionComboBox->addItem(QString::fromStdString(dispersionDataPair.first));
    }


    if(objectiveLabel == "") {
        objectiveLabel = "Native";
    }
    int index = objectiveSelectionComboBox->findText(QString::fromStdString(objectiveLabel));
    if (index != -1) {
        objectiveSelectionComboBox->setCurrentIndex(index);
    }
    return;
}

void SideBar::redoPreprocessingCheckBoxToggled() {
    /*
    This function sets the settings to the redo processing checkmark state.
    */

    settings->alwaysRedoPreprocessing = redoPreprocessingCheckBox->isChecked();

    return;
}

void SideBar::selectAlgorithmRadioButtonClicked(QAbstractButton* button) {
    /*
    button: FFT and RIAA selection button

    This function checks the state of the algorithm selection button and sets the settings accordingly. This function is triggered by changing the state of the select algorithm radio button.
    */
    settings->useRIAA = selectRIAARadioButton->isChecked();

    if(settings->useRIAA) {
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    } else {
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }


    return;
}


void SideBar::objectiveSelectionComboBoxChanged(int index) {
    /*
    index: new index of the to be selected objective lens

    This function sets the new objective lens based on the state of the combobox. This function is triggered by changing objectiveSelectionComboBox.
    */
    if(index < 0) {
        return;
    }
    settings->objectiveLabel =  objectiveSelectionComboBox->itemText(index).toStdString();


    //double* dispersionCoefficients = nullptr;
    //map<string, vector<double>> objectiveDispersionData;

    settings->dispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].data();
    settings->numberOfDispersionCoefficients = settings->objectiveDispersionData[settings->objectiveLabel].size();
    return;
}

void SideBar::floorPixelValueSpinBoxChanged(double newValue) {
    /*
    newValue: The new value that the user gave to the lower pixel value spinbox.

    This function sets either the decibel floor or the percentage floor based on the state of the useDecibelColorScale setting.
    If the log scale is on, it changes the decibel floor, and if the linear scaling is used, it uses a percentage as a lower pixel cutoff for the color map in the displayed image.
    */

    if(settings->useDecibelColorScale) {
        settings->decibelFloor = newValue;

        ceilPixelValueSpinBox->setRange(settings->decibelFloor, 0.0);
    } else {
        settings->percentageFloor = newValue;
        ceilPixelValueSpinBox->setRange(settings->percentageFloor, 100);
    }
    if(scan && !scan->imageFFT.empty()) {
        if(settings->useRIAA) {
            parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
        } else {
            parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
        }
    }
    return;
}

void SideBar::ceilPixelValueSpinBoxChanged(double newValue) {
    /*
    newValue: new maximum value for the color map to display

    This function sets either the decibel ceiling or the percentage floor based on the state of the useDecibelColorScale setting.
    If the log scale is on, it changes the decibel ceiling, and if the linear scaling is used, it uses a percentage as a upper pixel cutoff for the color map in the displayed image.
    */

    if(settings->useDecibelColorScale) {
        settings->decibelCeil = newValue;
        floorPixelValueSpinBox->setRange(-120, settings->decibelCeil);
    } else {
        settings->percentageCeil= newValue;
        floorPixelValueSpinBox->setRange(0, settings->percentageCeil);
    }
    if(scan) {
        if(settings->useRIAA && !scan->imageRIAA.empty()) {
            parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
        } else if(!scan->imageFFT.empty()) {
            parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
        }
    }
    return;
}


void SideBar::setColorScalingSpinboxes() {
    /*
    This function sets the color spinboxes.

    An update of the spinbox due to the change of the range also counts
    as a change in the value. Therefore, before changing to the new range,
    the spinBox range needs to be widened to encompass both the decibel
    and the percentage values;
    */
    floorPixelValueSpinBox->setRange(-9999, 9999);
    ceilPixelValueSpinBox->setRange(-9999, 9999);


    if(settings->useDecibelColorScale) {
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
    } else {
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



void SideBar::upscalingFactorComboBoxChanged(int index) {
    /*
    index: new index for the upscaling factor

    This function triggers when the upscaling factor combobox is triggered. This is a combobox, such that the user
    cannot set the upscaling to very large values, and that it is always a power of 2. An upscaling of 4 or 8 usually
    is fine.
    */


    settings->upscalingFactor =  upscalingFactorComboBox->itemText(index).toInt();


    return;
}



void SideBar::selectColorMapComboBoxToggled(int index) {
    /*
    index: the index of the new color map

    This function sets the color map based on the combobox input and is triggered by changing selectColorMapComboBox.
    */

    settings->colorMap = index;
    if(settings->useRIAA) {
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    } else {
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }
    return;
}


void SideBar::useDecibelColorScaleCheckBoxToggled() {
    /*
    This function sets useDecibelColorScale based on the useDecibelColorScaleCheckBox state, and sets the image accordingly.
    */
    settings->useDecibelColorScale = useDecibelColorScaleCheckBox->isChecked();
    setColorScalingSpinboxes();



    if(settings->useRIAA) {
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    } else {
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }



    return;
}


void SideBar::startButtonClicked() {
    /*
    After the start button is clicked, either a new FFT and RIAA calculation is made based on a new upscaling factor, or, when
    the settings->alwaysRedoPreprocessing is set to true, it reloads the data and reprocesses the spectra. the latter is needed if a
    new dispersion setting is chosen.
    */
    if(settings->alwaysRedoPreprocessing) {
        parentWindow->loadFile(settings->pathToData);
    }
    if(settings->useRIAA) {
        scan->processBScan();
        parentWindow->setImage(scan->imageRIAA,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum*settings->upscalingFactor);
    } else {
        scan->fftBScan();
        parentWindow->setImage(scan->imageFFT,scan->BScanSettings.sizeXSpectrum, scan->BScanSettings.sizeZSpectrum);
    }


    return;
}


SideBar::~SideBar() {
    //dtor
}


