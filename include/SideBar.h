#ifndef SIDEBAR_H
#define SIDEBAR_H


#include <QVBoxLayout>
#include <QApplication>
#include <QtWidgets>
#include <QPushButton>

#include "Settings.h"
#include "BScan.h"

class Window;


using namespace std;

class Settings;

extern shared_ptr<Settings> settings;



class SideBar : public QVBoxLayout  {
    Q_OBJECT
    public:

        SideBar(Window* parent);


        ~SideBar();
    private slots:
        void startButtonClicked();
        void selectAlgorithmRadioButtonClicked(QAbstractButton*);
        void useDecibelColorScaleCheckBoxToggled();
        void selectColorMapComboBoxToggled(int index);
        void upscalingFactorComboBoxChanged(int index);
        void floorPixelValueSpinBoxChanged(double newValue);
        void ceilPixelValueSpinBoxChanged(double newValue);
    private:
        Window* parentWindow = nullptr;
        void setColorScalingSpinboxes();

        QLabel* selectNumberComboBoxLabel;
        QComboBox* selectNumberComboBox;

        QLabel* selectColorMapComboBoxLabel;
        QComboBox* selectColorMapComboBox;



        QRadioButton* selectFFTRadioButton;
        QRadioButton* selectRIAARadioButton;

        QLabel* upscalingFactorComboBoxLabel;
        QComboBox* upscalingFactorComboBox;
        QVBoxLayout* upscalingFactorComboBoxLayout;

        QLabel* floorPixelValueSpinBoxLabel;
        QDoubleSpinBox* floorPixelValueSpinBox;
        QVBoxLayout* floorPixelValueSpinBoxLayout;
        QLabel* ceilPixelValueSpinBoxLabel;
        QDoubleSpinBox* ceilPixelValueSpinBox;
        QVBoxLayout* ceilPixelValueSpinBoxLayout;


        QPushButton* startButton;
        QDoubleSpinBox* decibelFloorSpinBox;
        QCheckBox* useDecibelColorScaleCheckBox;

};

#endif // SIDEBAR_H
