#ifndef WINDOW_H
#define WINDOW_H

#include <iostream>
#include <QApplication>
#include <QtWidgets>
#include <QPushButton>
#include <QFileDialog>
#include <cmath>
#include "Settings.h"
#include "UtilityMathFunctions.h"
#include "IO.h"
#include "BScan.h"
#include "globals.h"
#include "ColorMaps.h"



#include "SideBar.h"

using namespace std;

class ObjectiveSettingsWindow : public QWindow{
    Q_OBJECT
    public:
        ObjectiveSettingsWindow();

};


class Window : public QMainWindow {
    Q_OBJECT
public:
    Window();
    void displayImage(unsigned char*,size_t width,size_t height);
    void setImage(float** floatImage,const size_t floatImagewidth, const size_t floatImageheight_);
    virtual ~Window();

protected:

private slots:
    void confirmExit();
    void loadFileWithDialog();
    void saveFileWithDialog();
    void openchangeObjectiveSettingsWindow();

private:
    QWidget* centralWidget ;
    QVBoxLayout* mainLayout;
    void buildMenuBar();
    void buildSideBar();
    QHBoxLayout* topToolBarLayout;
    QLayout* mainCanvasLayout;
    QHBoxLayout* sideBarAndCanvasLayout;
    QLabel* imageLabel = nullptr;
    QMenu* fileMenu;
    QMenu* settingsMenu;
    SideBar* sideBar;
    unsigned char* bitMapImage = nullptr;
    QWindow changeObjectiveSettingsWindow();
    ObjectiveSettingsWindow* objectiveSettingsWindow;

    tuple<unsigned char,unsigned char,unsigned char> applyColorMap(double);
    size_t imageWidth = 0;
    size_t imageHeight = 0;
    size_t maxImageHeight = 512;
    size_t maxImageWidth = 1024;

};










#endif // WINDOW_H

