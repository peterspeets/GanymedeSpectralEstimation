#ifndef WINDOW_H
#define WINDOW_H

#include <iostream>
#include <QApplication>
#include <QtWidgets>
#include <QPushButton>
#include <QFileDialog>
#include "Settings.h"
#include "UtilityMathFunctions.h"
#include "IO.h"
#include "BScan.h"
#include "globals.h"


#include "SideBar.h"

using namespace std;

class Window : public QMainWindow {
    Q_OBJECT
public:
    Window();
    void displayImage(unsigned char*,size_t width,size_t height);
    virtual ~Window();

protected:

private slots:
    void confirmExit();
    void loadFileWithDialog();

private:
    QWidget* centralWidget ;
    QVBoxLayout* mainLayout;
    void buildMenuBar();
    void buildSideBar();
    QHBoxLayout* topToolBarLayout;
    QLayout* mainCanvasLayout;
    QHBoxLayout* sideBarAndCanvasLayout;

    QMenu* fileMenu;
    SideBar* sideBar;

};






#endif // WINDOW_H

