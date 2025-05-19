#ifndef SIDEBAR_H
#define SIDEBAR_H


#include <QVBoxLayout>
#include <QApplication>
#include <QtWidgets>
#include <QPushButton>

class SideBar : public QVBoxLayout  {
    Q_OBJECT
    public:

        SideBar();


        ~SideBar();
    private slots:

    private:
        QLabel* selectNumberComboBoxLabel;
        QComboBox* selectNumberComboBox;


};

#endif // SIDEBAR_H
