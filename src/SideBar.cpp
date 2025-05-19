#include "SideBar.h"

SideBar::SideBar()
{
    QLabel* selectNumberComboBoxLabel = new QLabel("Select a number:");
    QComboBox* selectNumberComboBox = new QComboBox();
    selectNumberComboBox->addItem("1");
    selectNumberComboBox->addItem("2");
    selectNumberComboBox->addItem("4");
    selectNumberComboBox->addItem("8");
    selectNumberComboBox->addItem("16");

    selectNumberComboBox->setMinimumWidth(100);
    selectNumberComboBox->setMaximumWidth(200);


    addWidget(selectNumberComboBoxLabel);
    addWidget(selectNumberComboBox);
}

SideBar::~SideBar()
{
    //dtor
}


