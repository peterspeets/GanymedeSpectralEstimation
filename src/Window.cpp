#include "window.h"

Window::Window() {
    resize(960, 500);
    setWindowTitle("Window title");

    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    mainLayout = new QVBoxLayout(centralWidget);
    centralWidget->setLayout(mainLayout);




    buildMenuBar();
    buildSideBar();

    sideBar = new SideBar();






    topToolBarLayout = new QHBoxLayout();
    sideBarAndCanvasLayout = new QHBoxLayout();
    mainCanvasLayout = new QGridLayout();

    mainLayout->addLayout(topToolBarLayout);
    mainLayout->addLayout(sideBarAndCanvasLayout);

    sideBarAndCanvasLayout->addLayout(sideBar);
    sideBarAndCanvasLayout->addLayout(mainCanvasLayout);




    show();
}

void Window::displayImage(unsigned char* bitMapImage,size_t width, size_t height){

    QImage image(bitMapImage,  width ,  height , QImage::Format_RGB888);


    QLabel* imageLabel = new QLabel();
    imageLabel->setAlignment(Qt::AlignCenter);
    QPixmap pixmap = QPixmap::fromImage(image);
    imageLabel->setPixmap(pixmap);
    mainCanvasLayout->addWidget(imageLabel);


}

void Window::buildSideBar() {

}

void Window::buildMenuBar() {
    fileMenu = menuBar()->addMenu(tr("&File"));
    QAction* loadFileAction = fileMenu->addAction(tr("&Load File"));
    connect(loadFileAction, &QAction::triggered, this, &Window::loadFileWithDialog);
    QAction* exitAction = fileMenu->addAction(tr("&Quit"));
    connect(exitAction, &QAction::triggered, this, &Window::confirmExit);


}



void Window::confirmExit() {
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "Exit Confirmation", "Are you sure you want to exit?",
                                  QMessageBox::Yes | QMessageBox::No);
    if (reply == QMessageBox::Yes) {
        QApplication::quit();
    }
}

void Window::loadFileWithDialog() {
    QString filePath = QFileDialog::getOpenFileName(nullptr,"Open File","C:\\","OCT Files (*.oct);;Text Files (*.txt)");

    if (filePath.isEmpty()) {
        return;
    }else {
        cout << "Loaded " << filePath.toStdString() << endl;
        //scan =
        BScan::BScan(filePath.toStdString());
    }

    return;

}

Window::~Window() {
    //dtor
}


