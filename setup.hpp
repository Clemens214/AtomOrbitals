#ifndef SETUP_H
#define SETUP_H

#include <cctype>
#include <cmath>

#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QOrbitCameraController>

#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>

#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

// Create and return the Qt3D window
Qt3DExtras::Qt3DWindow *createView()
{
    auto *view = new Qt3DExtras::Qt3DWindow();
    // set the background color
    view->defaultFrameGraph()->setClearColor(QColor(245, 245, 250));
    return view;
}

// Set up the camera
void setupCamera(Qt3DExtras::Qt3DWindow *view, Qt3DCore::QEntity *rootEntity)
{
    Qt3DRender::QCamera *camera = view->camera();
    camera->lens()->setPerspectiveProjection(45.f, 16.f/9.f, 0.1f, 1000.f);
    camera->setPosition(QVector3D(0, 20, 0));
    camera->setViewCenter(QVector3D(0, 0, 0));
    camera->setUpVector(QVector3D(0, 0, 1));

    auto *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);
    camController->setCamera(camera);
    camController->setLinearSpeed(0.f);
    camController->setLookSpeed(180.f);
    camController->setZoomInLimit(2.f);
}

void applyZoom(Qt3DRender::QCamera *camera, float distance)
{
    QVector3D dir = (camera->position() - camera->viewCenter()).normalized();
    camera->setPosition(camera->viewCenter() + dir * distance);
}

// Create the main widget
QWidget *createMainWidget(QWidget *container, QLabel *infoLabel,
                           QGroupBox *groupBox)
{
    auto *mainWidget = new QWidget();
    mainWidget->setWindowTitle("Qt3D Orbital Viewer");

    auto *layout = new QVBoxLayout(mainWidget);
    layout->addWidget(container);
    layout->addWidget(infoLabel);
    layout->addWidget(groupBox);

    mainWidget->resize(640, 620);
    return mainWidget;
}

#endif // SETUP_H