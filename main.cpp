#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>
#include <QGroupBox>
#include <QTimer>

#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QOrbitCameraController>

#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>

#include <QtConcurrent>
#include <QFutureWatcher>

#include "setup.hpp"
#include "update.hpp"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Scene setup
    auto *view       = createView();
    auto *rootEntity = new Qt3DCore::QEntity();
    setupCamera(view, rootEntity);

    // Initial orbital
    QuantumState qs = { 3, 2, 0};
    const float dotSize  = 0.1;
    const int   dotCount = 500;

    float count = dotCount * qs.n * qs.n;
    std::vector<Point> points = computePoints(qs.n, qs.l, qs.m, count);
    Qt3DCore::QEntity *cloudEntity = buildCloudFromPoints(points, dotSize, rootEntity);
    applyZoom(view->camera(), 5.0f * qs.n * qs.n);
    view->setRootEntity(rootEntity);

    // Embed 3D window
    QWidget *container = QWidget::createWindowContainer(view);
    container->setMinimumSize(600, 500);

    // Slider panel
    SliderPanel panel = createSliderPanel(qs.n, qs.l, qs.m);

    // Top-level widget
    QWidget *mainWidget = createMainWidget(container, panel.infoLabel, panel.groupBox);

    // Wire everything together
    connectSliders(panel, qs, cloudEntity, view->camera(), rootEntity, mainWidget, dotSize, dotCount);

    mainWidget->show();
    return app.exec();
}