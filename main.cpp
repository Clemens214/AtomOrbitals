#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>
#include <QGroupBox>

#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QOrbitCameraController>

#include <Qt3DCore/QEntity>
#include <Qt3DRender/QCamera>

#include <QtConcurrent>
#include <QFutureWatcher>

#include <iostream>
#include <cctype>
#include <random>
#include <cmath>

#include "slider.hpp"
#include "render.hpp"
#include "sample.hpp"
#include "Atom.hpp"
#include "Coordinates.hpp"
#include "Orbital.hpp"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // --- 3D Window ---
    auto *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(QColor(30, 30, 40));

    auto *rootEntity = new Qt3DCore::QEntity();

    Qt3DRender::QCamera *camera = view->camera();
    camera->lens()->setPerspectiveProjection(45.0f, 16.f/9.f, 0.1f, 1000.f);
    camera->setPosition(QVector3D(0, 0, 20));
    camera->setViewCenter(QVector3D(0, 0, 0));

    auto *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);
    camController->setCamera(camera);
    camController->setLinearSpeed(20.0f);
    camController->setLookSpeed(180.0f);

    // Initial quantum numbers
    int n = 1, l = 0, m = 0;
    const float dotSize  = 0.12f;
    const int   dotCount = 300;

    Qt3DCore::QEntity *cloudEntity = buildCloud(n, l, m, dotSize, dotCount, rootEntity);
    view->setRootEntity(rootEntity);

    // --- Embed 3D window ---
    QWidget *container = QWidget::createWindowContainer(view);
    container->setMinimumSize(600, 500);

    // --- Slider panel ---
    auto *groupBox  = new QGroupBox("Quantum Numbers");
    auto *grid      = new QGridLayout(groupBox);
    grid->setColumnStretch(1, 1);   // let sliders expand

    // n : 1 .. 7  (reasonable display range)
    // l : 0 .. n-1
    // m : -l .. +l  (we map slider range [0, 2l] → value - l)
    auto sN = makeLabelledSlider("n", 1,  7,         n,     grid, 0);
    auto sL = makeLabelledSlider("l", 0,  n - 1,     l,     grid, 1);
    auto sM = makeLabelledSlider("m", -l, l,         m,     grid, 2);

    // Orbital info label
    auto *infoLabel = new QLabel(QString("Orbital: n=%1  l=%2  m=%3").arg(n).arg(l).arg(m));
    infoLabel->setAlignment(Qt::AlignCenter);

    // --- Rebuild helper ---
    auto rebuild = [&]() {
        infoLabel->setText(QString("Computing: n=%1  l=%2  m=%3...").arg(n).arg(l).arg(m));

        // Disable sliders while computing to prevent queued rebuilds
        sN.slider->setEnabled(false);
        sL.slider->setEnabled(false);
        sM.slider->setEnabled(false);

        // Capture by value for thread safety
        int cn = n, cl = l, cm = m;

        auto *watcher = new QFutureWatcher<std::vector<Point>>(mainWidget);
        QFuture<std::vector<Point>> future = QtConcurrent::run([cn, cl, cm]() {
            const float boxHalf = 20.0f;
            Orbital orbit(cn, cl, cm);
            float pMax = computePMax(orbit, boxHalf);
            return rejectionSample(300, orbit, pMax, boxHalf);
        });

    // --- n changed: clamp l, clamp m, update all ranges ---
    QObject::connect(sN.slider, &QSlider::valueChanged, [&](int val) {
        n = val;
        sN.valueLabel->setText(QString::number(n));

        // Block signals to prevent cascade while we fix up l and m
        sL.slider->blockSignals(true);
        sM.slider->blockSignals(true);

        // Clamp l to [0, n-1] and update its range
        sL.slider->setMaximum(n - 1);
        l = std::min(l, n - 1);
        sL.slider->setValue(l);
        sL.valueLabel->setText(QString::number(l));

        // Clamp m to [-l, +l] and update its range
        sM.slider->setMinimum(-l);
        sM.slider->setMaximum( l);
        m = std::clamp(m, -l, l);
        sM.slider->setValue(m);
        sM.valueLabel->setText(QString::number(m));

        sL.slider->blockSignals(false);
        sM.slider->blockSignals(false);

        rebuild();
    });

    // --- l changed: clamp m, update m range ---
    QObject::connect(sL.slider, &QSlider::valueChanged, [&](int val) {
        l = val;
        sL.valueLabel->setText(QString::number(l));

        // Block signals to prevent cascade while we fix up m
        sM.slider->blockSignals(true);

        sM.slider->setMinimum(-l);
        sM.slider->setMaximum( l);
        m = std::clamp(m, -l, l);
        sM.slider->setValue(m);
        sM.valueLabel->setText(QString::number(m));

        sM.slider->blockSignals(false);

        rebuild();
    });

    // --- m changed: straightforward ---
    QObject::connect(sM.slider, &QSlider::valueChanged, [&](int val) {
        m = val;
        sM.valueLabel->setText(QString::number(m));
        rebuild();
    });

    // --- Main layout ---
    QWidget *mainWidget = new QWidget();
    mainWidget->setWindowTitle("Qt3D Orbital Viewer");

    auto *layout = new QVBoxLayout(mainWidget);
    layout->addWidget(container);
    layout->addWidget(infoLabel);
    layout->addWidget(groupBox);

    mainWidget->resize(640, 620);
    mainWidget->show();

    return app.exec();
}