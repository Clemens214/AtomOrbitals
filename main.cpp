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

#include <iostream>
#include <algorithm>

#include "setup.hpp"
#include "slider.hpp"
#include "render.hpp"
#include "sample.hpp"
#include "Atom.hpp"
#include "Coordinates.hpp"
#include "Orbital.hpp"

// -----------------------------------------------------------------------
// Quantum state — shared across all lambdas in main()
// -----------------------------------------------------------------------
struct QuantumState {
    int n = 1;
    int l = 0;
    int m = 0;
};

// -----------------------------------------------------------------------
// Wire up the rebuild logic and slider signal/slot connections.
//
// 'cloudEntity' is passed by reference so the lambda can swap it when
// a new cloud is ready. 'mainWidget' is used as the watcher's parent
// so it is cleaned up with the window.
// -----------------------------------------------------------------------
void connectSliders(SliderPanel       &panel,
                    QuantumState      &qs,
                    Qt3DCore::QEntity *&cloudEntity,
                    Qt3DCore::QEntity  *rootEntity,
                    QWidget            *mainWidget,
                    float               dotSize,
                    int                 dotCount)
{
    // --- Rebuild: runs sampling on a worker thread, builds entities on
    //             the main thread when done ---
    auto rebuild = [&panel, &qs, &cloudEntity, rootEntity,
                    mainWidget, dotSize, dotCount]()
    {
        panel.infoLabel->setText(
            QString("Computing: n=%1  l=%2  m=%3...")
                .arg(qs.n).arg(qs.l).arg(qs.m));

        panel.sN.slider->setEnabled(false);
        panel.sL.slider->setEnabled(false);
        panel.sM.slider->setEnabled(false);

        // Capture quantum numbers by value for thread safety
        const int cn = qs.n, cl = qs.l, cm = qs.m;

        QFuture<std::vector<Point>> future = QtConcurrent::run(
            [cn, cl, cm, dotCount]() {
                const float boxHalf = 5.0f * cn * cn;
                Orbital orbit(cn, cl, cm);
                float pMax = computePMax(orbit, boxHalf);
                return rejectionSample(dotCount, orbit, pMax, boxHalf);
            });

        auto *watcher = new QFutureWatcher<std::vector<Point>>(mainWidget);
        watcher->setFuture(future);

        // Back on the main thread once the worker finishes
        QObject::connect(watcher, &QFutureWatcher<std::vector<Point>>::finished,
                         [&cloudEntity, rootEntity, &panel, &qs, watcher, dotSize, cn, cl, cm]()
                         {
                             std::vector<Point> points = watcher->result();
                             watcher->deleteLater();

                             delete cloudEntity;
                             cloudEntity = buildCloudFromPoints(points, dotSize, rootEntity);

                             panel.infoLabel->setText(
                                 QString("Orbital: n=%1  l=%2  m=%3").arg(cn).arg(cl).arg(cm));

                             panel.sN.slider->setEnabled(true);
                             panel.sL.slider->setEnabled(true);
                             panel.sM.slider->setEnabled(true);
                         });
    };

    // --- Debounce timer: only rebuilds 300 ms after the last slider move ---
    auto *debounceTimer = new QTimer(mainWidget);
    debounceTimer->setSingleShot(true);
    debounceTimer->setInterval(300);
    QObject::connect(debounceTimer, &QTimer::timeout, rebuild);

    // --- n changed: clamp l and m, update all ranges ---
    QObject::connect(panel.sN.slider, &QSlider::valueChanged, [&, debounceTimer](int val) {
        qs.n = val;
        panel.sN.valueLabel->setText(QString::number(qs.n));

        panel.sL.slider->blockSignals(true);
        panel.sM.slider->blockSignals(true);

        panel.sL.slider->setMaximum(qs.n - 1);
        qs.l = std::min(qs.l, qs.n - 1);
        panel.sL.slider->setValue(qs.l);
        panel.sL.valueLabel->setText(QString::number(qs.l));

        panel.sM.slider->setMinimum(-qs.l);
        panel.sM.slider->setMaximum( qs.l);
        qs.m = std::clamp(qs.m, -qs.l, qs.l);
        panel.sM.slider->setValue(qs.m);
        panel.sM.valueLabel->setText(QString::number(qs.m));

        panel.sL.slider->blockSignals(false);
        panel.sM.slider->blockSignals(false);

        debounceTimer->start();
    });

    // --- l changed: clamp m, update m range ---
    QObject::connect(panel.sL.slider, &QSlider::valueChanged, [&, debounceTimer](int val) {
        qs.l = val;
        panel.sL.valueLabel->setText(QString::number(qs.l));

        panel.sM.slider->blockSignals(true);

        panel.sM.slider->setMinimum(-qs.l);
        panel.sM.slider->setMaximum( qs.l);
        qs.m = std::clamp(qs.m, -qs.l, qs.l);
        panel.sM.slider->setValue(qs.m);
        panel.sM.valueLabel->setText(QString::number(qs.m));

        panel.sM.slider->blockSignals(false);

        debounceTimer->start();
    });

    // --- m changed: no downstream effects ---
    QObject::connect(panel.sM.slider, &QSlider::valueChanged, [&, debounceTimer](int val) {
        qs.m = val;
        panel.sM.valueLabel->setText(QString::number(qs.m));
        debounceTimer->start();
    });
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Scene setup
    auto *view       = createView();
    auto *rootEntity = new Qt3DCore::QEntity();
    setupCamera(view, rootEntity);

    // Initial orbital
    QuantumState qs;                    // n=1, l=0, m=0
    const float dotSize  = 0.12f;
    const int   dotCount = 300;

    Qt3DCore::QEntity *cloudEntity =
        buildCloud(qs.n, qs.l, qs.m, dotSize, dotCount, rootEntity);
    view->setRootEntity(rootEntity);

    // Embed 3D window
    QWidget *container = QWidget::createWindowContainer(view);
    container->setMinimumSize(600, 500);

    // Slider panel
    SliderPanel panel = createSliderPanel(qs.n, qs.l, qs.m);

    // Top-level widget
    QWidget *mainWidget =
        createMainWidget(container, panel.infoLabel, panel.groupBox);

    // Wire everything together
    connectSliders(panel, qs, cloudEntity, rootEntity,
                   mainWidget, dotSize, dotCount);

    mainWidget->show();
    return app.exec();
}