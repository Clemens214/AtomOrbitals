#ifndef UPDATE_H
#define UPDATE_H

#include <QWidget>
#include <Qt3DCore/QEntity>
#include <QSlider>
#include <QLabel>

#include <cctype>
#include <cmath>

#include <QtConcurrent>
#include <QFutureWatcher>

#include "slider.hpp"
#include "render.hpp"

// -----------------------------------------------------------------------
// Wire up the inputs from the slots
// Create the rebuild logic to update the point cloud
// -----------------------------------------------------------------------
// 'cloudEntity' is passed by reference so the lambda can swap it when a new cloud is ready. 
// 'mainWidget' is used as the watcher's parent so it is cleaned up with the window.
// -----------------------------------------------------------------------
void connectSliders(SliderPanel       &panel,
                    QuantumState      &qs,
                    Qt3DCore::QEntity *&cloudEntity,
                    Qt3DCore::QEntity  *rootEntity,
                    QWidget            *mainWidget,
                    float               dotSize,
                    int                 dotCount)
{
    // -----------------------------------------------------------------------
    //      Rebuilding function for the point cloud, using updated values
    // -----------------------------------------------------------------------
    auto rebuild = [&panel, &qs, &cloudEntity, rootEntity, mainWidget, dotSize, dotCount]() 
    {
        // set the text of the label e.g.: n=1, l=0, m=0
        panel.infoLabel->setText(QString("Computing: n=%1  l=%2  m=%3...").arg(qs.n).arg(qs.l).arg(qs.m));
        
        // disable all the sliders
        panel.sN.slider->setEnabled(false);
        panel.sL.slider->setEnabled(false);
        panel.sM.slider->setEnabled(false);

        // Capture quantum numbers by value for thread safety
        const int cn = qs.n, cl = qs.l, cm = qs.m;

        // create a worker to calculate the new points
        QFuture<std::vector<Point>> future = QtConcurrent::run([cn, cl, cm, dotCount]() {
                return computePoints(cn, cl, cm, dotCount);
        });

        auto *watcher = new QFutureWatcher<std::vector<Point>>(mainWidget);
        watcher->setFuture(future);

        // return to main thread to render the points
        QObject::connect(watcher, &QFutureWatcher<std::vector<Point>>::finished,
                        [&cloudEntity, rootEntity, &panel, &qs, watcher, dotSize, cn, cl, cm]()
        {
            std::vector<Point> points = watcher->result();
            watcher->deleteLater();

            delete cloudEntity;
            cloudEntity = buildCloudFromPoints(points, dotSize, rootEntity);

            panel.infoLabel->setText(QString("Orbital: n=%1  l=%2  m=%3").arg(cn).arg(cl).arg(cm));

            panel.sN.slider->setEnabled(true);
            panel.sL.slider->setEnabled(true);
            panel.sM.slider->setEnabled(true);
        });
    };

    // -----------------------------------------------------------------------
    //                  Updater functions for the sliders
    // -----------------------------------------------------------------------

    //  add a debounce timer to only rebuild 0.5s after the slider move
    auto *debounceTimer = new QTimer(mainWidget);
    debounceTimer->setSingleShot(true);
    debounceTimer->setInterval(500);
    QObject::connect(debounceTimer, &QTimer::timeout, rebuild);

    // -----------------------------------------------------------------------
    //                  Updater functions for the n slider
    // -----------------------------------------------------------------------
    QObject::connect(panel.sN.slider, &QSlider::valueChanged, [&, debounceTimer](int val) 
    {
        // update label of n slider
        qs.n = val;
        panel.sN.valueLabel->setText(QString::number(qs.n));

        // block other updates of the l and m sliders
        panel.sL.slider->blockSignals(true);
        panel.sM.slider->blockSignals(true);

        // update l slider
        panel.sL.slider->setMaximum(qs.n - 1);
        qs.l = std::min(qs.l, qs.n - 1);
        panel.sL.slider->setValue(qs.l);
        panel.sL.valueLabel->setText(QString::number(qs.l));

        // update m slider
        panel.sM.slider->setMinimum(-qs.l);
        panel.sM.slider->setMaximum( qs.l);
        qs.m = std::clamp(qs.m, -qs.l, qs.l);
        panel.sM.slider->setValue(qs.m);
        panel.sM.valueLabel->setText(QString::number(qs.m));

        // unblock the other sliders
        panel.sL.slider->blockSignals(false);
        panel.sM.slider->blockSignals(false);

        debounceTimer->start();
    });

    // -----------------------------------------------------------------------
    //                  Updater functions for the l slider
    // -----------------------------------------------------------------------
    QObject::connect(panel.sL.slider, &QSlider::valueChanged, [&, debounceTimer](int val) 
    {
        // update label of l slider
        qs.l = val;
        panel.sL.valueLabel->setText(QString::number(qs.l));

        // block other updates of the m sliders
        panel.sM.slider->blockSignals(true);

        // update m slider
        panel.sM.slider->setMinimum(-qs.l);
        panel.sM.slider->setMaximum( qs.l);
        qs.m = std::clamp(qs.m, -qs.l, qs.l);
        panel.sM.slider->setValue(qs.m);
        panel.sM.valueLabel->setText(QString::number(qs.m));

        // unblock the m slider
        panel.sM.slider->blockSignals(false);

        debounceTimer->start();
    });

    // -----------------------------------------------------------------------
    //                  Updater functions for the m slider
    // -----------------------------------------------------------------------
    QObject::connect(panel.sM.slider, &QSlider::valueChanged, [&, debounceTimer](int val) 
    {
        // update label of m slider
        qs.m = val;
        panel.sM.valueLabel->setText(QString::number(qs.m));

        debounceTimer->start();
    });
}

#endif // UPDATE_H