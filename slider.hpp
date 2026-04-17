#ifndef SLIDER_H
#define SLIDER_H

#include <QWidget>
#include <Qt3DCore/QEntity>
#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>

#include <cctype>
#include <cmath>

struct LabelledSlider {
    QLabel  *nameLabel;
    QSlider *slider;
    QLabel  *valueLabel;
};

// The bundle of slider widgets and the info label
struct SliderPanel {
    LabelledSlider sN;
    LabelledSlider sL;
    LabelledSlider sM;
    QLabel        *infoLabel = nullptr;
    QGroupBox     *groupBox  = nullptr;
};

// Create a slider
LabelledSlider makeLabelledSlider(const QString &name,
                                  int min, int max, int value,
                                  QGridLayout *grid, int row)
{
    auto *nameLabel     = new QLabel(name);
    auto *slider        = new QSlider(Qt::Horizontal);
    auto *valueLabel    = new QLabel(QString::number(value));

    nameLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    valueLabel->setAlignment(Qt::AlignLeft  | Qt::AlignVCenter);
    valueLabel->setMinimumWidth(24);

    // set the values
    slider->setMinimum(min);
    slider->setMaximum(max);
    slider->setValue(value);

    grid->addWidget(nameLabel   , row   , 0);
    grid->addWidget(slider      , row   , 1);
    grid->addWidget(valueLabel  , row   , 2);

    return { nameLabel, slider, valueLabel };
}

// Create the quantum-number slider panel and return it
SliderPanel createSliderPanel(const int n, const int l, const int m)
{
    SliderPanel panel;

    panel.groupBox = new QGroupBox("Quantum Numbers");
    auto *grid     = new QGridLayout(panel.groupBox);
    grid->setColumnStretch(1, 1);

    panel.sN = makeLabelledSlider("n", 1    , 5     , n, grid, 0);
    panel.sL = makeLabelledSlider("l", 0    , n-1   , l, grid, 1);
    panel.sM = makeLabelledSlider("m", -l   , l     , m, grid, 2);

    panel.infoLabel = new QLabel(QString("Orbital: n=%1  l=%2  m=%3").arg(n).arg(l).arg(m));
    panel.infoLabel->setAlignment(Qt::AlignCenter);

    return panel;
}


#endif // SLIDER_H