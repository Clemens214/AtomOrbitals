#ifndef SLIDER_H
#define SLIDER_H

#include <QSlider>
#include <QLabel>

#include <cctype>
#include <cmath>

// -----------------------------------------------------------------------
// Helper: create a labelled slider row inside a grid layout
// -----------------------------------------------------------------------
struct LabelledSlider {
    QLabel  *nameLabel;   // e.g. "n"
    QSlider *slider;
    QLabel  *valueLabel;  // shows current value
};

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

#endif // SLIDER_H