#ifndef RENDER_H
#define RENDER_H

#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/QPhongMaterial>

#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>

#include <QVector3D>
#include <random>

#include "Orbital.hpp"

// --- Generate 'count' points from a 3D isotropic Gaussian(0, sigma) ---
QList<QVector3D> gaussianPoints(int count, float sigma)
{
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::normal_distribution<float> dist(0.0f, sigma);

    QList<QVector3D> points;
    points.reserve(count);
    for (int i = 0; i < count; ++i)
        points.append(QVector3D(dist(rng), dist(rng), dist(rng)));
    return points;
}

// --- Builds the point cloud entity under a given parent ---
Qt3DCore::QEntity *buildCloud(float sigma, float dotSize, int dotCount,
                              Qt3DCore::QEntity *parent)
{
    auto *cloudEntity = new Qt3DCore::QEntity(parent);

    // Shared mesh and material for all dots
    auto *mesh = new Qt3DExtras::QSphereMesh();
    mesh->setRadius(dotSize);
    mesh->setRings(8);
    mesh->setSlices(8);

    auto *material = new Qt3DExtras::QPhongMaterial();
    material->setDiffuse(QColor(70, 130, 210));
    material->setSpecular(QColor(255, 255, 255));
    material->setShininess(80.0f);

    const QList<QVector3D> points = gaussianPoints(dotCount, sigma);
    for (const QVector3D &pos : points) {
        auto *dotEntity = new Qt3DCore::QEntity(cloudEntity);

        auto *transform = new Qt3DCore::QTransform();
        transform->setTranslation(pos);

        dotEntity->addComponent(mesh);
        dotEntity->addComponent(material);
        dotEntity->addComponent(transform);
    }

    return cloudEntity;
}

#endif // RENDER_H