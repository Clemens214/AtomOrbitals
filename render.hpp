#ifndef RENDER_H
#define RENDER_H

#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/QPhongMaterial>

#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>

#include <QVector3D>

// include the rejection sampler to find the points
#include "sample.hpp"

// Build the point cloud for a given orbital (n, l, m)
std::vector<Point> computePoints(int n, int l, int m, int dotCount)
{
    Orbital orbit(n, l, m);
    float boxHalf = 5 * n*n;
    float pMax = computePMax(orbit, boxHalf);
    return rejectionSample(dotCount, orbit, pMax, boxHalf);
}

// Run this on the main thread — creates Qt3D entities
Qt3DCore::QEntity *buildCloudFromPoints(const std::vector<Point> &points,
                                         float dotSize,
                                         Qt3DCore::QEntity *parent)
{
    auto *cloudEntity = new Qt3DCore::QEntity(parent);

    auto *mesh = new Qt3DExtras::QSphereMesh();
    mesh->setRadius(dotSize);
    mesh->setRings(4);
    mesh->setSlices(4);

    auto *material = new Qt3DExtras::QPhongMaterial();
    material->setDiffuse(QColor(70, 130, 210));
    material->setSpecular(QColor(255, 255, 255));
    material->setShininess(80.0f);

    for (const Point &p : points) {
        auto *dotEntity = new Qt3DCore::QEntity(cloudEntity);
        auto *transform = new Qt3DCore::QTransform();
        transform->setTranslation(QVector3D(p.x, p.y, p.z));
        dotEntity->addComponent(mesh);
        dotEntity->addComponent(material);
        dotEntity->addComponent(transform);
    }
    return cloudEntity;
}

#endif // RENDER_H