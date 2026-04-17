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
    float rMax = computeRMax(orbit, boxHalf, pMax);
    return rejectionSample(dotCount, orbit, pMax, rMax);
}

// Map a normalized value [0, 1] to a heatmap (blue -> cyan -> green -> yellow -> red)
QColor heatmapColor(float t)
{
    t = std::clamp(t, 0.0f, 1.0f);
    float r, g, b;
    if (t < 0.25f) {
        float s = t / 0.25f;
        r = 0; g = s; b = 1;
    } else if (t < 0.5f) {
        float s = (t - 0.25f) / 0.25f;
        r = 0; g = 1; b = 1.0f - s;
    } else if (t < 0.75f) {
        float s = (t - 0.5f) / 0.25f;
        r = s; g = 1; b = 0;
    } else {
        float s = (t - 0.75f) / 0.25f;
        r = 1; g = 1.0f - s; b = 0;
    }
    return QColor::fromRgbF(r, g, b);
}

// Run this on the main thread — creates Qt3D entities
Qt3DCore::QEntity *buildCloudFromPoints(const std::vector<Point> &points, float dotSize, Qt3DCore::QEntity *parent)
{
    auto *cloudEntity = new Qt3DCore::QEntity(parent);

    auto *mesh = new Qt3DExtras::QSphereMesh();
    mesh->setRadius(dotSize);
    mesh->setRings(4);
    mesh->setSlices(4);

    // Find the value range for normalization
    double vMin = std::numeric_limits<double>::max();
    double vMax = std::numeric_limits<double>::lowest();
    for (const Point &p : points) {
        vMin = std::min(vMin, p.value);
        vMax = std::max(vMax, p.value);
    }
    double vRange = (vMax > vMin) ? (vMax - vMin) : 1.0;

    // Generate the color buckets
    constexpr int NUM_BUCKETS = 16;
    std::array<Qt3DExtras::QPhongMaterial*, NUM_BUCKETS> palette;
    for (int i = 0; i < NUM_BUCKETS; ++i) {
        float t = (i + 0.5f) / NUM_BUCKETS;
        auto *mat = new Qt3DExtras::QPhongMaterial(cloudEntity);
        mat->setDiffuse(heatmapColor(t));
        mat->setSpecular(QColor(0, 0, 0));
        mat->setShininess(0.0f);
        palette[i] = mat;
    }

    for (const Point &p : points) {
        float t = static_cast<float>((p.value - vMin) / vRange);
        int bucket = std::clamp((int)(t * NUM_BUCKETS), 0, NUM_BUCKETS - 1);

        auto *dotEntity = new Qt3DCore::QEntity(cloudEntity);
        auto *transform = new Qt3DCore::QTransform();
        transform->setTranslation(QVector3D(p.x, p.y, p.z));
        dotEntity->addComponent(mesh);
        dotEntity->addComponent(palette[bucket]);
        dotEntity->addComponent(transform);
    }
    return cloudEntity;
}

#endif // RENDER_H