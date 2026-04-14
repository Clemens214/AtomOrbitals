#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QSlider>
#include <QLabel>

#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QForwardRenderer>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QOrbitCameraController>

#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DRender/QCamera>

#include <QVector3D>
#include <random>
#include <cmath>

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

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // --- 3D Window ---
    auto *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(QColor(30, 30, 40));

    // Root entity
    auto *rootEntity = new Qt3DCore::QEntity();

    // Camera
    Qt3DRender::QCamera *camera = view->camera();
    camera->lens()->setPerspectiveProjection(45.0f, 16.f/9.f, 0.1f, 1000.f);
    camera->setPosition(QVector3D(0, 0, 20));
    camera->setViewCenter(QVector3D(0, 0, 0));

    // Orbit camera controller
    auto *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);
    camController->setCamera(camera);
    camController->setLinearSpeed(20.0f);
    camController->setLookSpeed(180.0f);

    // Initial parameters
    float sigma          = 3.0f;
    const float dotSize  = 0.12f;
    const int   dotCount = 600;

    // Build initial cloud
    Qt3DCore::QEntity *cloudEntity = buildCloud(sigma, dotSize, dotCount, rootEntity);

    view->setRootEntity(rootEntity);

    // --- Embed 3D window into a QWidget ---
    QWidget *container = QWidget::createWindowContainer(view);
    container->setMinimumSize(600, 500);

    // --- Sigma Slider ---
    // Slider range 1..100 maps to sigma 0.1..10.0
    QSlider *slider = new QSlider(Qt::Horizontal);
    slider->setMinimum(1);
    slider->setMaximum(100);
    slider->setValue(static_cast<int>(sigma * 10));

    QLabel *label = new QLabel(QString("σ (std dev): %1").arg(sigma, 0, 'f', 1));
    label->setAlignment(Qt::AlignCenter);

    // Rebuild cloud whenever sigma changes
    QObject::connect(slider, &QSlider::valueChanged,
                     [&cloudEntity, rootEntity, dotSize, dotCount, label](int value) {
                         float newSigma = value / 10.0f;
                         label->setText(QString("σ (std dev): %1").arg(newSigma, 0, 'f', 1));
                         delete cloudEntity;
                         cloudEntity = buildCloud(newSigma, dotSize, dotCount, rootEntity);
                     });

    // --- Layout ---
    QWidget *mainWidget = new QWidget();
    mainWidget->setWindowTitle("Qt3D Gaussian Point Cloud Viewer");

    auto *layout = new QVBoxLayout(mainWidget);
    layout->addWidget(container);
    layout->addWidget(label);
    layout->addWidget(slider);

    mainWidget->resize(640, 580);
    mainWidget->show();

    return app.exec();
}