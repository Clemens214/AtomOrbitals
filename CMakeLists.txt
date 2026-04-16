cmake_minimum_required(VERSION 3.16)

project(SphereViewer VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets LinguistTools)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets LinguistTools)
find_package(Qt6 REQUIRED COMPONENTS Widgets 3DCore 3DRender 3DExtras)

set(TS_FILES Project_de_DE.ts)

add_executable(SphereViewer
    main.cpp
    Atom.hpp
    Coordinates.hpp
    Orbital.hpp
    render.hpp
    sample.hpp
    slider.hpp
    setup.hpp
)

target_link_libraries(SphereViewer PRIVATE
    Qt6::Widgets
    Qt6::3DCore
    Qt6::3DRender
    Qt6::3DExtras
    Qt6::Concurrent
)