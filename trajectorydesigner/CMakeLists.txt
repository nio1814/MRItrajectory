cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake ${CMAKE_MODULE_PATH})
message(STATUS ${CMAKE_MODULE_PATH})

find_package(Qt5 REQUIRED COMPONENTS Widgets)
set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Widgets)
find_package(Qwt REQUIRED)
include_directories(${QWT_INCLUDE_DIR})
#include_directories(${Qwt_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/phantom)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/gridding)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory)

set(SOURCE
	main.cpp
	mainwindow.cpp
	generator.cpp
	plot.cpp
	timeseriesplot.cpp
	plot2d.cpp
	variabledensitydesigner.cpp
	variabledensityplot.cpp
	phantomreconstruction.cpp
)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(CMAKE_VERSION VERSION_LESS "3.1")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	set (CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")
if(CMAKE_VERSION VERSION_LESS "2.8.12")
	QT5_WRAP_UI(trajectorydesigner_UI mainwindow.ui)
endif()


add_executable(trajectorydesigner ${SOURCE} ${trajectorydesigner_UI})
target_link_libraries(trajectorydesigner ${QT_LIBRARIES} ${QWT_LIBRARIES} trajectory phantom gridding)
