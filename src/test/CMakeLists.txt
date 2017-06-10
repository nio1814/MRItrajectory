cmake_minimum_required(VERSION 3.5)

find_package(Qt5 REQUIRED COMPONENTS Test)
set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Test)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../phantom)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../gridding)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

set(CMAKE_AUTOMOC ON)

macro(set_test TEST_NAME TEST_SRC)
	add_executable(${TEST_NAME} ${TEST_SRC})
	target_link_libraries(${TEST_NAME} ${QT_LIBRARIES} trajectory gridding phantom)
	add_test(${TEST_NAME} ${TEST_NAME})
endmacro()

set_test(spiral_test spiraltest.cpp)
set_test(radial_test radialtest.cpp)
