cmake_minimum_required(VERSION 3.0.0)

project(generatecommandline VERSION 0.1 LANGUAGES CXX)
set(OUTPUT_NAME generate)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5Core)

include_directories(
  ${CMAKE_SOURCE_DIR}/src/trajectory
  ${CMAKE_SOURCE_DIR}/src
)

add_executable(${OUTPUT_NAME} "main.cpp"
  ${CMAKE_SOURCE_DIR}/src/trajectory/trajectorygenerator.cpp
)

target_link_libraries(${OUTPUT_NAME} Qt5::Core trajectory)
