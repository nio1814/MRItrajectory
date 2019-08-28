cmake_minimum_required(VERSION 3.10)

project(cmri)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(pybind11)

set(SOURCE_TRAJECTORY
  trajectory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory/trajectorygenerator.cpp
)
pybind11_add_module(_trajectory ${SOURCE_TRAJECTORY})
target_link_libraries(_trajectory PRIVATE trajectory)
#set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
set_target_properties(_trajectory PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}/trajectory/)
#$<${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}/trajectory:_trajectory>

configure_file(setup.py ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cmri DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
