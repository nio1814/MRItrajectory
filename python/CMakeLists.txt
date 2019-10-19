project(cmri)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE
  cmri.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory/trajectorygenerator.cpp
)

add_subdirectory(pybind11)

pybind11_add_module(${PROJECT_NAME} ${SOURCE})
target_link_libraries(${PROJECT_NAME} trajectory)
#set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
