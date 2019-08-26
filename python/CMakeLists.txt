project(cmri)

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")

add_subdirectory(pybind11)

set(SOURCE_TRAJECTORY
  trajectory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/trajectory/trajectorygenerator.cpp
)
pybind11_add_module(_trajectory ${SOURCE_TRAJECTORY})
target_link_libraries(_trajectory trajectory)
#set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
