set(SOURCE phantom.cpp ellipsoid.cpp)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

set(CMAKE_CXX_FLAGS "-Wall -Werror")

add_library(phantom ${SOURCE})
