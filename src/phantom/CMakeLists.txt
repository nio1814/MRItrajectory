set(SOURCE
	phantom.cpp
	ellipsoid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../mathops.c
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

set(CMAKE_CXX_FLAGS "-Wall -Werror")

add_library(phantom ${SOURCE})

add_subdirectory(test)
