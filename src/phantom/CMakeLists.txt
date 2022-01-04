set(SOURCE
	phantom.cpp
	shape.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../mathops.c
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

setCXX11()
add_library(phantom ${SOURCE})
set_target_properties(phantom PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_subdirectory(test)
