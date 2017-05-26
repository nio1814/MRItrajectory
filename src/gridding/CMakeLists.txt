find_package(FFTW3 REQUIRED)

set(GRIDDING_SOURCE
	mrdata.cpp
	gridding.cpp
	numericalrecipes.c
	${CMAKE_CURRENT_SOURCE_DIR}/../arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../trajectory.c
	${CMAKE_CURRENT_SOURCE_DIR}/../mathops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../convertendian.c
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "-Wall -Werror")

add_library(gridding ${GRIDDING_SOURCE})
target_link_libraries(gridding ${FFTW3_LIBRARIES})

add_subdirectory(test)
