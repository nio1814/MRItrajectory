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

add_library(gridding ${GRIDDING_SOURCE})

add_subdirectory(test)
