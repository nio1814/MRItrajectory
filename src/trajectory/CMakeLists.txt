INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

set(TRAJECTORY_SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/cones.c
	${CMAKE_CURRENT_SOURCE_DIR}/trajectory.c
	${CMAKE_CURRENT_SOURCE_DIR}/angles.c
	${CMAKE_CURRENT_SOURCE_DIR}/variabledensity.c
	${CMAKE_CURRENT_SOURCE_DIR}/../arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../mathops.c
	${CMAKE_CURRENT_SOURCE_DIR}/../convertendian.c
	${CMAKE_CURRENT_SOURCE_DIR}/spiral.c
	${CMAKE_CURRENT_SOURCE_DIR}/radial.c
	${CMAKE_CURRENT_SOURCE_DIR}/mrgradient.c
  ${CMAKE_CURRENT_SOURCE_DIR}/rings.c
  ${CMAKE_CURRENT_SOURCE_DIR}/phaseencode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/spinwarp.c
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(trajectory ${TRAJECTORY_SOURCE})

add_subdirectory(generate)
