set(TRAJECTORY_SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/cones.c
	${CMAKE_CURRENT_SOURCE_DIR}/trajectory.c
	${CMAKE_CURRENT_SOURCE_DIR}/angles.c
	${CMAKE_CURRENT_SOURCE_DIR}/variabledensity.c
	${CMAKE_CURRENT_SOURCE_DIR}/arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/mathops.c
	${CMAKE_CURRENT_SOURCE_DIR}/convertendian.c
	${CMAKE_CURRENT_SOURCE_DIR}/spiral.c
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall -Wextra")

add_library(trajectory ${TRAJECTORY_SOURCE})

add_subdirectory(gridding)
add_subdirectory(phantom)
add_subdirectory(test)
