set(CONES_SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/cones.c
	${CMAKE_CURRENT_SOURCE_DIR}/trajectory.c
	${CMAKE_CURRENT_SOURCE_DIR}/angles.c
	${CMAKE_CURRENT_SOURCE_DIR}/variabledensity.c
	${CMAKE_CURRENT_SOURCE_DIR}/arrayops.c
	${CMAKE_CURRENT_SOURCE_DIR}/mathops.c
	${CMAKE_CURRENT_SOURCE_DIR}/convertendian.c
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall -Wextra")

add_library(cones ${CONES_SOURCE})
