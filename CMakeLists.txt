PROJECT(mritrajectory)

set(CMAKE_MODULE_PATH
	${CMAKE_MODULE_PATH}
	${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

macro(setCXX11)
	if(CMAKE_VERSION VERSION_LESS "3.1")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	else()
		set (CMAKE_CXX_STANDARD 11)
	endif()
endmacro(setCXX11)

add_subdirectory(src)
if(ENABLE_MATLAB)
	add_subdirectory(matlab)
endif(ENABLE_MATLAB)
add_subdirectory(trajectorydesigner)
