PROJECT(mritrajectory)

option(ENABLE_PYTHON "Compile python module" ON)

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

if(MSVC)
  add_definitions(-D_USE_MATH_DEFINES)
#  add_compile_options(/W4 /WX)
  add_compile_options(/W1 /WX)
else()
  add_compile_options(-Wall -Wextra -Werror)
endif()

add_subdirectory(src)
if(ENABLE_MATLAB)
	add_subdirectory(matlab)
endif(ENABLE_MATLAB)
add_subdirectory(trajectorydesigner)
if(ENABLE_PYTHON)
  add_subdirectory(python)
endif(ENABLE_PYTHON)