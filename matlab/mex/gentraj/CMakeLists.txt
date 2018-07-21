#cmake_minimum_required(VERSION 2.8.12.2)

set_matlab_env()
SET(CMAKE_SKIP_BUILD_RPATH TRUE)

set(C_DIR ${PROJECT_SOURCE_DIR}/src)
set(SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/gentraj.cpp
	${C_DIR}/trajectory/trajectorygenerator.cpp
)

include_directories(${C_DIR})
include_directories(${C_DIR}/trajectory)
include_directories(${MATLAB_INCLUDE_DIR})

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "CXXFLAGS='$$CXXFLAGS -std=c++11 -Werror -Wall -Wextra'")

if (MATLAB_FOUND)
	add_executable(gentraj ${SOURCE})
	target_link_libraries(gentraj trajectory)
endif()
