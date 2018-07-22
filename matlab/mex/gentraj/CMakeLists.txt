#cmake_minimum_required(VERSION 2.8.12.2)

#set_matlab_env()
SET(CMAKE_SKIP_BUILD_RPATH TRUE)

set(SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/gentraj.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/../../../src/trajectory/trajectorygenerator.cpp
	trajectorygenerator.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../src/trajectory)
include_directories(${PROJECT_SOURCE_DIR}/src/trajectory)
include_directories(${MATLAB_INCLUDE_DIR})

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "CXXFLAGS='$$CXXFLAGS -std=c++11 -Werror -Wall -Wextra'")

if (MATLAB_FOUND)
	add_executable(gentraj ${SOURCE})
	target_link_libraries(gentraj trajectory)
endif()
