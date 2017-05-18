find_package(Qt5 REQUIRED COMPONENTS Test)

set(SOURCE griddingtest.cpp)

set(CMAKE_AUTOMOC ON)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(gridding_test ${SOURCE})
target_link_libraries(gridding_test gridding Qt5::Test)
