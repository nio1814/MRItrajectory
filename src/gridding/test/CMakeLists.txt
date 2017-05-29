find_package(Qt5 REQUIRED COMPONENTS Test)

set(SOURCE griddingtest.cpp)

set(CMAKE_AUTOMOC ON)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_executable(gridding_test ${SOURCE})
target_link_libraries(gridding_test gridding Qt5::Test)

add_executable(mrdata_test mrdatatest.cpp)
target_link_libraries(mrdata_test gridding Qt5::Test)

