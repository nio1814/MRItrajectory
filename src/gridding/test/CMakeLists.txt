find_package(Qt5 REQUIRED COMPONENTS Test)

set(SOURCE griddingtest.cpp)

set(CMAKE_AUTOMOC ON)

add_executable(gridding_test ${SOURCE})
target_link_libraries(gridding_test Qt5::Test)
