find_package(Qt5 REQUIRED COMPONENTS Test)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gridding)

set(CMAKE_AUTOMOC ON)

add_executable(phantom_test phantomtest.cpp)
target_link_libraries(phantom_test Qt5::Test phantom)
