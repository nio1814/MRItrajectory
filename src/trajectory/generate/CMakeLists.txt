include_directories("../../../src/trajectory")

add_executable(generatecones generatecones.cpp)
target_link_libraries(generatecones trajectory)
