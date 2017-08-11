find_package(Qt5 REQUIRED COMPONENTS Widgets)
set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Widgets)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(SOURCE main.cpp mainwindow.cpp generator.cpp)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall -Wextra")

add_executable(trajectorydesigner ${SOURCE})
target_link_libraries(trajectorydesigner ${QT_LIBRARIES})
