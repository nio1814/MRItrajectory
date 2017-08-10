find_package(Qt5 REQUIRED COMPONENTS Widgets)
set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Widgets)

set(SOURCE main.cpp mainwindow.cpp generator.cpp)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_CXX_STANDARD 11)

add_executable(trajectorydesigner ${SOURCE})
target_link_libraries(trajectorydesigner ${QT_LIBRARIES})
