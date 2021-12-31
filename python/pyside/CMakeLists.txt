project(pycmri)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
find_package(PySide2)
find_package(Qt5 REQUIRED COMPONENTS Core)
find_package(FFTW3 REQUIRED)

if(PySide2_FOUND)
  get_target_property(PYSIDE2_INCLUDE_DIR PySide2::pyside2 INTERFACE_INCLUDE_DIRECTORIES)

  find_package(Shiboken2 REQUIRED)
  set(SHIBOKEN2 Shiboken2::shiboken2)
  set(SHIBOKEN2_LIBRARY Shiboken2::libshiboken)
  get_target_property(SHIBOKEN_INCLUDE_DIR ${SHIBOKEN2_LIBRARY} INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "shiboken include ${SHIBOKEN_INCLUDE_DIR}")
else(PySide2_FOUND)
  execute_process(
    COMMAND python -c "import PySide2; print(PySide2.__path__[0])"
    OUTPUT_VARIABLE PYSIDE2_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  find_path(PYSIDE_TYPESYSTEMS
    PATHS ${PYSIDE2_DIR}/typesystems
    NAMES "typesystem_core.xml"
  )

  message(STATUS ${PYSIDE_TYPESYSTEMS})

  macro(pyside2_config option output_var)
    if(${ARGC} GREATER 2)
      set(is_list ${ARGV2})
    else()
      set(is_list "")
    endif()

    execute_process(
      COMMAND ${Python3_EXECUTABLE} "${PYSIDE2_DIR}/examples/utils/pyside2_config.py" ${option}
      OUTPUT_VARIABLE ${output_var}
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if ("${${output_var}}" STREQUAL "")
      message(FATAL_ERROR "Error: Calling pyside2_config.py ${option} returned no output.")
    endif()
    if(is_list)
      string (REPLACE " " ";" ${output_var} "${${output_var}}")
    endif()
  endmacro()

  pyside2_config(--shiboken2-generator-path SHIBOKEN2_GENERATOR_DIR)
  set(SHIBOKEN2 ${SHIBOKEN2_GENERATOR_DIR}/shiboken2${CMAKE_EXECUTABLE_SUFFIX})
  pyside2_config(--shiboken2-generator-include-path SHIBOKEN_INCLUDE_DIR)
  pyside2_config(--pyside2-include-path PYSIDE2_INCLUDE_DIR)
  pyside2_config(--shiboken2-module-shared-libraries-cmake SHIBOKEN2_LIBRARY)
endif(PySide2_FOUND)

message(STATUS "shiboken library ${SHIBOKEN2_LIBRARY}")

set(SHIBOKEN_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src/phantom
  ${CMAKE_SOURCE_DIR}/src/trajectory
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/gridding
  ${Python3_INCLUDE_DIRS}
  ${FFTW3_INCLUDE_DIR}
  ${Python3_NumPy_INCLUDE_DIRS}
)

# colon separated format for shiboken

set(SHIBOKEN_INCLUDES "${SHIBOKEN_INCLUDE_DIRS}")
list(TRANSFORM SHIBOKEN_INCLUDES PREPEND -I)
message(STATUS ${SHIBOKEN_INCLUDES})

set(SHIBOKEN_OPTIONS
  --generator-set=shiboken
  --enable-parent-ctor-heuristic
  --enable-return-value-heuristic
  --use-isnull-as-nb_nonzero
  --avoid-protected-hack
  --no-suppress-warnings
  -T${CMAKE_SOURCE_DIR}
  -T${PYSIDE_TYPESYSTEMS}
  --output-directory=${CMAKE_CURRENT_BINARY_DIR}
)

list(APPEND SHIBOKEN_OPTIONS "${SHIBOKEN_INCLUDES}")
message(STATUS ${SHIBOKEN_OPTIONS})

set(MODULE cmri)

set(PHANTOM_CLASSES
  Shape
  Phantom
)

set(TRAJECTORY_GENERATOR_CLASSES
  TrajectoryGenerator
  Trajectory
)

set(GRIDDING_CLASSES
  Gridding
  MRdata
)

set(CLASSES_WRAPPED
  ${PHANTOM_CLASSES}
  ${TRAJECTORY_GENERATOR_CLASSES}
  ${GRIDDING_CLASSES}
)

list(APPEND CLASSES_WRAPPED "${MODULE}_module")

list(TRANSFORM CLASSES_WRAPPED TOLOWER)
list(TRANSFORM CLASSES_WRAPPED APPEND _wrapper.cpp)
list(TRANSFORM CLASSES_WRAPPED PREPEND ${CMAKE_CURRENT_BINARY_DIR}/${MODULE}/)

set(TYPE_SYSTEM_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/pycmri.h)
set(TYPE_SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/pycmri.xml)

add_definitions(-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH)
add_custom_target(
  wrapcmri
  COMMAND ${SHIBOKEN2} ${SHIBOKEN_OPTIONS} ${TYPE_SYSTEM_HEADER} ${TYPE_SYSTEM}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  BYPRODUCTS ${CLASSES_WRAPPED}
  DEPENDS gridding.xml
)

add_library(${PROJECT_NAME} MODULE ${CLASSES_WRAPPED}
  ${CMAKE_SOURCE_DIR}/src/trajectory/trajectorygenerator.cpp
)

add_dependencies(${PROJECT_NAME} wrapcmri)
target_include_directories(${PROJECT_NAME} PRIVATE
  ${SHIBOKEN_INCLUDE_DIR}
  ${SHIBOKEN_INCLUDE_DIRS}
  ${PYSIDE2_INCLUDE_DIR}
  ${PYSIDE2_INCLUDE_DIR}/QtCore
  ${Qt5Core_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(${PROJECT_NAME} PRIVATE
  ${SHIBOKEN2_LIBRARY}
  Qt5::Core
  ${Python3_LIBRARIES}
  phantom
  trajectory
  gridding
)

if(NOT WIN32)
  target_compile_options(${PROJECT_NAME} PRIVATE
    -Wno-deprecated-copy
    -Wno-missing-field-initializers
    -Wno-unused-parameter
    -Wno-cast-function-type
  )
endif(NOT WIN32)

set_target_properties(
  ${PROJECT_NAME} PROPERTIES
  OUTPUT_NAME "${MODULE}${PYTHON_EXTENSION_SUFFIX}"
  PREFIX ""
)

set(SOURCE_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/out)
set(PYSIDE_INSTALL_DIR cmri-pyside)
set(PACKAGE_INSTALL_DIR ${PYSIDE_INSTALL_DIR}/${MODULE})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
  DESTINATION ${PYSIDE_INSTALL_DIR}
)

install(DIRECTORY ${SOURCE_PACKAGE_DIR}/
  DESTINATION ${PACKAGE_INSTALL_DIR}
)

install(FILES $<TARGET_FILE:${PROJECT_NAME}>
  DESTINATION ${PACKAGE_INSTALL_DIR}
)
