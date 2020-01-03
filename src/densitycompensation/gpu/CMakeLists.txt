if(NOT WIN32 OR (WIN32 AND MSVC))
  project(dcfgpu LANGUAGES CXX CUDA)

  find_package(CUDA)

  include_directories(../../trajectory)
  include_directories(../..)

  set(SOURCE
    dcfj.cu
	mp1-util.cu
  )

  cuda_add_executable(${PROJECT_NAME} ${SOURCE})
  target_link_libraries(${PROJECT_NAME} trajectory)
endif()