if(NOT WIN32 OR (WIN32 AND MSVC))
  find_package(CUDA)

  if(CUDA_FOUND)
    project(dcfgpu LANGUAGES CXX CUDA)
    
    include_directories(../../trajectory)
    include_directories(../..)
    
    set(SOURCE
    dcfj.cu
    mp1-util.cu
    )
    
    cuda_add_executable(${PROJECT_NAME} ${SOURCE})
    target_link_libraries(${PROJECT_NAME} trajectory)
  endif()
endif()
