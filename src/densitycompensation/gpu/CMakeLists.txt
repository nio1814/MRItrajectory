project(dcfgpu LANGUAGES CXX CUDA)

find_package(CUDA)

set(SOURCE
  dcfj.cu
)

cuda_add_executable(${PROJECT_NAME} ${SOURCE})
