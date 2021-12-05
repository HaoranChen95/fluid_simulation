cmake_minimum_required(VERSION 3.10)

project(fluid_simulation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER icpc)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_compile_options(-msse2 -qopt-report1 -qopt-report-phase=vec -qopt-report-phase=openmp -shared-intel -qopenmp -march=native -fp-model precise -std=c++11 -O2 -mcmodel=medium)

FILE(GLOB_RECURSE SRC_FILES "src/*.cpp")
include_directories(
	"${CMAKE_SOURCE_DIR}/include"
	"${CMAKE_SOURCE_DIR}/include/initialization"
)

ADD_EXECUTABLE(main "${SRC_FILES}")
