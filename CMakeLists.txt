cmake_minimum_required(VERSION 2.8.7)
project(ekf_core)


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_MODULE_PATH   ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system date_time)


# SSE flags
add_definitions("-DENABLE_SSE")

# Also add some useful compiler flag
set(CMAKE_CXX_FLAGS
   "${CMAKE_CXX_FLAGS} ${SSE_FLAGS} -O2 -march=native  -std=c++0x"
)

# NEON flags
#add_definitions("-DENABLE_NEON")

# Also add some useful compiler flag
#set(CMAKE_CXX_FLAGS
#   "${CMAKE_CXX_FLAGS} -O2 -march=native -mfpu=neon -std=c++0x"
#) 

# Set source files
set(ekf_SOURCE_FILES
  ${PROJECT_SOURCE_DIR}/src/EstimatorBase.cpp
  ${PROJECT_SOURCE_DIR}/src/EstimatorDelayHider.cpp
)
set(SOURCE_FILES
  ${ekf_SOURCE_FILES}
)

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src
)


# build shared library.
add_library(ekf_lib SHARED ${SOURCE_FILES})

target_link_libraries(ekf_lib ${Boost_LIBRARIES})


# build test program
add_executable(ekf_test src/test/main_test.cpp)
target_link_libraries(ekf_test ekf_lib)
