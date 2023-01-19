cmake_minimum_required(VERSION 3.5)
project (sshe-worker)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread -Wall -Wextra -Wl,--no-as-needed -funroll-loops -march=native -flto -pedantic -std=c++14 -fopenmp")

find_package(SCI REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../MPC/EzPC/SCI/build/install" NO_DEFAULT_PATH)
if (NOT SCI_FOUND)
    message(STATUS "SCI was not found")
endif()
include_directories()

# Paillier encryption dependency
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake/OPHE/cmake")
find_package(NTL REQUIRED)
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)
find_package(OPHELib REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP)

# emp two-party GC dependency
find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)
include(${CMAKE_FOLDER}/cmake/enable_rdseed.cmake)
# include(${CMAKE_FOLDER}/cmake/enable_float.cmake)

find_package(emp-sh2pc REQUIRED)
include_directories(${EMP-SH2PC_INCLUDE_DIRS})

# include_directories("${CMAKE_SOURCE_DIR}/FedGraph/include/")

set(Boost_USE_STATIC_LIBS        ON) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)

find_package(Boost  COMPONENTS system thread)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR}) 
  #message( "Found Boost at ${Boost_LIBRARIES}")
else()
   message(FATAL_ERROR "Failed to find boost at " ${Boost_Lib_Dirs} " Need system thread")
endif()

include_directories(include)
FILE(GLOB SRCS src/*.cpp)
add_library(TaskHandler ${SRCS})
target_compile_definitions(TaskHandler PRIVATE -DSSHEBACKEND)
target_link_libraries(TaskHandler cryptoTools ${EMP-OT_LIBRARIES} OPHELib::ophelib SCI::SCI-FloatingPoint)

# set(BOOST_ROOT "${CMAKE_SOURCE_DIR}/thirdparty/linux/boost/")
# find_package (Boost COMPONENTS thread system REQUIRED)

# FILE(GLOB_RECURSE TEST_COMM_SRCS ${CMAKE_SOURCE_DIR}/SSHE-Worker/test/*.cpp)
# add_executable(test-comm ${TEST_COMM_SRCS})
# target_compile_definitions(test-comm PRIVATE -DSSHEBACKEND)
# target_include_directories(test-comm PRIVATE "${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")
# target_include_directories(test-comm PRIVATE "${Boost_INCLUDE_DIRS}")
# target_link_libraries(test-comm PRIVATE ${Boost_LIBRARIES} cryptoTools TaskHandler ${EMP-OT_LIBRARIES} OPHELib::ophelib)

FILE(GLOB_RECURSE TEST_COMM_SRCS ${CMAKE_SOURCE_DIR}/SSHE-Worker/test/test.cpp)
add_executable(test-comm ${TEST_COMM_SRCS})
target_compile_definitions(test-comm PUBLIC -DSSHEBACKEND)
include_directories("${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")
target_link_libraries(test-comm TaskHandler)

FILE(GLOB_RECURSE TEST_CLIENT_SRC ${CMAKE_SOURCE_DIR}/SSHE-Worker/test/client_test.cpp)
add_executable(test-client ${TEST_CLIENT_SRC})
target_compile_definitions(test-client PUBLIC -DSSHEBACKEND)
include_directories("${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")
target_link_libraries(test-client TaskHandler)

FILE(GLOB_RECURSE TEST_SERVER_SRC ${CMAKE_SOURCE_DIR}/SSHE-Worker/test/server_test.cpp)
add_executable(test-server ${TEST_SERVER_SRC})
target_compile_definitions(test-server PUBLIC -DSSHEBACKEND)
include_directories("${CMAKE_SOURCE_DIR}/SSHE-Worker/include/")
target_link_libraries(test-server TaskHandler)

