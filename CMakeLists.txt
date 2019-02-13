set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 2.6)

set(SRC_DIR src)
set(TEST_DIR test)
set(BIN_DIR build)

set(DEPENDENCIES ${SRC_DIR}/mis.c ${SRC_DIR}/read_file.c)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(mis)

find_package(OpenMP COMPONENTS C REQUIRED)

set(LAPACKE_DIR /usr/include/lapacke)
find_package(LAPACKE REQUIRED)

find_package(MPI REQUIRED)
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# Add flag for profile tool gprof
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -lm -std=c99 ")

#add_executable(tests test/test.c ${DEPENDENCIES})
add_executable(bench src/bench.c ${DEPENDENCIES})

target_link_libraries(bench ${LAPACKE_LIBRARIES})
target_include_directories(bench PUBLIC ${LAPACKE_INCLUDE_DIRS})

# Include MPI library
target_link_libraries(bench ${MPI_C_LIBRARIES})
target_include_directories(bench PUBLIC ${MPI_C_INCLUDE_PATH})
