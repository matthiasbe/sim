set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 2.6)

set(SRC_DIR src)
set(TEST_DIR test)
set(BIN_DIR build)

set(DEPENDENCIES ${SRC_DIR}/mis.c ${SRC_DIR}/read_file.c)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(mis)

find_package(OpenMP COMPONENTS C REQUIRED)
find_package(GSL REQUIRED)
set(LAPACKE_DIR /usr/include/lapacke)
find_package(LAPACKE REQUIRED)

set(CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -pg")

add_executable(tests test/test.c ${DEPENDENCIES})
add_executable(bench src/bench.c ${DEPENDENCIES})

target_link_libraries(tests m ${GSL_LIBRARIES})
target_link_libraries(tests m ${LAPACKE_LIBRARIES})

target_include_directories(tests PUBLIC ${GSL_INCLUDE_DIRS})
target_include_directories(tests PUBLIC ${LAPACKE_INCLUDE_DIRS})
target_link_libraries(bench m ${GSL_LIBRARIES})
target_link_libraries(bench ${LAPACKE_LIBRARIES})
target_include_directories(bench PUBLIC ${GSL_INCLUDE_DIRS})
target_include_directories(bench PUBLIC ${LAPACKE_INCLUDE_DIRS})


