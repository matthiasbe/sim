set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 2.6)

set(SRC_DIR src)
set(TEST_DIR test)
set(BIN_DIR build)

set(DEPENDENCIES ${SRC_DIR}/mis.c ${SRC_DIR}/read_file.c)

project(mis)

find_package(OpenMP COMPONENTS C REQUIRED)
find_package(GSL REQUIRED)
set(CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(tests test/test.c ${DEPENDENCIES})
add_executable(bench src/bench.c ${DEPENDENCIES})

target_link_libraries(tests m ${GSL_LIBRARIES})
target_include_directories(tests PUBLIC ${GSL_INCLUDE_DIRS})
target_link_libraries(bench m ${GSL_LIBRARIES})
target_include_directories(bench PUBLIC ${GSL_INCLUDE_DIRS})
