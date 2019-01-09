set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 2.6)

project(mis)

find_package(OpenMP COMPONENTS C REQUIRED)
find_package(GSL REQUIRED)
set(CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(mis mis.c)
target_link_libraries(mis m ${GSL_LIBRARIES})
target_include_directories(mis PUBLIC ${GSL_INCLUDE_DIRS})