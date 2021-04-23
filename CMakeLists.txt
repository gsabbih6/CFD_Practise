cmake_minimum_required(VERSION 3.17)
project(HW)

set(CMAKE_CXX_STANDARD 20)

#add_executable(HW main.cpp)
#cmake_minimum_required(VERSION 3.16)
#project(Hw2_CFD)
#
#set(CMAKE_CXX_STANDARD 17)

add_executable(HW main.cpp ShockTubeLaxWendroff.cpp ThomasAlgorithm.cpp FiniteDifference.cpp CombinedMethod.cpp ViscousBurgers.cpp )
include(get_cpm.cmake)
CPMAddPackage(
        NAME matplot
        GITHUB_REPOSITORY alandefreitas/matplotplusplus GIT_TAG origin/master
)
#CPMAddPackage(
#        NAME Kokkos
#        GITHUB_REPOSITORY kokkos/kokkos GIT_TAG origin/master
#)

find_package(Boost 1.57.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(BoostTest main.cpp)
target_link_libraries(BoostTest ${Boost_LIBRARIES})

# ...
target_link_libraries(HW PUBLIC matplot)
#target_link_libraries(HW PUBLIC Kokkos::kokkos)