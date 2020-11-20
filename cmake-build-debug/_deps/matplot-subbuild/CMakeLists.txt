# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.17.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(matplot-populate NONE)

include(ExternalProject)
ExternalProject_Add(matplot-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/alandefreitas/matplotplusplus.git" "GIT_TAG" "origin/master"
                    SOURCE_DIR          "/Users/admin/CLionProjects/HW/cmake-build-debug/_deps/matplot-src"
                    BINARY_DIR          "/Users/admin/CLionProjects/HW/cmake-build-debug/_deps/matplot-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
