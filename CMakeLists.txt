# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman

cmake_minimum_required(VERSION 3.20)
project(TRIVENI VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

add_library(triveni
  src/Grid.cxx
  src/PPCollision.cxx
  src/Proton.cxx
  src/Random.cxx
  src/Settings.cxx
  src/Steiner.cxx)
target_include_directories(triveni PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(triveni PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(triveni PUBLIC Threads::Threads)
target_compile_definitions(triveni PUBLIC TRIVENI_VERSION="${PROJECT_VERSION}")

foreach(app triveni-proton triveni-pp)
  add_executable(${app} apps/${app}.cxx)
  target_link_libraries(${app} PRIVATE triveni)
  target_compile_options(${app} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()

enable_testing()
foreach(test test_random test_steiner test_proton test_pp)
  add_executable(${test} tests/${test}.cxx)
  target_link_libraries(${test} PRIVATE triveni)
  target_compile_options(${test} PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME ${test} COMMAND ${test})
endforeach()
