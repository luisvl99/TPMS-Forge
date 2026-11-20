cmake_minimum_required(VERSION 3.20)
project(tpms LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(tpms src/main.cpp src/mesher.cpp src/export.cpp src/smooth.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(tpms PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MSVC)
  target_compile_options(tpms PRIVATE /W4)
else()
  target_compile_options(tpms PRIVATE -Wall -Wextra)
endif()

# Smoke test: the tool exits non-zero when its own geometry checks fail
# (closed shell, fit error, volume against the field).
enable_testing()
foreach(t gyroid schwarz-p diamond iwp frd fischer-koch-s neovius)
  add_test(NAME smooth-${t} COMMAND tpms --type ${t} --cells 1 -o ${CMAKE_BINARY_DIR}/test-${t}.step)
endforeach()

# Standalone .exe without MinGW runtime DLLs.
if(MINGW)
  target_link_options(tpms PRIVATE -static)
endif()
