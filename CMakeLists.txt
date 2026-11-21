cmake_minimum_required(VERSION 3.20)
project(tpms LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TPMS_GUI "Build tpms-gui (downloads Dear ImGui and GLFW)" ON)

# Core shared by the CLI and the GUI.
add_library(tpmscore STATIC src/run.cpp src/mesher.cpp src/export.cpp src/smooth.cpp)
target_include_directories(tpmscore PUBLIC src)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(tpmscore PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(tpms src/main.cpp)
target_link_libraries(tpms PRIVATE tpmscore)
set(tpms_targets tpmscore tpms)

if(TPMS_GUI)
  include(FetchContent)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(glfw
    URL https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.zip
    URL_HASH SHA256=b5ec004b2712fd08e8861dc271428f048775200a2df719ccf575143ba749a3e9
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_Declare(imgui
    URL https://github.com/ocornut/imgui/archive/refs/tags/v1.91.5.tar.gz
    URL_HASH SHA256=2aa2d169c569368439e5d5667e0796d09ca5cc6432965ce082e516937d7db254
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(glfw imgui)

  add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp)
  target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)
  find_package(OpenGL REQUIRED)
  target_link_libraries(imgui PUBLIC glfw OpenGL::GL)

  add_executable(tpms-gui WIN32 src/app.cpp)
  target_link_libraries(tpms-gui PRIVATE tpmscore imgui)
  if(WIN32)
    target_link_libraries(tpms-gui PRIVATE comdlg32)
  endif()
  list(APPEND tpms_targets tpms-gui)
endif()

foreach(t ${tpms_targets})
  if(MSVC)
    target_compile_options(${t} PRIVATE /W4)
  else()
    target_compile_options(${t} PRIVATE -Wall -Wextra)
  endif()
  # Standalone .exe without MinGW runtime DLLs.
  if(MINGW AND NOT t STREQUAL tpmscore)
    target_link_options(${t} PRIVATE -static)
  endif()
endforeach()

# Smoke test: the tool exits non-zero when its own geometry checks fail
# (closed shell, fit error, volume against the field).
enable_testing()
foreach(t gyroid schwarz-p diamond iwp frd fischer-koch-s neovius)
  add_test(NAME smooth-${t} COMMAND tpms --type ${t} --cells 1 -o ${CMAKE_BINARY_DIR}/test-${t}.step)
endforeach()
