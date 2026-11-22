cmake_minimum_required(VERSION 3.20)
project(tpms LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TPMS_GUI "Build tpms-gui (downloads Dear ImGui and GLFW)" ON)

# Core shared by the CLI and the GUI.
add_library(tpmscore STATIC src/run.cpp src/report.cpp src/mesher.cpp src/export.cpp src/smooth.cpp)
target_include_directories(tpmscore PUBLIC src)
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND MINGW)
  # The imported target links libgomp.dll.a. With -static (below) these link
  # libgomp.a and the libdl.a it needs, so the .exe needs no MinGW DLL.
  target_compile_options(tpmscore PUBLIC -fopenmp)
  target_link_libraries(tpmscore PUBLIC gomp dl pthread)
elseif(OpenMP_CXX_FOUND)
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

# STEP tests are also read back by OpenCASCADE (tools/check_step.py: valid
# closed solids) when a Python with cadquery-ocp is found: the project's
# .venv, or set TPMS_PYTHON.
if(NOT TPMS_PYTHON)
  foreach(p ${CMAKE_SOURCE_DIR}/.venv/Scripts/python.exe ${CMAKE_SOURCE_DIR}/.venv/bin/python)
    if(EXISTS ${p})
      set(TPMS_PYTHON ${p})
    endif()
  endforeach()
endif()
set(occ OFF)
if(TPMS_PYTHON)
  execute_process(COMMAND ${TPMS_PYTHON} -c "import OCP" RESULT_VARIABLE no_occ OUTPUT_QUIET ERROR_QUIET)
  if(NOT no_occ)
    set(occ ON)
  endif()
endif()
message(STATUS "OpenCASCADE checks: ${occ} (${TPMS_PYTHON})")

# A tpms run writing test-NAME.step, then the OpenCASCADE check of it.
function(step_test name)
  set(out ${CMAKE_BINARY_DIR}/test-${name}.step)
  add_test(NAME ${name} COMMAND tpms ${ARGN} -o ${out})
  if(occ)
    set(check ${CMAKE_SOURCE_DIR}/tools/check_step.py)
    if(name STREQUAL "unit-cell")
      list(APPEND check --periodic)
    endif()
    set_tests_properties(${name} PROPERTIES FIXTURES_SETUP ${name})
    add_test(NAME occ-${name} COMMAND ${TPMS_PYTHON} ${check} ${out})
    set_tests_properties(occ-${name} PROPERTIES FIXTURES_REQUIRED ${name})
  endif()
endfunction()

foreach(t gyroid schwarz-p diamond iwp frd fischer-koch-s neovius)
  step_test(${t} --type ${t} --cells 1)
endforeach()
# Smooth STEP of the curved shapes: seams welded, curved caps, merged faces.
step_test(tube --shape tube --radius 10 20 --height 10 --cells 1 6 1)
step_test(torus-shell --shape torus-shell --radius 4 10 --major 25 --cells 1 12 4)
# Sphere shell: cone ends, and a dome that ends at the equator plane.
step_test(sphere-shell --shape sphere-shell --radius 15 25 --cells 1 12 5)
step_test(sphere-dome --shape sphere-shell --radius 15 25 --cells 1 12 3 --polar 10 90)
# Rotated lattice (exact edges), and a hollow trimmed sphere.
step_test(rotated --type diamond --rotate 30 20 10 --cells 1)
step_test(sphere-hollow --shape sphere --radius 6 12 --type schwarz-p --form network --density 0.3)
# Wall thickness in mm (the field divided by its gradient), also mapped.
step_test(wall --wall 0.6 --cells 1)
step_test(wall-tube --shape tube --radius 10 20 --height 10 --cells 1 6 1 --wall 0.6)
# Graded density and cell size (a mapped box).
step_test(graded --cell 5 --cells 3 2 2 --cell-grade x 10 --density 0.3 0.15 --grade x)
# Unit cell: exit 1 if opposite faces do not match.
step_test(unit-cell --unit-cell --type diamond --cell 8 10 12 --shift 0.25 0 0.1)
# Trimmed shapes: built around the axis, center line or center (poles).
step_test(cylinder --shape cylinder --radius 10 --height 20)
step_test(sphere --shape sphere --radius 10)
step_test(torus --shape torus --radius 5 --major 12)
# Inner cavities as voids: isolated network blobs in the fluid.
step_test(voids --type neovius --form network --density 0.3 --part void)
# A sheet's fluid: one file per channel, test-fluid_A.step and test-fluid_B.step.
add_test(NAME fluid-split COMMAND tpms --part void -o ${CMAKE_BINARY_DIR}/test-fluid.step)
if(occ)
  set_tests_properties(fluid-split PROPERTIES FIXTURES_SETUP fluid-split)
  add_test(NAME occ-fluid-split COMMAND ${TPMS_PYTHON} ${CMAKE_SOURCE_DIR}/tools/check_step.py
           ${CMAKE_BINARY_DIR}/test-fluid_A.step ${CMAKE_BINARY_DIR}/test-fluid_B.step)
  set_tests_properties(occ-fluid-split PROPERTIES FIXTURES_REQUIRED fluid-split)
endif()
# Plates and buffers: exit 1 if the volume (with the slabs) does not match.
step_test(plates --plates z 1 2)
step_test(buffers --part void --buffers x 5)
# Channel buffers (heat exchanger): the solid with its header plates, and each
# channel's fluid in a file of its own, test-channels_A.step and _B.step.
step_test(channel-plates --channel-buffers x y 5 1)
add_test(NAME channel-fluids COMMAND tpms --part void --channel-buffers z x 4 0.5 -o ${CMAKE_BINARY_DIR}/test-channels.step)
if(occ)
  set_tests_properties(channel-fluids PROPERTIES FIXTURES_SETUP channel-fluids)
  add_test(NAME occ-channel-fluids COMMAND ${TPMS_PYTHON} ${CMAKE_SOURCE_DIR}/tools/check_step.py
           ${CMAKE_BINARY_DIR}/test-channels_A.step ${CMAKE_BINARY_DIR}/test-channels_B.step)
  set_tests_properties(occ-channel-fluids PROPERTIES FIXTURES_REQUIRED channel-fluids)
endif()
add_test(NAME channel-fluids-stl COMMAND tpms --part void --channel-buffers x y 5 1 --res 12 -o ${CMAKE_BINARY_DIR}/test-channels.stl)
# The report with its inputs, as CSV.
add_test(NAME report-csv COMMAND tpms --material 4430 114 --flow 0.1 -o ${CMAKE_BINARY_DIR}/test-report.csv)
add_test(NAME report COMMAND tpms --shape tube --radius 10 20 --height 10 --cells 1 6 1 --report -o ${CMAKE_BINARY_DIR}/test-report.json)
# The surfaces' gradients (used for --wall) against finite differences.
add_executable(check-gradients tests/check_gradients.cpp)
target_include_directories(check-gradients PRIVATE src)
add_test(NAME gradients COMMAND check-gradients)
# Batch run: files per value, one joined CSV.
add_test(NAME sweep COMMAND ${CMAKE_COMMAND} -DTPMS=$<TARGET_FILE:tpms> -DOUT=${CMAKE_BINARY_DIR}
         -P ${CMAKE_SOURCE_DIR}/tests/sweep.cmake)
# Conflicting options must be refused (exit 2), not silently resolved.
add_test(NAME usage-density-and-level COMMAND tpms --density 0.3 --level 0.5 -o ${CMAKE_BINARY_DIR}/test-usage.stl)
set_tests_properties(usage-density-and-level PROPERTIES WILL_FAIL TRUE)
add_test(NAME usage-channel-buffers-network COMMAND tpms --form network --channel-buffers x y 5 -o ${CMAKE_BINARY_DIR}/test-usage.stl)
set_tests_properties(usage-channel-buffers-network PROPERTIES WILL_FAIL TRUE)
add_test(NAME usage-rotate-tube COMMAND tpms --shape tube --radius 10 20 --height 10 --rotate 10 0 0 -o ${CMAKE_BINARY_DIR}/test-usage.stl)
set_tests_properties(usage-rotate-tube PROPERTIES WILL_FAIL TRUE)
add_test(NAME usage-wall-network COMMAND tpms --form network --wall 0.5 -o ${CMAKE_BINARY_DIR}/test-usage.stl)
set_tests_properties(usage-wall-network PROPERTIES WILL_FAIL TRUE)
# Every shape as a closed STL mesh (exit 1 if not closed, 2 if empty).
add_test(NAME shape-tube COMMAND tpms --shape tube --radius 10 20 --height 20 --cells 1 8 2 --res 12 -o ${CMAKE_BINARY_DIR}/test-tube.stl)
add_test(NAME shape-torus-shell COMMAND tpms --shape torus-shell --radius 4 10 --major 25 --cells 1 12 4 --res 10 -o ${CMAKE_BINARY_DIR}/test-torus-shell.stl)
add_test(NAME shape-cylinder COMMAND tpms --shape cylinder --radius 10 --height 20 --res 12 -o ${CMAKE_BINARY_DIR}/test-cylinder.stl)
add_test(NAME shape-sphere COMMAND tpms --shape sphere --radius 10 --form network --res 12 -o ${CMAKE_BINARY_DIR}/test-sphere.stl)
add_test(NAME shape-torus COMMAND tpms --shape torus --radius 6 --major 20 --res 12 -o ${CMAKE_BINARY_DIR}/test-torus.stl)
