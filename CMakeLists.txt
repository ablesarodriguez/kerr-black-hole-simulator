cmake_minimum_required(VERSION 3.16)
project(relativistic_raytracer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# First CPU prototype (no raylib needed)
add_executable(prototype prototype/main.cpp)

# Real-time simulator (raylib + GLSL shader)
add_executable(simulator simulator/main.cpp)
if(WIN32 AND EXISTS ${CMAKE_SOURCE_DIR}/simulator/lib/libraylib.a)
  # raylib 6.0 static library bundled in the repo (MinGW)
  target_include_directories(simulator PRIVATE simulator/include)
  target_link_libraries(simulator PRIVATE ${CMAKE_SOURCE_DIR}/simulator/lib/libraylib.a opengl32 gdi32 winmm)
  # UCRT-based MinGW toolchains lack stat64i32, which that raylib build needs
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles("extern \"C\" int stat64i32(const char*, void*); int main() { return stat64i32(\"\", 0); }" HAVE_STAT64I32)
  if(NOT HAVE_STAT64I32)
    target_sources(simulator PRIVATE simulator/compat_ucrt.cpp)
  endif()
else()
  find_package(raylib QUIET)
  if(NOT raylib_FOUND)
    include(FetchContent)
    FetchContent_Declare(raylib GIT_REPOSITORY https://github.com/raysan5/raylib.git GIT_TAG 5.5)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(raylib)
  endif()
  target_link_libraries(simulator PRIVATE raylib)
endif()
# main.cpp loads the shader from the executable's folder; recopied on every build so shader edits are picked up
add_custom_target(copy_shader ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/simulator/blackhole.fs $<TARGET_FILE_DIR:simulator>/blackhole.fs
  DEPENDS simulator)

# Numerical tests of the Kerr physics (common/KerrPhysics.h)
add_executable(physics_tests tests/physics_tests.cpp)
enable_testing()
add_test(NAME physics_tests COMMAND physics_tests)

# CPU reference renderer to compare the shader against (tools/validate.py)
add_executable(reference_render tools/reference_render.cpp)
