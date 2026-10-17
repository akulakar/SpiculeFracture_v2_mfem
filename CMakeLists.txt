# ---------------------------------------------------------------------------------------------------
# SpiculeFracture_fem: MPI-parallel regularized variational fracture solver built on MFEM + mfemplus.
#
# Layout:
#   include/spicule/           header-only classes: PhaseFieldSolver, BoundaryConditions, OutputWriter
#   drivers/                   one executable per .cpp file, each linked against `spicule`
#
# Build:   cmake -S . -B build && cmake --build build
# Run:     cd build && mpirun -np 4 ./SpiculeFracture ../input/SpiculeFracture.json
# ---------------------------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.12...4.0)
project(spicule_fem LANGUAGES CXX)

# MFEM requires C++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to an optimized build; pass -DCMAKE_BUILD_TYPE=Debug for debugging.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose Debug or Release" FORCE)
endif()

# ---------------------------------------------------------------------------------------------------
# Dependencies (override paths with -DMFEM_DIR=..., -DMFEMPLUS_DIR=... etc.)
# ---------------------------------------------------------------------------------------------------
find_package(MPI REQUIRED)
# MFEM must be built with MPI, HYPRE and METIS (and ADIOS2 for the BP5 output).
find_package(MFEM REQUIRED)
# JSON input parsing. QUIET because the include path is also added explicitly below.
find_package(nlohmann_json QUIET HINTS /opt/homebrew/opt/nlohmann-json)

# ADIOS2 is needed for the default "adios2" output format (adios2stream in MFEM).
set(SPICULE_USE_ADIOS2 ON CACHE BOOL "Write results with ADIOS2")
if(SPICULE_USE_ADIOS2)
  find_package(ADIOS2 REQUIRED HINTS /opt/homebrew/opt/adios2/lib/cmake)
endif()

# Root of the mfemplus repository (provides the damaged-elasticity / strain-energy integrators).
# It is added as a subdirectory, which builds the static library target mfemplus::fracture.
set(MFEMPLUS_DIR "/Users/Aashutosh/Academic/Software/MFEM_v2/mfemplus" CACHE PATH "mfemplus repository root")
add_subdirectory(${MFEMPLUS_DIR} mfemplus_build)
# METIS is used by MFEM's mesh partitioner and must be linked explicitly.
set(METIS_INCLUDE_DIR /opt/homebrew/opt/metis/include CACHE PATH "")
set(METIS_LIBRARY /opt/homebrew/opt/metis/lib/libmetis.dylib CACHE FILEPATH "")

# ---------------------------------------------------------------------------------------------------
# Drivers: one executable per file in drivers/, named after the file.
# The solver, boundary conditions and output classes are header-only (include/spicule/PhaseFieldFracture.hpp).
# ---------------------------------------------------------------------------------------------------
file(GLOB DRIVERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/drivers/*.cpp")
foreach(drv ${DRIVERS})
  get_filename_component(name ${drv} NAME_WE)
  add_executable(${name} ${drv})

  target_include_directories(${name} PRIVATE
    ${CMAKE_SOURCE_DIR}/include   # spicule/PhaseFieldFracture.hpp
    ${METIS_INCLUDE_DIR}
    /opt/homebrew/opt/nlohmann-json/include)  # <nlohmann/json.hpp>

  # mfemplus::fracture carries its own include root (<mfemplus/SolidMechanics/...>) and MFEM/MPI deps.
  target_link_libraries(${name} PRIVATE mfemplus::fracture mfem ${METIS_LIBRARY} ${MPI_CXX_LIBRARIES})
  if(SPICULE_USE_ADIOS2)
    target_link_libraries(${name} PRIVATE adios2::adios2)
  endif()

  message(STATUS "Configured driver: ${name}")
endforeach()
