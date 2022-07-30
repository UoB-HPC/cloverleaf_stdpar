cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

project(cloverleaf_stdpar)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FetchContent)

if (USE_TBB)
    FetchContent_Declare(
            TBB
            GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
            GIT_TAG faaf43c4ab22cb4b4267d65d5e218fa58800eea8
    )
    # Don't fail builds on waring (TBB has -Wall while not being free of warnings from unused symbols...)
    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    set(TBB_STRICT OFF)
    # Not using FetchContent_MakeAvailable because we need EXCLUDE_FROM_ALL
    FetchContent_GetProperties(TBB)
    if (NOT TBB_POPULATED)
        FetchContent_Populate(TBB)
        add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif ()
    list(APPEND EXTRA_LIBS TBB::tbb)
    list(APPEND EXTRA_HEADERS ${TBB_SOURCE_DIR})
endif ()


if (USE_ONEDPL)
    FetchContent_Declare(
            oneDPL
            GIT_REPOSITORY https://github.com/oneapi-src/oneDPL.git
            GIT_TAG oneDPL-2021.7.0-release
    )
    string(TOLOWER ${USE_ONEDPL} ONEDPL_BACKEND)
    # XXX oneDPL looks for omp instead of openmp, which mismatches(!) with ONEDPL_PAR_BACKEND if using find_package
    if (ONEDPL_BACKEND STREQUAL "openmp")
        set(ONEDPL_BACKEND omp)
    endif ()
    # Not using FetchContent_MakeAvailable (CMake>= 3.14) because we need EXCLUDE_FROM_ALL
    FetchContent_GetProperties(oneDPL)
    if (NOT oneDPL_POPULATED)
        FetchContent_Populate(oneDPL)
        if (USE_TBB)
            macro(find_package NAME)
                if ("${NAME}" STREQUAL "TBB")
                    message(STATUS "Discarding oneDPL's call to find_package(${NAME} ${ARGN})")
                else ()
                    _find_package(${NAME} ${ARGN})
                endif ()
            endmacro()
        endif ()
        add_subdirectory(${onedpl_SOURCE_DIR} ${onedpl_BINARY_DIR} EXCLUDE_FROM_ALL)

        # Fixup oneDPL's omission on setting DPCPP definitions.
        # We do this after the creation of the oneDPL target.
        if (ONEDPL_BACKEND MATCHES "^(dpcpp|dpcpp_only)$")
            target_compile_definitions(oneDPL INTERFACE ONEDPL_USE_DPCPP_BACKEND=1)
        endif ()

    endif ()
    list(APPEND EXTRA_LIBS oneDPL)
endif ()

# Enable offloading support (via the non-standard `-stdpar`) for the new NVHPC SDK.
# The values are Nvidia architectures in ccXY format will be passed in via `-gpu=` (e.g `cc70`)
# Possible values are:
#   cc35  - Compile for compute capability 3.5
#   cc50  - Compile for compute capability 5.0
#   cc60  - Compile for compute capability 6.0
#   cc62  - Compile for compute capability 6.2
#   cc70  - Compile for compute capability 7.0
#   cc72  - Compile for compute capability 7.2
#   cc75  - Compile for compute capability 7.5
#   cc80  - Compile for compute capability 8.0
#   ccall - Compile for all supported compute capabilities"
if (NVHPC_OFFLOAD)
    set(OFFLOAD_FLAGS
            -stdpar
            -gpu=${NVHPC_OFFLOAD}
            --restrict
            -Mllvm-fast
            -Ktrap=none
            -Minfo=accel
            -Minfo=stdpar)
endif ()


if (MPI_AS_LIBRARY)

    if (NOT DEFINED MPI_C_LIB_DIR)
        message(FATAL_ERROR "MPI_C_LIB_DIR must be specified, typically <mpi_root_dir>/lib")
    endif ()

    if (NOT DEFINED MPI_C_INCLUDE_DIR)
        message(FATAL_ERROR "MPI_C_INCLUDE_DIR must be specified, typically <mpi_root_dir>/include")
    endif ()

    if (NOT DEFINED MPI_C_LIB)
        message(FATAL_ERROR "MPI_C_LIB must be specified, for example: mpich for libmpich.so in MPI_C_LIB_DIR")
    endif ()

    message(STATUS "Using MPI as a library (${MPI_C_LIB})")
    message(STATUS "MPI include dir: ${MPI_C_INCLUDE_DIR}")
    message(STATUS "MPI library dir: ${MPI_C_LIB_DIR}")

    include_directories(${MPI_C_INCLUDE_DIR})
    link_directories(${MPI_C_LIB_DIR})
else ()
    find_package(MPI REQUIRED)
    set(MPI_C_LIB MPI::MPI_CXX)
endif ()

set(SOURCES
        src/accelerate.cpp
        src/advec_cell.cpp
        src/advec_mom.cpp
        src/advection.cpp
        src/build_field.cpp
        src/calc_dt.cpp
        src/clover_leaf.cpp
        src/comms.cpp
        src/field_summary.cpp
        src/flux_calc.cpp
        src/generate_chunk.cpp
        src/hydro.cpp
        src/ideal_gas.cpp
        src/initialise_chunk.cpp
        src/initialise.cpp
        src/pack_kernel.cpp
        src/PdV.cpp
        src/read_input.cpp
        src/report.cpp
        src/reset_field.cpp
        src/revert.cpp
        src/start.cpp
        src/timer.cpp
        src/timestep.cpp
        src/update_halo.cpp
        src/update_tile_halo.cpp
        src/update_tile_halo_kernel.cpp
        src/viscosity.cpp
        src/visit.cpp)


add_executable(clover_leaf ${SOURCES})

target_include_directories(clover_leaf PUBLIC BEFORE src)

if (USE_VECTOR)
    target_compile_definitions(clover_leaf PUBLIC USE_VECTOR)
endif ()
if (SERIAL_COPY_CTOR)
    target_compile_definitions(clover_leaf PUBLIC SERIAL_COPY_CTOR)
endif ()
if (USE_ONEDPL)
    target_compile_definitions(clover_leaf PUBLIC USE_ONEDPL)
endif ()

target_compile_options(clover_leaf
        PUBLIC
        -Wall
        -Wextra
        -Wcast-align
        #        -Wfatal-errors
        -Werror=return-type
        -Wno-unused-parameter
        -Wno-unused-variable
        #        -Wno-ignored-attributes
        ${EXTRA_FLAGS}
        )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    # advec kernel has unused values, past attempts at removing them caused a regression (although that commit also changes how the input deck is parsed)
    set_source_files_properties(src/advec_cell.cpp PROPERTIES COMPILE_FLAGS -Wno-unused-but-set-variable)
endif()

if (USE_ONEDPL)
    # oneDPL's transform_reduce has sign comparisons, see https://github.com/oneapi-src/oneDPL/issues/472
    set_source_files_properties(src/calc_dt.cpp src/field_summary.cpp PROPERTIES COMPILE_FLAGS -Wno-sign-compare)
endif ()

target_include_directories(clover_leaf PUBLIC BEFORE ${EXTRA_HEADERS})

set(DEBUG_OPTIONS -O2 ${OFFLOAD_FLAGS} ${CXX_EXTRA_FLAGS})
set(RELEASE_OPTIONS -O3 ${OFFLOAD_FLAGS} ${CXX_EXTRA_FLAGS})

target_link_libraries(clover_leaf PUBLIC ${MPI_C_LIB} ${EXTRA_LIBS})
#target_link_libraries(clover_leaf PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C)

target_compile_options(clover_leaf PUBLIC "$<$<CONFIG:RelWithDebInfo>:${RELEASE_OPTIONS}>")
target_compile_options(clover_leaf PUBLIC "$<$<CONFIG:Release>:${RELEASE_OPTIONS}>")
target_compile_options(clover_leaf PUBLIC "$<$<CONFIG:Debug>:${DEBUG_OPTIONS}>")


target_link_options(clover_leaf PUBLIC ${OFFLOAD_FLAGS})
target_link_options(clover_leaf PUBLIC LINKER:${CXX_EXTRA_LINKER_FLAGS})

