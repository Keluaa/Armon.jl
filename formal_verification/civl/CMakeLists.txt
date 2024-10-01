cmake_minimum_required(VERSION 3.10)

project(Armon_CIVL C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP)
find_package(MPI)

set(NUM_THREADS      "4" CACHE STRING "number of threads to use")
set(NUM_CYCLES       "3" CACHE STRING "number of solver cycles to do")
set(MAX_SWEEPS       "2" CACHE STRING "number of domain sweeps to do per cycle")
set(GRID_SIZE_X      "3" CACHE STRING "number of blocks domain along the X axis")
set(GRID_SIZE_Y      "3" CACHE STRING "number of blocks domain along the Y axis")
set(DO_TIME_STEP     "1" CACHE STRING "enable the time step reduction (one per cycle)")
set(DO_HALO_EXCHANGE "1" CACHE STRING "enable the MPI halo exchange (one per sweep). Also enables MPI reduction if DO_TIME_STEP=1")
set(SIMPLE_XCHG      "0" CACHE STRING "simplifies block exchanges for CIVL")

set(SOURCES
        blocks.c blocks.h
        solver.c solver.h
        block_interface.c block_interface.h
        mpi_comms.c mpi_comms.h
        time_step.c time_step.h
        debug.c debug.h
        utils.c utils.h
        vars.h
)

add_executable(Armon_CIVL_main ${SOURCES} main.c)
target_link_libraries(Armon_CIVL_main PUBLIC OpenMP::OpenMP_C MPI::MPI_C)
set_property(TARGET Armon_CIVL_main PROPERTY C_STANDARD 11)  # CIVL-C extends C11
target_compile_options(Armon_CIVL_main PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(Armon_CIVL_main PUBLIC
        NUM_THREADS=${NUM_THREADS}
        NUM_CYCLES=${NUM_CYCLES}
        MAX_SWEEPS=${MAX_SWEEPS}
        GRID_SIZE_X=${GRID_SIZE_X}
        GRID_SIZE_Y=${GRID_SIZE_Y}
        DO_TIME_STEP=${DO_TIME_STEP}
        DO_HALO_EXCHANGE=${DO_HALO_EXCHANGE}
        SIMPLE_XCHG=${SIMPLE_XCHG}
)

add_custom_target(Armon_CIVL_run
        COMMAND civl run ${SOURCES} main.cvl
        SOURCES ${SOURCES} main.cvl
)

add_custom_target(Armon_CIVL_verify
        COMMAND civl verify ${SOURCES} main.cvl
        SOURCES ${SOURCES} main.cvl
)

add_custom_target(Armon_NO_CIVL
        COMMAND civl run -_CIVL ${SOURCES} main.cvl
        SOURCES ${SOURCES} main.cvl
)
