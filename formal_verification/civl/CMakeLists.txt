cmake_minimum_required(VERSION 3.10)

project(Armon_CIVL C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP)

set(SOURCES
        blocks.c blocks.h
        solver.c solver.h
        block_interface.c block_interface.h
        time_step.c time_step.h
        debug.c debug.h
        utils.c utils.h
        vars.h
)

add_executable(Armon_CIVL_main ${SOURCES} main.c)
target_link_libraries(Armon_CIVL_main PUBLIC OpenMP::OpenMP_C)
set_property(TARGET Armon_CIVL_main PROPERTY C_STANDARD 11)  # CIVL-C extends C11

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
