cmake_minimum_required(VERSION 3.29)
project(Grasp2_Exact)

set(CMAKE_CXX_STANDARD 20)

if(DEFINED ENV{GUROBI_HOME})
    set(GUROBI_HOME "$ENV{GUROBI_HOME}")
    message(STATUS "GUROBI_HOME found: ${GUROBI_HOME}")
else()
    message(FATAL_ERROR "GUROBI_HOME not defined")
endif()

include_directories(${GUROBI_HOME}/include)

link_directories(${GUROBI_HOME}/lib)


if(WIN32)
    include_directories(${GUROBI_HOME}/include)
    link_directories(${GUROBI_HOME}/lib)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(GUROBI_CXX_LIBRARY gurobi_c++mdd2017)
        #message(STATUS "Using libraries Gurobi C++: gurobi_c++mdd2017")
    else()
        set(GUROBI_CXX_LIBRARY gurobi_c++md2017)
        #message(STATUS "Using librerías Gurobi C++: gurobi_c++md2017")
    endif()

    set(GUROBI_C_LIBRARY gurobi120)

else()
    # In Linux/macOS
    set(GUROBI_CXX_LIBRARY gurobi_c++)
    set(GUROBI_C_LIBRARY gurobi120)
endif()


add_executable(Grasp2_Exact main.cpp)

target_link_libraries(Grasp2_Exact PRIVATE ${GUROBI_CXX_LIBRARY} ${GUROBI_C_LIBRARY})

if(UNIX AND NOT APPLE)
    set_target_properties(Grasp2_Exact PROPERTIES
            INSTALL_RPATH "${GUROBI_HOME}/lib"
            BUILD_WITH_INSTALL_RPATH TRUE)
endif()

