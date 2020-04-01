option(BSSN_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(BSSN_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(BSSN_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(BSSN_COMPUTE_CONSTRAINTS "Enable constraint equation computations" ON)
option(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(BSSN_ETA_FUNCTION "Use function as ETA damping" ON)
option(BSSN_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(BSSN_ENABLE_CUDA "enable rhs computation with GPUs" OFF)


if(BSSN_PROFILE_HUMAN_READABLE)
    add_definitions(-DBSSN_PROFILE_HUMAN_READABLE)
endif()

if(BSSN_CONSEC_COMM_SELECT)
    add_definitions(-DBSSN_CONSEC_COMM_SELECT)
endif()

if(BSSN_ENABLE_VTU_OUTPUT)
    add_definitions(-DBSSN_ENABLE_VTU_OUTPUT)
endif()

if(BSSN_COMPUTE_CONSTRAINTS)
    add_definitions(-DBSSN_COMPUTE_CONSTRAINTS)
endif()

if(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
    add_definitions(-DBSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(BSSN_ETA_FUNCTION)
    add_definitions(-DUSE_ETA_FUNC)
endif()

if(BSSN_KERR_SCHILD_TEST)
    add_definitions(-DBSSN_KERR_SCHILD_TEST)
endif()

if(BSSN_ENABLE_CUDA)
    add_definitions(-DBSSN_ENABLE_CUDA)
endif()

# some debug flags. uncomment to print the debug information.
#add_definitions(-DDEBUG_RK_SOLVER)
#add_definitions(-DDEBUG_UPWIND_INTERP)
#add_definitions(-DDEBUG_DOWNWIND_INTERP)
#add_definitions(-DDEBUG_DERIVS_COMP)
#add_definitions(-DMATVEC_PROFILE)
#add_definitions(-DFEM_ACCUMILATE_ONES_TEST)
#add_definitions(-DNLSM_DIRICHLET_BDY)


set(CCZ4_INC include/parameters.h
        include/rk45CCZ4.h
        include/rk4CCZ4.h
        include/rk3CCZ4.h
        include/bh.h
        include/grUtils.h
        include/grUtils.tcc
        include/rhs.h
        include/derivs.h
        include/physcon.h
        include/profile_params.h
        include/ccz4_constraints.h
        include/TPUtilities.h
        include/TwoPunctures.h
        )

set(CCZ4_SRC src/rk45CCZ4.cpp
        src/rk4CCZ4.cpp
        src/rk3CCZ4.cpp
        src/parameters.cpp
        src/grUtils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/physcon.cpp
        src/profile_params.cpp
        src/ccz4_constraints.cpp
        src/TPCoordTransf.cpp
        src/TPEquations.cpp
        src/TPFuncAndJacobian.cpp
        src/TPNewton.cpp
        src/TPUtilities.cpp
        src/TwoPunctures.cpp
        )

set(SOURCE_FILES src/gr.cpp  include/gr.h ${CCZ4_INC} ${CCZ4_SRC})
add_executable(ccz4Solver ${SOURCE_FILES})
install(TARGETS ccz4Solver RUNTIME DESTINATION bin)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(ccz4Solver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(ccz4Solver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(ccz4Solver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(bssnSolver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(ccz4Solver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)
