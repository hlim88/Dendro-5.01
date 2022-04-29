# only QUADGRAV (GR) options
option(QUADGRAV_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(QUADGRAV_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(QUADGRAV_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(QUADGRAV_COMPUTE_CONSTRAINTS "Enable constraint equation computations" ON)
option(QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(QUADGRAV_ETA_FUNCTION "Use function as ETA damping" OFF)
option(QUADGRAV_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(QUADGRAV_ENABLE_CUDA "enable rhs computation with GPUs" OFF)
#option(QUADGRAV_RHS_STAGED_COMP "enable rhs staged computation" OFF)
option(QUADGRAV_GAUGE_ROCHESTER "enable rochester gauge for QUADGRAV equations" OFF)
option(QUADGRAV_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(QUADGRAV_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" ON)
option(QUADGRAV_USE_4TH_ORDER_DERIVS "use 4th order derivs" OFF)
option(QUADGRAV_USE_6TH_ORDER_DERIVS "use 4th order derivs" ON)

#09/22/20 Milinda.
# uncomment below lines to use with adv derivs, adv dervs for 6th order derivs are not implemented. (Not recommended to enable)
#add_definitions(-DQUADGRAV_USE_ADVECTIVE_DERIVS)

if(QUADGRAV_USE_4TH_ORDER_DERIVS)
    set(QUADGRAV_USE_6TH_ORDER_DERIVS OFF)
elseif(QUADGRAV_USE_6TH_ORDER_DERIVS)
    set(QUADGRAV_USE_4TH_ORDER_DERIVS OFF)
endif()

if(QUADGRAV_USE_4TH_ORDER_DERIVS)
    add_definitions(-DQUADGRAV_USE_4TH_ORDER_DERIVS)
endif()

if(QUADGRAV_USE_6TH_ORDER_DERIVS)
    add_definitions(-DQUADGRAV_USE_6TH_ORDER_DERIVS)
endif()

if(QUADGRAV_PROFILE_HUMAN_READABLE)
    add_definitions(-DQUADGRAV_PROFILE_HUMAN_READABLE)
endif()

if(QUADGRAV_CONSEC_COMM_SELECT)
    add_definitions(-DQUADGRAV_CONSEC_COMM_SELECT)
endif()

if(QUADGRAV_ENABLE_VTU_OUTPUT)
    add_definitions(-DQUADGRAV_ENABLE_VTU_OUTPUT)
endif()

if(QUADGRAV_COMPUTE_CONSTRAINTS)
    add_definitions(-DQUADGRAV_COMPUTE_CONSTRAINTS)
endif()

if(QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
    add_definitions(-DQUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(QUADGRAV_ETA_FUNCTION)
    add_definitions(-DUSE_ETA_FUNC)
endif()

if(QUADGRAV_KERR_SCHILD_TEST)
    add_definitions(-DQUADGRAV_KERR_SCHILD_TEST)
endif()

if(QUADGRAV_ENABLE_CUDA)
    add_definitions(-DQUADGRAV_ENABLE_CUDA)
endif()

if(QUADGRAV_RHS_STAGED_COMP)
    add_definitions(-DQUADGRAV_RHS_STAGED_COMP)
endif()

if(QUADGRAV_EXTRACT_BH_LOCATIONS)
    add_definitions(-DQUADGRAV_EXTRACT_BH_LOCATIONS)
endif()


if(QUADGRAV_REFINE_BASE_EH)
    add_definitions(-DQUADGRAV_REFINE_BASE_EH)
endif()

if(QUADGRAV_GAUGE_ROCHESTER)
    add_definitions(-DUSE_ROCHESTER_GAUGE)
endif()


#message(${KWAY})

# some debug flags. uncomment to print the debug information.
#add_definitions(-DDEBUG_RK_SOLVER)
#add_definitions(-DDEBUG_UPWIND_INTERP)
#add_definitions(-DDEBUG_DOWNWIND_INTERP)
#add_definitions(-DDEBUG_DERIVS_COMP)
#add_definitions(-DMATVEC_PROFILE)
#add_definitions(-DFEM_ACCUMILATE_ONES_TEST)
#add_definitions(-DNLSM_DIRICHLET_BDY)

# generates the code to extract GW waves,

#add_custom_command(
#        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
#        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in
#        COMMAND python3 -W ignore ${CMAKE_SOURCE_DIR}/CodeGen/cog.py ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in > ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
#        COMMENT "Generating gwExtract.h file "
#        VERBATIM USES_TERMINAL
#)

#add_custom_target(GW_EXTRACT
#        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
#)


#add_custom_command(
#    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_standard_gauge.cpp
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py 
#    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py unstaged standard const . > ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_standard_gauge.cpp
#    COMMENT "generating QUADGRAV code unstaged standard gauge with constant eta damping (during intial compilation this might take a while. )"
#    VERBATIM USES_TERMINAL
#)

#add_custom_command(
#    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_standard_gauge.cpp
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py 
#    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py unstaged standard func . > ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_standard_gauge.cpp
#    COMMENT "generating QUADGRAV code unstaged standard gauge with varying (large mass ratios) eta damping (during intial compilation this might take a while. )"
#    VERBATIM USES_TERMINAL
#)


#add_custom_command(
#    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_rochester_gauge.cpp
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py 
#    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py unstaged rochester const . > ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_rochester_gauge.cpp
#    COMMENT "generating QUADGRAV code unstaged rochester gauge with constant eta damping (during intial compilation this might take a while. )"
#    VERBATIM USES_TERMINAL
#)

#add_custom_command(
#    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_rochester_gauge.cpp
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
#    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py 
#    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/QUADGRAV.py unstaged rochester func . > ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_rochester_gauge.cpp
#    COMMENT "generating QUADGRAV code unstaged rochester gauge with constant eta damping (during intial compilation this might take a while. )"
#    VERBATIM USES_TERMINAL
#)


#add_custom_target(QUADGRAV_CODEGEN
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_standard_gauge.cpp
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_standard_gauge.cpp
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_const_rochester_gauge.cpp
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/QUADGRAVeqs_eta_func_rochester_gauge.cpp
#    )


set(QUADGRAV_INC include/parameters.h
        include/bh.h
        include/grUtils.h
        include/grUtils.tcc
        include/rhs.h
        include/derivs.h
        include/physcon.h
        include/profile_params.h
        include/bssn_constraints.h
        include/TPUtilities.h
        include/TwoPunctures.h
        include/dataUtils.h
        include/aeh_solver.h
        include/bssnCtx.h
        )

set(QUADGRAV_SRC src/rkBSSN.cpp
           src/parameters.cpp
           src/grUtils.cpp
           src/rhs.cpp
           src/derivs.cpp
           src/physcon.cpp
           src/profile_params.cpp
           src/bssn_constraints.cpp
           src/TPCoordTransf.cpp
           src/TPEquations.cpp
           src/TPFuncAndJacobian.cpp
           src/TPNewton.cpp
           src/TPUtilities.cpp
           src/TwoPunctures.cpp
           src/dataUtils.cpp
           src/hadrhs.F90
           src/bssnCtx.cpp
           
        )

if(WITH_CUDA)

    set(CUDA_GR_INCLUDE_H
            cuda_gr/include/bssn_rhs_deriv_mem_cuda.h
            cuda_gr/include/params_cu.h
            cuda_gr/include/profile_gpu.h

            )



    set(CUDA_GR_SRC_CPP
            cuda_gr/src/bssn_rhs_deriv_mem_cuda.cpp
            cuda_gr/src/params_cu.cpp
            cuda_gr/src/profile_gpu.cpp

        )
endif()

set(CUDA_GR_INC
        cuda_gr/include/cudaUtils.cuh
        cuda_gr/include/derivs.cuh
        cuda_gr/include/rhs_bssn.cuh
        cuda_gr/include/rhs_cuda.cuh
        cuda_gr/include/gpuTest.cuh
        )

set(CUDA_GR_SRC
        cuda_gr/src/cudaUtils.cu
        cuda_gr/src/derivs.cu
        cuda_gr/src/rhs_bssn.cu
        cuda_gr/src/rhs_cuda.cu
        cuda_gr/src/gpuTest.cu
        )



if(WITH_CUDA)
    cuda_add_library(QUADGRAV_cuda ${CUDA_GR_INC} ${CUDA_GR_SRC})
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
    target_include_directories(QUADGRAV_cuda PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(QUADGRAV_cuda PRIVATE ${GSL_INCLUDE_DIRS})

    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(QUADGRAV_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)

endif()

#add_library(had_rhs STATIC src/hadrhs.F90)

set(SOURCE_FILES src/gr.cpp include/gr.h ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(quadgravSolver ${SOURCE_FILES})

target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(quadgravSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(quadgravSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(quadgravSolver PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(quadgravSolver QUADGRAV_cuda)
    target_link_libraries(quadgravSolver ${CUDA_LIBS})
endif()

target_link_libraries(quadgravSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)



#set(SOURCE_FILES src/gr_nuts.cpp include/gr.h ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#add_executable(quadgravSolverNUTS ${SOURCE_FILES})

#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#target_include_directories(quadgravSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#target_include_directories(quadgravSolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
#target_include_directories(quadgravSolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
#if(WITH_CUDA)
#    target_include_directories(quadgravSolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(quadgravSolverNUTS QUADGRAV_cuda)
#    target_link_libraries(quadgravSolverNUTS ${CUDA_LIBS})
#endif()

#target_link_libraries(quadgravSolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


#set(SOURCE_FILES src/checkDerivatives.cpp ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#add_executable(checkDerivs ${SOURCE_FILES})
#target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#target_include_directories(checkDerivs PRIVATE ${MPI_INCLUDE_PATH})
#target_include_directories(checkDerivs PRIVATE ${GSL_INCLUDE_DIRS})
#if(WITH_CUDA)
#    target_include_directories(checkDerivs PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(checkDerivs QUADGRAV_cuda)
#    target_link_libraries(checkDerivs ${CUDA_LIBS})
#endif()



#target_link_libraries(checkDerivs dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

#if(WITH_CUDA)
#    set(SOURCE_FILES cuda_gr/src/gpuQUADGRAVExample.cpp ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#    add_executable(cudaQUADGRAV ${SOURCE_FILES})
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#    target_include_directories(cudaQUADGRAV PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#    target_include_directories(cudaQUADGRAV PRIVATE ${MPI_INCLUDE_PATH})
#    target_include_directories(cudaQUADGRAV PRIVATE ${GSL_INCLUDE_DIRS})
#    target_include_directories(cudaQUADGRAV PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(cudaQUADGRAV QUADGRAV_cuda)
#    target_link_libraries(cudaQUADGRAV ${CUDA_LIBS})
#    target_link_libraries(cudaQUADGRAV dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES}  m)
#    add_dependencies(cudaQUADGRAV QUADGRAV_CODEGEN)
#    #add_dependencies(cudaQUADGRAV GW_EXTRACT)
#endif()



#if(WITH_CUDA)
#set(SOURCE_FILES cuda_gr/src/gpuTestBed.cpp cuda_gr/include/gpuTestBed.h ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#add_executable(gpuTestBed ${SOURCE_FILES})
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#    target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#    target_include_directories(gpuTestBed PRIVATE ${MPI_INCLUDE_PATH})
#    target_include_directories(gpuTestBed PRIVATE ${GSL_INCLUDE_DIRS})
#    target_include_directories(gpuTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(gpuTestBed QUADGRAV_cuda)
#    target_link_libraries(gpuTestBed ${CUDA_LIBS})
#    target_link_libraries(gpuTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)
#    add_dependencies(gpuTestBed QUADGRAV_CODEGEN)
#    #add_dependencies(gpuTestBed GW_EXTRACT)
#endif()


#set(SOURCE_FILES src/rhsTestBed.cpp include/rhsTestBed.h ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#add_executable(rhsTestBed ${SOURCE_FILES})
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#target_include_directories(rhsTestBed PRIVATE ${MPI_INCLUDE_PATH})
#target_include_directories(rhsTestBed PRIVATE ${GSL_INCLUDE_DIRS})

#if(WITH_CUDA)
#    target_include_directories(rhsTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(rhsTestBed QUADGRAV_cuda)
#    target_link_libraries(rhsTestBed ${CUDA_LIBS})
#endif()

#target_link_libraries(rhsTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


#set(SOURCE_FILES src/zipUnzipTest.cpp include/zipUnzipTest.h ${QUADGRAV_INC} ${QUADGRAV_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#add_executable(zipUnzipTest ${SOURCE_FILES})
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/test/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
#target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
#target_include_directories(zipUnzipTest PRIVATE ${MPI_INCLUDE_PATH})
#target_include_directories(zipUnzipTest PRIVATE ${GSL_INCLUDE_DIRS})

#if(WITH_CUDA)
#    target_include_directories(zipUnzipTest PRIVATE ${CUDA_INCLUDE_DIRS})
#    target_link_libraries(zipUnzipTest QUADGRAV_cuda)
#    target_link_libraries(zipUnzipTest ${CUDA_LIBS})
#endif()

#target_link_libraries(zipUnzipTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


#if(QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES)
#    add_definitions(-DQUADGRAV_EXTRACT_GRAVITATIONAL_WAVES)
    #add_dependencies(quadgravSolver GW_EXTRACT)
    #add_dependencies(checkDerivs GW_EXTRACT)
    #add_dependencies(rhsTestBed GW_EXTRACT)
    #add_dependencies(zipUnzipTest GW_EXTRACT)
#endif(QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES)


#add_dependencies(quadgravSolver QUADGRAV_CODEGEN)
#add_dependencies(rhsTestBed QUADGRAV_CODEGEN)


#set(SOURCE_FILES    src/TPID.cpp 
#                    ${QUADGRAV_INC}
#                    ${QUADGRAV_SRC}
#                    )

#add_executable(tpid ${SOURCE_FILES})
#target_link_libraries(tpid dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)
#target_include_directories(tpid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(tpid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#target_include_directories(tpid PRIVATE ${MPI_INCLUDE_PATH})
#target_include_directories(tpid PRIVATE ${GSL_INCLUDE_DIRS})
