# only BSSN (GR) options
option(BSSN_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(BSSN_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(BSSN_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(BSSN_COMPUTE_CONSTRAINTS "Enable constraint equation computations" ON)
option(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(BSSN_ETA_FUNCTION "Use function as ETA damping" OFF)
option(BSSN_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(BSSN_ENABLE_CUDA "enable rhs computation with GPUs" OFF)
option(BSSN_RHS_STAGED_COMP "enable rhs staged computation" OFF)
option(BSSN_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(BSSN_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(BSSN_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" OFF)


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

if(BSSN_RHS_STAGED_COMP)
    add_definitions(-DBSSN_RHS_STAGED_COMP)
endif()

if(BSSN_EXTRACT_BH_LOCATIONS)
    add_definitions(-DBSSN_EXTRACT_BH_LOCATIONS)
endif()


if(BSSN_REFINE_BASE_EH)
    add_definitions(-DBSSN_REFINE_BASE_EH)
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
add_custom_target(GW_EXTRACT
        COMMAND ${CMAKE_SOURCE_DIR}/CodeGen/cog.py ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in > ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in
)



set(BSSN_INC include/parameters.h
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
        )

set(BSSN_SRC src/rkBSSN.cpp
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
    cuda_add_library(bssn_cuda ${CUDA_GR_INC} ${CUDA_GR_SRC})
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
    target_include_directories(bssn_cuda PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(bssn_cuda PRIVATE ${GSL_INCLUDE_DIRS})

    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)

endif()

set(SOURCE_FILES src/gr.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(bssnSolver ${SOURCE_FILES})

target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(bssnSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bssnSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(bssnSolver PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(bssnSolver bssn_cuda)
    target_link_libraries(bssnSolver ${CUDA_LIBS})
endif()

target_link_libraries(bssnSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


set(SOURCE_FILES src/checkDerivatives.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(checkDerivs ${SOURCE_FILES})
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(checkDerivs PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(checkDerivs PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(checkDerivs PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(checkDerivs bssn_cuda)
    target_link_libraries(checkDerivs ${CUDA_LIBS})
endif()



target_link_libraries(checkDerivs dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

set(SOURCE_FILES cuda_gr/src/gpuBSSNExample.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(cudaBSSN ${SOURCE_FILES})
target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(cudaBSSN PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(cudaBSSN PRIVATE ${GSL_INCLUDE_DIRS})

if(WITH_CUDA)
    target_include_directories(cudaBSSN PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(cudaBSSN bssn_cuda)
    target_link_libraries(cudaBSSN ${CUDA_LIBS})
endif()
target_link_libraries(cudaBSSN dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES}  m)




set(SOURCE_FILES cuda_gr/src/gpuTestBed.cpp cuda_gr/include/gpuTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(gpuTestBed ${SOURCE_FILES})
target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(gpuTestBed PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(gpuTestBed PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(gpuTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(gpuTestBed bssn_cuda)
    target_link_libraries(gpuTestBed ${CUDA_LIBS})
endif()
target_link_libraries(gpuTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


set(SOURCE_FILES src/rhsTestBed.cpp include/rhsTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(rhsTestBed ${SOURCE_FILES})
target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(rhsTestBed PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(rhsTestBed PRIVATE ${GSL_INCLUDE_DIRS})

if(WITH_CUDA)
    target_include_directories(rhsTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(rhsTestBed bssn_cuda)
    target_link_libraries(rhsTestBed ${CUDA_LIBS})
endif()

target_link_libraries(rhsTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


set(SOURCE_FILES src/zipUnzipTest.cpp include/zipUnzipTest.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(zipUnzipTest ${SOURCE_FILES})
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/include/cuda)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(zipUnzipTest PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(zipUnzipTest PRIVATE ${GSL_INCLUDE_DIRS})

if(WITH_CUDA)
    target_include_directories(zipUnzipTest PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(zipUnzipTest bssn_cuda)
    target_link_libraries(zipUnzipTest ${CUDA_LIBS})
endif()

target_link_libraries(zipUnzipTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


if(BSSN_EXTRACT_GRAVITATIONAL_WAVES)
    add_definitions(-DBSSN_EXTRACT_GRAVITATIONAL_WAVES)
    add_dependencies(bssnSolver GW_EXTRACT)
    add_dependencies(checkDerivs GW_EXTRACT)
    add_dependencies(cudaBSSN GW_EXTRACT)
    add_dependencies(gpuTestBed GW_EXTRACT)
    add_dependencies(rhsTestBed GW_EXTRACT)
    add_dependencies(zipUnzipTest GW_EXTRACT)
endif(BSSN_EXTRACT_GRAVITATIONAL_WAVES)