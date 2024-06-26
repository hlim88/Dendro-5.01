//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for CCZ4 simulation.
*/


#ifndef SFCSORTBENCH_GRUTILS_H
#define SFCSORTBENCH_GRUTILS_H

#include "point.h"
#include "parameters.h"
#include "mesh.h"
#include "block.h"
#include "parUtils.h"
#include "json.hpp"
#include "dendroProfileParams.h"
#include "profile_params.h"

#define Rx (ccz4::CCZ4_COMPD_MAX[0]-ccz4::CCZ4_COMPD_MIN[0])
#define Ry (ccz4::CCZ4_COMPD_MAX[1]-ccz4::CCZ4_COMPD_MIN[1])
#define Rz (ccz4::CCZ4_COMPD_MAX[2]-ccz4::CCZ4_COMPD_MIN[2])

#define RgX (ccz4::CCZ4_OCTREE_MAX[0]-ccz4::CCZ4_OCTREE_MIN[0])
#define RgY (ccz4::CCZ4_OCTREE_MAX[1]-ccz4::CCZ4_OCTREE_MIN[1])
#define RgZ (ccz4::CCZ4_OCTREE_MAX[2]-ccz4::CCZ4_OCTREE_MIN[2])

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-ccz4::CCZ4_OCTREE_MIN[0]))+ccz4::CCZ4_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-ccz4::CCZ4_OCTREE_MIN[1]))+ccz4::CCZ4_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-ccz4::CCZ4_OCTREE_MIN[2]))+ccz4::CCZ4_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-ccz4::CCZ4_COMPD_MIN[0]))+ccz4::CCZ4_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-ccz4::CCZ4_COMPD_MIN[1]))+ccz4::CCZ4_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-ccz4::CCZ4_COMPD_MIN[2]))+ccz4::CCZ4_OCTREE_MIN[2])


using json = nlohmann::json;
namespace ccz4
{
/**
 * @brief These variable indexes are based on the variables defined in rkCCZ4.h
 * */

enum VAR {U_ALPHA=0,U_PSI,U_K,U_GH0,U_GH1,U_GH2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5,U_THETA_Z4};

enum VAR_CONSTRAINT {C_HAM=0, C_MOM0, C_MOM1, C_MOM2, C_PSI4_REAL, C_PSI4_IMG};

static const char * CCZ4_VAR_NAMES[]={"U_ALPHA","U_PSI","U_K","U_GH0","U_GH1","U_GH2","U_BETA0","U_BETA1","U_BETA2","U_B0","U_B1","U_B2","U_SYMGT0","U_SYMGT1","U_SYMGT2","U_SYMGT3","U_SYMGT4","U_SYMGT5","U_SYMAT0","U_SYMAT1","U_SYMAT2","U_SYMAT3","U_SYMAT4","U_SYMAT5","U_THETA_Z4"};

static const char * CCZ4_CONSTRAINT_VAR_NAMES[]={"C_HAM","C_MOM0","C_MOM1","C_MOM2","C_PSI4_REAL","C_PSI4_IMG"};

/**
 * @brief internal variables needed for rk update.
 * */


 /**
  * @brief: Read the parameter file and initialize the variables in parameters.h file.
  * @param[in] fName: file name
  * @param[in] comm: MPI communicator.
  * */
  void readParamFile(const char * fName,MPI_Comm comm);


/**
 * @brief Initialize all the variables for a given point in space.
 * @param [in] coord: coordinates of the point.
 * @param [out] var: pointer to the list of variables, computed. var size should be (VAR::U_SYMAT5+1)
 * @note This function is taken from the old single core ccz4 version.
 **/
 void punctureData(const double xx1,const double yy1,const double zz1, double *var);
 // Static Kerr-Schild
 void KerrSchildData(const double xx1,const double yy1,const double zz1, double *var);

 /**
  * @brief: Generates block adaptive octree for the given binary blockhole problem.
  *
  * */

  void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);


   /**
   * @brief wavelet tolerance as a function of space.
   * */
  double computeWTol(double x,double y,double z,double tol_min);

  template<typename T>
    double extractGravitationlWaves(const ot::Mesh* mesh,const T** constraintVar, double r,unsigned int timestep);

    template <typename T>
    double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm);

    template <typename T>
    double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec);

    template<typename T>
    double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,unsigned int timestep);




}// end of namespace ccz4



namespace ccz4
{

    namespace timer
    {

        /**@brief initialize all the flop counters. */
        void initFlops();

        /**@brief clears the snapshot counter for time profiler variables*/
        void resetSnapshot();


       /**@brief reduce min mean max.
        * @param [in] stat: local time
        * @param [out] stat_g 0-min, 1-mean 2-max
       * */
       template<typename T>
       void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm)
       {
           int rank,npes;
           MPI_Comm_size(comm,&npes);
           MPI_Comm_rank(comm,&rank);

           par::Mpi_Reduce(stat,stat_g,1,MPI_MIN,0,comm);
           par::Mpi_Reduce(stat,stat_g+1,1,MPI_SUM,0,comm);
           par::Mpi_Reduce(stat,stat_g+2,1,MPI_MAX,0,comm);
           stat_g[1]/=(npes);

       }


        /** @breif : printout the profile parameters. */
        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh);

        /** @breif : printout the profile parameters (intermediate profile information). */
        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep);


    }


}

#include "grUtils.tcc"



#endif //SFCSORTBENCH_GRUTILS_H
