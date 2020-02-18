//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for BSSN simulation.
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
#include "swsh.h"
#include "lebedev.h"
#include "grDef.h"


using json = nlohmann::json;
namespace bssn
{


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
 * @note This function is taken from the old single core bssn version.
 **/

 // Puncture BH
 void punctureData(const double xx1,const double yy1,const double zz1, double *var);
 // Static Kerr-Schild
 void KerrSchildData(const double xx1,const double yy1,const double zz1, double *var);
 void noiseData(const double xx1,const double yy1,const double zz1, double *var);

 void fake_initial_data(double x, double y, double z, double *u);

 /**
  * @brief: Generates block adaptive octree for the given binary blockhole problem.
  * @param[out] tmpNodes: created octree tmpNodes
  * @param[in] pt_min: block min point
  * @param[in] pt_max: block max point
  * @param[in] regLev: regular grid level
  * @param[in] maxDepth: maximum refinement level. 
  * @param[in] comm: MPI communicator. 
  * */

  void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);

  /**
   * @brief wavelet tolerance as a function of space.
   * */
  double computeWTol(double x,double y,double z,double tol_min);

   /**
    * @breif: Compute L2 constraint norms. 
    */
   template <typename T>
   double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm);

   /**
    * @breif: Compute L2 constraint norms. 
    */
   template <typename T>
   double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,const T* maskVector,T maskthreshoold);

   /**
    * @breif write constraints to a file. 
    */
   template<typename T>
   double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,const T* maskVec, double maskthreshoold,unsigned int timestep);

   void writeBLockToBinary(const double **unzipVarsRHS,unsigned int offset,const double *pmin, const double *pmax,double* bxMin,double * bxMax, const unsigned int *sz,unsigned int blkSz,double dxFactor,const char* fprefix);


}// end of namespace bssn



namespace bssn
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
