//
// Created by milinda on 05/02/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for CCZ4 simulation.
*/

#include "mesh.h"
#include <math.h>
#include "daUtils.h"
//#include "lebedev.h"
//#include "swsh.h"

namespace ccz4
{

    template<typename T>
    double extractGravitationlWaves(const ot::Mesh* mesh,const T** constraintVar, double r,unsigned int timestep)
    {


        double integral_real=0.0;
        double integral_imag=0.0;

        double integral_real_g=0.0;
        double integral_imag_g=0.0;
        
	//HL : Commented out for ccz4. We need a version for GW radiation about this
	//Also, need to check the integration part
	#if 0
        unsigned int rankGlobal=mesh->getMPIRankGlobal();
        unsigned int npesGlobal=mesh->getMPICommSizeGlobal();
        MPI_Comm commGlobal=mesh->getMPIGlobalCommunicator();

        if(mesh->isActive())
        {

            const unsigned int rankActive=mesh->getMPIRank();
            const unsigned int npesActive=mesh->getMPICommSize();

            const unsigned int numPts=LEBEDEV_025_NUM_PTS;

            std::vector<double> coords;
            coords.resize(3*numPts);


            for(unsigned int pts=0;pts<numPts;pts++)
            {
                coords[3*pts + 0]=X_TO_GRIDX(r*sin(LEBEDEV_025_THETA[pts])*cos(LEBEDEV_025_PHI[pts]));
                coords[3*pts + 1]=Y_TO_GRIDY(r*sin(LEBEDEV_025_THETA[pts])*sin(LEBEDEV_025_PHI[pts]));
                coords[3*pts + 2]=Z_TO_GRIDZ(r*cos(LEBEDEV_025_THETA[pts]));
            }

            std::vector<double> psi4_real;
            psi4_real.resize(numPts);

            std::vector<double> psi4_imag;
            psi4_imag.resize(numPts);


            std::vector<unsigned int > validIndex;
            ot::da::interpolateToCoords(mesh,constraintVar[VAR_CONSTRAINT::C_PSI4_REAL],&(*(coords.begin())),coords.size(),&(*(psi4_real.begin())),validIndex);
            ot::da::interpolateToCoords(mesh,constraintVar[VAR_CONSTRAINT::C_PSI4_IMG],&(*(coords.begin())),coords.size(),&(*(psi4_imag.begin())),validIndex);

            for(unsigned int index=0;index<validIndex.size();index++)
            {
                integral_real+=(psi4_real[validIndex[index]]*swsh::m2Y2_0_REAL[validIndex[index]]-psi4_imag[validIndex[index]]*swsh::m2Y2_0_IMAG[validIndex[index]])*sin(LEBEDEV_025_THETA[validIndex[index]])*LEBEDEV_025_WEIGHT[validIndex[index]];
                integral_imag+=(psi4_real[validIndex[index]]*swsh::m2Y2_0_IMAG[validIndex[index]]+psi4_imag[validIndex[index]]*swsh::m2Y2_0_REAL[validIndex[index]])*sin(LEBEDEV_025_THETA[validIndex[index]])*LEBEDEV_025_WEIGHT[validIndex[index]];
            }



        }


        par::Mpi_Reduce(&integral_real,&integral_real_g,1,MPI_SUM,0,commGlobal);
        par::Mpi_Reduce(&integral_imag,&integral_imag_g,1,MPI_SUM,0,commGlobal);


        if(!rankGlobal)
        {

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_GW_%d.dat",ccz4::CCZ4_PROFILE_FILE_PREFIX.c_str(),(int)r);
            fileGW.open (fName,std::ofstream::app);

            // writes the header
            if(timestep==0)
                fileGW<<"TimeStep\t"<<" r\t"<<" rxpsi_real_lm\t"<<" rxpsi_img_lm"<<std::endl;

            fileGW<<timestep<<"\t"<<r<<"\t"<<(r*integral_real_g)<<"\t"<<(r*integral_imag_g)<<"\t"<<std::endl;
            //std::cout<<"real: "<<integral_real_g<<" imag: "<<integral_imag_g<<std::endl;
            fileGW.close();
            return 0;

        }


	#endif


    }



    template <typename T>
    double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm)
    {

        double l2=0.0;
        double l2_g=0.0;
        const double MASK_THRESHOLD=1.0;
        for(unsigned int i=lbegin;i<lend;i++)
        {
            if(maskVec[i]<MASK_THRESHOLD)
                l2+=(constraintVec[i]*constraintVec[i]*maskVec[i]*maskVec[i]*maskVec[i]*maskVec[i]);
            else
                l2+=(constraintVec[i]*constraintVec[i]);
        }


        par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,comm);

        return (sqrt(l2_g));


    }

    
    template <typename T>
    double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec)
    {
        double l2_g=0.0;
        if(mesh->isActive())
        {
            
            unsigned int nodeLookUp_CG;
            unsigned int nodeLookUp_DG;
            unsigned int len,x,y,z;
            MPI_Comm comm=mesh->getMPICommunicator();
            const unsigned int eleLocalBegin=mesh->getElementLocalBegin();
            const unsigned int eleLocalEnd=mesh->getElementLocalEnd();
            
            const ot::TreeNode* pNodes=&(*(mesh->getAllElements().begin()));
            const unsigned int eleOrder=mesh->getElementOrder();
            unsigned int ownerID,ii_x,jj_y,kk_z;
            const unsigned int * e2n_cg=&(*(mesh->getE2NMapping().begin()));
            const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));
            const unsigned int nPe=mesh->getNumNodesPerElement();
            
            const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
            
            double origin[3];
            origin[0]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
            origin[1]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
            origin[2]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
            const double R0 = sqrt(
                 (ccz4::BH1.getBHCoordX()-ccz4::BH2.getBHCoordX())*(ccz4::BH1.getBHCoordX()-ccz4::BH2.getBHCoordX())
                +(ccz4::BH1.getBHCoordY()-ccz4::BH2.getBHCoordY())*(ccz4::BH1.getBHCoordY()-ccz4::BH2.getBHCoordY())
                +(ccz4::BH1.getBHCoordZ()-ccz4::BH2.getBHCoordZ())*(ccz4::BH1.getBHCoordZ()-ccz4::BH2.getBHCoordZ())
            )+3;
            //std::cout<<"R0: "<<R0<<std::endl;
            //ccz4::CCZ4_WAVELET_TOL_FUNCTION_R0;
            double l2=0.0;
            
            for(unsigned int elem=eleLocalBegin;elem<eleLocalEnd;elem++)
            {
                for(unsigned int k=0;k<(eleOrder+1);k++)
                  for(unsigned int j=0;j<(eleOrder+1);j++)
                    for(unsigned int i=0;i<(eleOrder+1);i++)
                    {
                     
                        nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            mesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            if(ownerID==elem)
                            {
                                len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                                x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                                y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                                z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));    
                                
                                double r = sqrt(
                                              GRIDX_TO_X(x)*GRIDX_TO_X(x)
                                            + GRIDY_TO_Y(y)*GRIDY_TO_Y(y)
                                            + GRIDZ_TO_Z(z)*GRIDZ_TO_Z(z)
                                            );
                                
                                if(r>R0)
                                {
                                    l2+=(constraintVec[nodeLookUp_CG]*constraintVec[nodeLookUp_CG]);
                                }
                                    
                                
                            }
                            
                        }
                        
                    }  
            
            }
            
            
            par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,mesh->getMPICommunicator());
            
        }
        
        return sqrt(l2_g);
        
        
    }

    template<typename T>
    double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,unsigned int timestep)
    {
        const unsigned int numConstraints=4;
        double constraintMaskedL2[numConstraints]; // remove the psi4

        unsigned int rankGlobal=mesh->getMPIRankGlobal();
        unsigned int npesGlobal=mesh->getMPICommSizeGlobal();
        MPI_Comm commGlobal=mesh->getMPIGlobalCommunicator();

        if(mesh->isActive())
        {

            unsigned int rankActive=mesh->getMPIRank();
            unsigned int npesActive=mesh->getMPICommSize();
            MPI_Comm commActive=mesh->getMPICommunicator();


            for(unsigned int index=0;index<numConstraints;index++)
            {
                constraintMaskedL2[index]=ccz4::computeConstraintL2Norm(mesh,constraintVar[index]);
                if(!rankActive)
                {
                    std::cout<<YLW<<"\tConstraint " <<ccz4::CCZ4_CONSTRAINT_VAR_NAMES[index]<< " L2 : ("<<constraintMaskedL2[index]<<" )"<<NRM<<std::endl;
                }

            }


            if(!rankActive)
            {
                std::ofstream fileGW;
                char fName[256];
                sprintf(fName,"%s_Constraints.dat",ccz4::CCZ4_PROFILE_FILE_PREFIX.c_str());
                fileGW.open (fName,std::ofstream::app);
                // writes the header
                if(timestep==0)
                    fileGW<<"TimeStep\t"<<" C_HAM\t"<<" C_MOM0\t"<<" C_MOM1\t"<<" C_MOM2\t"<<std::endl;

                fileGW<<timestep<<"\t"<<constraintMaskedL2[0]<<"\t"<<constraintMaskedL2[1]<<"\t"<<constraintMaskedL2[2]<<"\t"<<constraintMaskedL2[3]<<std::endl;
                fileGW.close();
                return 0;

            }



        }

    }

} // end of namespace ccz4

