//
// Created by milinda on 1/16/19.
//

#include "dataUtils.h"

namespace bssn
{

    void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,double tolerance, const Point* ptIn, unsigned int numPt,Point* ptOut)
    {
        if((pMesh->isActive()))
        {

            unsigned int rankActive=pMesh->getMPIRank();
            unsigned int npesActive=pMesh->getMPICommSize();

            MPI_Comm commActive=pMesh->getMPICommunicator();

            const ot::TreeNode* allElements=&(*(pMesh->getAllElements().begin()));
            const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));
            const unsigned int * e2n_dg=&(*(pMesh->getE2NMapping_DG().begin()));
            const unsigned int * cgToDg=&(*(pMesh->getCG2DGMap().begin()));

            unsigned int lookup=0;
            unsigned int ownerID, ii_x, jj_y, kk_z;

            ot::TreeNode tmpOct;
            const unsigned int eleOrder=pMesh->getElementOrder();
            double hx,x,y,z;


            std::vector<Point> ptList;
            for(unsigned int node=pMesh->getNodeLocalBegin();node<pMesh->getNodeLocalEnd();node++)
            {

                if(var[node]<tolerance)
                {
                    lookup=cgToDg[node];
                    pMesh->dg2eijk(lookup,ownerID,ii_x,jj_y,kk_z);
                    tmpOct=allElements[ownerID];
                    hx=(tmpOct.maxX()-tmpOct.minX())/((double) eleOrder);
                    x=tmpOct.minX() + ii_x*(hx);
                    y=tmpOct.minY() + jj_y*(hx);
                    z=tmpOct.minZ() + kk_z*(hx);

                    ptList.push_back(Point(GRIDX_TO_X(x),GRIDY_TO_Y(y),GRIDZ_TO_Z(z)));
                    
                    //std::cout<<" x: "<<GRIDX_TO_X(x)<<" y: "<<GRIDY_TO_Y(y)<<" z: "<<GRIDZ_TO_Z(z)<<std::endl;
                    

                }

            }

            
            
            

            std::vector<Point> * ptCluster=new std::vector<Point>[numPt];

            double min;
            unsigned int cID;
            for(unsigned int pt=0;pt<ptList.size();pt++)
            {
                cID=0;
                min=(ptIn[0]-ptList[pt]).abs();
                for(unsigned int c=1;c<numPt;c++)
                {
                    if(min>=(ptIn[c]-ptList[pt]).abs())
                    {
                        min=(ptIn[c]-ptList[pt]).abs();
                        cID=c;
                    }
                }

                ptCluster[cID].push_back(ptList[pt]);
            }


            ptList.clear();
            Point * ptMean=new Point[numPt];
            DendroIntL* ptCounts=new DendroIntL[numPt];
            DendroIntL* ptCounts_g=new DendroIntL[numPt];

            for(unsigned int c=0;c<numPt;c++)
            {
                ptMean[c]=Point(0,0,0);
                ptOut[c]=Point(0,0,0);

            }

            for(unsigned int c=0;c<numPt;c++)
            {
                ptCounts[c]=ptCluster[c].size();
                for(unsigned int pt=0;pt<ptCluster[c].size();pt++)
                    ptMean[c]+=ptCluster[c][pt];
            }

            
            
            par::Mpi_Allreduce(ptCounts,ptCounts_g,numPt,MPI_SUM,commActive);
            par::Mpi_Allreduce(ptMean,ptOut,numPt,par::Mpi_datatype<Point>::_SUM(),commActive);

            
            for(unsigned int c=0;c<numPt;c++)
               ptOut[c]/=(double)ptCounts_g[c];
            
            
            /*if(pMesh->getMPIRank()==0)std::cout<<"bh1 in : "<<ptIn[0].x()<<", "<<ptIn[0].y()<<", "<<ptIn[0].z()<<std::endl;
            if(pMesh->getMPIRank()==0)std::cout<<"bh2 in : "<<ptIn[1].x()<<", "<<ptIn[1].y()<<", "<<ptIn[1].z()<<std::endl;
            
            if(pMesh->getMPIRank()==0)std::cout<<"bh1 out: "<<ptOut[0].x()<<", "<<ptOut[0].y()<<", "<<ptOut[0].z()<<std::endl;
            if(pMesh->getMPIRank()==0)std::cout<<"bh2 out: "<<ptOut[1].x()<<", "<<ptOut[1].y()<<", "<<ptOut[1].z()<<std::endl;*/
            
            


            delete [] ptCounts;
            delete [] ptCounts_g;
            delete [] ptMean;
            delete [] ptCluster;

        }

    }

    void writeBHCoordinates(const ot::Mesh* pMesh,const Point* ptLocs, unsigned int numPt,unsigned int timestep)
    {
        
        unsigned int rankGlobal=pMesh->getMPIRankGlobal();
        if(!rankGlobal)
        {

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_BHLocations.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str());
            fileGW.open (fName,std::ofstream::app);

            // writes the header
            if(timestep==0)
                fileGW<<"TimeStep\t"<<" bh1_x\t"<<" bh1_y\t"<<" bh1_z\t"<<" bh2_x\t"<<" bh2_y\t"<<" bh2_z\t"<<std::endl;

            fileGW<<timestep<<"\t"<<ptLocs[0].x()<<"\t"<<ptLocs[0].y()<<"\t"<<ptLocs[0].z()<<"\t"<<ptLocs[1].x()<<"\t"<<ptLocs[1].y()<<"\t"<<ptLocs[1].z()<<std::endl;
            fileGW.close();
            return;

        }
    }

}// end of namespace bssn
