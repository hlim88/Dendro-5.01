//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for CCZ4 simulation.
*/
//

#include "grUtils.h"

namespace ccz4
{

    void readParamFile(const char * fName,MPI_Comm comm)
    {


        json parFile;
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int vtu_len;
        unsigned int chp_len;
        unsigned int prf_len;

        if(!rank)
        {
            std::ifstream infile(fName);
            if(!infile) {std::cout<<fName<<" parameter file open failed "<<std::endl;}
            infile>>parFile;

            ccz4::CCZ4_IO_OUTPUT_FREQ=parFile["CCZ4_IO_OUTPUT_FREQ"];
            ccz4::CCZ4_REMESH_TEST_FREQ=parFile["CCZ4_REMESH_TEST_FREQ"];
            ccz4::CCZ4_CHECKPT_FREQ=parFile["CCZ4_CHECKPT_FREQ"];
            ccz4::CCZ4_IO_OUTPUT_GAP=parFile["CCZ4_IO_OUTPUT_GAP"];
            ccz4::CCZ4_VTU_FILE_PREFIX=parFile["CCZ4_VTU_FILE_PREFIX"].get<std::string>();
            ccz4::CCZ4_CHKPT_FILE_PREFIX=parFile["CCZ4_CHKPT_FILE_PREFIX"].get<std::string>();
            ccz4::CCZ4_PROFILE_FILE_PREFIX=parFile["CCZ4_PROFILE_FILE_PREFIX"].get<std::string>();
            ccz4::CCZ4_RESTORE_SOLVER=parFile["CCZ4_RESTORE_SOLVER"];
            ccz4::CCZ4_ENABLE_BLOCK_ADAPTIVITY=parFile["CCZ4_ENABLE_BLOCK_ADAPTIVITY"];
            ccz4::CCZ4_ID_TYPE=parFile["CCZ4_ID_TYPE"];
            ccz4::CCZ4_BLK_MIN_X=parFile["CCZ4_BLK_MIN_X"];
            ccz4::CCZ4_BLK_MIN_Y=parFile["CCZ4_BLK_MIN_Y"];
            ccz4::CCZ4_BLK_MIN_Z=parFile["CCZ4_BLK_MIN_Z"];
            ccz4::CCZ4_BLK_MAX_X=parFile["CCZ4_BLK_MAX_X"];
            ccz4::CCZ4_BLK_MAX_Y=parFile["CCZ4_BLK_MAX_Y"];
            ccz4::CCZ4_BLK_MAX_Z=parFile["CCZ4_BLK_MAX_Z"];
            ccz4::CCZ4_DENDRO_GRAIN_SZ=parFile["CCZ4_DENDRO_GRAIN_SZ"];
            ccz4::CCZ4_ASYNC_COMM_K=parFile["CCZ4_ASYNC_COMM_K"];
            ccz4::CCZ4_DENDRO_AMR_FAC=parFile["CCZ4_DENDRO_AMR_FAC"];
            ccz4::CCZ4_WAVELET_TOL=parFile["CCZ4_WAVELET_TOL"];
            ccz4::CCZ4_LOAD_IMB_TOL=parFile["CCZ4_LOAD_IMB_TOL"];
            ccz4::CCZ4_RK45_TIME_BEGIN=parFile["CCZ4_RK45_TIME_BEGIN"];
            ccz4::CCZ4_RK45_TIME_END=parFile["CCZ4_RK45_TIME_END"];
            ccz4::CCZ4_RK45_TIME_STEP_SIZE=parFile["CCZ4_RK45_TIME_STEP_SIZE"];
            ccz4::CCZ4_RK45_DESIRED_TOL=parFile["CCZ4_RK45_DESIRED_TOL"];
            ccz4::CCZ4_DIM=parFile["CCZ4_DIM"];
            ccz4::CCZ4_MAXDEPTH=parFile["CCZ4_MAXDEPTH"];
            ccz4::BH1=BH((double)parFile["CCZ4_BH1"]["MASS"],(double)parFile["CCZ4_BH1"]["X"],(double)parFile["CCZ4_BH1"]["Y"],(double)parFile["CCZ4_BH1"]["Z"],(double)parFile["CCZ4_BH1"]["V_X"],(double)parFile["CCZ4_BH1"]["V_Y"],(double)parFile["CCZ4_BH1"]["V_Z"],(double)parFile["CCZ4_BH1"]["SPIN"],(double)parFile["CCZ4_BH1"]["SPIN_THETA"],(double)parFile["CCZ4_BH1"]["SPIN_PHI"]);
            ccz4::BH2=BH((double)parFile["CCZ4_BH2"]["MASS"],(double)parFile["CCZ4_BH2"]["X"],(double)parFile["CCZ4_BH2"]["Y"],(double)parFile["CCZ4_BH2"]["Z"],(double)parFile["CCZ4_BH2"]["V_X"],(double)parFile["CCZ4_BH2"]["V_Y"],(double)parFile["CCZ4_BH2"]["V_Z"],(double)parFile["CCZ4_BH2"]["SPIN"],(double)parFile["CCZ4_BH2"]["SPIN_THETA"],(double)parFile["CCZ4_BH2"]["SPIN_PHI"]);
            ccz4::CCZ4_GRID_MIN_X=parFile["CCZ4_GRID_MIN_X"];
            ccz4::CCZ4_GRID_MAX_X=parFile["CCZ4_GRID_MAX_X"];
            ccz4::CCZ4_GRID_MIN_Y=parFile["CCZ4_GRID_MIN_Y"];
            ccz4::CCZ4_GRID_MAX_Y=parFile["CCZ4_GRID_MAX_Y"];
            ccz4::CCZ4_GRID_MIN_Z=parFile["CCZ4_GRID_MIN_Z"];
            ccz4::CCZ4_GRID_MAX_Z=parFile["CCZ4_GRID_MAX_Z"];
            ccz4::ETA_CONST=parFile["ETA_CONST"];
            ccz4::ETA_R0=parFile["ETA_R0"];
            ccz4::ETA_DAMPING=parFile["ETA_DAMPING"];
            ccz4::ETA_DAMPING_EXP=parFile["ETA_DAMPING_EXP"];
            ccz4::CCZ4_LAMBDA[0]=(unsigned int) parFile["CCZ4_LAMBDA"]["CCZ4_LAMBDA_1"];
            ccz4::CCZ4_LAMBDA[1]=(unsigned int) parFile["CCZ4_LAMBDA"]["CCZ4_LAMBDA_2"];
            ccz4::CCZ4_LAMBDA[2]=(unsigned int) parFile["CCZ4_LAMBDA"]["CCZ4_LAMBDA_3"];
            ccz4::CCZ4_LAMBDA[3]=(unsigned int) parFile["CCZ4_LAMBDA"]["CCZ4_LAMBDA_4"];
            ccz4::CCZ4_LAMBDA_F[0]=parFile["CCZ4_LAMBDA_F"]["CCZ4_LAMBDA_F0"];
            ccz4::CCZ4_LAMBDA_F[1]=parFile["CCZ4_LAMBDA_F"]["CCZ4_LAMBDA_F1"];
            ccz4::CCZ4_KAPPA[0]=parFile["CCZ4_KAPPA"]["CCZ4_KAPPA_0"];
            ccz4::CCZ4_KAPPA[1]=parFile["CCZ4_KAPPA"]["CCZ4_KAPPA_1"];
            ccz4::PSI_FLOOR=parFile["PSI_FLOOR"];
            ccz4::CCZ4_TRK0=parFile["CCZ4_TRK0"];
            ccz4::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];
            ccz4::P_EXPO=parFile["P_EXPO"];
            ccz4::ANG_PAR=parFile["ANG_PAR"];

	    ccz4::CCZ4_NUM_REFINE_VARS=parFile["CCZ4_NUM_REFINE_VARS"];
            for(unsigned int i=0;i<ccz4::CCZ4_NUM_REFINE_VARS;i++)
                ccz4::CCZ4_REFINE_VARIABLE_INDICES[i]=parFile["CCZ4_REFINE_VARIABLE_INDICES"][i];

            ccz4::CCZ4_NUM_EVOL_VARS_VTU_OUTPUT=parFile["CCZ4_NUM_EVOL_VARS_VTU_OUTPUT"];
            ccz4::CCZ4_NUM_CONST_VARS_VTU_OUTPUT=parFile["CCZ4_NUM_CONST_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<ccz4::CCZ4_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                ccz4::CCZ4_VTU_OUTPUT_EVOL_INDICES[i]=parFile["CCZ4_VTU_OUTPUT_EVOL_INDICES"][i];

            for(unsigned int i=0;i<ccz4::CCZ4_NUM_CONST_VARS_VTU_OUTPUT;i++)
                ccz4::CCZ4_VTU_OUTPUT_CONST_INDICES[i]=parFile["CCZ4_VTU_OUTPUT_CONST_INDICES"][i];
	    /* Parameter for TwoPuncure */
    	    TPID::target_M_plus=parFile["TPID_TARGET_M_PLUS"];
            TPID::target_M_minus=parFile["TPID_TARGET_M_MINUS"];
            TPID::par_m_plus=TPID::target_M_plus;
            TPID::par_m_minus=TPID::target_M_minus;
            TPID::par_b=parFile["TPID_PAR_B"];

            TPID::par_P_plus[0]=parFile["TPID_PAR_P_PLUS"]["X"];
            TPID::par_P_plus[1]=parFile["TPID_PAR_P_PLUS"]["Y"];
            TPID::par_P_plus[2]=parFile["TPID_PAR_P_PLUS"]["Z"];

            TPID::par_P_minus[0]=parFile["TPID_PAR_P_MINUS"]["X"];
            TPID::par_P_minus[1]=parFile["TPID_PAR_P_MINUS"]["Y"];
            TPID::par_P_minus[2]=parFile["TPID_PAR_P_MINUS"]["Z"];

            TPID::par_S_plus[0]=parFile["TPID_PAR_S_PLUS"]["X"];
            TPID::par_S_plus[1]=parFile["TPID_PAR_S_PLUS"]["Y"];
            TPID::par_S_plus[2]=parFile["TPID_PAR_S_PLUS"]["Z"];

            TPID::par_S_minus[0]=parFile["TPID_PAR_S_MINUS"]["X"];
            TPID::par_S_minus[1]=parFile["TPID_PAR_S_MINUS"]["Y"];
            TPID::par_S_minus[2]=parFile["TPID_PAR_S_MINUS"]["Z"];

            TPID::center_offset[0]=parFile["TPID_CENTER_OFFSET"]["X"];
            TPID::center_offset[1]=parFile["TPID_CENTER_OFFSET"]["Y"];
            TPID::center_offset[2]=parFile["TPID_CENTER_OFFSET"]["Z"];

            TPID::initial_lapse_psi_exponent=parFile["TPID_INITIAL_LAPSE_PSI_EXPONENT"];
            TPID::npoints_A=parFile["TPID_NPOINTS_A"];
            TPID::npoints_B=parFile["TPID_NPOINTS_B"];
            TPID::npoints_phi=parFile["TPID_NPOINTS_PHI"];

            TPID::give_bare_mass=parFile["TPID_GIVE_BARE_MASS"];
            TPID::initial_lapse=parFile["INITIAL_LAPSE"];
            TPID::solve_momentum_constraint=parFile["TPID_SOLVE_MOMENTUM_CONSTRAINT"];
            TPID::grid_setup_method=parFile["TPID_GRID_SETUP_METHOD"];
            TPID::verbose=parFile["TPID_VERBOSE"];
            TPID::adm_tol=parFile["TPID_ADM_TOL"];
            TPID::Newton_tol=parFile["TPID_NEWTON_TOL"];

            vtu_len=CCZ4_VTU_FILE_PREFIX.size();
            chp_len=CCZ4_CHKPT_FILE_PREFIX.size();
            prf_len=CCZ4_PROFILE_FILE_PREFIX.size();

        }


        par::Mpi_Bcast(&CCZ4_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&CCZ4_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&CCZ4_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&CCZ4_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&CCZ4_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&CCZ4_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&CCZ4_ASYNC_COMM_K,1,0,comm);

        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=CCZ4_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=CCZ4_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=CCZ4_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        CCZ4_VTU_FILE_PREFIX=std::string(vtu_name);
        CCZ4_CHKPT_FILE_PREFIX=std::string(chp_name);
        CCZ4_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&CCZ4_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&CCZ4_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&CCZ4_WAVELET_TOL,1,0,comm);
        par::Mpi_Bcast(&CCZ4_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&CCZ4_RK45_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&CCZ4_RK45_TIME_END,1,0,comm);
        par::Mpi_Bcast(&CCZ4_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&CCZ4_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&CCZ4_DIM,1,0,comm);
        par::Mpi_Bcast(&CCZ4_MAXDEPTH,1,0,comm);

        MPI_Bcast(&(ccz4::BH1),sizeof(double)*10,MPI_BYTE,0,comm);
        MPI_Bcast(&(ccz4::BH2),sizeof(double)*10,MPI_BYTE,0,comm);

	par::Mpi_Bcast(&CCZ4_ID_TYPE,1,0,comm);

        par::Mpi_Bcast(&CCZ4_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&CCZ4_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&CCZ4_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&CCZ4_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&CCZ4_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&CCZ4_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&CCZ4_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&CCZ4_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&CCZ4_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&CCZ4_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&CCZ4_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&CCZ4_BLK_MAX_Z,1,0,comm);

        CCZ4_OCTREE_MAX[0]=(double )(1u<<ccz4::CCZ4_MAXDEPTH);
        CCZ4_OCTREE_MAX[1]=(double )(1u<<ccz4::CCZ4_MAXDEPTH);
        CCZ4_OCTREE_MAX[2]=(double )(1u<<ccz4::CCZ4_MAXDEPTH);

        CCZ4_COMPD_MIN[0]=CCZ4_GRID_MIN_X;
        CCZ4_COMPD_MIN[1]=CCZ4_GRID_MIN_Y;
        CCZ4_COMPD_MIN[2]=CCZ4_GRID_MIN_Z;

        CCZ4_COMPD_MAX[0]=CCZ4_GRID_MAX_X;
        CCZ4_COMPD_MAX[1]=CCZ4_GRID_MAX_Y;
        CCZ4_COMPD_MAX[2]=CCZ4_GRID_MAX_Z;


        par::Mpi_Bcast(&ETA_CONST,1,0,comm);
        par::Mpi_Bcast(&ETA_R0,1,0,comm);
        par::Mpi_Bcast(&ETA_DAMPING,1,0,comm);
        par::Mpi_Bcast(&ETA_DAMPING_EXP,1,0,comm);

        par::Mpi_Bcast(&PSI_FLOOR,1,0,comm);
        par::Mpi_Bcast(&CCZ4_TRK0,1,0,comm);
        par::Mpi_Bcast(&KO_DISS_SIGMA, 1, 0, comm);
        par::Mpi_Bcast(&P_EXPO, 1, 0, comm);
        par::Mpi_Bcast(&ANG_PAR, 1, 0, comm);

        MPI_Bcast(&(ccz4::CCZ4_LAMBDA),4,MPI_UNSIGNED,0,comm);
        MPI_Bcast(&(ccz4::CCZ4_LAMBDA_F),2,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(ccz4::CCZ4_KAPPA),2,MPI_DOUBLE,0,comm);

	par::Mpi_Bcast(&CCZ4_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&CCZ4_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);
        par::Mpi_Bcast(&CCZ4_NUM_CONST_VARS_VTU_OUTPUT,1,0,comm);

  	if(CCZ4_NUM_REFINE_VARS>CCZ4_NUM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of CCZ4_NUM_VARS"<<std::endl; exit(0);}
        if(CCZ4_NUM_EVOL_VARS_VTU_OUTPUT>CCZ4_NUM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of CCZ4_NUM_VARS"<<std::endl; exit(0);}
        if(CCZ4_NUM_CONST_VARS_VTU_OUTPUT>CCZ4_CONSTRAINT_NUM_VARS){std::cout<<"Error[parameter file]: Number of constraint VTU variables should be less than number of CCZ4_CONSTRAINT_NUM_VARS"<<std::endl; exit(0);}  

        par::Mpi_Bcast(CCZ4_REFINE_VARIABLE_INDICES,CCZ4_NUM_VARS,0,comm);
        par::Mpi_Bcast(CCZ4_VTU_OUTPUT_EVOL_INDICES,CCZ4_NUM_VARS,0,comm);
        par::Mpi_Bcast(CCZ4_VTU_OUTPUT_CONST_INDICES,CCZ4_CONSTRAINT_NUM_VARS,0,comm);

	par::Mpi_Bcast(CCZ4_REFINE_VARIABLE_INDICES,CCZ4_NUM_VARS,0,comm);
        par::Mpi_Bcast(CCZ4_VTU_OUTPUT_EVOL_INDICES,CCZ4_NUM_VARS,0,comm);
        par::Mpi_Bcast(CCZ4_VTU_OUTPUT_CONST_INDICES,CCZ4_CONSTRAINT_NUM_VARS,0,comm);

        par::Mpi_Bcast(&TPID::target_M_plus,1,0,comm);
        par::Mpi_Bcast(&TPID::target_M_minus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_m_plus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_m_minus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_b,1,0,comm);

        MPI_Bcast(&(TPID::par_P_plus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_P_minus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_S_plus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_S_minus),3,MPI_DOUBLE,0,comm);

        MPI_Bcast(&(TPID::center_offset),3,MPI_DOUBLE,0,comm);

        par::Mpi_Bcast(&TPID::initial_lapse_psi_exponent,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_A,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_B,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_phi,1,0,comm);

        par::Mpi_Bcast(&TPID::give_bare_mass,1,0,comm);
        par::Mpi_Bcast(&TPID::initial_lapse,1,0,comm);
        par::Mpi_Bcast(&TPID::solve_momentum_constraint,1,0,comm);
        par::Mpi_Bcast(&TPID::grid_setup_method,1,0,comm);
        par::Mpi_Bcast(&TPID::verbose,1,0,comm);
        par::Mpi_Bcast(&TPID::adm_tol,1,0,comm);
        par::Mpi_Bcast(&TPID::Newton_tol,1,0,comm);


    }



    void punctureData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double xx=GRIDX_TO_X(xx1);
        const double yy=GRIDY_TO_Y(yy1);
        const double zz=GRIDZ_TO_Z(zz1);

        /* Define the Levi-Cevita pseudo-tensor and Kroneckar delta */
        double epijk[3][3][3];
        int i,j,k;
        for (k=0;k<3;k++) {
            for (j=0;j<3;j++) {
                for (i=0;i<3;i++) {
                    epijk[k][j][i] = 0.0;
                }
            }
        }
        epijk[0][1][2] = 1.0;epijk[1][2][0] = 1.0;epijk[2][0][1] = 1.0;
        epijk[0][2][1] = -1.0;epijk[2][1][0] = -1.0;epijk[1][0][2] = -1.0;

        double deltaij[3][3];
        for (j=0;j<3;j++) {
            for (i=0;i<3;i++) {
                deltaij[j][i] = 0.0;
            }
        }

        deltaij[0][0] = 1.0;deltaij[1][1] = 1.0;deltaij[2][2] = 1.0;

        double x1,y1,z1,rv1;
        double x2,y2,z2,rv2;
        double vn1[3],vn2[3];

        double vpsibl;
        double v_u_corr,amp_capj,amp_capr,l_r,u0_j,u2_j,mu_j,p2_mu_j,v_u_j1;
        double v1,v2,v3,v4,vt1,vt2;

        int i1,i2,i3,i4;
        double amp_capp,u0_p,u2_p,mu_p,p2_mu_p;
        double v_u_p1,v_u_c1,v_u_j2,v_u_p2;
        double v_u_c2,vpsibl_u,vpsibl_u2;


        // bh 1
        double mass1 = BH1.getBHMass();
        double bh1x = BH1.getBHCoordX();
        double bh1y = BH1.getBHCoordY();
        double bh1z = BH1.getBHCoordZ();

        double vp1[3];
        vp1[0] = BH1.getVx();
        vp1[1] = BH1.getVy();
        vp1[2] = BH1.getVz();

        double vp1tot = sqrt( vp1[0]*vp1[0] + vp1[1]*vp1[1] + vp1[2]*vp1[2] );
        double spin1 = BH1.getBHSpin();
        double spin1_th = BH1.getBHSpinTheta();
        double spin1_phi = BH1.getBHSpinPhi();
        double vs1[3];

        vs1[0] = spin1*sin(spin1_th)*cos(spin1_phi);
        vs1[1] = spin1*sin(spin1_th)*sin(spin1_phi);
        vs1[2] = spin1*cos(spin1_th);

        // bh 2
        double mass2 = BH2.getBHMass();
        double bh2x =  BH2.getBHCoordX();
        double bh2y =  BH2.getBHCoordY();
        double bh2z =  BH2.getBHCoordZ();

        double vp2[3];
        vp2[0] = BH2.getVx();
        vp2[1] = BH2.getVy();
        vp2[2] = BH2.getVz();

        double vp2tot = sqrt( vp2[0]*vp2[0] + vp2[1]*vp2[1] + vp2[2]*vp2[2] );
        double spin2 = BH2.getBHSpin();
        double spin2_th = BH2.getBHSpinTheta();
        double spin2_phi = BH2.getBHSpinPhi();

        double vs2[3];
        vs2[0] = spin2*sin(spin2_th)*cos(spin2_phi);
        vs2[1] = spin2*sin(spin2_th)*sin(spin2_phi);
        vs2[2] = spin2*cos(spin2_th);


        // coordinates with respect to center of bh1
        x1 = xx - bh1x;
        y1 = yy - bh1y;
        z1 = zz - bh1z;

        //locating as a radial form
        rv1 = sqrt(x1*x1 + y1*y1 + z1*z1);
        vn1[0] = x1/rv1;
        vn1[1] = y1/rv1;
        vn1[2] = z1/rv1;

        //same as BH2
        x2 = xx - bh2x;
        y2 = yy - bh2y;
        z2 = zz - bh2z;

        rv2 = sqrt(x2*x2 + y2*y2 + z2*z2);
        vn2[0] = x2/rv2;
        vn2[1] = y2/rv2;
        vn2[2] = z2/rv2;

        //Initial data is related with the paper: http://arxiv.org/abs/0711.1165
        //Brill-Lindquist conformal factor
        vpsibl = 1.0 + mass1/(2.0*rv1);
        vpsibl = vpsibl + mass2/(2.0*rv2);

        v_u_corr = 0.0;
        // bh 1

        //For spinning puncture
        if ( fabs(spin1) > 1.e-6 ) {
            amp_capj = 4.0*spin1/(mass1*mass1);
            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0+amp_capr);
            u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
            u2_j = -pow(l_r,5)/20.0;
            mu_j = vn1[0]*vs1[0];
            mu_j = mu_j + vn1[1]*vs1[1];
            mu_j = (mu_j + vn1[2]*vs1[2])/fabs(spin1);
            p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
            v_u_j1 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
            v_u_corr = v_u_corr + v_u_j1;
        }
        //For boosting puncture
        if (vp1tot > 1.e-6) {
            amp_capp = 2.0*vp1tot/mass1;
            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0 + amp_capr);
            u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
            u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
            u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
            u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
            u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
            u2_p = (u2_p)/(80.0*amp_capr);
            mu_p =        vn1[0]*vp1[0]/vp1tot;
            mu_p = mu_p + vn1[1]*vp1[1]/vp1tot;
            mu_p = mu_p + vn1[2]*vp1[2]/vp1tot;
            p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
            v_u_p1 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
            v_u_corr = v_u_corr + v_u_p1;
        }
        //For spinning boosted pucture
        if ( vp1tot > 1.e-6 && fabs(spin1) > 1.e-6 ) {
            v1 =      (vp1[1]*vs1[2]-vp1[2]*vs1[1])*vn1[0];
            v1 = v1 + (vp1[2]*vs1[0]-vp1[0]*vs1[2])*vn1[1];
            v1 = v1 + (vp1[0]*vs1[1]-vp1[1]*vs1[0])*vn1[2];
            v1 = v1*(16.0/pow(mass1,4))*rv1;

            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0 + amp_capr);

            v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);

            v_u_c1 = (v1*v2*pow(l_r,5))/80.0;
            v_u_corr = v_u_corr + v_u_c1;
        }
        // bh 2 same puncture as bh 1
        if ( fabs(spin2) > 1.e-6 ) {
            amp_capj = 4.0*spin2/(mass2*mass2);
            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0+amp_capr);
            u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
            u2_j = -pow(l_r,5)/20.0;
            mu_j = vn2[0]*vs2[0];
            mu_j = mu_j + vn2[1]*vs2[1];
            mu_j = (mu_j + vn2[2]*vs2[2])/fabs(spin2);
            p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
            v_u_j2 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
            v_u_corr = v_u_corr + v_u_j2;
        }

        if (vp2tot > 1.e-6) {
            amp_capp = 2.0*vp2tot/mass2;
            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0 + amp_capr);
            u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
            u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
            u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
            u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
            u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
            u2_p = (u2_p)/(80.0*amp_capr);
            mu_p =        vn2[0]*vp2[0]/vp2tot;
            mu_p = mu_p + vn2[1]*vp2[1]/vp2tot;
            mu_p = mu_p + vn2[2]*vp2[2]/vp2tot;
            p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
            v_u_p2 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
            v_u_corr = v_u_corr + v_u_p2;
        }

        if ( vp2tot > 1.e-6 && fabs(spin2) > 1.e-6 ) {
            v1 =      (vp2[1]*vs2[2]-vp2[2]*vs2[1])*vn2[0];
            v1 = v1 + (vp2[2]*vs2[0]-vp2[0]*vs2[2])*vn2[1];
            v1 = v1 + (vp2[0]*vs2[1]-vp2[1]*vs2[0])*vn2[2];
            v1 = v1*(16.0/pow(mass2,4))*rv2;

            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0 + amp_capr);

            v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);

            v_u_c2 = (v1*v2*pow(l_r,5))/80.0;
            v_u_corr = v_u_corr + v_u_c2;
        }

        // vpsibl_u will be used for the conformal factor,
        vpsibl_u  = vpsibl + v_u_corr;
        // vpsibl_u2 is for the Aij terms...
        // ! since the corrections are first order...
        // ! adding half of the correction seems to give the best results...
        // ! update - do a fit for spin = 0.6...
        vpsibl_u2 = vpsibl + v_u_corr;

        var[VAR::U_ALPHA] = 1.0/(vpsibl_u*vpsibl_u);
        //std::cout<<"Alpha: "<<u[U_ALPHA]<<" vpsibl_u: "<< vpsibl_u<<std::endl;

        v2 = 1.0/pow(vpsibl_u,4);
        //var[VAR::U_CHI] = v2;
        var[VAR::U_PSI] = v2;
        var[VAR::U_K] = 0.0;
        var[VAR::U_THETA_Z4] = 0.0;

        var[VAR::U_BETA0] = 0.0;
        var[VAR::U_BETA1] = 0.0;
        var[VAR::U_BETA2] = 0.0;

        var[VAR::U_GH0] = 0.0;
        var[VAR::U_GH1] = 0.0;
        var[VAR::U_GH2] = 0.0;

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;

        var[VAR::U_SYMGT0] = 1.0; //XX
        var[VAR::U_SYMGT1] = 0.0; //XY
        var[VAR::U_SYMGT2] = 0.0; //XZ
        var[VAR::U_SYMGT3] = 1.0; //YY
        var[VAR::U_SYMGT4] = 0.0; //YZ
        var[VAR::U_SYMGT5] = 1.0; //ZZ

        for (i1=0;i1<3;i1++) {
            for (i2=0;i2<3;i2++) {
                // first BH
                v2 = 0.0;
                for (i3=0;i3<3;i3++) {
                    for (i4=0;i4<3;i4++) {
                        vt1 = epijk[i1][i3][i4]*vs1[i3]*vn1[i4]*vn1[i2];
                        vt2 = epijk[i2][i3][i4]*vs1[i3]*vn1[i4]*vn1[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp1[i1]*vn1[i2] + vp1[i2]*vn1[i1];
                vt1 = 0.0;
                for (i3=0;i3<3;i3++) {
                    vt1 = vt1 + vp1[i3]*vn1[i3];
                }
                vt1 = vt1*(vn1[i1]*vn1[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv1,3));
                v4 = v1*(v2+(rv1/2.0)*v3);

                // second BH
                v2 = 0.0;
                for (i3=0;i3<3;i3++) {
                    for (i4=0;i4<3;i4++) {
                        vt1 = epijk[i1][i3][i4]*vs2[i3]*vn2[i4]*vn2[i2];
                        vt2 = epijk[i2][i3][i4]*vs2[i3]*vn2[i4]*vn2[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp2[i1]*vn2[i2] + vp2[i2]*vn2[i1];
                vt1 = 0.0;
                for (i3=0;i3<3;i3++) {
                    vt1 = vt1 + vp2[i3]*vn2[i3];
                }
                vt1 = vt1*(vn2[i1]*vn2[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv2,3));
                v4 = v4 + v1*(v2+(rv2/2.0)*v3);

                if ( i1 == 0 && i2 == 0 ) {
                    var[VAR::U_SYMAT0] = v4; //XX
                } else if ( i1 == 0 && i2 == 1 ) {
                    var[VAR::U_SYMAT1] = v4; //XY
                } else if ( i1 == 0 && i2 == 2 ) {
                    var[VAR::U_SYMAT2] = v4; //XZ
                } else if ( i1 == 1 && i2 == 1 ) {
                    var[VAR::U_SYMAT3] = v4; //YY
                } else if ( i1 == 1 && i2 == 2 ) {
                    var[VAR::U_SYMAT4] = v4; //YZ
                } else if ( i1 == 2 && i2 == 2 ) {
                    var[VAR::U_SYMAT5] = v4; //ZZ
                }

            }
        }

    }


    void KerrSchildData(const double xx1,const double yy1,const double zz1, double *var)
    {
        
        const double xx=GRIDX_TO_X(xx1);
        const double yy=GRIDY_TO_Y(yy1);
        const double zz=GRIDZ_TO_Z(zz1);

        double x,y,z,r;

        double M = BH1.getBHMass();
        double bh1x = BH1.getBHCoordX();
        double bh1y = BH1.getBHCoordY();
        double bh1z = BH1.getBHCoordZ();

        double vp1[3];
        vp1[0] = BH1.getVx();
        vp1[1] = BH1.getVy();
        vp1[2] = BH1.getVz();

        double vp1tot = sqrt( vp1[0]*vp1[0] + vp1[1]*vp1[1] + vp1[2]*vp1[2] );
        double spin1 = BH1.getBHSpin();
        double spin1_th = BH1.getBHSpinTheta();
        double spin1_phi = BH1.getBHSpinPhi();
        double vs1[3];
        vs1[0] = spin1*sin(spin1_th)*cos(spin1_phi);
        vs1[1] = spin1*sin(spin1_th)*sin(spin1_phi);
        vs1[2] = spin1*cos(spin1_th);
 
        // coordinates with respect to center of bh1
        x = xx - bh1x;
        y = yy - bh1y;
        z = zz - bh1z;

        //locating as a radial form
        r = sqrt(x*x + y*y + z*z);
        double vn1[3];
        vn1[0] = x/r;
        vn1[1] = y/r;
        vn1[2] = z/r;

        //HL : Angular momentum parameter will be added as param file after testing
        //     value for ANG_PAR should be between 0 and 1
        double a = ccz4::ANG_PAR;

        double gtd[3][3], Atd[3][3];
        double alpha, Gamh[3];
        double Psi, TrK, Betau[3];

        #include "ks_vars.cpp"
        #include "ksinit.cpp"

        var[VAR::U_ALPHA] = alpha;
        var[VAR::U_PSI] = Psi;
        var[VAR::U_K] = TrK;

        var[VAR::U_BETA0] = Betau[0];
        var[VAR::U_BETA1] = Betau[1];
        var[VAR::U_BETA2] = Betau[2];

        var[VAR::U_GH0] = Gamh[0];
        var[VAR::U_GH1] = Gamh[1];
        var[VAR::U_GH2] = Gamh[2];

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;

        var[VAR::U_SYMGT0] = gtd[0][0];
        var[VAR::U_SYMGT1] = gtd[0][1];
        var[VAR::U_SYMGT2] = gtd[0][2];
        var[VAR::U_SYMGT3] = gtd[1][1];
        var[VAR::U_SYMGT4] = gtd[1][2];
        var[VAR::U_SYMGT5] = gtd[2][2];

        var[VAR::U_SYMAT0] = Atd[0][0];
        var[VAR::U_SYMAT1] = Atd[0][1];
        var[VAR::U_SYMAT2] = Atd[0][2];
        var[VAR::U_SYMAT3] = Atd[1][1];
        var[VAR::U_SYMAT4] = Atd[1][2];
        var[VAR::U_SYMAT5] = Atd[2][2];

    }


    void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_size(comm,&npes);
        MPI_Comm_rank(comm,&rank);

        double pt_g_min[3];
        double pt_g_max[3];

        pt_g_min[0]=X_TO_GRIDX(pt_min.x());
        pt_g_min[1]=Y_TO_GRIDY(pt_min.y());
        pt_g_min[2]=Z_TO_GRIDZ(pt_min.z());

        pt_g_max[0]=X_TO_GRIDX(pt_max.x());
        pt_g_max[1]=Y_TO_GRIDY(pt_max.y());
        pt_g_max[2]=Z_TO_GRIDZ(pt_max.z());

        assert(pt_g_min[0]>=0 && pt_g_min[0]<(1u<<maxDepth));
        assert(pt_g_min[1]>=0 && pt_g_min[1]<(1u<<maxDepth));
        assert(pt_g_min[2]>=0 && pt_g_min[2]<(1u<<maxDepth));

        assert(pt_g_max[0]>=0 && pt_g_max[0]<(1u<<maxDepth));
        assert(pt_g_max[1]>=0 && pt_g_max[1]<(1u<<maxDepth));
        assert(pt_g_max[2]>=0 && pt_g_max[2]<(1u<<maxDepth));


        unsigned int xRange_b,xRange_e;
        unsigned int yRange_b=pt_g_min[1],yRange_e=pt_g_max[1];
        unsigned int zRange_b=pt_g_min[2],zRange_e=pt_g_max[2];

        xRange_b=pt_g_min[0];//(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
        xRange_e=pt_g_max[1];//((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

        unsigned int stepSz=1u<<(maxDepth-regLev);

       /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
        std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
        std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/


        for(unsigned int x=xRange_b;x<xRange_e;x+=stepSz)
            for(unsigned int y=yRange_b;y<yRange_e;y+=stepSz)
              for(unsigned int z=zRange_b;z<zRange_e;z+=stepSz)
                  tmpNodes.push_back(ot::TreeNode(x,y,z,regLev,m_uiDim,maxDepth));


        return ;


    }

      double computeWTol(double x,double y,double z,double tol_min)
      {
         double origin[3];
         origin[0]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
         origin[1]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
         origin[2]=(double)(1u<<ccz4::CCZ4_MAXDEPTH-1);
 
         double r2 = (x-origin[0])*(x-origin[0]) + (y-origin[1])*(y-origin[1]) + (z-origin[2])*(z-origin[2]);
           return (r2+1)*tol_min;
 
      }







}// end of namespace ccz4





namespace ccz4
{

    namespace timer
    {
        void initFlops()
        {
            total_runtime.start();
            t_f2o.start();
            t_cons.start();
            t_bal.start();
            t_mesh.start();
            t_rkSolve.start();
            t_ghostEx_sync.start();
            t_unzip_sync.start();

            for(unsigned int i=0;i<NUM_FACES;i++)
                dendro::timer::t_unzip_sync_face[i].start();

            dendro::timer::t_unzip_async_internal.start();
            dendro::timer::t_unzip_sync_edge.start();
            dendro::timer::t_unzip_sync_vtex.start();
            dendro::timer::t_unzip_p2c.start();
            dendro::timer::t_unzip_sync_nodalval.start();
            dendro::timer::t_unzip_sync_cpy.start();
            dendro::timer::t_unzip_sync_f_c1.start();
            dendro::timer::t_unzip_sync_f_c2.start();
            dendro::timer::t_unzip_sync_f_c3.start();

            t_unzip_async.start();
            dendro::timer::t_unzip_async_comm.start();

            dendro::timer::t_unzip_async_internal.start();
            dendro::timer::t_unzip_async_external.start();
            dendro::timer::t_unzip_async_comm.start();
            t_deriv.start();
            t_rhs.start();

            t_rhs_a.start();
            t_rhs_b.start();
            t_rhs_gt.start();
            t_rhs_psi.start();
            t_rhs_At.start();
            t_rhs_K.start();
            t_rhs_Gh.start();
            t_rhs_B.start();
            t_rhs_theta_z4.start();

            t_bdyc.start();

            t_zip.start();
            t_rkStep.start();
            t_isReMesh.start();
            t_gridTransfer.start();
            t_ioVtu.start();
            t_ioCheckPoint.start();
        }

        void resetSnapshot()
        {

            total_runtime.snapreset();
            t_f2o.snapreset();
            t_cons.snapreset();
            t_bal.snapreset();
            t_mesh.snapreset();
            t_rkSolve.snapreset();
            t_ghostEx_sync.snapreset();
            t_unzip_sync.snapreset();

            for(unsigned int i=0;i<NUM_FACES;i++)
                dendro::timer::t_unzip_sync_face[i].snapreset();

            dendro::timer::t_unzip_sync_internal.snapreset();
            dendro::timer::t_unzip_sync_edge.snapreset();
            dendro::timer::t_unzip_sync_vtex.snapreset();
            dendro::timer::t_unzip_p2c.snapreset();
            dendro::timer::t_unzip_sync_nodalval.snapreset();
            dendro::timer::t_unzip_sync_cpy.snapreset();

            dendro::timer::t_unzip_sync_f_c1.snapreset();
            dendro::timer::t_unzip_sync_f_c2.snapreset();
            dendro::timer::t_unzip_sync_f_c3.snapreset();

            t_unzip_async.snapreset();
            dendro::timer::t_unzip_async_internal.snapreset();
            dendro::timer::t_unzip_async_external.snapreset();
            dendro::timer::t_unzip_async_comm.snapreset();

            t_deriv.snapreset();
            t_rhs.snapreset();

            t_rhs_a.snapreset();
            t_rhs_b.snapreset();
            t_rhs_gt.snapreset();
            t_rhs_psi.snapreset();
            t_rhs_At.snapreset();
            t_rhs_K.snapreset();
            t_rhs_Gh.snapreset();
            t_rhs_B.snapreset();

            t_rhs_theta_z4.snapreset();

            t_bdyc.snapreset();

            t_zip.snapreset();
            t_rkStep.snapreset();
            t_isReMesh.snapreset();
            t_gridTransfer.snapreset();
            t_ioVtu.snapreset();
            t_ioCheckPoint.snapreset();

        }

        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh)
        {


            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            if(!activeRank)
            {
                sprintf(fName,"%s_final.prof",filePrefix);
                outfile.open (fName);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"partition tol : "<<ccz4::CCZ4_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<ccz4::CCZ4_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<ccz4::CCZ4_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;


            t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_bal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=dendro::timer::t_unzip_async_internal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_external.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_comm.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
#endif

            t_stat=t_deriv.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_bdyc.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();




        }


        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep)
        {

            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            DendroIntL ghostElements;
            DendroIntL localElements;

            DendroIntL ghostNodes;
            DendroIntL localNodes;

            DendroIntL totalSendNode;
            DendroIntL totalRecvNode;

            DendroIntL numCalls;


#ifdef BSSN_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<ccz4::CCZ4_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<ccz4::CCZ4_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<ccz4::CCZ4_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            if(!rank)outfile<<"========================= MESH ======================================================================= "<<std::endl;

            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(#)"<<std::endl;

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"send Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"recv Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank)outfile<<"========================= RUNTIME =================================================================== "<<std::endl;
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;




            /* t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/


            t_stat=t_bal.snap;
            //numCalls=t_bal.num_calls;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm_wait (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_nodalval.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_nodalVal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c1.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c1";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c2.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c2";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c3.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c3";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_cpy.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_cpy";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[0].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_left";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[1].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_right";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[2].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_down";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[3].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_up";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[4].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_back";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[5].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_front";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_edge.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_edge";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_vtex.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_vtex";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_p2c.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_p2c";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif

            /*
            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=t_unzip_async_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async_external.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif
            */
            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++isReMesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++gridTransfer";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_rhs_a.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_a ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_b.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_b ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_gt.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_gt ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rhs_psi.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_psi ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_theta_z4.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_theta_z4 ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rhs_At.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_At ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_K.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_K ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_Gh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_Gh ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_B.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_B ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_bdyc.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();
#else

            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                //writes the header
                if(currentStep==0)
                 outfile<<"step\t act_npes\t glb_npes\t part_tol\t wave_tol\t maxdepth\t numOcts\t dof_zip\t dof_unzip\t"<<\
                 "element_ghost_min\t element_ghost_mean\t element_ghost_max\t"<<\
                 "element_local_min\t element_local_mean\t element_local_max\t"<<\
                 "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"<<\
                 "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"<<\
                 "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"<<\
                 "bal_min\t bal_mean\t bal_max\t"<<\
                 "mesh_min\t mesh_mean\t mesh_max\t"<<\
                 "rkstep_min\t rkstep_mean\t rkstep_max\t"<<\
                 "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"<<\
                 "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"<<\
                 "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"<<\
                 "unzip_async_wait_min\t unzip_async_wait_mean\t unzip_async_wait_max\t"<<\
                 "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"<<\
                 "GT_min\t GT_mean\t GT_max\t"<<\
                 "deriv_min\t deriv_mean\t deriv_max\t"<<\
                 "rhs_min\t rhs_mean\t rhs_max\t"<<std::endl;

            }
	    
            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            if(!rank) outfile<<currentStep<<"\t ";
            if(!rank) outfile<<activeNpes<<"\t ";
            if(!rank) outfile<<globalNpes<<"\t ";
            if(!rank) outfile<<ccz4::CCZ4_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<ccz4::CCZ4_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<ccz4::CCZ4_MAXDEPTH<<"\t ";

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            /*t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";*/

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_bal.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            if(!rank) outfile<<std::endl;
            if(!rank) outfile.close();
#endif




        }


    }

}
