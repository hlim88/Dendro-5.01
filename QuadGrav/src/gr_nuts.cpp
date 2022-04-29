/**
 * @file gr_nuts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Solving BSSN equations using non uniform time stepping. 
 * @version 0.1
 * @date 2020-03-04
 * 
 * @copyright Copyright (c) 2020, School of Computing, University of Utah. 
 * 
 */


#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rkBSSN.h"
#include "octUtils.h"
#include "meshUtils.h"


int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    unsigned int ts_mode=0;     
    
    if(argc<2)
    {
        std::cout<<"Usage: "<<argv[0]<<"paramFile TSMode(0){0-Spatially Adaptive Time Stepping(SATS, "<<GRN<<"default"<<NRM<<") , 1- Uniform Time Stepping.  }"<<std::endl;
        return 0;
    }
        
    if(argc>2)
        ts_mode = std::atoi(argv[2]);


    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);


    // Print out CMAKE options
    if (!rank) {
        #ifdef BSSN_COMPUTE_CONSTRAINTS
          std::cout<<GRN<<"  Compiled with BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ETA_FUNCTION 
          std::cout<<GRN<<"  Compiled with  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_BH_LOCATIONS 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_GAUGE_ROCHESTER 
          std::cout<<GRN<<"  Compiled with  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_KERR_SCHILD_TEST 
          std::cout<<GRN<<"  Compiled with  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #endif

        #ifdef BSSN_REFINE_BASE_EH 
          std::cout<<GRN<<"  Compiled with  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #endif

        #ifdef USE_FD_INTERP_FOR_UNZIP 
          std::cout<<GRN<<"  Compiled with  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #endif

    }

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    bssn::readParamFile(argv[1],comm);


    int root = std::min(1,npes-1);
    bssn::dumpParamFile(std::cout,root,comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth=bssn::BSSN_MAXDEPTH;
    
    if(bssn::BSSN_NUM_VARS%bssn::BSSN_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total BSSN_NUM_VARS: "<<bssn::BSSN_NUM_VARS<<" is not divisable by BSSN_ASYNC_COMM_K: "<<bssn::BSSN_ASYNC_COMM_K<<std::endl;
        MPI_Abort(comm,0);
    }

    if(bssn::BSSN_GW_EXTRACT_FREQ> bssn::BSSN_IO_OUTPUT_FREQ)
    {
      if(!rank) std::cout<<" BSSN_GW_EXTRACT_FREQ  should be less BSSN_IO_OUTPUT_FREQ "<<std::endl;
      MPI_Abort(comm,0);
    }


    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::punctureData(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; bssn::punctureData(x,y,z,var); return var[0];};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<bssn::BSSN_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    bssn::timer::t_f2o.start();

    if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(bssn::BSSN_BLK_MIN_X,bssn::BSSN_BLK_MIN_Y,bssn::BSSN_BLK_MIN_Z);
        const Point pt_max(bssn::BSSN_BLK_MAX_X,bssn::BSSN_BLK_MAX_Y,bssn::BSSN_BLK_MAX_Z);

        bssn::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        const unsigned int f2olmin = std::min(bssn::BSSN_BH1_MAX_LEV,bssn::BSSN_BH2_MAX_LEV);
        if(f2olmin < MAXDEAPTH_LEVEL_DIFF + 2)
        {
          if(!rank)
            std::cout<<"BH min level should be larger than "<<(MAXDEAPTH_LEVEL_DIFF+2)<<std::endl;

          MPI_Abort(comm,0);
          
        }
        function2Octree(f_init,bssn::BSSN_NUM_VARS,varIndex,interpVars,tmpNodes,(f2olmin-MAXDEAPTH_LEVEL_DIFF-2),bssn::BSSN_WAVELET_TOL,bssn::BSSN_ELE_ORDER,comm);
        
    }

    if(ts_mode == 0)
    { 

        //std::vector<ot::TreeNode> f2Octants(tmpNodes);
        ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
        mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
        unsigned int lmin, lmax;
        mesh->computeMinMaxLevel(lmin,lmax);
        if(!rank)
        {
          std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
          std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
          std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
          std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
          std::cout<<"ts mode: "<<ts_mode<<std::endl;
          std::cout<<"==============================================================================================================="<<std::endl;
        }    
        bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        tmpNodes.clear();

        //NUTS
        assert(ot::test::isBlkFlagsValid(mesh));
        
        bssn::BSSNCtx *  bssnCtx = new bssn::BSSNCtx(mesh); 
        ts::ExplicitNUTS<DendroScalar,bssn::BSSNCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,bssn::BSSNCtx>(bssnCtx);

        std::vector<double> ld_stat_g;
        enuts->set_evolve_vars(bssnCtx->get_evolution_vars());
        
        if((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
            enuts->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
            enuts->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
            enuts->set_ets_coefficients(ts::ETSType::RK5);
        
        const unsigned int rank_global = enuts->get_global_rank();
        const unsigned int lts_remesh_freq=2;

        //for(enuts->init(); enuts->curr_time() < bssn::BSSN_RK_TIME_END ; enuts->evolve_with_remesh(lts_remesh_freq))
        for(enuts->init(); enuts->curr_time() < bssn::BSSN_RK_TIME_END ; enuts->evolve())
        {
            
          const DendroIntL step = enuts->curr_step();
          const DendroScalar time = enuts->curr_time();    
          const bool isActive = enuts->is_active();

          if(step==0 || bssn::BSSN_RESTORE_SOLVER)
            bssnCtx->compute_lts_ts_offset();
           
          if(!rank_global)
          {
            std::cout<<GRN<<"[Explicit LTS] : lts offset: "<<bssn::BSSN_LTS_TS_OFFSET<<NRM<<std::endl;
            std::cout<<GRN<<"[Explicit LTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;
          }
              

          if(step>0)
            bssnCtx->evolve_bh_loc(bssnCtx->get_evolution_vars(),enuts->get_dt_max());
          
          bssnCtx->terminal_output();
          enuts->dump_load_statistics(std::cout);

          const bool is_merged = bssnCtx->is_bh_merged(0.1);
          if(is_merged)
            bssn::BSSN_REMESH_TEST_FREQ=bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;  

          bool isRemesh = false;    
          if( (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 )
              isRemesh = bssnCtx->is_remesh();

          if(isRemesh)
          {
              if(!rank_global)
                  std::cout<<"[Explicit LTS]: Remesh triggered "<<std::endl;;

              bssnCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX,true,false,false);
              bssnCtx->compute_lts_ts_offset();
          }
            
          enuts->sync_with_mesh();

          if((step % bssn::BSSN_GW_EXTRACT_FREQ) == 0 )
            bssnCtx -> write_vtu();   

          if( (step % bssn::BSSN_CHECKPT_FREQ) == 0 )
            bssnCtx -> write_checkpt();

        }

        delete bssnCtx->get_mesh();    
        delete bssnCtx;
        delete enuts;


    }else if(ts_mode==1)
    { 


        //std::vector<ot::TreeNode> f2Octants(tmpNodes);
        ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
        mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
        unsigned int lmin, lmax;
        mesh->computeMinMaxLevel(lmin,lmax);
        if(!rank)
        {
          std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
          std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
          std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
          std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
          std::cout<<"ts mode: "<<ts_mode<<std::endl;
          std::cout<<"==============================================================================================================="<<std::endl;
        }    
        bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        tmpNodes.clear();

        //UTS
        bssn::BSSNCtx *  bssnCtx = new bssn::BSSNCtx(mesh);
        ts::ETS<DendroScalar,bssn::BSSNCtx>* ets = new ts::ETS<DendroScalar,bssn::BSSNCtx>(bssnCtx);
        ets->set_evolve_vars(bssnCtx->get_evolution_vars());
        
        if((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
            ets->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
            ets->set_ets_coefficients(ts::ETSType::RK5);


        for(ets->init(); ets->curr_time() < bssn::BSSN_RK_TIME_END ; ets->evolve())
        {
          const DendroIntL   step = ets->curr_step();
          const DendroScalar time = ets->curr_time();    

          const bool isActive = ets->is_active();
          const unsigned int rank_global = ets->get_global_rank();

          if(!rank_global)
          std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

          if( (step% bssn::BSSN_TIME_STEP_OUTPUT_FREQ) == 0)
            bssnCtx->terminal_output();  

          if(step>0)
            bssnCtx->evolve_bh_loc(bssnCtx->get_evolution_vars(),ets->ts_size());

          const bool is_merged = bssnCtx->is_bh_merged(0.1);
          if(is_merged)
            bssn::BSSN_REMESH_TEST_FREQ=bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;  

          bool isRemesh = false;    
          if( (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 )
            isRemesh = bssnCtx->is_remesh();

          if(isRemesh)
          {
              if(!rank_global)
                  std::cout<<"[ETS] : Remesh is triggered.  \n";

              bssnCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX,true,false,false);
          }
          
          ets->sync_with_mesh();

          if((step % bssn::BSSN_GW_EXTRACT_FREQ) == 0 )
          bssnCtx -> write_vtu();   

          if( (step % bssn::BSSN_CHECKPT_FREQ) == 0 )
          bssnCtx -> write_checkpt();
            
        }

        delete bssnCtx->get_mesh();    
        delete bssnCtx;
        delete ets;

    }else if(ts_mode ==2)
    {
      profiler_t t_rt;
      t_rt.clear();


      ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
      mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
      unsigned int lmin, lmax;
      mesh->computeMinMaxLevel(lmin,lmax);
      if(!rank)
      {
        std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
        std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
        std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
        std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
        std::cout<<"ts mode: "<<ts_mode<<std::endl;
        std::cout<<"==============================================================================================================="<<std::endl;
      }    
      bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));

      // perform a comparison test between ets and enuts. 
      bssn::BSSNCtx *  bssnCtx_enuts = new bssn::BSSNCtx(mesh); 
      
      ts::ExplicitNUTS<DendroScalar,bssn::BSSNCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,bssn::BSSNCtx>(bssnCtx_enuts);
      enuts -> set_evolve_vars(bssnCtx_enuts->get_evolution_vars());


      if((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
      {
        enuts -> set_ets_coefficients(ts::ETSType::RK3);

      }else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
      {
        enuts -> set_ets_coefficients(ts::ETSType::RK4);

      }else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
      {
        enuts -> set_ets_coefficients(ts::ETSType::RK5);
      }
      
      
      const unsigned int rank_global = enuts->get_global_rank();

      if(!rank_global)
        std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

      const unsigned int num_lts_steps=1;

      
      
      enuts->init();
      
      bssnCtx_enuts->compute_lts_ts_offset();
      enuts->dump_load_statistics(std::cout);
      
      std::ofstream fout_est_speedup;
      fout_est_speedup.open ("lts_est_speedup.txt");
      enuts->dump_est_speedup(fout_est_speedup,false);
      fout_est_speedup.close();
      
      if(!rank_global)
        std::cout<<" LTS warm up run begin"<<std::endl;

      t_rt.snapreset();

      t_rt.start();
        enuts->evolve();
      t_rt.stop();

      if(!rank_global)
        std::cout<<" LTS warm up run end"<<std::endl;
      
      t_stat=t_rt.snap;
      bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
      if(!rank_global)
        printf("[LTS] warm up time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

      t_rt.snapreset();

      enuts->get_mesh()->waitAll();
      
      while (enuts->curr_step() < num_lts_steps+1)
      {
        t_rt.start();
          
          enuts->evolve();

        t_rt.stop();

        if(!rank_global)
          {
            std::cout<<GRN<<"[Explicit LTS] : lts offset: "<<bssn::BSSN_LTS_TS_OFFSET<<NRM<<std::endl;
            std::cout<<GRN<<"[Explicit LTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;
          }


      }
  
      t_stat=t_rt.snap;
      bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
      
      if(!rank_global)
        printf("[LTS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

      //=========================================== Now do GTS ===============================================================================================================================

      mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
      mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
      mesh->computeMinMaxLevel(lmin,lmax);
      if(!rank)
      {
        std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
        std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
        std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
        std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
        std::cout<<"ts mode: "<<ts_mode<<std::endl;
        std::cout<<"==============================================================================================================="<<std::endl;
      }    
      bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));

      bssn::BSSNCtx *  bssnCtx_ets = new bssn::BSSNCtx(mesh); 
      ts::ETS<DendroScalar,bssn::BSSNCtx>*           ets   = new ts::ETS<DendroScalar,bssn::BSSNCtx>(bssnCtx_ets);
      ets   -> set_evolve_vars(bssnCtx_ets->get_evolution_vars());
      
      if((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
      {
        ets -> set_ets_coefficients(ts::ETSType::RK3);

      }else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
      {
        ets -> set_ets_coefficients(ts::ETSType::RK4);

      }else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
      {
        ets -> set_ets_coefficients(ts::ETSType::RK5);
      }

      t_rt.snapreset();

      unsigned int num_steps = ( enuts->get_dt_max() / enuts->get_dt_min() );
      if(!rank_global) 
        std::cout<<"GTS to LTS step ratio: "<<num_steps<<std::endl;
      const unsigned int num_gts_steps = num_lts_steps*num_steps;
      
      

      ets->init();
      ets->dump_load_statistics(std::cout);
      
      t_rt.snapreset();
      if(!rank_global)
        std::cout<<" GTS warm up run begin "<<std::endl;

      t_rt.start();
        ets->evolve();
      t_rt.stop();

      if(!rank_global)
        std::cout<<" GTS warm up run end "<<std::endl;

      t_stat=t_rt.snap;
      bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
      
      if(!rank_global)
        printf("[LTS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);
      

      // while(ets->curr_step() < num_gts_steps+1)
      // {
        
      //   if(!rank_global)
      //       std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

      //   t_rt.start();
      //     ets->evolve();
      //   t_rt.stop();

      //   t_stat=t_rt.snap;
      //   bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
      //   if(!rank_global)
      //     printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

      // }

      t_rt.snapreset();
      const unsigned int num_avg_steps=5;
      while(ets->curr_step() < num_avg_steps+1)
      {
        
        if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

        t_rt.start();
          ets->evolve();
        t_rt.stop();

        t_stat=t_rt.snap;
        bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
        if(!rank_global)
          printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

      }

      t_stat=t_rt.snap;
      bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);

      if(!rank_global)
              std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

      const double tscaling = num_gts_steps/num_avg_steps;
      if(!rank_global)
        printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0]*tscaling, t_stat_g[1]*tscaling, t_stat_g[2]*tscaling);

      
      delete bssnCtx_enuts->get_mesh();    
      delete bssnCtx_enuts;

      delete bssnCtx_ets->get_mesh();    
      delete bssnCtx_ets;
      
      delete enuts;
      delete ets;
        


    }


    MPI_Finalize();
    return 0; 


}
