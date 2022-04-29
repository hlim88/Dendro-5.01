     program      main
      implicit     none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program takes clms (up to l=4) from the already
! recombined (from desired batches) run with harmonic3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
            
      integer :: i, j, k, N, itheta, iphi, output, initial, itime
      integer :: correction, integration, ISTART, imax
      real*8  :: pi, theta, phi, integral, h, h2, x, y, localtime
      real*8  :: rp, ip, norm, dt, factor, temp1, time_tort,integralclm(2:4,-4:4,2)
      real*8  :: massI, rext, test, smallnumber, ra, rr, fpR, fpI
      real*8  :: surface_integral,surface_integral2, givefile,cleanhdot, hitwithr
      integer :: cphase, l, m
            
! ------- this parameter changes depending on the size of u------------
      integer Ndim
!      parameter ( Ndim = 80 )


      character*20    file_psi(2:4,-4:4)
      character*20    file_psi_inf(2:4,-4:4)
      integer         read_shape(2), shape(2)
      character*5     read_cnames
      integer         read_rank, myi
      real*8          read_coords(1000), bbox(4) 
      real*8          read_time, timeold, time_min
      real*8          r8arg
      external        r8arg
      
      integer         iargc, gf3_rc, gft_read_full,index
      character*5     name, sarg
      external        sarg, gft_read_full
      character*32    c2p2file,c2p2_inffile,massname,jzname




      integer     MAXMEM
!      parameter ( MAXMEM = 1000000 )
      parameter ( MAXMEM = 50000 )
  

      integer     MAXTIME
      parameter ( MAXTIME = 50000 )

      real*8  :: C0, factorl2,factorl3,factorl4, cr2p2,ci2p2,dtcI2p2,dtcR2p2

      real*8  :: clm(2:4,-4:4,2,MAXTIME),dtclm(2:4,-4:4,2,MAXTIME)
      real*8  :: clm_inf(2:4,-4:4,2,MAXTIME),dtclm_inf(2:4,-4:4,2,MAXTIME)
      real*8  :: d2tclm(2:4,-4:4,2,MAXTIME), d2tclm_inf(2:4,-4:4,2,MAXTIME)
      real*8  :: hlm(2:4,-4:4,2,MAXTIME),hlm_inf(2:4,-4:4,2,MAXTIME)

      real*8  :: omega1(MAXTIME),omega2(MAXTIME),LGW(MAXTIME),LGW_inf(MAXTIME)
      real*8  :: JGW(MAXTIME),JGW_inf(MAXTIME)
      real*8  :: time(MAXTIME)
      real*8  :: tortutime(MAXTIME)
      real*8  :: phase_inf(MAXTIME), phase(MAXTIME), amp_inf
      real*8  :: HHphase_inf(MAXTIME), HHamp_inf, HHphase(MAXTIME), HHamp
      real*8  :: hlm_ave(2:4,-4:4,2), hlm_inf_ave(2:4,-4:4,2), time_ave
      real*8  :: xxtime_ave, xyhlm_ave(2:4,-4:4,2), xyhlm_inf_ave(2:4,-4:4,2)
      real*8 ::  count_ave
      real*8  :: slopehlm(2:4,-4:4,2),interhlm(2:4,-4:4,2), &
     &           slopehlm_inf(2:4,-4:4,2),interhlm_inf(2:4,-4:4,2)

!for quadratic regretion. things have become too messy, so defining new variables, sorry about that
     real*8 ::s11,s12,s22
     real*8 ::sy1(2:4,-4:4,2),sy2(2:4,-4:4,2)
     real*8 ::s11I,s12I,s22I
     real*8 ::sy1I(2:4,-4:4,2),sy2I(2:4,-4:4,2)
     real*8 ::x1ave,x2ave,yave(2:4,-4:4,2),yaveI(2:4,-4:4,2)
     real*8 ::bet1(2:4,-4:4,2),bet2(2:4,-4:4,2),bet3(2:4,-4:4,2)
     real*8 ::bet1I(2:4,-4:4,2),bet2I(2:4,-4:4,2),bet3I(2:4,-4:4,2)
      
!for output with letter so as to keep the file names
      character*32 :: Fhlm(2:4,-4:4),Fhlm_inf(2:4,-4:4),Fdlm(2:4,-4:4),Fdlm_inf(2:4,-4:4),&
     &                Fomega,Fphase, Fampphase,Fampphase_inf,Fhhampphase,Fhhampphase_inf, &
     &                Fenergy, Fangmom




       if (iargc() .lt. 1) then
          write(*,*) '*************************************************'
          write(*,*) '* proccessc2p2: get strains   *'
          write(*,*) '* Usage:                                        *'
          write(*,*) '*************************************************'
          write(*,*) '  proccessc2p2  <surface letter> <mass> <radius> <timeforh22> <hitwithr> <cleanhdot?>'
          write(*,*) ' cleanhdot>0 will clean hdot'
          write(*,*) '   E.g. processc2p2 B 2.7 140.d0 300 1 -1'
          write(*,*) '                *****************'
          stop
       end if
    
       name = sarg( 1,'c')
       massI= r8arg(2,1.d0)
       rext = r8arg(3,140.d0)
       time_min = r8arg(4,1.d0)
       hitwithr = r8arg(5,1.d0)
       cleanhdot = r8arg(6,1.d0)
       
!names for input files from clm (that must exist from harmonic)
!and outputfiles for their infinity versions

       file_psi(2,2) = 'psi4_c2p2'//name(:1)//'.dat'
       file_psi_inf(2,2)  = 'P4c2p2_inf'//name(:1)//'.dat'
       file_psi(2,1) = 'psi4_c2p1'//name(:1)//'.dat'
       file_psi_inf(2,1)  = 'P4c2p1_inf'//name(:1)//'.dat'
       file_psi(2,0) = 'psi4_c20'//name(:1)//'.dat'
       file_psi_inf(2,0)  = 'P4c20_inf'//name(:1)//'.dat'
       file_psi(2,-2) = 'psi4_c2m2'//name(:1)//'.dat'
       file_psi_inf(2,-2)  = 'P4c2m2_inf'//name(:1)//'.dat'
       file_psi(2,-1) = 'psi4_c2m1'//name(:1)//'.dat'
       file_psi_inf(2,-1)  = 'P4c2m1_inf'//name(:1)//'.dat'


       file_psi(3,3) = 'psi4_c3p3'//name(:1)//'.dat'
       file_psi_inf(3,3)  = 'P4c3p3_inf'//name(:1)//'.dat'
       file_psi(3,2) = 'psi4_c3p2'//name(:1)//'.dat'
       file_psi_inf(3,2)  = 'P4c3p2_inf'//name(:1)//'.dat'
       file_psi(3,1) = 'psi4_c3p1'//name(:1)//'.dat'
       file_psi_inf(3,1)  = 'P4c3p1_inf'//name(:1)//'.dat'
       file_psi(3,0) = 'psi4_c30'//name(:1)//'.dat'
       file_psi_inf(3,0)  = 'P4c30_inf'//name(:1)//'.dat'
       file_psi(3,-3) = 'psi4_c3m3'//name(:1)//'.dat'
       file_psi_inf(3,-3)  = 'P4c3m3_inf'//name(:1)//'.dat'
       file_psi(3,-2) = 'psi4_c3m2'//name(:1)//'.dat'
       file_psi_inf(3,-2)  = 'P4c3m2_inf'//name(:1)//'.dat'
       file_psi(3,-1) = 'psi4_c3m1'//name(:1)//'.dat'
       file_psi_inf(3,-1)  = 'P4c3m1_inf'//name(:1)//'.dat'


       file_psi(4,4) = 'psi4_c4p4'//name(:1)//'.dat'
       file_psi_inf(4,4)  = 'P4c4p4_inf'//name(:1)//'.dat'
       file_psi(4,3) = 'psi4_c4p3'//name(:1)//'.dat'
       file_psi_inf(4,3)  = 'P4c4p3_inf'//name(:1)//'.dat'
       file_psi(4,2) = 'psi4_c4p2'//name(:1)//'.dat'
       file_psi_inf(4,2)  = 'P4c4p2_inf'//name(:1)//'.dat'
       file_psi(4,1) = 'psi4_c4p1'//name(:1)//'.dat'
       file_psi_inf(4,1)  = 'P4c4p1_inf'//name(:1)//'.dat'
       file_psi(4,0) = 'psi4_c40'//name(:1)//'.dat'
       file_psi_inf(4,0)  = 'P4c40_inf'//name(:1)//'.dat'
       file_psi(4,-4) = 'psi4_c4m4'//name(:1)//'.dat'
       file_psi_inf(4,-4)  = 'P4c4m4_inf'//name(:1)//'.dat'
       file_psi(4,-3) = 'psi4_c4m3'//name(:1)//'.dat'
       file_psi_inf(4,-3)  = 'P4c4m3_inf'//name(:1)//'.dat'
       file_psi(4,-2) = 'psi4_c4m2'//name(:1)//'.dat'
       file_psi_inf(4,-2)  = 'P4c4m2_inf'//name(:1)//'.dat'
       file_psi(4,-1) = 'psi4_c4m1'//name(:1)//'.dat'
       file_psi_inf(4,-1)  = 'P4c4m1_inf'//name(:1)//'.dat'



!for output
       Fphase = 'phase'//name(:1)//'.dat'
       Fomega = 'omega'//name(:1)//'.dat'

       Fenergy = 'Energy'//name(:1)//'.dat'
       Fangmom = 'Angmom'//name(:1)//'.dat'

       Fampphase = 'amp_and_phase'//name(:1)//'.dat'
       Fampphase_inf = 'amp_and_phase_inf'//name(:1)//'.dat'

       Fhhampphase = 'HHamp_and_phase'//name(:1)//'.dat'
       Fhhampphase_inf = 'HHamp_and_phase_inf'//name(:1)//'.dat'

       Fhlm(2,2) = 'h_2p2'//name(:1)//'.dat'
       Fhlm_inf(2,2)  = 'h_2p2_inf'//name(:1)//'.dat'
       Fhlm(2,1) = 'h_2p1'//name(:1)//'.dat'
       Fhlm_inf(2,1)  = 'h_2p1_inf'//name(:1)//'.dat'
       Fhlm(2,0) = 'h_20'//name(:1)//'.dat'
       Fhlm_inf(2,0)  = 'h_20_inf'//name(:1)//'.dat'
       Fhlm(2,-2) = 'h_2m2'//name(:1)//'.dat'
       Fhlm_inf(2,-2)  = 'h_2m2_inf'//name(:1)//'.dat'
       Fhlm(2,-1) = 'h_2m1'//name(:1)//'.dat'
       Fhlm_inf(2,-1)  = 'h_2m1_inf'//name(:1)//'.dat'


       Fhlm(3,3) = 'h_3p3'//name(:1)//'.dat'
       Fhlm_inf(3,3)  = 'h_3p3_inf'//name(:1)//'.dat'
       Fhlm(3,2) = 'h_3p2'//name(:1)//'.dat'
       Fhlm_inf(3,2)  = 'h_3p2_inf'//name(:1)//'.dat'
       Fhlm(3,1) = 'h_3p1'//name(:1)//'.dat'
       Fhlm_inf(3,1)  = 'h_3p1_inf'//name(:1)//'.dat'
       Fhlm(3,0) = 'h_30'//name(:1)//'.dat'
       Fhlm_inf(3,0)  = 'h_30_inf'//name(:1)//'.dat'
       Fhlm(3,-3) = 'h_3m3'//name(:1)//'.dat'
       Fhlm_inf(3,-3)  = 'h_3m3_inf'//name(:1)//'.dat'
       Fhlm(3,-2) = 'h_3m2'//name(:1)//'.dat'
       Fhlm_inf(3,-2)  = 'h_3m2_inf'//name(:1)//'.dat'
       Fhlm(3,-1) = 'h_3m1'//name(:1)//'.dat'
       Fhlm_inf(3,-1)  = 'h_3m1_inf'//name(:1)//'.dat'


       Fhlm(4,4) = 'h_4p4'//name(:1)//'.dat'
       Fhlm_inf(4,4)  = 'h_4p4_inf'//name(:1)//'.dat'
       Fhlm(4,3) = 'h_4p3'//name(:1)//'.dat'
       Fhlm_inf(4,3)  = 'h_4p3_inf'//name(:1)//'.dat'
       Fhlm(4,2) = 'h_4p2'//name(:1)//'.dat'
       Fhlm_inf(4,2)  = 'h_4p2_inf'//name(:1)//'.dat'
       Fhlm(4,1) = 'h_4p1'//name(:1)//'.dat'
       Fhlm_inf(4,1)  = 'h_4p1_inf'//name(:1)//'.dat'
       Fhlm(4,0) = 'h_40'//name(:1)//'.dat'
       Fhlm_inf(4,0)  = 'h_40_inf'//name(:1)//'.dat'
       Fhlm(4,-4) = 'h_4m4'//name(:1)//'.dat'
       Fhlm_inf(4,-4)  = 'h_4m4_inf'//name(:1)//'.dat'
       Fhlm(4,-3) = 'h_4m3'//name(:1)//'.dat'
       Fhlm_inf(4,-3)  = 'h_4m3_inf'//name(:1)//'.dat'
       Fhlm(4,-2) = 'h_4m2'//name(:1)//'.dat'
       Fhlm_inf(4,-2)  = 'h_4m2_inf'//name(:1)//'.dat'
       Fhlm(4,-1) = 'h_4m1'//name(:1)//'.dat'
       Fhlm_inf(4,-1)  = 'h_4m1_inf'//name(:1)//'.dat'

!!! this is getting long and uggly... but i am lazzy

       Fdlm(2,2) = 'dth_2p2'//name(:1)//'.dat'
       Fdlm_inf(2,2)  = 'dth_2p2_inf'//name(:1)//'.dat'
       Fdlm(2,1) = 'dth_2p1'//name(:1)//'.dat'
       Fdlm_inf(2,1)  = 'dth_2p1_inf'//name(:1)//'.dat'
       Fdlm(2,0) = 'dth_20'//name(:1)//'.dat'
       Fdlm_inf(2,0)  = 'dth_20_inf'//name(:1)//'.dat'
       Fdlm(2,-2) = 'dth_2m2'//name(:1)//'.dat'
       Fdlm_inf(2,-2)  = 'dth_2m2_inf'//name(:1)//'.dat'
       Fdlm(2,-1) = 'dth_2m1'//name(:1)//'.dat'
       Fdlm_inf(2,-1)  = 'dth_2m1_inf'//name(:1)//'.dat'


       Fdlm(3,3) = 'dth_3p3'//name(:1)//'.dat'
       Fdlm_inf(3,3)  = 'dth_3p3_inf'//name(:1)//'.dat'
       Fdlm(3,2) = 'dth_3p2'//name(:1)//'.dat'
       Fdlm_inf(3,2)  = 'dth_3p2_inf'//name(:1)//'.dat'
       Fdlm(3,1) = 'dth_3p1'//name(:1)//'.dat'
       Fdlm_inf(3,1)  = 'dth_3p1_inf'//name(:1)//'.dat'
       Fdlm(3,0) = 'dth_30'//name(:1)//'.dat'
       Fdlm_inf(3,0)  = 'dth_30_inf'//name(:1)//'.dat'
       Fdlm(3,-3) = 'dth_3m3'//name(:1)//'.dat'
       Fdlm_inf(3,-3)  = 'dth_3m3_inf'//name(:1)//'.dat'
       Fdlm(3,-2) = 'dth_3m2'//name(:1)//'.dat'
       Fdlm_inf(3,-2)  = 'dth_3m2_inf'//name(:1)//'.dat'
       Fdlm(3,-1) = 'dth_3m1'//name(:1)//'.dat'
       Fdlm_inf(3,-1)  = 'dth_3m1_inf'//name(:1)//'.dat'


       Fdlm(4,4) = 'dth_4p4'//name(:1)//'.dat'
       Fdlm_inf(4,4)  = 'dth_4p4_inf'//name(:1)//'.dat'
       Fdlm(4,3) = 'dth_4p3'//name(:1)//'.dat'
       Fdlm_inf(4,3)  = 'dth_4p3_inf'//name(:1)//'.dat'
       Fdlm(4,2) = 'dth_4p2'//name(:1)//'.dat'
       Fdlm_inf(4,2)  = 'dth_4p2_inf'//name(:1)//'.dat'
       Fdlm(4,1) = 'dth_4p1'//name(:1)//'.dat'
       Fdlm_inf(4,1)  = 'dth_4p1_inf'//name(:1)//'.dat'
       Fdlm(4,0) = 'dth_40'//name(:1)//'.dat'
       Fdlm_inf(4,0)  = 'dth_40_inf'//name(:1)//'.dat'
       Fdlm(4,-4) = 'dth_4m4'//name(:1)//'.dat'
       Fdlm_inf(4,-4)  = 'dth_4m4_inf'//name(:1)//'.dat'
       Fdlm(4,-3) = 'dth_4m3'//name(:1)//'.dat'
       Fdlm_inf(4,-3)  = 'dth_4m3_inf'//name(:1)//'.dat'
       Fdlm(4,-2) = 'dth_4m2'//name(:1)//'.dat'
       Fdlm_inf(4,-2)  = 'dth_4m2_inf'//name(:1)//'.dat'
       Fdlm(4,-1) = 'dth_4m1'//name(:1)//'.dat'
       Fdlm_inf(4,-1)  = 'dth_4m1_inf'//name(:1)//'.dat'

!	do l=2,4
!	 do m=-l,l
!	  print*,l,m,file_psi(l,m)
!	  end do
!	end do
!	STOP

       if (.true.) then
          write(*,*) '   name = ',name
          write(*,*) 'Reading in the following files:'
          write(*,*) '   c2p2file = ',file_psi(2,2)
          write(*,*) '   c2p2_inffile  = ',file_psi_inf(2,2)
          write(*,*) '   massI       = ',massI
          write(*,*) '   rext        = ',rext
       end if

	if(givefile.lt.0) then
	print*,'read files dont assume names: real and im'
	read*,c2p2file
	read*,c2p2_inffile
	print*, 'this code no longer accepts this'
	STOP
	end if

!-----------------------------------------------
!-----------------------------------------------
       output = 0        !-- 0, only .dat, 1 also output on the screen
!       massI  = 1.0      !-- the initial (m1 + m2) of the isolated bodies
!       rext   = 140.0     !-- the radius of the surface extraction
!       time_min = 300.0    !-- time to start the integration of Ds and Es
       factor  = 1.0     ! arbitrary normalization factor to match other results
      
       N      = Ndim
       pi     = dacos(-1.0d0)
       h      = pi/N
       h2     = h*h
!------------------------------------------------------



      do l=2, 4
       do m = -l, l
         open(unit=200+10*l+m, file=file_psi(l,m))
         open(unit=400+10*l+m, file=file_psi_inf(l,m))
         open(unit=500+10*l+m, file=Fhlm(l,m))
         open(unit=600+10*l+m, file=Fhlm_inf(l,m))
         open(unit=700+10*l+m, file=Fdlm(l,m))
         open(unit=800+10*l+m, file=Fdlm_inf(l,m))
       end do
      end do

 
      open(unit=122, file=Fphase)            
                 
      open(unit=130, file=Fomega)                        
      open(unit=131, file=Fampphase)
      open(unit=137, file=Fampphase_inf)                                                

      open(unit=138, file=Fhhampphase)
      open(unit=139, file=Fhhampphase_inf)         

      open(unit=132, file=Fenergy)
      open(unit=133, file=Fangmom)          


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j =0
	 do i=1, 500000
	 j = j+1
          do l=2,4
           do m= -l, l
	    read (200+10*l+m,*,end=10) tortutime(j), clm(l,m,1,j),clm(l,m,2,j), time(j)
	   end do
	  end do
	  if(j.gt.1.and.time(j).eq.time(j-1)) j=j-1
	end do
10 	continue
	 imax = j - 1

	print*, imax, tortutime(j), clm(2,2,1,j), clm(2,2,2,j), time(j)
!	STOP


        cphase = 0
	do itime = 1, imax
        smallnumber =1e-10
        x = clm(2,2,1,itime)
        y = clm(2,2,2,itime)
        if (x .gt. smallnumber.and.y.ge.0.) then
          phase(itime) = atan(y/x) 
        else if ((x .lt. -smallnumber) .AND. (y .ge. 0.0)) then 
          phase(itime) = atan(y/x)+pi
        else if ((x .lt. -smallnumber) .AND. (y .le. 0.0)) then 
          phase(itime) = atan(y/x)+pi
        else if ((x .gt. smallnumber) .AND. (y .le. 0.0)) then 
          phase(itime) = atan(y/x)+2.*pi
        end if

	if(itime.gt.2.and.x.gt.smallnumber) then
	  if (clm(2,2,2,itime-1).gt.0.and.clm(2,2,2,itime).lt.0.) then
	     print*,phase(itime-1),phase(itime)
             cphase=cphase-1

	  else if(clm(2,2,2,itime-1).lt.0.and.clm(2,2,2,itime).gt.0.) then
	     print*,phase(itime-1),phase(itime)
             cphase=cphase+1

	  end if
	end if

	phase(itime) = phase(itime) + (cphase)*2.*pi

	write(122,*), time(itime), phase(itime)




        !----output on the screen ---------------------------
        if (output .EQ. 1) then
  	  print*, "-------------------------"      !		
	  print*, "-------------------------"      
          print*, "-----L=2 MODES ----------"
  	  print*, "-------------------------"
	  print*, "c2p2=",clm(2,2,1,itime),"c2p2=",clm(2,2,2,itime)
        end if	


      end do

	
      !-computing the frequency from the dominant mode l=m=2
      !-either by omega=dphase/dt
      ! or omega = -1/m*Im [ (d C_{lm}/dt)/ C_{lm}] 
      ! (a+i*b)/(c+i*d) = (a*c + b*d + i*(b*c - a*d))/(c^2 + d^2)


      !-nakano extrapolation to infinity
      ! first define some constant
      ra        = rext * ( 1.0 + massI/(2.0*rext) )**2
      rr        = rext * ( 1.0 + massI/(2.0*rext) )**2
      C0        = (1.0 - 2.0*massI/ra)
      factorl2  = (1.0*4.0)/(2.0*ra) 
      factorl3  = (2.0*5.0)/(2.0*ra) 
      factorl4  = (3.0*6.0)/(2.0*ra) 

      ISTART = 2

      print*, "ra=",ra,"factorl2", factorl2


!to define time all the way till the end
!notice itime was dangerously used first for a varying index
!and from now on as a maxed index. So i will make it more explicit
	imax = itime - 1

	i = imax
        localtime = time(i)
        time_tort = tortutime(i)

      do i=ISTART, imax-1
	dt = time(i)-time(i-1)
        if (time(i)+3*dt .GT. time_min) then 
	  do l=2,4
	   do m = -l,l
	    integralclm(l,m,:) = integralclm(l,m,:) + clm(l,m,:,i)*dt
	   end do
	  end do
        end if

!this is hard coded to give c2p2 at infinity using nakano's extrapolation
!notice however that what was there before, was getting c2p2_inf but using
!c20 in the correction. I fixed it 
	do m=-2,2
         clm_inf(2,m,:,i) = C0 * (clm(2,m,:,i) - factorl2*integralclm(2,m,:))
	end do 

	do m=-3,3
         clm_inf(3,m,:,i) = C0 * (clm(3,m,:,i) - factorl3*integralclm(3,m,:))
	end do 

	do m=-4,4
         clm_inf(4,m,:,i) = C0 * (clm(4,m,:,i) - factorl4*integralclm(4,m,:))
	end do 

	do l = 2, 4
	 do m = -l, l
           write(400+10*l+m,100) tortutime(i), clm_inf(l,m,1,i), clm_inf(l,m,2,i), time(i)
	 end do
	end do 

         omega1(i) = 0.5*( phase(i+1)-phase(i-1) )/(2.0*dt)
         cR2p2   = clm(2,2,1,i)
         cI2p2   = clm(2,2,2,i)

	if(i.gt.1.and.i.lt.imax) then
         dtcR2p2 = (clm(2,2,1,i+1) - clm(2,2,1,i-1))/(2.0*dt)
         dtcI2p2 = (clm(2,2,2,i+1) - clm(2,2,2,i-1))/(2.0*dt)
         temp1 = (dtcI2p2*cR2p2 - dtcR2p2*cI2p2)/(cR2p2**2 + cI2p2**2)    
         omega2(i) = abs((1.0d0/2.0d0)*temp1)
         write(130,100) tortutime(i), omega2(i),omega1(i), time(i)
	end if

	 write(131,100) tortutime(i), &
     &                  sqrt( cR2p2**2 + cI2p2**2 ), &
     &                             phase(i), time(i)

      end do 

!perhaps c_extrapolated picks a drift, see what happens if we take it off
!!!!!!!!!!!!!!!!!!!!! the infinity part is picking a linear drift, 
! let us remove to see  if we can fix it somehow
! we will reuse variables defined for the h2p2 part


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hlm = 0.0
      hlm_inf = 0.0
      hlm_ave = 0.0
      hlm_inf_ave = 0.0
      xyhlm_ave = 0.0
      xyhlm_inf_ave = 0.0
      count_ave = 0
      slopehlm_inf = 0.0
      interhlm_inf = 0.0

	do l=2,4
	  do m=-l,l
      do i=2, imax-1
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+3.*dt .GT. time_min) then


	hlm_inf_ave(l,m,:) = hlm_inf_ave(l,m,:) + clm_inf(l,m,:,i)

	xyhlm_inf_ave(l,m,:) = xyhlm_inf_ave(l,m,:) + &
     &                              clm_inf(l,m,:,i)*time(i)

	time_ave = time_ave + time(i)
	xxtime_ave = xxtime_ave + time(i)*time(i)
	count_ave = count_ave + 1

	else
	time_ave =  0
	xxtime_ave = 0
	count_ave = 0
	end if
	  end do
	 end do
	end do


	do l=2,4
	  do m=-l,l
           do i=2, imax-1
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+3.*dt .GT. time_min) then

	  slopehlm_inf(l,m,:) =  (xyhlm_inf_ave(l,m,:) - hlm_inf_ave(l,m,:) &
     &                        *time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)


	  interhlm_inf(l,m,:) = hlm_inf_ave(l,m,:)/count_ave - &
     &                         slopehlm_inf(l,m,:)*time_ave/count_ave

!what if i just get the average out?
	interhlm_inf(l,m,:) = hlm_inf_ave(l,m,:)/count_ave

	  clm_inf(l,m,:,i) =  clm_inf(l,m,:,i) - ( slopehlm_inf(l,m,:)*time(i) *0 + interhlm_inf(l,m,:))

!	write(1000+l*10+m,*) tortutime(i),clm_inf(l,m,1,i),clm_inf(l,m,2,i)

	end if

	  end do
	  end do
	end do

!STOP

!!!!!!! OK the above did a bunch of things wrt c22, but we have
!! the exptrapolated to infinity wavefor, so let us redo with that
!! as well... this file is getting big and uggly... but i saw beaty and the beast...

	cphase = 0

	do i = ISTART, imax - 1
        localtime = time(i)
        x = clm_inf(2,2,1,i)
        y = clm_inf(2,2,2,i)
        if (x .gt. smallnumber.and.y.ge.0.) then
          phase_inf(i) = atan(y/x) 
        else if ((x .lt. -smallnumber) .AND. (y .ge. 0.0)) then 
          phase_inf(i) = atan(y/x)+pi
        else if ((x .lt. -smallnumber) .AND. (y .le. 0.0)) then 
          phase_inf(i) = atan(y/x)+pi
        else if ((x .gt. smallnumber) .AND. (y .le. 0.0)) then 
          phase_inf(i) = atan(y/x)+2.*pi
        end if

	if(i.gt.2.and.x.gt.smallnumber) then
	  if (clm_inf(2,2,2,i-1).gt.0.and.clm_inf(2,2,2,i).lt.0.) then
	     print*,phase_inf(i-1),phase_inf(i)
             cphase=cphase-1

	  else if(clm_inf(2,2,2,i-1).lt.0.and.clm_inf(2,2,2,i).gt.0.) then
	     print*,phase_inf(i-1),phase_inf(i)
             cphase=cphase+1

	  end if
	end if
	phase_inf(i) = phase_inf(i) + (cphase)*2.*pi

	 write(137,100) tortutime(i), &
     &                  sqrt( x**2 + y**2 ), &
     &                             phase_inf(i), localtime

	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!

      ! integrate twice in time to get the strain
!first get the first integral, which for some odd reason
!will be named d2tclm

      d2tclm = 0.0
      d2tclm_inf = 0.0

      ! first time integral 
      do i=2, imax-1
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+2*dt .GT. time_min) then
        
	do l=2,4
	 do m = -l, l
          d2tclm(l,m,:,i) = d2tclm(l,m,:,i-1) &
     &               + 0.5*dt*(clm(l,m,:,i-1) + clm(l,m,:,i))
          d2tclm_inf(l,m,:,i) = d2tclm_inf(l,m,:,i-1) &
     &               + 0.5*dt*(clm_inf(l,m,:,i-1) + clm_inf(l,m,:,i))
	
!	write(700+10*l+m,100) tortutime(i),d2tclm(l,m,1,i),    d2tclm(l,m,2,i),    time(i)
!	write(800+10*l+m,100) tortutime(i),d2tclm_inf(l,m,1,i),d2tclm_inf(l,m,2,i),time(i)

	  end do
	end do

        end if
      end do 




!!!!!!!!!!!!!!!!!!!!! the infinity part is picking a drift, 
! let us remove to see  if we can fix it somehow
! we will reuse variables defined for the h2p2 part


      hlm = 0.0
      hlm_inf = 0.0
      hlm_ave = 0.0
      hlm_inf_ave = 0.0
      xyhlm_ave = 0.0
      xyhlm_inf_ave = 0.0
      count_ave = 0
      s11=0.;s12=0;s22=0.;sy1=0.;sy2=0.;x1ave=0.;x2ave=0.;yave=0.
      s11I=0.;s12I=0;s22I=0.;sy1I=0.;sy2I=0.;x1ave=0.;x2ave=0.;yaveI=0.

     do l = 2, 4
      do m = -l, l

       do i=2, imax-1
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation

        if (localtime+2*dt .GT. time_min) then

	s11=s11+time(i)**2
	s12=s12+time(i)**3
	s22=s22+time(i)**4
	x1ave=x1ave+time(i)
	x2ave=x2ave+time(i)**2

!for both
        yave(l,m,:) =yave(l,m,:) +d2tclm(l,m,:,i)
        yaveI(l,m,:)=yaveI(l,m,:)+d2tclm_inf(l,m,:,i)

	sy1(l,m,:) =sy1(l,m,:) +time(i)*d2tclm(l,m,:,i)
	sy1I(l,m,:)=sy1I(l,m,:)+time(i)*d2tclm_inf(l,m,:,i)

	sy2(l,m,:) =sy2(l,m,:) +time(i)**2*d2tclm(l,m,:,i)
	sy2I(l,m,:)=sy2I(l,m,:)+time(i)**2*d2tclm_inf(l,m,:,i)

	time_ave = time_ave + time(i)
	xxtime_ave = xxtime_ave + time(i)*time(i)
	count_ave = count_ave + 1
	else
	count_ave = 0

	s11=0.
	s12=0.
	s22=0.
	s11I=0.
	s12I=0.
	s22I=0.

	x1ave=0.
	x2ave=0.
	end if

      end do 
!
	s11 =s11 -x1ave**2/count_ave
	s11I=s11I-x1ave**2/count_ave

	s12 =s12 -(x1ave*x2ave)/count_ave
	s12I=s12I-(x1ave*x2ave)/count_ave

	s22 =s22 -x2ave**2/count_ave
	s22I=s22I-x2ave**2/count_ave

	sy1(l,m,:) =sy1(l,m,:) -yave(l,m,:) *x1ave/count_ave
	sy1I(l,m,:)=sy1I(l,m,:)-yaveI(l,m,:)*x1ave/count_ave

	sy2(l,m,:) =sy2(l,m,:) -yave(l,m,:) *x2ave/count_ave
	sy2I(l,m,:)=sy2I(l,m,:)-yaveI(l,m,:)*x2ave/count_ave


	bet2(l,m,:) =(sy1(l,m,:)*s22 -sy2(l,m,:) *s12)/(s22*s11-s12**2)
	bet2I(l,m,:)=(sy1I(l,m,:)*s22-sy2I(l,m,:)*s12)/(s22*s11-s12**2)

	bet3(l,m,:) =(sy2(l,m,:) *s11-sy1(l,m,:) *s12)/(s22*s11-s12**2)
	bet3I(l,m,:)=(sy2I(l,m,:)*s11-sy1I(l,m,:)*s12)/(s22*s11-s12**2)

!try removing something but not all
	bet2(l,m,:) =0.0;bet3(l,m,:) =0.0
!	bet2I(l,m,:) =0.0;
        bet3I(l,m,:) =0.0

        bet1(l,m,:) =(yave(l,m,:) -bet2(l,m,:)*x1ave- bet3(l,m,:) *x2ave)/count_ave
        bet1I(l,m,:)=(yaveI(l,m,:)-bet2I(l,m,:)*x1ave-bet3I(l,m,:)*x2ave)/count_ave



!remove drift
	do i=2, imax
	 if(time(i)+2*dt.gt.time_min) then


!	d2tclm(l,m,:,i)     = d2tclm(l,m,:,i)     - (bet1(l,m,:) +bet2(l,m,:)*time(i) +  0*bet3(l,m,:)*time(i)**2)
!	d2tclm_inf(l,m,:,i) = d2tclm_inf(l,m,:,i) - (bet1I(l,m,:)+bet2I(l,m,:)*time(i)+  bet3I(l,m,:)*time(i)**2)


	bet2I(l,m,:) =  (sy1I(l,m,:))/(x2ave-x1ave**2/count_ave)
	bet1I(l,m,:) = yaveI(l,m,:)/count_ave - bet2I(l,m,:)*x1ave/count_ave

	bet2(l,m,:) = 0.! (sy1I(l,m,:))/(x2ave-x1ave**2/count_ave)
	bet1(l,m,:) = yave(l,m,:)/count_ave - bet2(l,m,:)*x1ave/count_ave

!get rid of average offset
	d2tclm_inf(l,m,:,i)     = d2tclm_inf(l,m,:,i) - (bet1I(l,m,:)+bet2I(l,m,:)*time(i) )
	d2tclm(l,m,:,i)     = d2tclm(l,m,:,i) - (bet1(l,m,:)+bet2(l,m,:)*time(i) )

	write(700+10*l+m,100) tortutime(i),d2tclm(l,m,1,i),    d2tclm(l,m,2,i),    time(i)
	write(800+10*l+m,100) tortutime(i),d2tclm_inf(l,m,1,i),d2tclm_inf(l,m,2,i),time(i)


	 end if
	end do

         end do
        end do 
	print*,'slope dth', bet2I(2,2,1),bet1I(2,2,1)


!	print*,count_ave,time_ave
!	STOP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hlm = 0.0
      hlm_inf = 0.0
      count_ave = 0
      s11=0.;s12=0;s22=0.;sy1=0.;sy2=0.;x1ave=0.;x2ave=0.;yave=0.
      s11I=0.;s12I=0;s22I=0.;sy1I=0.;sy2I=0.;x1ave=0.;x2ave=0.;yaveI=0.

     do l = 2, 4
      do m = -l, l

      do i=2, imax-1
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime .GT. time_min) then
          hlm(l,m,:,i) = hlm(l,m,:,i-1) &
     &                   + 0.5*dt*(d2tclm(l,m,:,i-1) + d2tclm(l,m,:,i))
          hlm_inf(l,m,:,i) = hlm_inf(l,m,:,i-1) &
     &                   + 0.5*dt*(d2tclm_inf(l,m,:,i-1) + d2tclm_inf(l,m,:,i))
!        end if

	s11=s11+time(i)**2
	s12=s12+time(i)**3
	s22=s22+time(i)**4
	x1ave=x1ave+time(i)
	x2ave=x2ave+time(i)**2

!for both
        yave(l,m,:) =yave(l,m,:) +hlm(l,m,:,i)
        yaveI(l,m,:)=yaveI(l,m,:)+hlm_inf(l,m,:,i)

	sy1(l,m,:) =sy1(l,m,:) +time(i)*hlm(l,m,:,i)
	sy1I(l,m,:)=sy1I(l,m,:)+time(i)*hlm_inf(l,m,:,i)

	sy2(l,m,:) =sy2(l,m,:) +time(i)**2*hlm(l,m,:,i)
	sy2I(l,m,:)=sy2I(l,m,:)+time(i)**2*hlm_inf(l,m,:,i)

	count_ave = count_ave + 1
	else
	count_ave = 0

	s11=0.;s12=0.;s22=0.
	s11I=0.;s12I=0.;s22I=0.
	x1ave=0.
	x2ave=0.


	end if


      end do 

	s11 =s11 -x1ave**2/count_ave
	s11I=s11I-x1ave**2/count_ave

	s12 =s12 -(x1ave*x2ave)/count_ave
	s12I=s12I-(x1ave*x2ave)/count_ave

	s22 =s22 -x2ave**2/count_ave
	s22I=s22I-x2ave**2/count_ave

	sy1(l,m,:) =sy1(l,m,:) -yave(l,m,:) *x1ave/count_ave
	sy1I(l,m,:)=sy1I(l,m,:)-yaveI(l,m,:)*x1ave/count_ave

	sy2(l,m,:) =sy2(l,m,:) -yave(l,m,:) *x2ave/count_ave
	sy2I(l,m,:)=sy2I(l,m,:)-yaveI(l,m,:)*x2ave/count_ave


	bet2(l,m,:) =(sy1(l,m,:)*s22 -sy2(l,m,:) *s12)/(s22*s11-s12**2)
	bet2I(l,m,:)=(sy1I(l,m,:)*s22-sy2I(l,m,:)*s12)/(s22*s11-s12**2)

	bet3(l,m,:) =(sy2(l,m,:) *s11-sy1(l,m,:) *s12)/(s22*s11-s12**2)
	bet3I(l,m,:)=(sy2I(l,m,:)*s11-sy1I(l,m,:)*s12)/(s22*s11-s12**2)


        bet1(l,m,:) =(yave(l,m,:) -bet2(l,m,:)*x1ave- bet3(l,m,:) *x2ave)/count_ave
        bet1I(l,m,:)=(yaveI(l,m,:)-bet2I(l,m,:)*x1ave-bet3I(l,m,:)*x2ave)/count_ave



!remove drift
	do i=2, imax-1
	 if(time(i).gt.time_min) then

!perhaps only linear for h

	bet2(l,m,:) = (sy1(l,m,:))/(x2ave-x1ave**2/count_ave)
	bet1(l,m,:) = yave(l,m,:)/count_ave - bet2(l,m,:)*x1ave/count_ave


	hlm(l,m,:,i)     = hlm(l,m,:,i)     - (bet1(l,m,:) +bet2(l,m,:)*time(i) +0*bet3(l,m,:)*time(i)**2)
	hlm_inf(l,m,:,i) = hlm_inf(l,m,:,i) - (bet1I(l,m,:)+bet2I(l,m,:)*time(i)+bet3I(l,m,:)*time(i)**2)

!	hlm(l,m,:,i)     = hlm(l,m,:,i)     - (yave(l,m,:)/count_ave)
!	hlm_inf(l,m,:,i) = hlm_inf(l,m,:,i) - (yaveI(l,m,:)/count_ave)


        write(500+10*l+m,100) tortutime(i), hlm(l,m,1,i), hlm(l,m,2,i), time(i)
        write(600+10*l+m,100) tortutime(i), hlm_inf(l,m,1,i) , &
     &                       hlm_inf(l,m,2,i), time(i)

	 end if
	end do

	end do
	end do

	print*,'facts',bet1I(2,2,1),bet2I(2,2,1),bet3I(2,2,1)
	print*,'factsBAS',bet1(2,2,1),bet2(2,2,1),bet3(2,2,1)


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!finally with all the above, can compute the energy/ang mom radiated... supposedly

     ! compute the GW luminosity and total energy radiated in GW
      ! first define some constant
      factor    = 1.0/(16.0*pi) 

      print*, "ra=",ra,"factor", factor 
      LGW = 0.
      LGW_inf = 0.
      JGW = 0.
      JGW_inf = 0.
          
! for weird reasons the first derivative is in d2tclm
      do i=2, imax-1
      LGW(i) = 0.
      LGW_inf(i) = 0.
      JGW(i) = 0.
      JGW_inf(i) = 0.
        !dt         = 0.5*(time(i+1)-time(i-1))*factor
        dt         = (time(i)-time(i-1))*factor
        ! with the raw data and the infinity extrapolation
        if (time(i) .GT. time_min) then
           do l=2,4
              do m = -l ,l
                 if(m.ne.0) then
                    LGW(i)     = LGW(i)     + (d2tclm(l,m,1,i)**2    +d2tclm(l,m,2,i)**2    )*dt
                    LGW_inf(i) = LGW_inf(i) + (d2tclm_inf(l,m,1,i)**2+d2tclm_inf(l,m,2,i)**2)*dt
                    JGW(i)     = JGW(i)     + dt* m* (-d2tclm(l,m,1,i)*hlm(l,m,2,i) &
     &                                                  +d2tclm(l,m,2,i)*hlm(l,m,1,i) )
                    JGW_inf(i) = JGW_inf(i) + dt* m* (-d2tclm_inf(l,m,1,i)*hlm_inf(l,m,2,i) &
     &                                                  +d2tclm_inf(l,m,2,i)*hlm_inf(l,m,1,i))
                 end if
             end do
            end do
!now do the cumulative sum
	            LGW(i) =     LGW(i-1)     + LGW(i)
	            LGW_inf(i) = LGW_inf(i-1) + LGW_inf(i)
	            JGW(i) =     JGW(i-1)     + JGW(i)
	            JGW_inf(i) = JGW_inf(i-1) + JGW_inf(i)

            !
            write(132,100) tortutime(i), LGW(i), LGW_inf(i), time(i)
            write(133,100) tortutime(i), JGW(i), JGW_inf(i), time(i)
        end if
      end do 


!!!!!! finally we note the phase and amplitudes above were obtained
!!! from c2p2, but that isnot the strain. so let us redo here with strains

!!! first from h2p2

	cphase = 0

	do i = 2, imax - 1
        localtime = time(i)
        x = hlm(2,2,1,i)
        y = hlm(2,2,2,i)

        if (localtime .GT. time_min) then

        if (x .gt. smallnumber.and.y.ge.0.) then
          HHphase(i) = atan(y/x) 

        else if ((x .lt. -smallnumber) .AND. (y .ge. 0.0)) then 
          HHphase(i) = atan(y/x)+pi
        else if ((x .lt. -smallnumber) .AND. (y .le. 0.0)) then 
          HHphase(i) = atan(y/x)+pi
        else if ((x .gt. smallnumber) .AND. (y .le. 0.0)) then 
          HHphase(i) = atan(y/x)+2.*pi
        end if

	if(i.gt.2.and.x.gt.smallnumber) then
	  if (hlm(2,2,2,i-1).gt.0.and.hlm(2,2,2,i).lt.0.) then
	     print*,HHphase(i-1),HHphase(i)
             cphase=cphase-1

	  else if(hlm(2,2,2,i-1).lt.0.and.hlm(2,2,2,i).gt.0.) then
	     print*,HHphase(i-1),HHphase(i)
             cphase=cphase+1

	  end if
	end if
	HHphase(i) = HHphase(i) + (cphase)*2.*pi

	 write(138,100) tortutime(i), &
     &                  sqrt( x**2 + y**2 ), &
     &                             HHphase(i), time(i)

	end if
	end do

!!!! then with h2p2_inf
	cphase = 0

	do i = ISTART, imax - 1
        localtime = time(i)
        x = hlm_inf(2,2,1,i)
        y = hlm_inf(2,2,2,i)
        if (localtime .GT. time_min) then

        if (x .gt. smallnumber.and.y.ge.0.) then
          HHphase_inf(i) = atan(y/x) 
        else if ((x .lt. -smallnumber) .AND. (y .ge. 0.0)) then 
          HHphase_inf(i) = atan(y/x)+pi
        else if ((x .lt. -smallnumber) .AND. (y .le. 0.0)) then 
          HHphase_inf(i) = atan(y/x)+pi
        else if ((x .gt. smallnumber) .AND. (y .le. 0.0)) then 
          HHphase_inf(i) = atan(y/x)+2.*pi
        end if

	if(i.gt.2.and.x.gt.smallnumber) then
	  if (hlm_inf(2,2,2,i-1).gt.0.and.hlm_inf(2,2,2,i).lt.0.) then
	     print*,HHphase_inf(i-1),HHphase_inf(i)
             cphase=cphase-1

	  else if(hlm_inf(2,2,2,i-1).lt.0.and.hlm_inf(2,2,2,i).gt.0.) then
	     print*,HHphase_inf(i-1),HHphase_inf(i)
             cphase=cphase+1

	  end if
	end if
	HHphase_inf(i) = HHphase_inf(i) + (cphase)*2.*pi

	 write(139,100) tortutime(i), &
     &                  sqrt( x**2 + y**2 ), &
     &                             HHphase_inf(i), time(i)

	end if
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      close(unit=114)
      close(unit=122)

      close(unit=130)
      close(unit=131)                                                              
      close(unit=137)   
      close(unit=138)                                              
      close(unit=139)                                              
      close(unit=133)                                              
      close(unit=134) 
      close(unit=135)
      close(unit=136)                                              
      close(unit=314)

	do l=2,4
	 do m=-l,l
         close(unit=200+10*l+m)
         close(unit=400+10*l+m)	
         close(unit=500+10*l+m)
         close(unit=600+10*l+m)
         close(unit=700+10*l+m)
         close(unit=800+10*l+m)
	 end do
	end do 


!-----------------------------------------
      
  100 format(4e16.7)      

      end program

