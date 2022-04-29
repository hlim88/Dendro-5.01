     program      main
      implicit     none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This program takes c2p2 and c2p2_inf of an already
! recombined (from desired batches) run with harmonic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
            
      integer :: i, j, k, N, itheta, iphi, output, initial, itime
      integer :: correction, integration, ISTART, imax
      real*8  :: pi, theta, phi, integral, h, h2, x, y, localtime
      real*8  :: rp, ip, norm, dt, factor, temp1, time_tort,integralc2p2(2)
      real*8  :: massI, rext, test, smallnumber, ra, rr, fpR, fpI
      real*8  :: surface_integral,surface_integral2, givefile
      integer :: cphase
            
! ------- this parameter changes depending on the size of u------------
      integer Ndim
!      parameter ( Ndim = 80 )


      character*20    file_psi
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

 
      real*8  :: cI2p2,cR2p2,dtcR2p2,dtCI2p2
      real*8  :: C0, factorl2,factorl3,factorl4
      real*8  :: d2p2(2),d2p1(2),d20(2),d2m1(2),d2m2(2)
      real*8  :: d2p2_inf(2),d2p1_inf(2),d20_inf(2),d2m1_inf(2)
      real*8  :: d2m2_inf(2)

      integer     MAXMEM
!      parameter ( MAXMEM = 1000000 )
      parameter ( MAXMEM = 50000 )
  

      integer     MAXTIME
      parameter ( MAXTIME = 50000 )

      real*8  :: c2p2(2,MAXTIME)
      real*8  :: c2p2_inf(2,MAXTIME)

      real*8  :: omega1(MAXTIME),omega2(MAXTIME),LGW(MAXTIME),LGW_inf(MAXTIME)
      real*8  :: d2p2T(2,MAXTIME),h2p2(2,MAXTIME), time(MAXTIME)
      real*8  :: d2p2T_inf(2,MAXTIME),h2p2_inf(2,MAXTIME), tortutime(MAXTIME)
      real*8  :: phase_inf(MAXTIME), phase(MAXTIME), amp_inf
      real*8  :: HHphase_inf(MAXTIME), HHamp_inf, HHphase(MAXTIME), HHamp
      real*8  :: h22_ave(2), h22_inf_ave(2), time_ave
      real*8  :: xxtime_ave, xyh22_ave(2), xyh22_inf_ave(2)
      real*8 ::  count_ave
      real*8  :: slopeh22(2),interh22(2),slopeh22_inf(2),interh22_inf(2)
      
!for output with letter so as to keep the file names
      character*32 :: Fh22,Fh22_inf,Fd22,Fd22_inf,Fomega,Fphase, &
     &                Fampphase,Fampphase_inf,Fhhampphase,Fhhampphase_inf




       if (iargc() .lt. 1) then
          write(*,*) '*************************************************'
          write(*,*) '* proccessc2p2: get strains   *'
          write(*,*) '* Usage:                                        *'
          write(*,*) '*************************************************'
          write(*,*) '  proccessc2p2  <surface letter> <mass> <radius> <timeforh22> <givefiles?>'
          write(*,*) ' givefiles>0 assume names as in had, <1 give them below'
          write(*,*) '   E.g. processc2p2 B 2.7 140.d0 300'
          write(*,*) '                *****************'
          stop
       end if
    
       name = sarg( 1,'c')
       massI= r8arg(2,1.d0)
       rext = r8arg(3,140.d0)
       time_min = r8arg(4,1.d0)
       givefile = r8arg(5,1.d0)
       
       c2p2file = 'psi4_c2p2'//name(:1)//'.dat'
       c2p2_inffile  = 'P4c2p2_inf'//name(:1)//'.dat'

!for output
       Fphase = 'phase'//name(:1)//'.dat'
       Fomega = 'omega'//name(:1)//'.dat'
       Fh22 = 'h2p2'//name(:1)//'.dat'
       Fh22_inf = 'h2p2_inf'//name(:1)//'.dat'
       Fd22 = 'd2p2'//name(:1)//'.dat'
       Fd22_inf = 'd2p2_inf'//name(:1)//'.dat'

       Fampphase = 'amp_and_phase'//name(:1)//'.dat'
       Fampphase_inf = 'amp_and_phase_inf'//name(:1)//'.dat'

       Fhhampphase = 'HHamp_and_phase'//name(:1)//'.dat'
       Fhhampphase_inf = 'HHamp_and_phase_inf'//name(:1)//'.dat'





       if (.true.) then
          write(*,*) '   name = ',name
          write(*,*) 'Reading in the following files:'
          write(*,*) '   c2p2file = ',c2p2file
          write(*,*) '   c2p2_inffile  = ',c2p2_inffile
          write(*,*) '   massI       = ',massI
          write(*,*) '   rext        = ',rext
       end if

	if(givefile.lt.0) then
	print*,'read files dont assume names: real and im'
	read*,c2p2file
	read*,c2p2_inffile
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

      open(unit=114, file=c2p2file)
      open(unit=314, file=c2p2_inffile)
 
      open(unit=122, file=Fphase)            
                 
      open(unit=130, file=Fomega)                        
      open(unit=131, file=Fampphase)
      open(unit=137, file=Fampphase_inf)                                                

      open(unit=138, file=Fhhampphase)
      open(unit=139, file=Fhhampphase_inf)                   

      open(unit=133, file=Fd22)
      open(unit=134, file=Fd22_inf)

      open(unit=135, file=Fh22)
      open(unit=136, file=Fh22_inf)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j =0
	 do i=1, 500000
	 j = j+1
	 read (114,*,end=10) tortutime(j), c2p2(1,j), c2p2(2,j), time(j)
	  if(j.gt.1.and.time(j).eq.time(j-1)) j=j-1
	end do
10 	continue
	 imax = j - 1

	print*, imax, tortutime(j), c2p2(1,j), c2p2(2,j), time(j)
!	STOP


        cphase = 0
	do itime = 1, imax
        smallnumber =1e-10
        x = c2p2(1,itime)
        y = c2p2(2,itime)
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
	  if (c2p2(2,itime-1).gt.0.and.c2p2(2,itime).lt.0.) then
	     print*,phase(itime-1),phase(itime)
             cphase=cphase-1

	  else if(c2p2(2,itime-1).lt.0.and.c2p2(2,itime).gt.0.) then
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
	  print*, "c2p2=",c2p2(1,itime),"c2p2=",c2p2(2,itime)
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
!to define time all the way till the end
!notice itime was dangerously used first for a varying index
!and from now on as a maxed index. So i will make it more explicit
	imax = itime - 1

	i = imax
        localtime = time(i)
        time_tort = tortutime(i)

      do i=ISTART, imax-1
	dt = time(i)-time(i-1)
        if (time(i) .GT. time_min) then 
          integralc2p2 = integralc2p2 + c2p2(:,i)*dt 
        end if

!this is hard coded to give c2p2 at infinity using nakano's extrapolation
!notice however that what was there before, was getting c2p2_inf but using
!c20 in the correction. I fixed it 

        c2p2_inf(:,i) = C0 * (c2p2(:,i) - factorl2*integralc2p2)

        write(314,100) tortutime(i), c2p2_inf(1,i), c2p2_inf(2,i), time(i)

         omega1(i) = 0.5*( phase(i+1)-phase(i-1) )/(2.0*dt)
         cR2p2   = c2p2(1,i)
         cI2p2   = c2p2(2,i)

	if(i.gt.1.and.i.lt.imax) then
         dtcR2p2 = (c2p2(1,i+1) - c2p2(1,i-1))/(2.0*dt)
         dtcI2p2 = (c2p2(2,i+1) - c2p2(2,i-1))/(2.0*dt)
         temp1 = (dtcI2p2*cR2p2 - dtcR2p2*cI2p2)/(cR2p2**2 + cI2p2**2)    
         omega2(i) = abs((1.0d0/2.0d0)*temp1)
         write(130,100) tortutime(i), omega2(i),omega1(i), time(i)
	end if

	 write(131,100) tortutime(i), &
     &                  sqrt( c2p2(1,i)**2 + c2p2(2,i)**2 ), &
     &                             phase(i), time(i)

      end do 


!!!!!!! OK the above did a bunch of things wrt c22, but we have
!! the exptrapolated to infinity wavefor, so let us redo with that
!! as well... this file is getting big and uggly... but i saw beaty and the beast...

	cphase = 0

	do i = ISTART, imax - 1
        localtime = time(i)
        x = c2p2_inf(1,i)
        y = c2p2_inf(2,i)
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
	  if (c2p2_inf(2,i-1).gt.0.and.c2p2_inf(2,i).lt.0.) then
	     print*,phase_inf(i-1),phase_inf(i)
             cphase=cphase-1

	  else if(c2p2_inf(2,i-1).lt.0.and.c2p2_inf(2,i).gt.0.) then
	     print*,phase_inf(i-1),phase_inf(i)
             cphase=cphase+1

	  end if
	end if
	phase_inf(i) = phase_inf(i) + (cphase)*2.*pi

	 write(137,100) tortutime(i), &
     &                  sqrt( c2p2_inf(1,i)**2 + c2p2_inf(2,i)**2 ), &
     &                             phase_inf(i), localtime

	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!

      ! integrate twice in time to get the strain
      d2p2T = 0.0
      d2p2T_inf = 0.0

      ! first time integral 
      do i=2, imax
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+2.*dt .GT. time_min) then
          d2p2T(:,i) = d2p2T(:,i-1) &
     &               + 0.5*dt*(c2p2(:,i-1) + c2p2(:,i))
          d2p2T_inf(:,i) = d2p2T_inf(:,i-1) &
     &               + 0.5*dt*(c2p2_inf(:,i-1) + c2p2_inf(:,i))


        write(133,100) tortutime(i), d2p2T(1,i), d2p2T(2,i), localtime
!        write(134,100) tortutime(i), d2p2T_inf(1,i), d2p2T_inf(2,i), localtime

        end if
      end do 

!!!!!!!!!!!!!!!!!!!!! the infinity part is picking a linear drift, 
! let us remove to see  if we can fix it somehow
! we will reuse variables defined for the h2p2 part


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      h2p2 = 0.0
      h2p2_inf = 0.0
      ! second time integral 
!!!! and we want to remove
!!! using Steve's suggestion to do so
      h22_ave = 0.0
      h22_inf_ave = 0.0
      xyh22_ave = 0.0
      xyh22_inf_ave = 0.0
      count_ave = 0

      do i=2, imax
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+2.*dt .GT. time_min) then

	h22_inf_ave = h22_inf_ave + d2p2T_inf(:,i)
	xyh22_inf_ave = xyh22_inf_ave + d2p2T_inf(:,i)*time(i)

	time_ave = time_ave + time(i)
	xxtime_ave = xxtime_ave + time(i)*time(i)
	count_ave = count_ave + 1
	else
	time_ave =  time(i)
	xxtime_ave = time(i)*time(i)
	count_ave = 1
	end if

      end do 

	slopeh22_inf =  (xyh22_inf_ave - h22_inf_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22_inf = h22_inf_ave/count_ave - slopeh22_inf*time_ave/count_ave


!remove drift
	do i=2, imax
	 if(time(i)+2.*dt.gt.time_min) then

	 d2p2T_inf(:,i) = d2p2T_inf(:,i) - ( slopeh22_inf*time(i) + interh22_inf)
        write(134,100) tortutime(i), d2p2T_inf(1,i), d2p2T_inf(2,i), time(i)

	 end if
	end do

	print*,'slope',interh22_inf, slopeh22_inf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      h2p2 = 0.0
      h2p2_inf = 0.0
      ! second time integral 
!!!! and we want to remove
!!! using Steve's suggestion to do so
      h22_ave = 0.0
      h22_inf_ave = 0.0
      xyh22_ave = 0.0
      xyh22_inf_ave = 0.0
      count_ave = 0

      do i=2, imax
        localtime = time(i)
	dt = time(i)-time(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime .GT. time_min) then
          h2p2(:,i) = h2p2(:,i-1) &
     &               + 0.5*dt*(d2p2T(:,i-1) + d2p2T(:,i))
          h2p2_inf(:,i) = h2p2_inf(:,i-1) &
     &                   + 0.5*dt*(d2p2T_inf(:,i-1) + d2p2T_inf(:,i))
!        end if


	h22_ave = h22_ave + h2p2(:,i)
	h22_inf_ave = h22_inf_ave + h2p2_inf(:,i)

	xyh22_ave = xyh22_ave + h2p2(:,i)*time(i)
	xyh22_inf_ave = xyh22_inf_ave + h2p2_inf(:,i)*time(i)

	time_ave = time_ave + time(i)
	xxtime_ave = xxtime_ave + time(i)*time(i)
	count_ave = count_ave + 1
	else
	time_ave =  time(i)
	xxtime_ave = time(i)*time(i)
	count_ave = 1
	end if
	write(101,*) i, time_ave, h22_ave, xyh22_ave

      end do 

	slopeh22 =  (xyh22_ave - h22_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22 = h22_ave/count_ave - slopeh22*time_ave/count_ave

	slopeh22_inf =  (xyh22_inf_ave - h22_inf_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22_inf = h22_inf_ave/count_ave - slopeh22_inf*time_ave/count_ave


!remove drift
	do i=2, imax
	 if(time(i).gt.time_min) then

	 h2p2(:,i) = h2p2(:,i) - ( slopeh22*time(i) + interh22)

	 h2p2_inf(:,i) = h2p2_inf(:,i) - ( slopeh22_inf*time(i) + interh22_inf)

        write(135,100) tortutime(i), h2p2(1,i), h2p2(2,i), time(i)
        write(136,100) tortutime(i), h2p2_inf(1,i) , &
     &                       h2p2_inf(2,i), time(i)

	 end if
	end do

	print*,'slope',interh22, slopeh22



!!!!!! finally we note the phase and amplitudes above were obtained
!!! from c2p2, but that isnot the strain. so let us redo here with strains

!!! first from h2p2

	cphase = 0

	do i = 2, imax - 1
        localtime = time(i)
        x = h2p2(1,i)
        y = h2p2(2,i)

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
	  if (h2p2(2,i-1).gt.0.and.h2p2(2,i).lt.0.) then
	     print*,HHphase(i-1),HHphase(i)
             cphase=cphase-1

	  else if(h2p2(2,i-1).lt.0.and.h2p2(2,i).gt.0.) then
	     print*,HHphase(i-1),HHphase(i)
             cphase=cphase+1

	  end if
	end if
	HHphase(i) = HHphase(i) + (cphase)*2.*pi

	 write(138,100) tortutime(i), &
     &                  sqrt( h2p2(1,i)**2 + h2p2(2,i)**2 ), &
     &                             HHphase(i), time(i)

	end if
	end do

!!!! then with h2p2_inf
	cphase = 0

	do i = ISTART, imax - 1
        localtime = time(i)
        x = h2p2_inf(1,i)
        y = h2p2_inf(2,i)
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
	  if (h2p2_inf(2,i-1).gt.0.and.h2p2_inf(2,i).lt.0.) then
	     print*,HHphase_inf(i-1),HHphase_inf(i)
             cphase=cphase-1

	  else if(h2p2_inf(2,i-1).lt.0.and.h2p2_inf(2,i).gt.0.) then
	     print*,HHphase_inf(i-1),HHphase_inf(i)
             cphase=cphase+1

	  end if
	end if
	HHphase_inf(i) = HHphase_inf(i) + (cphase)*2.*pi

	 write(139,100) tortutime(i), &
     &                  sqrt( h2p2_inf(1,i)**2 + h2p2_inf(2,i)**2 ), &
     &                             HHphase_inf(i), time(i)

	end if
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      open(unit=114, file="c2p2.dat")
      open(unit=314, file="c2p2_inf.dat")
 
      open(unit=122, file="phase.dat")            

      open(unit=130, file="omega.dat")                        
      open(unit=131, file="amp_and_phase.dat")
      open(unit=137, file="amp_and_phase_inf.dat")                                                

      open(unit=138, file="HHamp_and_phase.dat")
      open(unit=139, file="HHamp_and_phase_inf.dat")                   

      open(unit=133, file="d2p2.dat")
      open(unit=134, file="d2p2_inf.dat")

      open(unit=135, file="h2p2.dat")
      open(unit=136, file="h2p2_inf.dat")


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



!-----------------------------------------
      
  100 format(4e16.7)      

      end program

