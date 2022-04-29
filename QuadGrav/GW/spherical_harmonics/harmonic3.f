
      program      main
      implicit     none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     decomposition into spin-weighted spherical harmonics
!
!     cos(theta) = z/r,  tg(phi) = y/x
!
!     u = u(phi, theta)
!
!     phi = {0, 2*pi}, theta = {0, pi}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
            
      integer :: i, j, k, N, itheta, iphi, output, initial, itime
      integer :: correction, integration, ISTART, imax, hitwithr
      real*8  :: pi, theta, phi, integral, h, h2, x, y, localtime
      real*8  :: rp, ip, norm, dt, factor, temp1, time_tort
      real*8  :: massI, rext, test, smallnumber, ra, rr, fpR, fpI
      real*8  :: surface_integral,surface_integral2, givefile
      integer :: cphase
      real*8  :: rfactorhit
            
! ------- this parameter changes depending on the size of u------------
      integer Ndim
!      parameter ( Ndim = 80 )


      character*20    file_psi
      integer         read_shape(2), shape(2)
      character*5     read_cnames
      integer         read_rank, myi
      real*8          read_coords(1000), bbox(4) 
      real*8          read_time, time, timeold, time_min
      real*8          r8arg
      external        r8arg
      
      integer         iargc, gf3_rc, gft_read_full,index
      character*5     name, sarg
      external        sarg, gft_read_full
      character*32    realpsiname,imgpsiname,massname,jzname

      real*8  :: tY2m2,tY2m1,tY20,tY2p1
      real*8  :: tY2p2,tY3m3,tY3m2,tY3m1
      real*8  :: tY30,tY3p1,tY3p2,tY3p3
      real*8  :: tY4m4,tY4m3,tY4m2,tY4m1
      real*8  :: tY40,tY4p1,tY4p2,tY4p3,tY4p4
      real*8  :: cI2p2,cR2p2,dtcR2p2,dtCI2p2
      real*8  :: C0, factorl2,factorl3,factorl4
      real*8  :: d2p2(2),d2p1(2),d20(2),d2m1(2),d2m2(2)
      real*8  :: d2p2_inf(2),d2p1_inf(2),d20_inf(2),d2m1_inf(2)
      real*8  :: d2m2_inf(2)
      real*8  :: integralc2m2(2),integralc2m1(2),integralc20(2),integralc2p1(2),integralc2p2(2)
      real*8  :: integralc3m2(2),integralc3m1(2),integralc30(2),integralc3p1(2),integralc3p2(2)
      real*8  :: integralc3m3(2),integralc3p3(2)
      real*8  :: integralc4m2(2),integralc4m1(2),integralc40(2),integralc4p1(2),integralc4p2(2)
      real*8  :: integralc4m3(2),integralc4m4(2),integralc4p3(2),integralc4p4(2)

      integer     MAXMEM
!      parameter ( MAXMEM = 1000000 )
      parameter ( MAXMEM = 14000 )
      real*8  :: timein(MAXMEM), tortutime(MAXMEM)
      real*8  :: rpsi(MAXMEM),ipsi(MAXMEM),mass(MAXMEM),Jz(MAXMEM)
      real*8  :: Y2m2(2,MAXMEM),Y2m1(2,MAXMEM),Y20(2,MAXMEM),Y2p1(2,MAXMEM)
      real*8  :: Y2p2(2,MAXMEM),Y3m3(2,MAXMEM),Y3m2(2,MAXMEM),Y3m1(2,MAXMEM)
      real*8  :: Y30(2,MAXMEM),Y3p1(2,MAXMEM),Y3p2(2,MAXMEM),Y3p3(2,MAXMEM)
      real*8  :: Y4m4(2,MAXMEM),Y4m3(2,MAXMEM),Y4m2(2,MAXMEM),Y4m1(2,MAXMEM)
      real*8  :: Y40(2,MAXMEM),Y4p1(2,MAXMEM),Y4p2(2,MAXMEM),Y4p3(2,MAXMEM)
      real*8  :: Y4p4(2,MAXMEM), temp(2,MAXMEM)

      integer     MAXTIME
      parameter ( MAXTIME = 10000 )
      real*8  :: c2m2(2,MAXTIME),c2m1(2,MAXTIME),c20(2,MAXTIME)
      real*8  :: c2p2(2,MAXTIME),c2p1(2,MAXTIME)
      real*8  :: c2m2_inf(2,MAXTIME),c2m1_inf(2,MAXTIME),c20_inf(2,MAXTIME)
      real*8  :: c2p2_inf(2,MAXTIME),c2p1_inf(2,MAXTIME)
      real*8  :: c3m3(2,MAXTIME),c3m2(2,MAXTIME),c3m1(2,MAXTIME),c30(2,MAXTIME)
      real*8  :: c3p1(2,MAXTIME),c3p2(2,MAXTIME),c3p3(2,MAXTIME)
      real*8  :: c4m4(2,MAXTIME),c4m3(2,MAXTIME),c4m2(2,MAXTIME),c4m1(2,MAXTIME)
      real*8  :: c40(2,MAXTIME),c4p1(2,MAXTIME),c4p2(2,MAXTIME),c4p3(2,MAXTIME)
      real*8  :: c4p4(2,MAXTIME),tmass(MAXTIME),tJz(MAXTIME),phase(MAXTIME) 
      real*8  :: c3m3_inf(2,MAXTIME),c3m2_inf(2,MAXTIME),c3m1_inf(2,MAXTIME),c30_inf(2,MAXTIME)
      real*8  :: c3p1_inf(2,MAXTIME),c3p2_inf(2,MAXTIME),c3p3_inf(2,MAXTIME)
      real*8  :: c4m4_inf(2,MAXTIME),c4m3_inf(2,MAXTIME),c4m2_inf(2,MAXTIME),c4m1_inf(2,MAXTIME)
      real*8  :: c40_inf(2,MAXTIME),c4p1_inf(2,MAXTIME),c4p2_inf(2,MAXTIME),c4p3_inf(2,MAXTIME)
      real*8  :: c4p4_inf(2,MAXTIME)
      real*8  :: omega1(MAXTIME),omega2(MAXTIME),LGW(MAXTIME),LGW_inf(MAXTIME)
      real*8  :: d2p2T(2,MAXTIME),h2p2(2,MAXTIME)
      real*8  :: d2p2T_inf(2,MAXTIME),h2p2_inf(2,MAXTIME)
      real*8  :: phase_inf(MAXTIME), amp_inf
      real*8  :: HHphase_inf(MAXTIME), HHamp_inf, HHphase(MAXTIME), HHamp
      real*8  :: h22_ave(2), h22_inf_ave(2), time_ave
      real*8  :: xxtime_ave, xyh22_ave(2), xyh22_inf_ave(2)
      real*8 ::  count_ave
      real*8  :: slopeh22(2),interh22(2),slopeh22_inf(2),interh22_inf(2)

      if (iargc() .lt. 1) then
          write(*,*) '*************************************************'
          write(*,*) '* harmonic3: Produces a decomposition           *'
          write(*,*) '*            into -2 spin-weighted harmonics.   *'
          write(*,*) '* Usage:                                        *'
          write(*,*) '*************************************************'
          write(*,*) '  harmonic3  <surface letter> <mass> <radius> <timeforh22> <hitwithr> <givefiles?>'
          write(*,*) ' hitwithr>0 will multiply psi4 by rext'
          write(*,*) ' givefiles>0 assume names as in had, <1 give them below'
          write(*,*) '   E.g. harmonic3 b 2.7 140.d0 300 1 -1'
          write(*,*) '                *****************'
          stop
       end if

       name = sarg( 1,'c')
       massI= r8arg(2,1.d0)
       rext = r8arg(3,140.d0)
       time_min = r8arg(4,1.d0)
       hitwithr = r8arg(5,1.d0)
       givefile = r8arg(6,1.d0)

!for tortoise time    
      !rr        = rext * ( 1.0 + 2.0*massI/rext )**2
      rr        = rext * ( 1.0 + massI/(2.0*rext) )**2
    

       realpsiname = 'r'//name(:1)//'p4.sdf'
       imgpsiname  = 'i'//name(:1)//'p4.sdf'
       massname    = name(:1)//'mass.sdf'
       jzname      = name(:1)//'Jz.sdf'


       if (.true.) then
          write(*,*) '   name = ',name
          write(*,*) 'Reading in the following files:'
          write(*,*) '   realpsiname = ',realpsiname
          write(*,*) '   imgpsiname  = ',imgpsiname
          write(*,*) '   massname    = ',massname
          write(*,*) '   jzname      = ',jzname
          write(*,*) '   massI       = ',massI
          write(*,*) '   rext        = ',rext
       end if

	if(givefile.lt.0) then
	print*,'read files dont assume names: real and im'
	read*,realpsiname
	read*,imgpsiname
	end if

	rfactorhit = 1.0d0
	if(hitwithr.gt.0) then
	print*,'will hit psi4 with the extraction radius'
	rfactorhit = rext
	end if


       !----reading parameters from a file------------------
       itime = 2                                   
       gf3_rc =  gft_read_full(realpsiname, itime,  read_shape, &
      &     read_cnames, read_rank, timein(itime), read_coords, rpsi)      
       if (gf3_rc .ne.1) goto 99
       if (read_shape(1)*read_shape(2).gt.MAXMEM) then
          write(*,*)'read_shape(1)/(2) = ',read_shape(1),read_shape(2)
          write(*,*)'read_shape(1)*(2) = ',read_shape(1)*read_shape(2)
          write(*,*)'MAXMEM            = ',MAXMEM
          write(*,*)'Not enough memory, recompile. Quitting.'
          goto 99
       end if

       Ndim = read_shape(2)-1     

       print*, "read_shape=",read_shape
       print*, "Ndim=",Ndim
       print*,'MAXMEM = ',read_shape(1)*read_shape(2)
 
!       stop

	cphase = 0.
!-----------------------------------------------
!-----------------------------------------------
       output = 0        !-- 0, only .dat, 1 also output on the screen
       !massI  = 1.0      !-- the initial (m1 + m2) of the isolated bodies
       !rext   = 140.0     !-- the radius of the surface extraction
       !time_min = 300.0    !-- time to start the integration of Ds and Es
       factor  = 1.0     ! arbitrary normalization factor to match other results
      
       N      = Ndim
       pi     = dacos(-1.0d0)
       h      = pi/N
       h2     = h*h
!------------------------------------------------------
      open(unit=110, file="c2m2.dat")
      open(unit=111, file="c2m1.dat")
      open(unit=112, file="c20.dat")
      open(unit=113, file="c2p1.dat")
      open(unit=114, file="c2p2.dat")

      open(unit=210, file="c3m3.dat")
      open(unit=211, file="c3m2.dat")
      open(unit=212, file="c3m1.dat")
      open(unit=213, file="c30.dat")
      open(unit=214, file="c3p1.dat")
      open(unit=215, file="c3p2.dat")
      open(unit=216, file="c3p3.dat")

      open(unit=220, file="c4m4.dat")
      open(unit=221, file="c4m3.dat")
      open(unit=222, file="c4m2.dat")
      open(unit=223, file="c4m1.dat")
      open(unit=224, file="c40.dat")
      open(unit=225, file="c4p1.dat")
      open(unit=226, file="c4p2.dat")
      open(unit=227, file="c4p3.dat")
      open(unit=228, file="c4p4.dat")

      open(unit=310, file="c2m2_inf.dat")
      open(unit=311, file="c2m1_inf.dat")
      open(unit=312, file="c20_inf.dat")
      open(unit=313, file="c2p1_inf.dat")
      open(unit=314, file="c2p2_inf.dat")

      open(unit=320, file="c3m3_inf.dat")
      open(unit=321, file="c3m2_inf.dat")
      open(unit=322, file="c3m1_inf.dat")
      open(unit=323, file="c30_inf.dat")
      open(unit=324, file="c3p1_inf.dat")
      open(unit=325, file="c3p2_inf.dat")
      open(unit=326, file="c3p3_inf.dat")

      open(unit=330, file="c4m4_inf.dat")
      open(unit=331, file="c4m3_inf.dat")
      open(unit=332, file="c4m2_inf.dat")
      open(unit=333, file="c4m1_inf.dat")
      open(unit=334, file="c40_inf.dat")
      open(unit=335, file="c4p1_inf.dat")
      open(unit=336, file="c4p2_inf.dat")
      open(unit=337, file="c4p3_inf.dat")
      open(unit=338, file="c4p4_inf.dat")


      open(unit=120, file="mass.dat")      
      open(unit=121, file="Jz.dat")      
      open(unit=122, file="phase.dat")            

      open(unit=130, file="omega.dat")                        
      open(unit=131, file="amp_and_phase.dat")
      open(unit=137, file="amp_and_phase_inf.dat")                                                
      open(unit=132, file="lum_gw.dat")      

      open(unit=138, file="HHamp_and_phase.dat")
      open(unit=139, file="HHamp_and_phase_inf.dat")                   

      open(unit=133, file="d2p2.dat")
      open(unit=134, file="d2p2_inf.dat")

      open(unit=135, file="h2p2.dat")
      open(unit=136, file="h2p2_inf.dat")

!------------------------------------------------------
      itime = 1                                    

!------------------------------------------------------------            
!-----DEFINE SPIN WEIGHTED SPHERICAL HARMONIC DECOMPOSITION--
!------------------------------------------------------------

!!------computing the spin-weighted spherical harmonics Y'_l,m------------	    
!-----  Y_l,m = Y'_l,m * exp ( +i*m*phi) ---------------------------------
!-----  Y_l,m = Y'_l,m * [cos( +i*m*phi) + i*sin( +i*m*phi)] -------------
!----- *Y_l,m = Y'_l,m * [cos( +i*m*phi) - i*sin( +i*m*phi)] -------------
            
      do iphi = 1, 2*N+1
      phi = (iphi-1)*h
      do itheta = 1, N+1
        theta = (itheta-1)*h   
	x = dcos(theta)

        index = (itheta-1)*read_shape(1)+iphi
	    
!        tY2m2 = (x-1)**2*sqrt(5.D0)/sqrt(0.3141592653589793D1)/8
!        tY2m1 = sqrt(x+1)*sqrt(1-x)**3*sqrt(5.D0)/sqrt(0.3141592653589793D1)/4
!        tY20 = -(x-1)*(x+1)*sqrt(5.D0)*sqrt(6.D0)/sqrt(0.3141592653589793D1)/8
!        tY2p1 = sqrt(5.D0)*sqrt(1-x)*sqrt(x+1)**3/sqrt(0.3141592653589793D1)/4
!        tY2p2 = sqrt(5.D0)*(x+1)**2/sqrt(0.3141592653589793D1)/8

        !from Brugman
!        tY2m2 = sqrt(5.0/(64.0*pi)) * (1.0 - cos(theta))**2
!        tY2m1 =-sqrt(5.0/(16.0*pi)) * sin(theta) * (1.0 - cos(theta))
!        tY20  = sqrt(15.0/(32.0*pi)) * (sin(theta))**2
!        tY2p1 =-sqrt(5.0/(16.0*pi)) * sin(theta) * (1.0 + cos(theta)) 
!        tY2p2 = sqrt(5.0/(64.0*pi)) * (1.0 + cos(theta))**2
       	    		    	    
!	Y2m2(1,index) = tY2m2*cos(-2.0*phi)
!	Y2m1(1,index) = tY2m1*cos(-1.0*phi)
!	Y20(1,index)  = tY20*cos(0.0*phi)
!	Y2p1(1,index) = tY2p1*cos(+1.0*phi)
!	Y2p2(1,index) = tY2p2*cos(+2.0*phi)

!	Y2m2(2,index) =-tY2m2*sin(-2.0*phi)
!	Y2m1(2,index) =-tY2m1*sin(-1.0*phi)
!	Y20(2,index)  =-tY20*sin(0.0*phi)
!	Y2p1(2,index) =-tY2p1*sin(+1.0*phi)
!	Y2p2(2,index) =-tY2p2*sin(+2.0*phi)

        !from mathematica
        ! l=2
        tY2m2 = 0.5*dsqrt(5.0/pi) * (dsin(theta/2.0))**4
        tY2m1 = 0.5*dsqrt(5.0/pi) * dsin(theta) * (dsin(theta/2.0))**2
        tY20  = 0.25*dsqrt(15.0/(2.0*pi)) * (dsin(theta))**2
        tY2p1 = 0.25*dsqrt(5.0/pi) * dsin(theta) * (1.0 + dcos(theta)) 
        tY2p2 = 0.5*dsqrt(5.0/pi) * (dcos(theta/2.0))**4
       	    		    	    
	Y2m2(1,index) = tY2m2*dcos(-2.0*phi)
	Y2m1(1,index) = tY2m1*dcos(-1.0*phi)
	Y20(1,index)  = tY20
	Y2p1(1,index) = tY2p1*dcos(+1.0*phi)
	Y2p2(1,index) = tY2p2*dcos(+2.0*phi)

	Y2m2(2,index) =-tY2m2*dsin(-2.0*phi)
	Y2m1(2,index) =-tY2m1*dsin(-1.0*phi)
	Y20(2,index)  =-tY20*0.0
	Y2p1(2,index) =-tY2p1*dsin(+1.0*phi)
	Y2p2(2,index) =-tY2p2*dsin(+2.0*phi)

        ! l=3
        tY3m3 = 0.5*sqrt(21.0/(2.0*pi)) * sin(theta) * (sin(theta/2.0))**4 
        tY3m2 = 0.5*sqrt(7.0/pi) * (2.0 + 3.0*cos(theta)) * (sin(theta/2.0))**4
        tY3m1 = (1.0/32.0)*sqrt(35.0/(2.0*pi)) * (sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))
        tY30  = 0.25*sqrt(105.0/(2.0*pi)) * cos(theta) * (sin(theta))**2
        tY3p1 =-(1.0/32.0)*sqrt(35.0/(2.0*pi)) * (sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))
        tY3p2 = 0.5*sqrt(7.0/pi) * (-2.0 + 3.0*cos(theta)) * (cos(theta/2.0))**4
        tY3p3 =-sqrt(21.0/(2.0*pi)) * sin(theta/2.0) * (cos(theta/2.0))**5

	Y3m3(1,index) = tY3m3*cos(-3.0*phi)       	    		    	    
	Y3m2(1,index) = tY3m2*cos(-2.0*phi)
	Y3m1(1,index) = tY3m1*cos(-1.0*phi)
	Y30(1,index)  = tY30*cos(0.0*phi)
	Y3p1(1,index) = tY3p1*cos(+1.0*phi)
	Y3p2(1,index) = tY3p2*cos(+2.0*phi)
	Y3p3(1,index) = tY3p3*cos(+3.0*phi)

	Y3m3(2,index) =-tY3m3*sin(-3.0*phi)
	Y3m2(2,index) =-tY3m2*sin(-2.0*phi)
	Y3m1(2,index) =-tY3m1*sin(-1.0*phi)
	Y30(2,index)  =-tY30*sin(0.0*phi)
	Y3p1(2,index) =-tY3p1*sin(+1.0*phi)
	Y3p2(2,index) =-tY3p2*sin(+2.0*phi)
	Y3p3(2,index) =-tY3p3*sin(+3.0*phi)

        ! l=4
        tY4m4 = (3.0/4.0)*sqrt(7.0/pi) * (sin(theta))**2 * (sin(theta/2.0))**4 
        tY4m3 = 3.0*sqrt(7.0/(2.0*pi)) * (1.0 + 2.0*cos(theta)) * (sin(theta/2.0))**5
        tY4m2 = (3.0/4.0)*sqrt(1.0/pi) * (9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta)) * (sin(theta/2.0))**4
        tY4m1 = (3.0/32.0)*sqrt(1.0/(2.0*pi)) * (3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta))
        tY40  = (3.0/16.0)*sqrt(5.0/(2.0*pi)) * (5.0 + 7.0*cos(2.0*theta)) * (sin(theta))**2
        tY4p1 = (3.0/32.0)*sqrt(1.0/(2.0*pi)) * (3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + sin(4.0*theta))
        tY4p2 = (3.0/4.0)*sqrt(1.0/pi) * (9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)) * (cos(theta/2.0))**4
        tY4p3 =-3.0*sqrt(7.0/(2.0*pi)) * (-1.0 + 2.0*cos(theta)) * (sin(theta/2.0)) * (cos(theta/2.0))**5
        tY4p4 = (3.0/4.0)*sqrt(7.0/pi) * (sin(theta))**2 * (cos(theta/2.0))**(4) 

	Y4m4(1,index) = tY4m4*cos(-4.0*phi)
	Y4m3(1,index) = tY4m3*cos(-3.0*phi)
	Y4m2(1,index) = tY4m2*cos(-2.0*phi)
	Y4m1(1,index) = tY4m1*cos(-1.0*phi)
	Y40(1,index)  = tY40*cos(0.0*phi)
	Y4p1(1,index) = tY4p1*cos(+1.0*phi)
	Y4p2(1,index) = tY4p2*cos(+2.0*phi)
	Y4p3(1,index) = tY4p3*cos(+3.0*phi)
	Y4p4(1,index) = tY4p4*cos(+4.0*phi)

	Y4m4(2,index) =-tY4m4*sin(-4.0*phi)
	Y4m3(2,index) =-tY4m3*sin(-3.0*phi)
	Y4m2(2,index) =-tY4m2*sin(-2.0*phi)
	Y4m1(2,index) =-tY4m1*sin(-1.0*phi)
	Y40(2,index)  =-tY40*sin(0.0*phi)
	Y4p1(2,index) =-tY4p1*sin(+1.0*phi)
	Y4p2(2,index) =-tY4p2*sin(+2.0*phi)
	Y4p3(2,index) =-tY4p3*sin(+3.0*phi)
	Y4p4(2,index) =-tY4p4*sin(+4.0*phi)
	       	       
      end do
      end do


!-----------------------------------------
      ! for testing only 
      ! first that the integral of 1 is 4*pi
!       mass = 1.0
!       rpsi = 1.0
!       ipsi = 0.0


!       test = surface_integral(mass,MAXMEM,Ndim)                                        
!       print*, "4*pi is",test,"and it should be", 4.0*pi
!       test = surface_integral2(mass,MAXMEM,Ndim)
!       print*, "4*pi is nowd",test,"and it should be", 4.0*pi

!       temp(1,:) = rpsi(:)*Y2p2(1,:) - ipsi(:)*Y2p2(2,:)
!       temp(2,:) = rpsi(:)*Y2p2(2,:) + ipsi(:)*Y2p2(1,:)
!       print*, "N=", N
!       test = surface_integral2(temp(1,:),MAXMEM,Ndim)
!       print*, "Y22 ortonormalization",test
!       test = surface_integral2(temp(2,:),MAXMEM,Ndim)
!       print*, "Y22 ortonormalization",test
!       stop     

!-----------------------------------------
      itime = 0

      do while (.true.)

        itime = itime + 1

        !----reading from a file----------------------             
        gf3_rc =  gft_read_full(realpsiname, itime,  read_shape, &
        &     read_cnames, read_rank, timein(itime), read_coords, rpsi)      
        if (gf3_rc .ne.1) goto 99
        gf3_rc =  gft_read_full(imgpsiname, itime,  read_shape, &
        &     read_cnames, read_rank, timein(itime), read_coords, ipsi)            
        if (gf3_rc .ne.1) goto 99
!        gf3_rc = gft_read_full(massname,   itime,  read_shape, &
!        &     read_cnames, read_rank, timein(itime), read_coords, mass)            
!        if (gf3_rc .ne.1) goto 99
!        gf3_rc = gft_read_full(jzname,    itime,  read_shape, &
!        &     read_cnames, read_rank, timein(itime), read_coords, Jz)            
!       if (gf3_rc .ne.1) goto 99

	mass = massI
	Jz = 0.0
        print*, "time=", itime, timein(itime)
        time = timein(itime)
         
        ! Some strangeness sometimes appears at the boundaries
        ! and this will regularize them:
        do myi= 1, read_shape(1)*read_shape(2)
           if (abs(rpsi(myi)).gt.1e5)rpsi(myi)=0.d0
           if (abs(ipsi(myi)).gt.1e5)ipsi(myi)=0.d0
        end do

        rpsi = rpsi*factor*rfactorhit
        ipsi = ipsi*factor*rfactorhit

        ! l=2         	 	 
        temp(1,:) = rpsi(:)*Y2m2(1,:) - ipsi(:)*Y2m2(2,:)
        temp(2,:) = rpsi(:)*Y2m2(2,:) + ipsi(:)*Y2m2(1,:)
        c2m2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c2m2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y2m1(1,:) - ipsi(:)*Y2m1(2,:)
        temp(2,:) = rpsi(:)*Y2m1(2,:) + ipsi(:)*Y2m1(1,:)
        c2m1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c2m1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y20(1,:) - ipsi(:)*Y20(2,:)
        temp(2,:) = rpsi(:)*Y20(2,:) + ipsi(:)*Y20(1,:)
        c20(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c20(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y2p1(1,:) - ipsi(:)*Y2p1(2,:)
        temp(2,:) = rpsi(:)*Y2p1(2,:) + ipsi(:)*Y2p1(1,:)
        c2p1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c2p1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y2p2(1,:) - ipsi(:)*Y2p2(2,:)
        temp(2,:) = rpsi(:)*Y2p2(2,:) + ipsi(:)*Y2p2(1,:)
        c2p2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c2p2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        ! l=3        	 	 
        temp(1,:) = rpsi(:)*Y3m3(1,:) - ipsi(:)*Y3m3(2,:)
        temp(2,:) = rpsi(:)*Y3m3(2,:) + ipsi(:)*Y3m3(1,:)
        c3m3(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3m3(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y3m2(1,:) - ipsi(:)*Y3m2(2,:)
        temp(2,:) = rpsi(:)*Y3m2(2,:) + ipsi(:)*Y3m2(1,:)
        c3m2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3m2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y3m1(1,:) - ipsi(:)*Y3m1(2,:)
        temp(2,:) = rpsi(:)*Y3m1(2,:) + ipsi(:)*Y3m1(1,:)
        c3m1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3m1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y30(1,:) - ipsi(:)*Y30(2,:)
        temp(2,:) = rpsi(:)*Y30(2,:) + ipsi(:)*Y30(1,:)
        c30(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c30(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y3p1(1,:) - ipsi(:)*Y3p1(2,:)
        temp(2,:) = rpsi(:)*Y3p1(2,:) + ipsi(:)*Y3p1(1,:)
        c3p1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3p1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y3p2(1,:) - ipsi(:)*Y3p2(2,:)
        temp(2,:) = rpsi(:)*Y3p2(2,:) + ipsi(:)*Y3p2(1,:)
        c3p2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3p2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y3p3(1,:) - ipsi(:)*Y3p3(2,:)
        temp(2,:) = rpsi(:)*Y3p3(2,:) + ipsi(:)*Y3p3(1,:)
        c3p3(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c3p3(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        ! l=4
        temp(1,:) = rpsi(:)*Y4m4(1,:) - ipsi(:)*Y4m4(2,:)
        temp(2,:) = rpsi(:)*Y4m4(2,:) + ipsi(:)*Y4m4(1,:)
        c4m4(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4m4(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4m3(1,:) - ipsi(:)*Y4m3(2,:)
        temp(2,:) = rpsi(:)*Y4m3(2,:) + ipsi(:)*Y4m3(1,:)
        c4m3(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4m3(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4m2(1,:) - ipsi(:)*Y4m2(2,:)
        temp(2,:) = rpsi(:)*Y4m2(2,:) + ipsi(:)*Y4m2(1,:)
        c4m2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4m2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4m1(1,:) - ipsi(:)*Y4m1(2,:)
        temp(2,:) = rpsi(:)*Y4m1(2,:) + ipsi(:)*Y4m1(1,:)
        c4m1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4m1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y40(1,:) - ipsi(:)*Y40(2,:)
        temp(2,:) = rpsi(:)*Y40(2,:) + ipsi(:)*Y40(1,:)
        c40(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c40(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4p1(1,:) - ipsi(:)*Y4p1(2,:)
        temp(2,:) = rpsi(:)*Y4p1(2,:) + ipsi(:)*Y4p1(1,:)
        c4p1(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4p1(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4p2(1,:) - ipsi(:)*Y4p2(2,:)
        temp(2,:) = rpsi(:)*Y4p2(2,:) + ipsi(:)*Y4p2(1,:)
        c4p2(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4p2(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4p3(1,:) - ipsi(:)*Y4p3(2,:)
        temp(2,:) = rpsi(:)*Y4p3(2,:) + ipsi(:)*Y4p3(1,:)
        c4p3(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4p3(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        temp(1,:) = rpsi(:)*Y4p4(1,:) - ipsi(:)*Y4p4(2,:)
        temp(2,:) = rpsi(:)*Y4p4(2,:) + ipsi(:)*Y4p4(1,:)
        c4p4(1,itime) = surface_integral2(temp(1,:),MAXMEM,Ndim)
        c4p4(2,itime) = surface_integral2(temp(2,:),MAXMEM,Ndim)

        !ADM mass and Jz
        tmass(itime) = surface_integral2(mass(:),MAXMEM,Ndim)
        tJz(itime)   = surface_integral2(Jz(:),MAXMEM,Ndim)

        !if z = x + iy we define ϕ=Arg(z) as
        ! arctan(y/x) when x > 0
        ! arctan(y/x)+pi when x < 0 and y ≥ 0 
        ! arctan(y/x)-pi when x < 0 and y < 0 
        ! pi/2 when x = 0 and y > 0
        ! -pi/2 when x = 0 and y < 0 
        !indeterminate when x = 0 and y = 0.

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

!	print*,itime,phase(itime)



        !----output on the screen ---------------------------
        if (output .EQ. 1) then
  	  print*, "-------------------------"      !
          print*, "-----MASS--------------"
          print*, " mass=", tmass(itime)
	  print*, "-------------------------"      
          print*, "-----ANGULAR MOMENTUM----"
          print*, " Jz=", tJz(itime)		
	  print*, "-------------------------"      
          print*, "-----L=2 MODES ----------"
  	  print*, "-------------------------"
	  print*, "c2m2=",c2m2(1,itime),"c2m2=",c2m2(2,itime)
	  print*, "c2m1=",c2m1(1,itime),"c2m1=",c2m1(2,itime)
	  print*, "c20=",c20(1,itime),  "c20=",c20(2,itime)
	  print*, "c2p1=",c2p1(1,itime),"c2p1=",c2p1(2,itime)
	  print*, "c2p2=",c2p2(1,itime),"c2p2=",c2p2(2,itime)
        end if	

        time_tort = time - (rr + 2.0*massI*log(0.5*rr/massI -1.0))
	!tortutime(i) = time_tort

        write(110,100) time_tort, c2m2(1,itime), c2m2(2,itime), time 
        write(111,100) time_tort, c2m1(1,itime), c2m1(2,itime), time 
        write(112,100) time_tort, c20(1,itime), c20(2,itime), time 
        write(113,100) time_tort, c2p1(1,itime), c2p1(2,itime), time 
        write(114,100) time_tort, c2p2(1,itime), c2p2(2,itime), time 

        write(210,100) time_tort, c3m3(1,itime), c3m3(2,itime), time 
        write(211,100) time_tort, c3m2(1,itime), c3m2(2,itime), time 
        write(212,100) time_tort, c3m1(1,itime), c3m1(2,itime), time 
        write(213,100) time_tort, c30(1,itime), c30(2,itime), time 
        write(214,100) time_tort, c3p1(1,itime), c3p1(2,itime), time 
        write(215,100) time_tort, c3p2(1,itime), c3p2(2,itime), time 
        write(216,100) time_tort, c3p3(1,itime), c3p3(2,itime), time 

        write(220,100) time_tort, c4m4(1,itime), c4m4(2,itime), time 
        write(221,100) time_tort, c4m3(1,itime), c4m3(2,itime), time 
        write(222,100) time_tort, c4m2(1,itime), c4m2(2,itime), time 
        write(223,100) time_tort, c4m1(1,itime), c4m1(2,itime), time 
        write(224,100) time_tort, c40(1,itime), c40(2,itime), time 
        write(225,100) time_tort, c4p1(1,itime), c4p1(2,itime), time 
        write(226,100) time_tort, c4p2(1,itime), c4p2(2,itime), time 
        write(227,100) time_tort, c4p3(1,itime), c4p3(2,itime), time 
        write(228,100) time_tort, c4p4(1,itime), c4p4(2,itime), time

        write(120,100) time_tort, tmass(itime), time
        write(121,100) time_tort, tJz(itime), time
        write(122,100) time_tort, phase(itime), time


      end do

                  
 99   continue


       dt   = timein(2) - timein(1)
  
       print*, "timespaceing", dt

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
        localtime = timein(i)
        time_tort = localtime - (rr + 2.0*massI*log(0.5*rr/massI -1.0))
	tortutime(i) = time_tort

      do i=ISTART, imax
        time = timein(i)
        time_tort = time - (rr + 2.0*massI*log(0.5*rr/massI -1.0))
	tortutime(i) = time_tort
        if (time .GT. time_min) then
          integralc2m2 = integralc2m2 + c2m2(:,i)*dt 
          integralc2m1 = integralc2m1 + c2m1(:,i)*dt 
          integralc20  = integralc20  + c20(:,i)*dt 
          integralc2p1 = integralc2p1 + c2p1(:,i)*dt 
          integralc2p2 = integralc2p2 + c2p2(:,i)*dt 

          integralc3m3 = integralc3m3 + c3m3(:,i)*dt 
          integralc3m2 = integralc3m2 + c3m2(:,i)*dt 
          integralc3m1 = integralc3m1 + c3m1(:,i)*dt 
          integralc30  = integralc30  + c30(:,i)*dt 
          integralc3p1 = integralc3p1 + c3p1(:,i)*dt 
          integralc3p2 = integralc3p2 + c3p2(:,i)*dt 
          integralc3p3 = integralc3p3 + c3p3(:,i)*dt

          integralc4m4 = integralc4m4 + c4m4(:,i)*dt 
          integralc4m3 = integralc4m3 + c4m3(:,i)*dt 
          integralc4m2 = integralc4m2 + c4m2(:,i)*dt 
          integralc4m1 = integralc4m1 + c4m1(:,i)*dt 
          integralc40  = integralc40  + c40(:,i)*dt 
          integralc4p1 = integralc4p1 + c4p1(:,i)*dt 
          integralc4p2 = integralc4p2 + c4p2(:,i)*dt 
          integralc4p3 = integralc4p3 + c4p3(:,i)*dt
          integralc4p4 = integralc4p4 + c4p4(:,i)*dt
        end if

!this is hard coded to give c2p2 at infinity using nakano's extrapolation
!notice however that what was there before, was getting c2p2_inf but using
!c20 in the correction. I fixed it 

        c2m2_inf(:,i) = C0 * (c2m2(:,i) - factorl2*integralc2m2)
        c2m1_inf(:,i) = C0 * (c2m1(:,i) - factorl2*integralc2m1)
        c20_inf(:,i)  = C0 * (c20(:,i)  - factorl2*integralc20)
        c2p1_inf(:,i) = C0 * (c2p1(:,i) - factorl2*integralc2p1)
        c2p2_inf(:,i) = C0 * (c2p2(:,i) - factorl2*integralc2p2)

        c3m3_inf(:,i) = C0 * (c3m3(:,i) - factorl3*integralc3m3)
        c3m2_inf(:,i) = C0 * (c3m2(:,i) - factorl3*integralc3m2)
        c3m1_inf(:,i) = C0 * (c3m1(:,i) - factorl3*integralc3m1)
        c30_inf(:,i)  = C0 * (c30(:,i)  - factorl3*integralc30)
        c3p1_inf(:,i) = C0 * (c3p1(:,i) - factorl3*integralc3p1)
        c3p2_inf(:,i) = C0 * (c3p2(:,i) - factorl3*integralc3p2)
        c3p3_inf(:,i) = C0 * (c3p3(:,i) - factorl3*integralc3p3)

        c4m4_inf(:,i) = C0 * (c4m4(:,i) - factorl4*integralc4m4)
        c4m3_inf(:,i) = C0 * (c4m3(:,i) - factorl4*integralc4m3)
        c4m2_inf(:,i) = C0 * (c4m2(:,i) - factorl4*integralc4m2)
        c4m1_inf(:,i) = C0 * (c4m1(:,i) - factorl4*integralc4m1)
        c40_inf(:,i)  = C0 * (c40(:,i)  - factorl4*integralc40)
        c4p1_inf(:,i) = C0 * (c4p1(:,i) - factorl4*integralc4p1)
        c4p2_inf(:,i) = C0 * (c4p2(:,i) - factorl4*integralc4p2)
        c4p3_inf(:,i) = C0 * (c4p3(:,i) - factorl4*integralc4p3)
        c4p4_inf(:,i) = C0 * (c4p4(:,i) - factorl4*integralc4p4)

        write(310,100) time_tort, c2m2_inf(1,i), c2m2_inf(2,i), time
        write(311,100) time_tort, c2m1_inf(1,i), c2m1_inf(2,i), time
        write(312,100) time_tort, c20_inf(1,i), c20_inf(2,i), time
        write(313,100) time_tort, c2p1_inf(1,i), c2p1_inf(2,i), time
        write(314,100) time_tort, c2p2_inf(1,i), c2p2_inf(2,i), time

        write(320,100) time_tort, c3m3_inf(1,i), c3m3_inf(2,i), time
        write(321,100) time_tort, c3m2_inf(1,i), c3m2_inf(2,i), time
        write(322,100) time_tort, c3m1_inf(1,i), c3m1_inf(2,i), time
        write(323,100) time_tort, c30_inf(1,i), c30_inf(2,i), time
        write(324,100) time_tort, c3p1_inf(1,i), c3p1_inf(2,i), time
        write(325,100) time_tort, c3p2_inf(1,i), c3p2_inf(2,i), time
        write(326,100) time_tort, c3p3_inf(1,i), c3p3_inf(2,i), time

        write(330,100) time_tort, c4m4_inf(1,i), c4m4_inf(2,i), time
        write(331,100) time_tort, c4m3_inf(1,i), c4m3_inf(2,i), time
        write(332,100) time_tort, c4m2_inf(1,i), c4m2_inf(2,i), time
        write(333,100) time_tort, c4m1_inf(1,i), c4m1_inf(2,i), time
        write(334,100) time_tort, c40_inf(1,i), c40_inf(2,i), time
        write(335,100) time_tort, c4p1_inf(1,i), c4p1_inf(2,i), time
        write(336,100) time_tort, c4p2_inf(1,i), c4p2_inf(2,i), time
        write(337,100) time_tort, c4p3_inf(1,i), c4p3_inf(2,i), time
        write(338,100) time_tort, c4p4_inf(1,i), c4p4_inf(2,i), time



         omega1(i) = 0.5*( phase(i+1)-phase(i-1) )/(2.0*dt)
         cR2p2   = c2p2(1,i)
         cI2p2   = c2p2(2,i)

	if(i.gt.1.and.i.lt.imax) then
         dtcR2p2 = (c2p2(1,i+1) - c2p2(1,i-1))/(2.0*dt)
         dtcI2p2 = (c2p2(2,i+1) - c2p2(2,i-1))/(2.0*dt)
         temp1 = (dtcI2p2*cR2p2 - dtcR2p2*cI2p2)/(cR2p2**2 + cI2p2**2)    
         omega2(i) = abs((1.0d0/2.0d0)*temp1)
         write(130,100) time_tort, omega2(i),omega1(i), time
	end if

	 write(131,100) time_tort, &
     &                  sqrt( c2p2(1,i)**2 + c2p2(2,i)**2 ), &
     &                             phase(i), time

      end do 


!!!!!!! OK the above did a bunch of things wrt c22, but we have
!! the exptrapolated to infinity wavefor, so let us redo with that
!! as well... this file is getting big and uggly... but i saw beauty and the beast...

	cphase = 0

	do i = ISTART, imax
        time = timein(i)
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
     &                             phase_inf(i), time

	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!


      ! compute the GW luminosity and total energy radiated in GW
      ! first define some constant
      factor    = 1.0/(16.0*pi) 
      d2p2 = 0.0
      d2p1 = 0.0
      d20  = 0.0
      d2m1 = 0.0
      d2m2 = 0.0
      d2p2_inf = 0.0
      d2p1_inf = 0.0
      d20_inf  = 0.0
      d2m1_inf = 0.0
      d2m2_inf = 0.0
      print*, "ra=",ra,"factor", factor 

      do i=2, imax
        time = timein(i)
        ! with the raw data and the infinity extrapolation
        if (time .GT. time_min) then
          d2p2 = d2p2 + c2p2(:,i)*dt 
          d2p1 = d2p1 + c2p1(:,i)*dt 
          d20  = d20  + c20(:,i)*dt 
          d2m1 = d2m1 + c2m1(:,i)*dt 
          d2m2 = d2m2 + c2m2(:,i)*dt 

          d2p2_inf = d2p2_inf + c2p2_inf(:,i)*dt 
          d2p1_inf = d2p1_inf + c2p1_inf(:,i)*dt 
          d20_inf  = d20_inf  + c20_inf(:,i)*dt 
          d2m1_inf = d2m1_inf + c2m1_inf(:,i)*dt 
          d2m2_inf = d2m2_inf + c2m2_inf(:,i)*dt 
        end if

        LGW(i) = factor*( d2p2(1)**2 + d2p2(2)**2 +&
     &           d2p1(1)**2 + d2p1(2)**2 +&
     &           d20(1)**2 + d20(2)**2 +&
     &           d2m1(1)**2 + d2m1(2)**2 +&
     &           d2m2(1)**2 + d2m2(2)**2 )
    
        LGW_inf(i) = factor * (d2p2_inf(1)**2 + d2p2_inf(2)**2 +&
     &               d2p1_inf(1)**2 + d2p1_inf(2)**2 +&
     &               d20_inf(1)**2 + d20_inf(2)**2 +&
     &               d2m1_inf(1)**2 + d2m1_inf(2)**2 +&
     &               d2m2_inf(1)**2 + d2m2_inf(2)**2 )

        write(132,100) tortutime(i), LGW(i), LGW_inf(i), time
!        write(132,100) time, d20(1), d20(2), LGW(i)
      end do 


      ! integrate twice in time to get the strain
      d2p2T = 0.0
      d2p2T_inf = 0.0

      ! first time integral 
      do i=2, imax
        localtime = timein(i)
	dt = timein(i)-timein(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+2.*dt .GT. time_min) then
          d2p2T(:,i) = d2p2T(:,i-1) &
     &               + 0.5*dt*(c2p2(:,i-1) + c2p2(:,i))
          d2p2T_inf(:,i) = d2p2T_inf(:,i-1) &
     &                   + 0.5*dt*(c2p2_inf(:,i-1) + c2p2_inf(:,i))
        write(133,100) tortutime(i), d2p2T(1,i), d2p2T(2,i), timein(i)
!        write(134,100) time_tort, d2p2T_inf(1,i), d2p2T_inf(2,i), time
        end if


      end do 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! the infinity part is picking a linear drift, 
! let us remove to see  if we can fix it somehow
! we will reuse variables defined for the h2p2 part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      h2p2 = 0.0
      h2p2_inf = 0.0
!!!! and we want to remove
!!! using Steve's suggestion to do so
      h22_ave = 0.0
      h22_inf_ave = 0.0
      xyh22_ave = 0.0
      xyh22_inf_ave = 0.0
      count_ave = 0

      do i=2, imax
        localtime = timein(i)
	dt = timein(i)-timein(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime+2.*dt .GT. time_min) then

	h22_inf_ave = h22_inf_ave + d2p2T_inf(:,i)
	xyh22_inf_ave = xyh22_inf_ave + d2p2T_inf(:,i)*timein(i)

	time_ave = time_ave + timein(i)
	xxtime_ave = xxtime_ave + timein(i)*timein(i)
	count_ave = count_ave + 1
	else
	time_ave =  timein(i)
	xxtime_ave = timein(i)*timein(i)
	count_ave = 1
	end if

      end do 

	slopeh22_inf =  (xyh22_inf_ave - h22_inf_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22_inf = h22_inf_ave/count_ave - slopeh22_inf*time_ave/count_ave

!remove drift
	do i=2, imax
	print*,timein(i)
	 if(timein(i)+2.*dt.gt.time_min) then
	 d2p2T_inf(:,i) = d2p2T_inf(:,i) - ( slopeh22_inf*timein(i) + interh22_inf)
        write(134,100) tortutime(i), d2p2T_inf(1,i), d2p2T_inf(2,i), timein(i)
	 end if

	end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Now 2nd time integral, and then remove drift

      h2p2 = 0.0
      h2p2_inf = 0.0
      h22_ave = 0.0
      h22_inf_ave = 0.0
      xyh22_ave = 0.0
      xyh22_inf_ave = 0.0
      count_ave = 0

      do i=2, imax
        localtime = timein(i)
	dt = timein(i)-timein(i-1)
        ! with the raw data and the infinity extrapolation
        if (localtime .GT. time_min) then
          h2p2(:,i) = h2p2(:,i-1) &
     &               + 0.5*dt*(d2p2T(:,i-1) + d2p2T(:,i))
          h2p2_inf(:,i) = h2p2_inf(:,i-1) &
     &                   + 0.5*dt*(d2p2T_inf(:,i-1) + d2p2T_inf(:,i))
 !        end if


	h22_ave = h22_ave + h2p2(:,i)
	h22_inf_ave = h22_inf_ave + h2p2_inf(:,i)

	xyh22_ave = xyh22_ave + h2p2(:,i)*timein(i)
	xyh22_inf_ave = xyh22_inf_ave + h2p2_inf(:,i)*timein(i)

	time_ave = time_ave + timein(i)
	xxtime_ave = xxtime_ave + timein(i)*timein(i)
	count_ave = count_ave + 1
	else
	h22_ave =  h2p2(:,i)
	h22_inf_ave =  h2p2_inf(:,i)
	xyh22_ave =  h2p2(:,i)*timein(i)
	xyh22_inf_ave =  h2p2_inf(:,i)*timein(i)
	time_ave =  timein(i)
	xxtime_ave = timein(i)*timein(i)
	count_ave = 1
	end if

      end do 

	slopeh22 =  (xyh22_ave - h22_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22 = h22_ave/count_ave - slopeh22*time_ave/count_ave

	slopeh22_inf =  (xyh22_inf_ave - h22_inf_ave*time_ave/count_ave) &
     &             /(xxtime_ave - time_ave*time_ave/count_ave)

	interh22_inf = h22_inf_ave/count_ave - slopeh22_inf*time_ave/count_ave


!remove drift
	do i=2, imax
	 if(timein(i).gt.time_min) then

	 h2p2(:,i) = h2p2(:,i) - (slopeh22*timein(i) + interh22)

	 h2p2_inf(:,i) = h2p2_inf(:,i) - (slopeh22_inf*timein(i) + interh22_inf)

        write(135,100) tortutime(i), h2p2(1,i), h2p2(2,i), timein(i)
        write(136,100) tortutime(i), h2p2_inf(1,i) , &
     &                       h2p2_inf(2,i), timein(i)

	 end if
	end do

!!!!!! finally we note the phase and amplitudes above were obtained
!!! from c2p2, but that isnot the strain. so let us redo here with strains

!!! first from h2p2

	cphase = 0

	do i = 2, imax 
        localtime = timein(i)
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
     &                             HHphase(i), localtime

	end if
	end do

!!!! then with h2p2_inf
	cphase = 0

	do i = ISTART, imax 
        localtime = timein(i)
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
     &                             HHphase_inf(i), localtime

	end if
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      close(unit=110)
      close(unit=111)
      close(unit=112)
      close(unit=113)
      close(unit=114)

      close(unit=120)
      close(unit=121)
      close(unit=122)

      close(unit=130)
      close(unit=131)                 
      close(unit=132)                 
      close(unit=133)                 
      close(unit=134)                 
      close(unit=135)                 
      close(unit=136)                                              
      close(unit=137)                                              

      close(unit=210)
      close(unit=211)
      close(unit=212)
      close(unit=213)
      close(unit=214)
      close(unit=215)
      close(unit=216)

      close(unit=220)
      close(unit=221)
      close(unit=222)
      close(unit=223)
      close(unit=224)
      close(unit=225)
      close(unit=226)
      close(unit=227)
      close(unit=228)

      close(unit=310)
      close(unit=311)
      close(unit=312)
      close(unit=313)
      close(unit=314)

      close(unit=320)
      close(unit=321)
      close(unit=322)
      close(unit=323)
      close(unit=324)
      close(unit=325)
      close(unit=326)

      close(unit=330)
      close(unit=331)
      close(unit=332)
      close(unit=333)
      close(unit=334)
      close(unit=335)
      close(unit=336)
      close(unit=337)
      close(unit=338)



!-----------------------------------------
      
  100 format(5e16.7)      

      end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----the integral on a sphere-------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------
!-----basic integral on a sphere--------------------
!---------------------------------------------------

      real*8 function surface_integral(u,MAXMEM,N)     
      
      integer :: N, iphi,itheta,MAXMEM
      real*8  :: u(MAXMEM)
      real*8  :: phi,theta,pi,integral

      pi     = dacos(-1.0d0)
      h      = pi/N
      h2     = h*h
           
      do iphi = 1, 2*N+1
        phi = (iphi-1)*h
        do itheta = 1, N+1
          theta = (itheta-1)*h   
          index = (itheta-1)*(2*N+1) + iphi

	  integral = integral + h2*dsin(theta)*u(index)

      end do
      end do      	

      surface_integral = integral

      return 
      end function surface_integral

!---------------------------------------------------
!-----better integral on a sphere-------------------
!---------------------------------------------------

      real*8 function surface_integral2(u,MAXMEM,N)     
      
      integer :: N, iphi,itheta,MAXMEM
      real*8  :: u(MAXMEM)
      real*8  :: phi,theta,pi,integral

      pi       = dacos(-1.0d0)
      h        = pi/N
      h2       = h*h
      integral = 0.0
  
      ! corner points
      iphi = 1
      itheta = 1
      phi = (iphi-1)*h
      theta = (itheta-1)*h   
      index = (itheta-1)*(2*N+1) + iphi
      integral = integral + dsin(theta)*u(index)
     
      iphi = 2*N+1
      itheta = N+1
      phi = (iphi-1)*h
      theta = (itheta-1)*h   
      index = (itheta-1)*(2*N+1) + iphi
      integral = integral + dsin(theta)*u(index)

      iphi = 1
      itheta = N+1
      phi = (iphi-1)*h
      theta = (itheta-1)*h   
      index = (itheta-1)*(2*N+1) + iphi
      integral = integral + dsin(theta)*u(index)

      iphi = 2*N+1
      itheta = 1
      phi = (iphi-1)*h
      theta = (itheta-1)*h   
      index = (itheta-1)*(2*N+1) + iphi
      integral = integral + dsin(theta)*u(index)
      
   	
      ! boundary points
      iphi = 1
      phi = (iphi-1)*h
      do itheta = 2, N
          theta = (itheta-1)*h   
          index = (itheta-1)*(2*N+1) + iphi
	  integral = integral + 2.0*dsin(theta)*u(index)
      end do

      iphi = 2*N+1
      phi = (iphi-1)*h
      do itheta = 2, N
        theta = (itheta-1)*h   
        index = (itheta-1)*(2*N+1) + iphi
	integral = integral + 2.0*dsin(theta)*u(index)
      end do

     itheta = 1
     theta = (itheta-1)*h    
      do iphi = 2, 2*N
        phi = (iphi-1)*h
        index = (itheta-1)*(2*N+1) + iphi
	integral = integral + 2.0*dsin(theta)*u(index)
      end do      	

     itheta = N+1
     theta = (itheta-1)*h    
      do iphi = 2, 2*N
        phi = (iphi-1)*h
        index = (itheta-1)*(2*N+1) + iphi
	integral = integral + 2.0*dsin(theta)*u(index)
      end do      	


      ! interior points
      do iphi = 2, 2*N
        phi = (iphi-1)*h
        do itheta = 2, N
          theta = (itheta-1)*h   
          index = (itheta-1)*(2*N+1) + iphi
	  integral = integral + 4.0*dsin(theta)*u(index)
      end do
      end do   

      surface_integral2 = 0.25*h2*integral
      return 

      end function surface_integral2


