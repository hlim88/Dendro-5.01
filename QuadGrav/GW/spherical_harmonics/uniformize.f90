! Data at different grid resolutions
! have non-coincident times, for richardson
! extrapolation, we would want things at the 
! same time. let's uniformize things then here
! we will call things rp4 and ip4 simply as 
! variables.

	program main	
	implicit none
	real*8 :: rp4(1:100000), ip4(1:100000), t1(1:100000) 
	integer :: imax, i, j, KK
	real*8 :: dtspace, tout,tini, rp4out, ip4out, trash
	real*8 :: fac1, fac2, fac3, fac4, tfinish
	character*32 :: filein(5), fileout(5)
	integer :: fileloop

      character*5     name, sarg
      external        sarg 

      real*8          r8arg
      external        r8arg

!	print*,'name of file to uniformize'
!	read*, filein	
!	open (unit = 8, file = filein)	

	print*,'what time space to use?, and when to start and finish'
	read*, dtspace, tini, tfinish


       if (iargc() .lt. 1) then
          write(*,*) '*************************************************'
          write(*,*) '* uniformize: get output on uniform time grid   *'
          write(*,*) '* from t=tini, to t=tfinish with spacing dt of  *'
          write(*,*) '* strain,straininf, amp/phase from c2p2, h,hinf *'
          write(*,*) '* Usage:                                        *'
          write(*,*) '*************************************************'
          write(*,*) '  uniformize <surface letter> <tini> <tfinish> <dt>'
          write(*,*) '   E.g. uniformize B 200 4000 1.5 '
          write(*,*) '                *****************'
          stop
       end if
    
       name = sarg( 1,'c')
       tini= r8arg(2,1.d0)
       tfinish = r8arg(3,140.d0)
       dtspace = r8arg(4,1.d0)


	filein(1) = 'h2p2'//name(:1)//'.dat'
	filein(2) = 'h2p2_inf'//name(:1)//'.dat'
	filein(3) = 'amp_and_phase'//name(:1)//'.dat'
	filein(4) = 'HHamp_and_phase'//name(:1)//'.dat'
	filein(5) = 'HHamp_and_phase_inf'//name(:1)//'.dat'


	fileout(1) = 'h2p2_unif'//name(:1)//'.dat'
	fileout(2) = 'h2p2_inf_unif'//name(:1)//'.dat'
	fileout(3) = 'amp_and_phase_unif'//name(:1)//'.dat'
	fileout(4) = 'HHamp_and_phase_unif'//name(:1)//'.dat'
	fileout(5) = 'HHamp_and_phase_inf_unif'//name(:1)//'.dat'


	DO KK=1,5

	open (unit = 4+KK, file = filein(KK))	
	open (unit = 14+KK, file = fileout(KK))
		

	j =0
	 do i=1, 100000
	 j = j+1
	 read (4+KK,*,end=10) t1(j), rp4(j), ip4(j),trash
	  if(j.gt.1.and.t1(j).eq.t1(j-1)) j=j-1
	end do
10 	continue
	 imax = j - 1
	
	print*, imax

	tout = -1000000.

!ok, now will try to get something come out uniform
!with a spacing of dt
!	 write(13,*) t1(1), rp4(1), ip4(1)

	do j = 1, 100000
	 if(tout.gt.t1(imax).or.tout.gt.tfinish) exit
	tout = tini + j*dtspace

	  do i=1, imax
	   if(t1(i).gt.tout) exit
	  end do


	   if(i.le.2.or.i.ge.imax-2) then
	   fac1 = (tout-t1(i))/(t1(i-1)-t1(i))
	   fac2 = (tout-t1(i-1))/(t1(i)-t1(i-1))

	   rp4out = rp4(i-1)* fac1 + rp4(i)* fac2
	   ip4out = ip4(i-1)* fac1 + ip4(i)* fac2

	   else

	   fac1 = (tout- t1(i-1))*(tout-   t1(i))*(tout-   t1(i+1)) &
     &        /((t1(i-2)-t1(i-1))*(t1(i-2)-t1(i))*(t1(i-2)-t1(i+1)))

	   fac2 = (tout- t1(i-2))*(tout-   t1(i))*(tout-   t1(i+1)) &
     &        /((t1(i-1)-t1(i-2))*(t1(i-1)-t1(i))*(t1(i-1)-t1(i+1))) 

	   fac3 = (tout- t1(i-2))*(tout- t1(i-1))*(tout-   t1(i+1)) &
     &          /((t1(i)-t1(i-2))*(t1(i)-t1(i-1))*(t1(i)  -t1(i+1)))

	   fac4 = (tout- t1(i-2))*(tout-   t1(i-1))*(tout-   t1(i)) &
     &        /((t1(i+1)-t1(i-2))*(t1(i+1)-t1(i-1))*(t1(i+1)-t1(i))) 

	   rp4out = rp4(i-2) * fac1 &
     &	          + rp4(i-1) * fac2 &
     &	          + rp4(i)   * fac3 &
     &	          + rp4(i+1) * fac4

	   ip4out = ip4(i-2) * fac1 &
     &	          + ip4(i-1) * fac2 &
     &	          + ip4(i)   * fac3 &
     &	          + ip4(i+1) * fac4  

	   end if
	 write(14+KK,100) tout, rp4out, ip4out
!	   print*,t1(i),i,tout, t1(i-1),rp4out, ip4out,rp4(i),ip4(i)
!	   STOP
	end do

	close(unit=14+KK)

	END DO

100 format(3e16.7)      
	end
