! Data at different grid resolutions
! have non-coincident times, for richardson
! extrapolation, we would want things at the 
! same time. let's uniformize things then here
! we will call things rp4 and ip4 simply as 
! variables.


	program main	
	implicit none
	real*8 :: rp4(3,1:100000), ip4(3,1:100000), t1(1:100000) 
	real*8 :: rp4UNI(3,1:100000), ip4UNI(3,1:100000), t1UNI(1:100000) 

	real*8 :: AA(1:100000), PHI(1:100000)

	integer :: imax(3), i, j, KK, jmax
	real*8 :: dtspace, tout,tend,tini, rp4out, ip4out, trash, h(3)
	real*8 :: fac1, fac2, fac3, fac4, tfinish
	character*32 :: filein(5), fileout(5)
	integer :: fileloop


	do i=1,3
	  print*,'give me spacing=',i
	  read*, h(i)
	end do
	
	print*,'tini,tout,dtspace'
	read*, tini, tfinish, dtspace


	filein(1) = 'file1.dat'
	filein(2) = 'file2.dat'
	filein(3) = 'file3.dat'

	fileout(1) = 'phase_RE.dat'
	fileout(2) = 'amp_RE.dat'
	fileout(3) = 'strain_RE.dat'


	DO KK=1,3
	print*,'doing file', KK

	open (unit = 20+KK, file = filein(KK))	
	open (unit = 14+KK, file = fileout(KK))
		
	print*,'opened files', KK

	j =0
	 do i=1, 100000
	 j = j+1
	 read (20+KK,*,end=10) t1(j), rp4(KK,j), ip4(KK,j),trash
	  if(j.gt.1.and.t1(j).eq.t1(j-1)) j=j-1
	end do
10 	continue
	 imax(KK) = j - 1
	
	print*, imax(KK)

!ok, now will try to get something come out uniform
!with a spacing of dt
!	 write(13,*) t1(1), rp4(1), ip4(1)
	tout = -10000

	jmax =int( (tfinish-tini)/dtspace )

	do j = 1, jmax
	tout = tini + j*dtspace
	t1uni(j) = tout

	  do i=1, imax(KK)
	   if(t1(i).gt.tout) exit
	  end do


	   if(i.le.2.or.i.ge.imax(KK)-2) then
	   fac1 = (tout-t1(i))/(t1(i-1)-t1(i))
	   fac2 = (tout-t1(i-1))/(t1(i)-t1(i-1))

	   rp4out = rp4(KK,i-1)* fac1 + rp4(KK,i)* fac2
	   ip4out = ip4(KK,i-1)* fac1 + ip4(KK,i)* fac2

	   else

	   fac1 = (tout- t1(i-1))*(tout-   t1(i))*(tout-   t1(i+1)) &
     &        /((t1(i-2)-t1(i-1))*(t1(i-2)-t1(i))*(t1(i-2)-t1(i+1)))

	   fac2 = (tout- t1(i-2))*(tout-   t1(i))*(tout-   t1(i+1)) &
     &        /((t1(i-1)-t1(i-2))*(t1(i-1)-t1(i))*(t1(i-1)-t1(i+1))) 

	   fac3 = (tout- t1(i-2))*(tout- t1(i-1))*(tout-   t1(i+1)) &
     &          /((t1(i)-t1(i-2))*(t1(i)-t1(i-1))*(t1(i)  -t1(i+1)))

	   fac4 = (tout- t1(i-2))*(tout-   t1(i-1))*(tout-   t1(i)) &
     &        /((t1(i+1)-t1(i-2))*(t1(i+1)-t1(i-1))*(t1(i+1)-t1(i))) 

	   rp4out = rp4(KK,i-2) * fac1 &
     &	          + rp4(KK,i-1) * fac2 &
     &	          + rp4(KK,i)   * fac3 &
     &	          + rp4(KK,i+1) * fac4

	   ip4out = ip4(KK,i-2) * fac1 &
     &	          + ip4(KK,i-1) * fac2 &
     &	          + ip4(KK,i)   * fac3 &
     &	          + ip4(KK,i+1) * fac4  

	  rp4UNI(KK,j) = rp4out
	  ip4UNI(KK,j) = ip4out
	  t1UNI(j) = tout

	   end if

!	   print*,t1(i),i,tout, t1(i-1),rp4out, ip4out,rp4(i),ip4(i)
!	   STOP
	end do

	END DO

	jmax=(tfinish-tini)/dtspace
	print*,'upto jmax', jmax,j


!OK the above gave me things at the right times, can now time by time
!do an extrapolation assuming a quadratic expression
!i need only output the asymptote to h->0 value (y=a+bx+c^2... ie, just a)

	do i = 1, jmax

!amplitude is in comumn rp4
	AA(i) = rp4UNI(3,i)*h(1)*h(2)*(h(1)-h(2))+h(3)*(rp4UNI(1,i)*h(2)*(h(2)-h(3))+ &
     &                              rp4UNI(2,i)*h(1)*(h(3)-h(1)) )
	AA(i) = AA(i)/( (h(1)-h(2)) * (h(1) - h(3)) * (h(2)-h(3)) )	

!amplitude is in comumn rp4
	PHI(i) = ip4UNI(3,i)*h(1)*h(2)*(h(1)-h(2))+h(3)*(ip4UNI(1,i)*h(2)*(h(2)-h(3))+ &
     &                               ip4UNI(2,i)*h(1)*(h(3)-h(1)) )
	PHI(i) = PHI(i)/( (h(1)-h(2)) * (h(1) - h(3)) * (h(2)-h(3)) )	


	 write(15,100) t1UNI(i), PHI(i), AA(i)
	 write(16,100) t1UNI(i), AA(i), PHI(i)
	 write(17,100) t1UNI(i), AA(i)*cos(PHI(i)), AA(i)*sin(PHI(i))
	print*, 'out',i, t1UNI(i), AA(i)*cos(PHI(i))
	write(31,*), t1UNI(i),rp4UNI(1,i),ip4UNI(1,i)
	write(32,*), t1UNI(i),rp4UNI(2,i),ip4UNI(2,i)
	write(33,*), t1UNI(i),rp4UNI(3,i),ip4UNI(3,i)
	end do


	close(unit=15)
	close(unit=16)
	close(unit=17)

100 format(3e16.7)      
	end
