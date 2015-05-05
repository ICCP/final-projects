   ! putting this aside momentarily
   program gridding
   implicit none
   real*8 parm(6),step,pzero,val(6),conversion(6,2)
   integer i,j,k,m,n,nstep
   
   ! initialize parm to zero
   parm = 0
   ! print *, (parm(k),k=1,6)
   
   ! define a step size
   step = 0.05
   ! number of steps to take
   nstep = int(1/step)
   
   ! define the conversion between grid and parameters
   conversion(1,1) = 170d0 - 50d0
   conversion(1,2) = 50d0
   conversion(2,1) = 1.4d0 - 0.8d0
   conversion(2,2) = 0.8d0
   conversion(3,1) = 0.8d0 - 0.4d0
   conversion(3,2) = 0.4d0
   conversion(4,1) = 60d0 - 10d0
   conversion(4,2) = 10d0
   conversion(5,1) = 1.4d0 - 0.8d0
   conversion(5,2) = 0.8d0
   conversion(6,1) = 0.5d0 - 0.1d0
   conversion(6,2) = 0.1d0
   
   !print *, (conversion(k,2),k=1,6)
   
   open(unit=4,file="parmgrid.txt")
   !write(4,20) (conversion(k,2),k=1,6)
   pzero = 0
      
   ! set the rest of the parameters to some fixed value
   ! grid the rest
   do m=0,10
      pzero = m*0.1d0
      
      ! create the grid
      ! keep everything fixed except pairs of parameters
      do i=1,6
         do j=i+1,6
	    write(4,*) i,j
	    parm(i) = 0.d0
	    do k=1,nstep
	       parm(j) = 0.d0
	       do while (parm(j) <= 1.d0) 
	          write(4,20) (parm(n)*conversion(n,1)+conversion(n,2),n=1,6)
		  parm(j) = parm(j) + step
	       enddo 
	       parm(i) = parm(i) + step
	       parm(j) = pzero
	    enddo 
	 enddo 
	 parm(i) = pzero
      enddo 
   
   enddo 
   
   close(4)
   
20 format(4x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)
   
   end program gridding