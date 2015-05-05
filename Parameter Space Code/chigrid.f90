   program chigrid
   implicit none
   integer nlines,io,i,j,a,b,k
   real*8 parm(6),chi
   real*8, allocatable :: data(:,:)
   
   ! read in the data points
   open(unit=2,file="xexp.txt")
   nlines = 0
   do 
      read(2,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
   enddo 
   
   rewind(2)
   
   allocate (data(3,nlines))
   
   do i=1,nlines
      read(2,*) data(1,i),data(2,i),data(3,i)
   enddo
   
   close(2)
   
   ! read in parameter grid
   open(unit=3,file="parmgrid.txt")
   open(unit=40,file="chigrid.txt")
   do i=1,165
      read(3,*) a,b
      write(40,*) a,b
      !print *, a,b
      do j=1,380
         read(3,*) parm(1),parm(2),parm(3),parm(4),parm(5),parm(6)
	 call inputfile(parm)
	 call system("fresco < newfit.in > newfit.out")
	 call chicalc(data,nlines,chi)
	 write(40,20) parm(1),parm(2),parm(3),parm(4),parm(5),parm(6),chi
20 format(4x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.6)
      enddo 
   enddo
   
   close(3)
   close(40)
   
   end program chigrid