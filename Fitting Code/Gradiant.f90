   program gradmin
   implicit none
   real*8 parm(12),initparm(12),chi,chinew,step,mean(12),sigma(12),w(12)
   integer nlines,io,i,nfit,j,k
   real*8, allocatable :: data(:,:)
   integer, allocatable :: fitparm(:)
   real*8, allocatable :: grad(:)
   
   ! define the step size
   step = 0.0001d0
   
   ! read in inital parameters
   nfit = 0
   open(unit=1,file="fittingparm.txt")
   do i=1,12
      read(1,*) parm(i)
      if (parm(i)==0) then
         nfit = nfit + 1
      end if 
   enddo 
   close(1)
   !print *, "nfit", nfit
   
   allocate (fitparm(nfit))
   
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
   
   ! read in the initial parameters
   j = 1
   open(unit=3,file="initialparm.txt")
   do i=1,12
      read(3,*) initparm(i)
      if (parm(i)==0) then
         parm(i) = initparm(i)
      else 
         fitparm(j) = i
	 !print *,fitparm(j),j
	 j = j + 1
      end if 
      !print *, parm(i)
   enddo 
   
   ! read in the parameter bounds
   open(unit=8,file="limits.txt")
   do i=1,12
      read(8,*) mean(i),sigma(i),w(i)
   enddo 

   ! count the number of iterations
   k=0
   do 
   
   ! replace values in the input file
   call inputfile(parm)
   
   ! run fresco
   call system("fresco < newfit.in > newfit.out")
   
   ! calculate chi square
   call chicalc(data,nlines,chi)
   
   ! find the gradiant of chi square function
   allocate (grad(nfit))
   call calcgrad(parm,fitparm,nfit,data,nlines,grad,mean,sigma,w)
   
   ! move the points - gradient descent
   do i=1,12
   do j=1,nfit
      if (fitparm(j)==i) then
         parm(i) = parm(i) - step*grad(j)
      end if
   enddo 
   enddo 
   
   deallocate (grad)
   
   ! remake input file and run fresco
   call inputfile(parm)
   call system("fresco < newfit.in > newfit.out")
   
   ! calculate chi square
   call chicalc(data,nlines,chinew)
   print *, k, chinew
   
   ! need a condition for ending - if |newchi - chi| < some value
   if (abs(chinew - chi) < 0.1d0) then
      exit
   end if 

   k = k + 1
   ! if a minimium isn't found, don't let the code run infinitely
   if (k>500) exit
   enddo
   
   ! print final values
   do i=1,12
      print *, "#",parm(i)
   enddo 
   
   end program gradmin