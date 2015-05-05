   subroutine chicalc(data,n,chi)
   implicit none
   integer, intent(in) :: n
   real*8, intent(in) :: data(3,n)
   real*8, intent(inout) :: chi
   real*8 res(n),angle,xsth
   integer i,j
   
   ! read in calculation from fresco
   open(unit=8,file="fort.16")
   ! first 13 lines are plotting options
   do i=1,13
      read(8,*)
   enddo 
   ! next 181 lines contain data
   j=1
   do i=0,180
      read(8,*) angle, xsth
      ! if the angles between theory and experiment/data are the same
      ! calculate the residual
      if (idnint(data(1,j))==i) then
         res(j) = (xsth-data(2,j))/data(3,j)
	 j = j + 1
      end if 
   enddo 
   close(8)
   
   ! calculate chisquare value
   chi = 0
   do i=1,n
      chi = chi + res(i)**2
   enddo 
   
   end subroutine chicalc