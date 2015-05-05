   real function gauss(m,s,x)
   real*8 m,s,x,const,pi
   
   pi = 4.d0*atan(1.d0)
   
   !const = 1.d0/(sqrt(2.d0*pi*s**2))
   const = 1.d0
   
   gauss = const*exp(-(x-m)**2/(2.d0*s**2))
   
   return
   end
      
   subroutine calcgrad(parm,fparm,nf,data,nl,grad,mean,sigma)
   use constants
   use F90library
   implicit none
   integer, intent(in) :: nf,nl,fparm(nf)
   real*8, intent(inout) :: parm(12),grad(nf)
   real*8, intent(in) :: data(nl),mean(12),sigma(12)
   integer i,j
   real*8 jac(nf,nf),chip,chim,step,delta,gaup,gaum,gauss
   real*8, allocatable :: d(:)
   real*8, allocatable :: e(:)
   
   !do i=1,nf
   !   print*, fparm(i)
   !enddo 
   
   ! for each parameter that's being fit, calculate chi^2 at \pm step
   ! weight the chi^2 by a gaussian to keep the parameters physical
   step = 0.1d0
   do i=1,nf
      ! origianal parm + step, calculate chi square value
      parm(fparm(i)) = step + parm(fparm(i))
      call inputfile(parm)
      call system("fresco < newfit.in > newfit.out")
      call chicalc(data,nl,chip)
      gaup = gauss(mean(fparm(i)),sigma(fparm(i)),parm(fparm(i)))
      !print *, gaup
      
      ! original parm - step, calculate chi square value
      parm(fparm(i)) = parm(fparm(i)) - 2*step
      call inputfile(parm)
      call system("fresco < newfit.in > newfit.out")
      call chicalc(data,nl,chim)
      gaum = gauss(mean(fparm(i)),sigma(fparm(i)),parm(fparm(i)))
      !print *, gaum
      
      ! put parameter back at it's original value and move to next
      parm(fparm(i)) = parm(fparm(i)) + step
      !print *, parm(fparm(i)),chip,chim
      
      ! calculate derivative
      grad(i) = (chip*gaup - chim*gaum)/(2*step)
      !print *, grad(i),gaup,gaum,chip,chim
   enddo 
   
   ! form the approximate jacobian
   !do i=1,nf
   !   do j=i,nf
   !      jac(i,j) = grad(i)*grad(j)
   !      jac(j,i) = jac(i,j)
	 !print *, jac(i,j),i,j
   !   enddo 
   !enddo 
   
   ! diagonalize the jacobian and return it in grad
   ! uses a diagonalization routine recieved from M. Hjorth-Jensen
   ! define step size
   !delta = 0.1
   !allocate (e(nf))
   !allocate (d(nf))
   !call tred2(jac,nf,d,e)
   !call tqli(d,e,nf,jac)
   !do i=1,nf
   !   grad(i) = d(i)
   !   print *, grad(i)
      !grad(i) = grad(i) + delta
   !   grad(i) = sqrt(grad(i))
   !enddo 
   
   !do i=1,nf
   !   print *, (jac(i,j),j=1,nf)
   !enddo 
   
   ! matrix inversion doesn't work - talk to Jenni about method
   ! take the inverse of jac after going through tqli
   !call matinv(jac,nf,1.d0)
   
   ! reconstruct new parameters
   !parm = matmul(jac,grad)
   !do i=1,nf
   !  print *, parm(i)
   !enddo 
   !deallocate(e)
   !deallocate(d)
   
   end subroutine calcgrad