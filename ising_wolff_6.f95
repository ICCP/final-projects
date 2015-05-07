program ising_wolff

IMPLICIT NONE

integer, parameter :: L=100
integer, parameter :: N=L*L
integer, parameter :: out_unit=20
integer, parameter :: out_unit2=10
integer, parameter :: out_unit3=30
integer, parameter :: out_unit4=40
integer, parameter :: cluster_flips=100000
integer, dimension(N,5) :: neighbors
integer, dimension(N) :: spins, queue
integer :: i,x,initial_spin,T_step,j,every_hundred
integer :: mCount,avgCount,k
real, dimension(20) :: m_bins
real, dimension(20) :: e_bins

real :: T,M,s,beta,avgM,diffSqSum,diff,stdev,diffsq
real :: m_avg,m_sq_avg,chi,E,cv,e_avg,avgE,e_sq_avg

character(len=10) ::  out_file
character(len=7) ::  out_file_base
character(len=3) ::  out_file_no

call seed_number_generator()
call get_random_number(s)
call initialize_positions_neighbors(neighbors)

out_file_base="Results"
open (unit=out_unit2,file="mvsT",action="write",status="replace")
open (unit=out_unit3,file="chi",action="write",status="replace")
open (unit=out_unit4,file="cv",action="write",status="replace")

do T_step=1,2 
    T=1.0*T_step
    call randomize_spins(spins)
    write(*,*) T
    write(out_file_no,"(F3.1)") T
    out_file=out_file_base//out_file_no
    open (unit=out_unit,file=out_file,action="write",status="replace")

   
    beta=(1/T)*1.0

   m_avg=0
   e_avg=0
   m_sq_avg=0
   e_sq_avg=0
   avgM=0
   every_hundred=0
   avgCount=1
   mCount=1 
   avgE=0   
   diffSqSum=0
	
    do i=1,cluster_flips
        queue(:)=0
        x=select_spin_to_flip()
        initial_spin=spins(x)
        spins(x)=(-1)*spins(x)
        queue(1)=x

        call grow_cluster(spins,neighbors,initial_spin,queue,beta)
        if (every_hundred==100) then
            every_hundred=0
            M=0
			E=0
            do j=1,N
                M=M+spins(j)
				do k=2,5
				    E=E+(spins(j)*spins(neighbors(j,k)))
				end do
            end do
			E=((-1)*E)/(2.0*N)
            M=M/N
			if (M<0) then
			    M=(-1)*M
			end if
            write(out_unit,*) i, M

! If the iteration is 10 thousand or above use to calculate average
	   if (cluster_flips>9999) then
	       m_avg=m_avg+m
		   e_avg=e_avg+E
! Calculate average m every 10,000 spin flip iterations
		   if (avgCount==50) then
			   m_avg=m_avg/avgCount
			   e_avg=e_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   e_bins(mCount)=e_avg
			   mCount=mCount+1
			   m_avg=0
			   e_avg=0
		   end if
		   avgCount=avgCount+1
	    end if
        end if
        every_hundred=every_hundred+1

		
    end do

    close(out_unit)
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
	   avgE=avgE+e_bins(i)
	   m_sq_avg=m_sq_avg+(m_bins(i)*m_bins(i)*1.0)
	   e_sq_avg=e_sq_avg+(e_bins(i)*e_bins(i)*1.0)
   end do
   avgM=avgM/mCount
   m_sq_avg=m_sq_avg/mCount
   chi=(1.0/T)*(m_sq_avg-(avgM*avgM))
   
   avgE=avgE/mCount
   e_sq_avg=e_sq_avg/mCount
   cv=(1.0/(T**2))*(e_sq_avg-(avgE*avgE))
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   write (out_unit3,*) T, chi
   write (out_unit4,*) T, cv
	
end do


do T_step=1,9 
    T=2.0+(T_step/10.0)
    call randomize_spins(spins)
    write(*,*) T
    write(out_file_no,"(F3.1)") T
    out_file=out_file_base//out_file_no
    open (unit=out_unit,file=out_file,action="write",status="replace")

   
    beta=(1/T)*1.0

   m_avg=0
   e_avg=0
   m_sq_avg=0
   e_sq_avg=0
   avgM=0
   every_hundred=0
   avgCount=1
   mCount=1 
   avgE=0   
   diffSqSum=0
	
    do i=1,cluster_flips
        queue(:)=0
        x=select_spin_to_flip()
        initial_spin=spins(x)
        spins(x)=(-1)*spins(x)
        queue(1)=x

        call grow_cluster(spins,neighbors,initial_spin,queue,beta)
        if (every_hundred==100) then
            every_hundred=0
            M=0
			E=0
            do j=1,N
                M=M+spins(j)
				do k=2,5
				    E=E+(spins(j)*spins(neighbors(j,k)))
				end do
            end do
			E=((-1)*E)/(2.0*N)
            M=M/N
			if (M<0) then
			    M=(-1)*M
			end if
            write(out_unit,*) i, M

! If the iteration is 10 thousand or above use to calculate average
	   if (cluster_flips>9999) then
	       m_avg=m_avg+m
		   e_avg=e_avg+E
! Calculate average m every 10,000 spin flip iterations
		   if (avgCount==50) then
			   m_avg=m_avg/avgCount
			   e_avg=e_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   e_bins(mCount)=e_avg
			   mCount=mCount+1
			   m_avg=0
			   e_avg=0
		   end if
		   avgCount=avgCount+1
	   end if
        end if
        every_hundred=every_hundred+1

		
    end do

    close(out_unit)
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
	   avgE=avgE+e_bins(i)
	   m_sq_avg=m_sq_avg+(m_bins(i)*m_bins(i)*1.0)
	   e_sq_avg=e_sq_avg+(e_bins(i)*e_bins(i)*1.0)
   end do
   avgM=avgM/mCount
   m_sq_avg=m_sq_avg/mCount
   chi=(1.0/T)*(m_sq_avg-(avgM*avgM))
   
   avgE=avgE/mCount
   e_sq_avg=e_sq_avg/mCount
   cv=(1.0/(T**2))*(e_sq_avg-(avgE*avgE))
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   write (out_unit3,*) T, chi
   write (out_unit4,*) T, cv
	
end do


do T_step=3,4 
    T=1.0*T_step
    call randomize_spins(spins)
    write(*,*) T
    write(out_file_no,"(F3.1)") T
    out_file=out_file_base//out_file_no
    open (unit=out_unit,file=out_file,action="write",status="replace")

   
    beta=(1/T)*1.0

   m_avg=0
   e_avg=0
   m_sq_avg=0
   e_sq_avg=0
   avgM=0
   every_hundred=0
   avgCount=1
   mCount=1 
   avgE=0   
   diffSqSum=0
	
    do i=1,cluster_flips
        queue(:)=0
        x=select_spin_to_flip()
        initial_spin=spins(x)
        spins(x)=(-1)*spins(x)
        queue(1)=x

        call grow_cluster(spins,neighbors,initial_spin,queue,beta)
        if (every_hundred==100) then
            every_hundred=0
            M=0
			E=0
            do j=1,N
                M=M+spins(j)
				do k=2,5
				    E=E+(spins(j)*spins(neighbors(j,k)))
				end do
            end do
			E=((-1)*E)/(2.0*N)
            M=M/N
			if (M<0) then
			    M=(-1)*M
			end if
            write(out_unit,*) i, M

! If the iteration is 10 thousand or above use to calculate average
	   if (cluster_flips>9999) then
	       m_avg=m_avg+m
		   e_avg=e_avg+E
! Calculate average m every 10,000 spin flip iterations
		   if (avgCount==50) then
			   m_avg=m_avg/avgCount
			   e_avg=e_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   e_bins(mCount)=e_avg
			   mCount=mCount+1
			   m_avg=0
			   e_avg=0
		   end if
		   avgCount=avgCount+1
	   end if
        end if
        every_hundred=every_hundred+1

		
    end do

    close(out_unit)
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
	   avgE=avgE+e_bins(i)
	   m_sq_avg=m_sq_avg+(m_bins(i)*m_bins(i)*1.0)
	   e_sq_avg=e_sq_avg+(e_bins(i)*e_bins(i)*1.0)
   end do
   avgM=avgM/mCount
   m_sq_avg=m_sq_avg/mCount
   chi=(1.0/T)*(m_sq_avg-(avgM*avgM))
   
   avgE=avgE/mCount
   e_sq_avg=e_sq_avg/mCount
   cv=(1.0/(T**2))*(e_sq_avg-(avgE*avgE))
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   write (out_unit3,*) T, chi
   write (out_unit4,*) T, cv
	
end do
    close(out_unit2)
    close(out_unit3)
	close(out_unit4)
contains

Subroutine initialize_positions_neighbors(neighbors)
    IMPLICIT NONE
    integer, dimension(N,5) :: neighbors
    integer :: i
    real :: x
	
    do i=1,N
	    
        neighbors(i,1)=i
        
        if ((((i*1.0)/(L*1.0))/floor((1.0*i)/(1.0*L)))==1.0) then
            neighbors(i,2)=i-(L-1)
        else
            neighbors(i,2)=i+1
        end if
		
        if ((floor((1.0*i)/(1.0*L))==(1.0*(L-1)) .OR. ((i*1.0)/(L*L*1.0))*1.0==1.0) .AND. (1.0*i)/(1.0*L)/=(1.0*(L-1))) then
            neighbors(i,3)=i-(L*(L-1))
        else
            neighbors(i,3)=i+L
        end if

        x=(((i-1)*1.0)/(L*1.0))/floor((1.0*(i-1))/(1.0*L))
        if ((x==1.0) .or. (i==1)) then
           neighbors(i,4)=i+(L-1)
        else
            neighbors(i,4)=i-1
        end if

		
        if (floor((1.0*i)/(1.0*L))==0 .OR. ((i*1.0)/(L*1.0))*1.0==1.0 ) then
            neighbors(i,5)=i+(L*(L-1))
        else
            neighbors(i,5)=i-L
        end if
		
    end do
end subroutine

Subroutine seed_number_generator()
    Implicit None
    integer :: i_seed
    integer, dimension(:), ALLOCATABLE :: Set_seed
    integer, dimension(1:8) :: dateSet_seed
    
    CALL RANDOM_SEED(size=i_seed)
    ALLOCATE(Set_seed(1:i_seed))
    CALL RANDOM_SEED(get=Set_seed)
    CALL DATE_AND_TIME(values=dateSet_seed)
    Set_seed(i_seed)=dateSet_seed(8)
    Set_seed(1)=dateSet_seed(8)*dateSet_seed(6)
    CALL RANDOM_SEED(put=Set_seed)
    DEALLOCATE(Set_seed)
End Subroutine

Subroutine get_random_number(s)
    IMPLICIT NONE
    real :: r,s

    CALL RANDOM_NUMBER(r)
    s = r
End Subroutine

function select_spin_to_flip()
   Implicit None
   integer :: select_spin_to_flip
   real :: r
   
   call get_random_number(r)
   select_spin_to_flip = (r*N)+1
   if (select_spin_to_flip==(N+1)) then
       select_spin_to_flip = (N*r)+1
   end if
end function

Subroutine initialize_spins(spins)
    Implicit None
	
    integer, dimension(N) :: spins
    integer :: i
	
    do i=1,N
        spins(i)=1
    end do
end subroutine

Subroutine randomize_spins(spins)
    Implicit None
	
    integer, dimension(N) :: spins
    real :: s
    integer :: i
	
    do i=1,N
        call get_random_number(s)
        if (s>0.5) then
            spins(i)=1
        else
            spins(i)=-1
        end if
    end do
end subroutine

subroutine grow_cluster(spins,neighbors,initial_spin,queue,beta)
    Implicit None
	
    integer, dimension(N) :: spins
    integer :: initial_spin,start_index,end_index,i
    integer, dimension(N) :: queue
    integer, dimension(N,5) :: neighbors
    integer :: spin_could_flip
	
    real :: r,prob,beta
	
    start_index=0
    end_index=1
	
	
100 start_index=start_index+1
    do i=2,5
        if (initial_spin==spins(neighbors(queue(start_index),i))) then

            spin_could_flip=neighbors(queue(start_index),i)
        
            call get_random_number(r)
            prob=1-exp(((-2)*beta))
!		    write(*,*) prob 	
            if(r<=prob) then
                spins(spin_could_flip)=initial_spin*(-1)
                end_index=end_index+1
                queue(end_index)=spin_could_flip
           end if
        end if
    end do
    if (start_index/=end_index) then
        go to 100
    end if
	    
end subroutine

end program ising_wolff