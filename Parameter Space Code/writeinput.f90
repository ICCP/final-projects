   subroutine inputfile(parm)
   implicit none
   real*8, intent(in) :: parm(12)
   character (len=100) line
   integer i
   
   ! make new input file - replace current parameters with fitting ones
   open(unit=4,file="pull_ca48-ddfit102.in")
   open(unit=7,file="newfit.in")
   ! first 15 lines are the same
   do i=1,15
      read(4,'(A)') line
      write(7,'(A)') line
   enddo 
   ! replace parameters in three lines
   read(4,'(A)') line
   write(line,' (" &pot kp= 1 type= 1  p(1:6)= " F7.3 "  " F8.5 "  " F9.6 "  " F7.3 "  " F8.5 "  " F9.6 " /")') parm(1),parm(2),parm(3),parm(4),parm(5),parm(6)
   write(7,'(A)') line
   read(4,'(A)') line
   write(line,' (" &pot kp=1  type= 2  p(1:6)= 0.0  0.0	 0.0  "F7.3"  "F8.5"  "F9.6" /")') parm(7),parm(8),parm(9)
   write(7,'(A)') line
   read(4,'(A)') line
   write(line,' (" &pot kp= 1 type= 3  p(1:6)=  "F7.3"  "F8.5"  "F9.6"  0.000  0.000  0.000 /")') parm(10),parm(11),parm(12)
   write(7,'(A)') line
   ! last three lines are the same
   do i=1,3
      read(4,'(A)') line
      write(7,'(A)') line
   enddo 
   close(4)
   close(7)
   
   end subroutine inputfile