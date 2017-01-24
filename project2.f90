program project1

implicit none

integer:: itr,b,b_delta, j,Nmax, n, m, t,i,it,k
real:: z,x,c,y,p_left,p_right,c_cell, max_particles
integer, dimension(:,:), allocatable :: a, out

open (7, file = 'rsvrbndry.txt')

p_left = 0.2
p_right = 0.8
Nmax = 50
itr = 16
n = 100
m = 1
t = 1000000
max_particles = 1.25
b_delta = 0

b = int(n*(abs(p_left + p_right))*Nmax)
write(*,*), 'b=', b



write(*,*), 'int(b*max_particles) = ', int((b*max_particles))
allocate(a(int((b*max_particles)),2))
allocate(out(n,m))



Mainloop: do it = 1,itr

initialize: do i = 1,b
            call random_number(x)
            x = x * n
            call random_number(y)
            if (y*Nmax<(abs((p_left-p_right)/n)*(x)+max(p_left,p_right)*Nmax)) then
            a(i,1) = int(x)+1
            a(i,2) = 1
         endif
         enddo initialize
            
         b_delta = 2*Nmax  
         a(b+1:b+Nmax,1)=0           ! set to left reservoir
         a(b+1:b+Nmax,2)=0           ! species = 0
         
         a(b+Nmax+1:b+b_delta,1)=n+1    ! set to right reservoir
         a(b+Nmax+1:b+b_delta,2)=0      ! species = 0
 
b_delta = int(p_left*Nmax+p_right*Nmax) !i.e. b_delta =Nmax(p_l+p_r)
write(*,*) , 'b_delta = ', b_delta
write(*,*) , b+int(p_left*Nmax)
a(b+1:b+int(p_left*Nmax),1) = 0
a(b+1:b+int(p_left*Nmax),2) = 0
a(b+int(p_left*Nmax)+1:b+b_delta,1) = n+1
a(b+int(p_left*Nmax)+1:b+b_delta,2) = 0
write(*,*) ,'a_1 = ', a(:,1)
write(*,*) ,'a_2 = ', a(:,2)

out = 0;

do i = 1, n                  ! write results to screen/output file
   write(*, *), i, out(i,:)
enddo

do i = 1, b  
   if ((a(i,1) /= 0) .and. (a(i,1) /= n+1)) then
      write(*,*), 'i=',i,'a_1', a(i,1)! write output array
      write(*,*), 'i=',i,'a_2', a(i,2)! write output array
      out(a(i,1), a(i,2)) = out(a(i,1), a(i,2)) + 1
   endif
enddo

do i = 1, n                  ! write results to screen/output file
   write(*, *) i, out(i,:)
enddo

do i = 1,t
   call random_number(x)
   call random_number(c)
   call random_number(y) !picking random numbers x,c,y
   x = x*(b+b_delta)
   j = min((b+b_delta), int(x)+1)

   if  (a(j,1) == 0) then
      write(*,*), 'Moving Particle!', a(j,1)
      if (c>0.5) then
         c_cell = (out(1,1)*1.0)/Nmax
         if (y<(1-c_cell)) then
            a(j,1) = 1
            a(j,2) = 1
            out(1,1) = out(1,1) + 1
            b_delta = b_delta + 1
            a(b+b_delta,1) = 0
            a(b+b_delta,2) = 0
         endif
      endif
   elseif(a(j,1)==n+1) then
      write(*,*), 'Moving Particle!', a(j,1) 
      if (c>0.5) then
         c_cell = (out(n,1)*1.0)/Nmax
         if(y<(1-c_cell)) then
            a(j,1) = n
            a(j,2) = 1
            out(n,1) = out(n,1) + 1
            b_delta = b_delta + 1
            a(b+b_delta,1) = n+1
            a(b+b_delta,2) = 0
         endif
      endif
   else
      write(*,*), 'Moving Particle!', a(j,1)
      chooseDirection: if(c>0.5) then
         if(a(j,1) == n) then
            c_cell = p_right
         else
            c_cell = (out(a(j,1)+1,1)*1.0)/Nmax
         endif
         if(y<(1-c_cell)) then
            out(a(j,1),1) = out(a(j,1),1) - 1
            a(j,1) = a(j,1) + 1
            overBoundsRight: if(a(j,1)==(n+1)) then
               a(j,1) = a(b+b_delta,1)
               a(j,2) = a(b+b_delta,2)
               b_delta = b_delta - 1
            else
               out(a(j,1),1) = out(a(j,1),1) + 1
            endif overBoundsRight
         endif
      else
         if (a(j,1)==1) then
            c_cell = p_left
         else
            c_cell = (out(a(j,1)-1,1)*1.0)/Nmax
         endif
         if(y<(1-c_cell)) then
            out(a(j,1),1) = out(a(j,1),1) - 1
            a(j,1) = a(j,1)-1
            overBoundsLeft: if(a(j,1)==0) then
               a(j,1) = a(b+b_delta,1)
               a(j,2) = a(b+b_delta,2)
               b_delta = b_delta - 1
            else
               out(a(j,1),1) = out(a(j,1),1) + 1
            endif overBoundsLeft
         endif
      endif chooseDirection
   endif
enddo
   out = 0
   do i = 1, b+b_delta
      if((a(i,1) /= 0) .and. (a(i,1) /= n+1)) then
         out( a(i,1) , a(i,2) ) = out( a(i,1) , a(i,2) ) + 1
      endif
   enddo

   do i = 1, n                  ! write results to output file
      write(7, *) i, out(i,:)
   enddo


enddo Mainloop

close(7)

contains 

subroutine Concentration_neighbor(concentration_cell, cell)

integer, intent(in) :: cell
real, intent(inout) :: concentration_cell

concentration_cell = 0

do k = 1, b+b_delta
   if (a(k,1) == cell) then
      concentration_cell = concentration_cell + 1.0
   endif
enddo

concentration_cell = concentration_cell/Nmax

end subroutine

end program
