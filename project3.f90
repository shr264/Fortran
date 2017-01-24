program project3
   implicit none

   integer :: t, i, dt, st                              ! t time steps
   integer   ::    seed, boundary, dist, skip, stp, tofiles, Nmax
                  ! seed, boundary type, distribution type
   integer   ::    n, m, b, b_delta, j, k, itr, it, center, steps, snapshots  
                  ! n cells, m particle types, b particles
                  ! dt time, i j k u indicies, stp stops
   
   real  ::    x, c, z, y  ! random variables, gradient concentrations
   real  ::    c_cell
   character :: filename
   character  :: itr_string
   
   integer, parameter :: max_species=10
   real, parameter :: max_particles=1.25             ! the proportion used to allocate the particle array, needed for
                                                         ! reservoir bounds   
   ! Arrays carrying concentration levels
   real :: p(max_species)=0.0
   real :: p_left(max_species)=0.0
   real :: p_right(max_species)=0.0
   
   ! Arrays carrying cumulative concentration levels
   real :: q(max_species)=0.0
   real :: q_left(max_species)=0.0
   real :: q_right(max_species)=0.0
   
   namelist /Diffusion/ boundary, filename, n, m, dt, p, dist, st, skip, p_left, p_right, itr, tofiles, Nmax, steps

   integer, dimension(:,:), allocatable    ::    a       ! array of particles
   integer, dimension(:,:), allocatable    ::    out     ! output array
   integer, dimension(:), allocatable      ::    stops   ! stop times for status/checks
   !real(wp), dimension(:), allocatable     ::    p,q     ! proportion arrays
   real, dimension(:,:), allocatable   ::   av, av_means, c_squared, mean_c_squared, r_squared, mean_r_squared, concentration, mean_c, av_r_prod, mean_av_r_prod
                                                ! running average, averages over iterations
   real, dimension(:), allocatable     ::   av_c_prod, av_prod, mean_av_c_prod, mean_av_prod, correlation, correlation_squared, mean_correlation, mean_correlation_squared
	
   ! Make sure all variables have default values

   open(11,file='periodic.txt')

   boundary=0
   n = 100
   m = 1
   Nmax = 50
   dt = 20000
   dist = 1
   st = 0
   itr = 16
   tofiles = 0
   snapshots = 0

   ! set number of particles
   b = n * Nmax
   
   ! set time steps, stop times
   t = b * dt
   
   ! calculate a center point for calculating total particle correlation

   if (mod(n,2)==1) then
      center = (n/2)+1;
   else
      center = (n/2);
   endif

   if (m>1) then
      p(m) = 1 - sum(p)
      p_left(m) = 1 - sum(p_left)
      p_right(m) = 1 - sum(p_right)
   endif
   
   if ( any( (p<0) .or. (p>1) ) ) then
      write(*,*) 'ERROR: Concentrations ', p, ' are not between 0 and 1'
      stop
   endif   
   if ( any( (p_left<0) .or. (p_left>1) ) ) then
      write(*,*) 'ERROR: Concentrations ', p_left, ' are not between 0 and 1'
      stop
   endif  
   if ( any( (p_right<0) .or. (p_right>1) ) ) then
      write(*,*) 'ERROR: Concentrations ', p_right, ' are not between 0 and 1'
      stop
   endif  
      
   ! Calculate cumulative probability arrays:
   
   q(1) = p(1)
   do i = 2,m
      q(i) = p(i) + q(i-1)
   enddo
   
   q_left(1) = p_left(1)
   do i = 2,m
      q_left(i) = p_left(i) + q_left(i-1)
   enddo
   
   q_right(1) = p_right(1)
   do i = 2,m
      q_right(i) = p_right(i) + q_right(i-1)
   enddo

   allocate(a(int(max_particles*b),2))
   allocate(out(n,m))

   


   
   Mainloop: do it = 1, itr
      write(*,*) 'Working through iteration: ', it, ' of ', itr, ', ', ((100.0*it)/(itr*1.0)), 'percentage complete.'
 
      
      !call check()

      ifUniform: if (dist == 0) then
         initialCenter: do i = 1, b             ! initialize particles to center cell
            a(i,1) = center
            call ChooseSpecies(a(i,2),q)
         enddo initialCenter
      endif ifUniform

      ifGradient: if (dist == 1) then
         initialUniform: do i = 1, b
            call random_number(x)
            x = x * n
            a(i,1) = min(n, int(x) + 1)
            if(boundary == 3) then
               z = (1.0*a(i,1))/(n+1.0)
            else
               z = (1.0*a(i,1)-0.5)/(n*1.0)
            endif
            a(i,2) = 1
            !call random_number(c)
            !chooseSpeciesGradient: do j = 1,m
               !if (c < (q_left(j) + z * (q_right(j) - q_left(j)))) then
                  !a(i,2) = j
                  !exit chooseSpeciesGradient
               !endif
            !enddo chooseSpeciesGradient
         enddo initialUniform
      endif ifGradient
      
      ifReservoir: if (boundary==3) then               ! add paricles for cells 0, n+1
         b_delta = 2*Nmax
         a(b+1:b+Nmax,1)=0           ! set to left reservoir
         a(b+1:b+Nmax,2)=0           ! species = 0
         
         a(b+Nmax+1:b+b_delta,1)=n+1    ! set to right reservoir
         a(b+Nmax+1:b+b_delta,2)=0      ! species = 0
         
      endif ifReservoir
      
      call check()
   
      do i = 1, b                  
         write(*,*), i,'a_1=',a(i,1)
         write(*,*), i,'a_2=',a(i,2)
      enddo

      out = 0
      do i = 1, b                  ! write output array
         if ((a(i,1) /= 0) .and. (a(i,1) /= n+1)) then
                  out(a(i,1), a(i,2)) = out(a(i,1), a(i,2)) + 1
         endif
      enddo

         do i = 1, n                  ! write results to output file
            write(*, *) i, out(i,:)
         enddo

      do i = 1, t
         call random_number(x)
         call random_number(c)
         call random_number(y)

         x = x * (b + b_delta)
         j = min((b + b_delta), int(x) + 1)
         
   if  (a(j,1) == 0) then
      write(*,*), 'Moving Particle!', a(j,1)
      if (c>0.5) then
         c_cell = (out(1,1)*1.0)/Nmax
         if (y<(1-c_cell)) then
               a(j,1) = 1
               call ChooseSpecies(a(j,2), q_left)
               b_delta = b_delta + 1               ! particle enters system
               a(b + b_delta,1) = 0
               a(b + b_delta,2) = 0
            endif
         endif
      elseif (a(j,1) == n+1) then
      write(*,*), 'Moving Particle!', a(j,1) 
      if (c>0.5) then
         c_cell = (out(n,1)*1.0)/Nmax
         if(y<(1-c_cell)) then
               a(j,1) = n
               call ChooseSpecies(a(j,2),q_right)
               b_delta = b_delta + 1               ! particle enters system
               a(b + b_delta,1) = n+1
               a(b + b_delta,2) = 0
            endif
         endif
         else         
            chooseDirection: if (c > 0.5) then
               if(a(j,1) == n) then
                  if (boundary == 0) then
                     c_cell = (out(a(1,1),1)*1.0)/Nmax
                  elseif (boundary == 3) then
                     c_cell = sum(p_right)
                  else
                     c_cell = (out(a(j,1)+1,1)*1.0)/Nmax
                  endif
               endif
               if(y<(1-c_cell)) then
               out(a(j,1),1) = out(a(j,1),1) - 1   
               a(j,1) = a(j,1)+1
               overBoundsRight: if (a(j,1) == (n+1)) then
                  if (boundary == 0) then 
                     a(j,1) = 1
                  elseif (boundary == 1) then
                     a(j,1) = a(j,1) - 1
                  elseif (boundary == 2) then
                     a(j,1) = a(j,1) - 1
                     call ChooseSpecies(a(j,2), q_right)
                  elseif (boundary == 3) then
                     a(j,1) = a(b+b_delta,1)
                     a(j,2) = a(b+b_delta,2)
                     b_delta = b_delta - 1
                  else
                     out(a(j,1),1) = out(a(j,1),1) + 1
                  endif
               endif overBoundsRight
            endif
            else
                if (a(j,1)==1) then
                   if (boundary == 0) then
                      c_cell = (out(a(1,1),1)*1.0)/Nmax
                   elseif(boundary == 3) then
                      c_cell = sum(p_left)
                   else
                      c_cell = (out(a(j,1)-1,1)*1.0)/Nmax
                   endif
                endif
         if(y<(1-c_cell)) then
            out(a(j,1),1) = out(a(j,1),1) - 1
               a(j,1) = a(j,1)-1
               overBoundsLeft: if (a(j,1) == 0) then
                  if (boundary == 0) then
                     a(j,1) = n
                  elseif (boundary == 1) then
                     a(j,1) = a(j,1) + 1
                  elseif (boundary == 2) then
                     a(j,1) = a(j,1) + 1
                     call ChooseSpecies(a(j,2), q_left)
                  elseif (boundary == 3) then
                     a(j,1) = a(b+b_delta,1)
                     a(j,2) = a(b+b_delta,2)
                     b_delta = b_delta - 1
                  else
                     out(a(j,1),1) = out(a(j,1),1) + 1
                  endif
               endif overBoundsLeft
            endif
         endif chooseDirection
      endif
   enddo

            out = 0
            
            do j = 1, b+b_delta                  ! write output array
               if ((a(j,1) /= 0) .and. (a(j,1) /= n+1)) then
                  out(a(j,1), a(j,2)) = out(a(j,1), a(j,2)) + 1
               endif
            enddo

   do i = 1, n                  ! write results to output file
      write(11, *) i, out(i,:)
   enddo
      
   enddo Mainloop

close(11)

contains

   subroutine Check()
      ! sanity check code.
      write(*,*) 'Working... ', 100*(i/(t*1.0)),'% complete. Iteration ', it, ' of', itr
      write(*,*) 'Total particles: ', (b + b_delta), 'in possible array size: ', int(max_particles*b)
      write(*,*) 'Step: i = ', i
      ! write(*, '(a)',advance='no') '-'
      ! write(*, *) 'Total particles: ', b+b_delta
   end subroutine
   
   subroutine ChooseSpecies(species, cumulative)
      ! ChooseSpecies takes an integer species number and a cumulative probability array
      ! and changes the 'species' variable to an index of the probability array.
      integer  :: species
      real, dimension(:), intent(in) :: cumulative
      
      call random_number(x)
      
      Choose: do k = 1,m         ! choose type based on given probabilities
         if (x < cumulative(k)) then
            species = k
            exit Choose
         endif
      enddo Choose
      
   end subroutine

end program
