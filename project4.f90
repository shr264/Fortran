program ColoredExclusion

implicit none

integer, dimension(:,:), allocatable :: particle_position_array, output
real, dimension(:), allocatable :: p_left, p_right, initialspecies
integer, dimension(:), allocatable :: distance
integer :: nspaces, avgparticlesperspace, maxparticlesperspace, nparticles, nspecies
integer :: total,  particletomove, species, nparticles_delta, positiontomove
integer :: repeats, directiontomove, end1, end2, timesteps, boundary
real :: extraspace, initialspecies_i

real ::rndm1, rndm2, rndm3, rndm4, rndm5, rndm6, rndm7
integer :: count1, count2, count3, count4, count5, count6, extcount

open (7, file = 'periodic.txt')
open (8, file = 'periodicdistances.txt')


boundary = 1
timesteps = 100000
repeats = 16
nspaces = 20
avgparticlesperspace = 20
maxparticlesperspace = 40
nparticles_delta = 0
nparticles = nspaces*avgparticlesperspace + nparticles_delta
extraspace = 1.5
nspecies = 1
total = int(extraspace*nparticles)


allocate(particle_position_array(total,2))
allocate(output(nspaces,nspecies+1))
allocate(initialspecies(nspecies))
allocate(p_left(nspecies))
allocate(p_right(nspecies))
allocate(distance(total))

particle_position_array = 0
p_left = 1
p_right = 1
initialspecies = 1

Mainloop: do count4 = 1,repeats

output = 0;
distance = 0;

Initialize: do count1 = 1,nparticles
   call random_number(rndm1)
   particletomove = int(rndm1*nparticles) + 1

   species = pickspecies(initialspecies,nspecies)
   
   extcount = 0
   call findposition

   if(boundary==3) then
      nparticles_delta = int(sum(p_left)+sum(p_right))*maxparticlesperspace

      end1 = int(sum(p_left)*maxparticlesperspace)
      particle_position_array(nparticles+1:nparticles+end1,1) = 0
      particle_position_array(nparticles+1:nparticles+end1,2) = 1
      end2 = int(sum(p_right)*maxparticlesperspace)
      particle_position_array(nparticles+1:nparticles+end2,1) = nspaces+1
      particle_position_array(nparticles+1:nparticles+end2,2) = 1
   endif

enddo Initialize

write(*,*), 'Initial Positions'
do count3 = 1,nspaces
   write(*,*), count3, output(count3,:)
enddo

Timestep: do count6 = 1,timesteps
   
   write(*,*), 'Timestep:', count6 
   call random_number(rndm4)
   particletomove = int(rndm4*nparticles) + 1
   write(*,*), 'Moving Particle: ', particletomove 
   write(*,*), 'In position: ', particle_position_array(particletomove,1) 
   write(*,*), 'Of type:', particle_position_array(particletomove,2)
   call random_number(rndm5)
   if(rndm5<0.5) then
      directiontomove = 0 !this is left
   else
      directiontomove = 1 !this is right
   endif
!left boundary going right
   if(particle_position_array(particletomove,1)==0) then
      if(directiontomove==1) then
         positiontomove = particle_position_array(particletomove,1) + 1
         write(*,*), 'Position Moving to: ', positiontomove
         if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
            particle_position_array(particletomove,1) = particle_position_array(particletomove,1) + 1
            species = pickspecies(p_left, nspecies)
            particle_position_array(particletomove,2) = species
            output(positiontomove,species) = output(positiontomove,species) + 1
            write(*,*), 'moving to',positiontomove, species, output(positiontomove,species)
            nparticles_delta = nparticles_delta+1
            distance(particletomove) = distance(particletomove) + 1
            write(*,*), 'Particle has entered system'
            write(*,*), '# of particles: ', nparticles+nparticles_delta
         endif
      endif
!right boundary going left
   elseif(particle_position_array(particletomove,1)==nspaces+1) then
      if(directiontomove==0) then
         positiontomove = particle_position_array(particletomove,1) - 1
         write(*,*), 'Position Moving to: ', positiontomove
         if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
            particle_position_array(particletomove,1) = particle_position_array(particletomove,1) - 1
            species = pickspecies(p_right, nspecies)
            particle_position_array(particletomove,2) = species
            output(positiontomove,species) = output(positiontomove,species) + 1
            write(*,*), 'moving to',positiontomove, species, output(positiontomove,species)
            nparticles_delta = nparticles_delta + 1
            distance(particletomove) = distance(particletomove) + 1
            write(*,*), 'Particle has entered the system'
            write(*,*), '# of particles: ', nparticles+nparticles_delta
         endif
      endif
!at postion nspaces
   elseif(particle_position_array(particletomove,1)==nspaces) then
      if(directiontomove==1) then
         if(boundary==3) then
            positiontomove = particle_position_array(particletomove,1) + 1
            write(*,*), 'Position Moving to: ', positiontomove
            particle_position_array(particletomove,1) = particle_position_array(particletomove,1) + 1
            species = particle_position_array(particletomove,2)
            particle_position_array(particletomove,2) = 1
            output(nspaces,species) = output(nspaces,species) - 1
            write(*,*), 'moving from',nspaces, species, output(nspaces,species)
            nparticles_delta = nparticles_delta - 1
            distance(particletomove) = distance(particletomove) + 1
            write(*,*), 'Particle has left the system'
            write(*,*), '# of particles: ', nparticles+nparticles_delta
         elseif(boundary==1) then
            positiontomove = 1
            if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
               write(*,*), 'Position Moving to: ', positiontomove
               particle_position_array(particletomove,1) = 1
               species = particle_position_array(particletomove,2)
               output(positiontomove,species) = output(positiontomove,species) + 1
               write(*,*),'moving to', positiontomove, species, output(positiontomove,species)
               output(nspaces,species) = output(nspaces,species) - 1
               write(*,*), 'moving from',nspaces, species, output(nspaces,species)
               distance(particletomove) = distance(particletomove) + 1
            endif
         elseif(boundary==2) then
            species = pickspecies(p_right, nspecies)
            particle_position_array(particletomove,2) = species
         endif
         elseif(directiontomove==0) then
            positiontomove = particle_position_array(particletomove,1) - 1
            write(*,*), 'Position Moving to: ', positiontomove
            if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
               particle_position_array(particletomove,1) = particle_position_array(particletomove,1) - 1
               species = particle_position_array(particletomove,2)
               output(positiontomove,species) = output(positiontomove,species) + 1
               write(*,*), 'moving to',positiontomove, species, output(positiontomove,species)
               output(positiontomove+1,species) = output(positiontomove+1,species) - 1
               write(*,*), 'moving from',positiontomove+1, species, output(positiontomove+1,species)
               distance(particletomove) = distance(particletomove) + 1
            endif
         endif
!at position 1
   elseif(particle_position_array(particletomove,1)==1) then
      if(directiontomove==0) then
         if (boundary==3) then
            positiontomove = particle_position_array(particletomove,1) - 1
            write(*,*), 'Position Moving to: ', positiontomove
            particle_position_array(particletomove,1) = particle_position_array(particletomove,1) - 1
            species = particle_position_array(particletomove,2)
            particle_position_array(particletomove,2) = 1
            output(1,species) = output(1,species) - 1
            write(*,*), 'moving from',1, species, output(1,species)
            distance(particletomove) = distance(particletomove) + 1
            nparticles_delta = nparticles_delta - 1
            write(*,*), 'Particle has left the system'
            write(*,*), '# of particles: ', nparticles+nparticles_delta
      elseif (boundary==1) then
            positiontomove = nspaces
            if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
               write(*,*), 'Position Moving to: ', positiontomove
               particle_position_array(particletomove,1) = nspaces
               species = particle_position_array(particletomove,2)
               output(positiontomove,species) = output(positiontomove,species) + 1
               write(*,*),'moving to', positiontomove, species, output(positiontomove,species)
               output(1,species) = output(1,species) - 1
               write(*,*), 'moving from',1, species, output(1,species)
               distance(particletomove) = distance(particletomove) + 1
            endif  
         elseif(boundary==2) then
            species = pickspecies(p_left, nspecies)
            particle_position_array(particletomove,2) = species
      endif
      elseif(directiontomove==1) then
         positiontomove = particle_position_array(particletomove,1) + 1
         write(*,*), 'Position Moving to: ', positiontomove
         if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
            particle_position_array(particletomove,1) = particle_position_array(particletomove,1) + 1
            species = particle_position_array(particletomove,2)
            output(positiontomove,species) = output(positiontomove,species) + 1
            write(*,*), 'moving to',positiontomove, species, output(positiontomove,species)
            output(positiontomove-1,species) = output(positiontomove-1,species) - 1
            write(*,*), 'moving from',positiontomove-1, species, output(positiontomove-1,species)
            distance(particletomove) = distance(particletomove) + 1
         endif
      endif
!all other posistions
   else
      if(particle_position_array(particletomove,1)/=1) then
         if(particle_position_array(particletomove,1)/=nspaces) then
            if(particle_position_array(particletomove,1)/=0) then
               if(particle_position_array(particletomove,1)/=nspaces+1) then
                  if(directiontomove==0) then
                     positiontomove = particle_position_array(particletomove,1) - 1
                     write(*,*), 'Position Moving to: ', positiontomove
                     if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
                        particle_position_array(particletomove,1) = particle_position_array(particletomove,1) - 1
                        species = particle_position_array(particletomove,2)
                        output(positiontomove,species) = output(positiontomove,species) + 1
                        write(*,*),'moving to', positiontomove, species, output(positiontomove,species)
                        output(positiontomove+1,species) = output(positiontomove+1,species) - 1
                        write(*,*), 'moving from',positiontomove+1, species, output(positiontomove+1,species)
                        distance(particletomove) = distance(particletomove) + 1
                     endif
                  elseif(directiontomove==1) then
                     positiontomove = particle_position_array(particletomove,1) + 1
                     write(*,*), 'Position Moving to: ', positiontomove
                     if(sum(output(positiontomove,:)) /= maxparticlesperspace) then
                        particle_position_array(particletomove,1) = particle_position_array(particletomove,1) + 1
                        species = particle_position_array(particletomove,2)
                        output(positiontomove,species) = output(positiontomove,species) + 1
                        write(*,*),'moving to', positiontomove, species, output(positiontomove,species)
                        output(positiontomove-1,species) = output(positiontomove-1,species) - 1
                        write(*,*), 'moving from',positiontomove-1, species, output(positiontomove-1,species)
                        distance(particletomove) = distance(particletomove) + 1
                     endif
                  endif
               endif
            endif
         endif
      endif
   endif
enddo Timestep

write(*,*), 'Positions'
do count3 = 1,nspaces
   write(*,*), count3, output(count3,:)
   write(7,*), count3, output(count3,:)
enddo

write(*,*), 'Distances'
do count3 = 1,nparticles+nparticles_delta
   write(*,*), count3, distance(count3)
   write(8,*), count3, distance(count3)
enddo

enddo Mainloop



contains

  recursive subroutine findposition
   call random_number(rndm3)
   positiontomove = int(rndm3*(nspaces))+1
   write (*,*), 'Particles in position', positiontomove ,'Type = ', species, 'Total = ' ,sum(output(positiontomove,:))
   call random_number(rndm7)
   if((sum(output(positiontomove,:))/maxparticlesperspace)<rndm7) then
      particle_position_array(particletomove,1) = positiontomove
      particle_position_array(particletomove,2) = species
      output(positiontomove,species) = output(positiontomove,species) + 1
   else
      write(*,*), 'Particles didnot find place in system. Looking for a replacement' 
      call findposition
   endif
 end subroutine findposition

 integer function pickspecies(init_species,nspecies)

   integer :: nspecies
   real, dimension(nspecies) :: init_species
   real :: rndm2

   call random_number(rndm2)
   rndm2 = rndm2*sum(init_species)
   write(*,*), 'RN = ', rndm2
   initialspecies_i = 0
   pickspecies = 0

   if(nspecies==1) then
      pickspecies = 2
   elseif(nspecies>1) then
      count2 = 1
      do while((initialspecies_i<rndm2).and.(count2<nspecies))
         initialspecies_i = initialspecies_i + init_species(count2)
         write(*,*), 'is=' ,initialspecies_i
         count2 = count2 + 1
         write(*,*), 'count2=' ,count2
      enddo
      if(rndm2<initialspecies_i) then
         pickspecies = count2
      elseif((rndm2>initialspecies_i).and.(rndm2<(initialspecies_i+init_species(count2)))) then
         pickspecies = count2 + 1
      endif
   endif
   write(*,*), 'Species = ',pickspecies
   return
 end function pickspecies

end program
