!!! woptic/src/obtain_dist.f90
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!
!!! $Id: obtain_dist.f90 408 2015-06-05 12:07:46Z assmann $

program obtain_dist

implicit none

character*80 seed,arg,dummy
integer :: unit1,unit2,nw,j1,j2,stat
real*8,allocatable :: centres(:,:),distmatx(:,:),distmaty(:,:),distmatz(:,:)

  unit1 = 1
  unit2 = 2

  call get_command_argument(1, VALUE=arg, STATUS=stat)
  if (stat /= 0) stop 'error getting seed'
  read(arg, *) seed
  write(*,*)seed
  
  open(UNIT=unit1, FILE=trim(seed)//"_centres.xyz", STATUS='old', FORM='formatted')
  read(unit1,*)nw
  read(unit1,*)
  allocate(centres(nw,3),distmatx(nw,nw),distmaty(nw,nw),distmatz(nw,nw))
  do j1=1,nw
     read(unit1,*)dummy,centres(j1,:)
     write(*,*)centres(j1,:)
  enddo
  close(unit1)

  do j1=1,nw
     do j2=1,nw
        distmatx(j1,j2) = 1.889725d0*(centres(j1,1)-centres(j2,1))
        distmaty(j1,j2) = 1.889725d0*(centres(j1,2)-centres(j2,2))
        distmatz(j1,j2) = 1.889725d0*(centres(j1,3)-centres(j2,3))
     enddo
  enddo

  open(UNIT=unit2, FILE=trim(seed)//".intrahop", STATUS='replace', FORM='formatted')
  do j1=1,nw
     write(unit2,"(100F12.4)")distmatx(j1,:)
  enddo
  do j1=1,nw
     write(unit2,"(100F12.4)")distmaty(j1,:)
  enddo
  do j1=1,nw
     write(unit2,"(100F12.4)")distmatz(j1,:)
  enddo
  close(unit2)
end program obtain_dist


!! Time-stamp: <2015-06-05 13:48:35 assman@faepop23.tu-graz.ac.at>
