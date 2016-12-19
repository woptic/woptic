!!! woptic/src/kanalysis.f90
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann
!!!

program kanalysis
  use util,      only:
  use const,     only: BUFSZ, DPk
  use structmod, only: struct_t, struct_read

  implicit none

  character(BUFSZ) :: seedname, argdummy
  integer nk,jk,iarg,unit1,unit2
  integer niter,nw,dumi
  integer jw,kmaxidx
  real(DPk), allocatable :: kp(:,:),kwdata(:,:,:)
  real(DPk) :: kmax,convfac,dumr,dw
  character :: dum1
  type(struct_t) :: stru
  logical matlabmode

  iarg=iargc() 
  if (iarg.lt.2) then
     write(*,*)"Usage: kanalysis startingfrequency case mode(optional)"
     stop
  endif
  matlabmode = .false.
   if (iarg.eq.3) then
       call getarg(3,argdummy)
       if (argdummy(1:1).eq."1") then
          matlabmode = .true.
          write(*,*)"Matlab mode"
       endif
   endif
  
  call getarg(2,argdummy)
  write(seedname,"(A80)")argdummy
  call getarg(1,argdummy)
  read(argdummy,*)niter

  unit1 = 1
  open(unit=unit1,file=trim(seedname)//'.inwop',status='old')
  read(unit1,*)
  read(unit1,*)dumr,dw
  close(unit1)
  unit2 =2
  open(unit=unit1,file=trim(seedname)//'.struct',status='old')
  call struct_read(unit1, stru)
  close(unit1)

   open(unit=unit1,file=trim(seedname)//'.kcontribw_band',status='old')
   open(unit=unit2,file=trim(seedname)//'.optanalysis_band',status='replace')
   read(unit1,*)dum1,nk,nw,dumi,convfac
   write(unit2,*)dum1,nk,nw+1-niter,convfac
   kmax = 0d0
   allocate(kwdata(nk,0:nw,6),kp(nk,3))
    do jk=1,nk
         read(unit1,*)dumi,kwdata(jk,0,:)
         do jw=1,nw
             read(unit1,*)kwdata(jk,jw,:)
             if (jw.ge.niter) then
                 write(unit2,*)jk,dw*jw,kwdata(jk,jw,1)*convfac/1000d0
             endif
         enddo
         if (kwdata(jk,150,1).gt.kmax) then
               kmax = kwdata(jk,150,1)
               kmaxidx = jk
        endif
        if  (.not.matlabmode) then
           write(unit2,*)
        endif
    enddo
    write(*,*)'maximumk',kmax,kmaxidx
!   stop
   close(unit1)
   open(unit=unit1,file=trim(seedname)//'.K1w_band',status='old')
   open(unit=unit2,file=trim(seedname)//'.thermanalysis_band',status='replace')
!    read(unit1,*)dum1,nk,nw,dumi,convfac
   write(unit2,*)dum1,nk,nw+1,convfac
   kmax = 0d0
    do jk=1,nk
!          read(unit1,*)dumi,kwdata(jk,0,:)
!          do jw=1,nw
             read(unit1,*)dumi,kwdata(jk,0,:)
              write(unit2,*)jk,kwdata(jk,0,1)
!          enddo
         if (kwdata(jk,0,1).gt.kmax) then
               kmax = kwdata(jk,0,1)
               kmaxidx = jk
        endif
!         if  (.not.matlabmode) then
!            write(unit2,*)
!         endif
    enddo
    write(*,*)'maximumk',kmax,kmaxidx
!   stop
   close(unit1)
 end program kanalysis
