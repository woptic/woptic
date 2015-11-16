!!! woptic/src/inwopcheck.f90
!!!
!!!    Parses ‘inwop‘ file, checks for consistency, and outputs
!!!    variables needed by ‘woptic’
!!!
!!! Copyright 2014-2015 Elias Assmann
!!!
!!! $Id: inwopcheck.f90 386 2015-06-01 13:11:11Z assmann $

program check_inwop
  use woptic, only: inwop_t, inwop_read
  use clio,   only: argstr, fetcharg

  type(inwop_t) :: inw
  type(argstr)  :: fname
  integer       :: i
  character(*), parameter :: names(6) = (/ &
       'runlapw1  ', &
       'runoptic  ', &
       'convert_vr', &
       'mixed_vr  ', &
       'peierls   ', &
       'need_ham  '  /)
  character(3) :: values(size(names,1)) = ''

  call fetcharg(1, fname)

  call inwop_read(fname, inw)

  if (inw%read_energy .or. inw%read_mommat) values(1)='yes'
  if (                     inw%read_mommat) values(2)='yes'
  if (                     inw%read_vk    ) values(3)='yes'
  if (                     inw%read_vvk   ) values(4)='yes'
  if (                     inw%peierls    ) values(5)='yes'
  if (                     inw%read_ham   ) values(6)='yes'
  if (inw%orig_umatrix) values((/1, 2/)) = ''

  print '("matelmode=", I0)', inw%matelmode
  do i=1,size(values, 1)
     print '(A, "=", A)', trim(names(i)), values(i)
  end do

  print '("# Emax ΔE δ β = ", F7.3, 2F6.3, F8.3)', &
       inw%Emax, inw%dE, inw%delterl, inw%beta
  print '("# chem. pot. = ", F6.3)', inw%chempot
  print '("# use self energy?; mixed transitions?", 2L2)', &
       inw%selfE, inw%mixed
  print '("# rotate WF?; use U(k) from Wannier90?", 2L2)', &
       inw%wfrot, inw%orig_umatrix
  print '("# compute opt. cond.?, DOS?, joint DOS?", 3L2)', &
       inw%optcond, inw%dos, inw%joint
end program check_inwop


!! Time-stamp: <2015-06-01 12:36:57 assman@faepop23.tu-graz.ac.at>
