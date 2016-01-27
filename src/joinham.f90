!!! woptic/src/joinham.f90
!!!
!!!    Joins two Hamiltonian files and two mommat files for the woptic
!!!    algorithm.  The header of the first file is kept and the second
!!!    one is appended to the first.
!!!
!!!    It's case.hk_old + case.hk     = case.hk_joined      and
!!!     case.mommat_old + case.mommat = case.mommat_joined
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann
!!!

program combine_hamfiles
  use const,     only: BUFSZ, DPk
  use clio,      only: fetcharg, croak
  use util,      only: lowercase
  use woptic_io, only: suf_ham, suf_mommat
  use woptic,    only: fmt_vk_head
  use maybebin,  only: maybin_open, maybin_read, maybin_write, formatted

  implicit none

  character(*), parameter :: rev_str = "$version: v0.1.0-27-ga77d9d3$"
  character(*), parameter :: woptic_version = rev_str(11 : len (rev_str)-1)

  integer,          parameter :: unit_ham1=11, unit_ham2=12, unit_mom1=13
  integer,          parameter :: unit_mom2=14, unit_ham3=15, unit_mom3=16
  character(len=*), parameter :: old='_old', cur='', jnd='_joined'

  integer          :: Nkb1(2), Nkb2(2), i,j
  logical          :: info, mation, mommat=.false.
  character(BUFSZ) :: case, ham1, ham2, ham3, mom1, mom2, mom3, buf

  ! The “HK” file may be binary (vvk!).  During the first woptic
  ! iteration, HK1 will be plain text even if HK2 is binary.
  logical :: binham1, binham2

  ! These will be used to copy binary ‘vvk’ files.  Caveat: size(kw)=4
  ! assumes that only vvk, not hk or vk, files will ever be binary (3
  ! k-coordinates and ω).
  real(DPk)                 :: kw(4)
  complex(DPk), allocatable :: vv(:)

  select case (command_argument_count())
  case (1)
     call fetcharg(1, case)
     option: select case (case)
     case ('-h','-H','-help','--help')
        print '(A)', 'USAGE: joinham CASE'
        print '(A)', '   or  joinham HK MOMMAT'
        print '(A)', '   or  joinham HK1 HK2 HKOUT'
        print '(A)', '   or  joinham HK1 HK2 HKOUT MOM1 MOM2 MOMOUT'
        print '(A)', '   or  joinham { OPTION }'
        print '(A)'
        print '(A)', 'Splice two ‘hk’ files together, and possibly two ‘mommat2’ files.'
        print '(A)'
        print '(A)', 'OPTIONS:'
        print '(A)', '  --help, -h'
        print '(A)', '  --version, -v'
        call exit(0)

     case ('-v', '-V', '-version', '--version')
        print '("woptic_main ", A)', WOPTIC_VERSION
        call exit(0)

     case default
        if (case(1:1) == '-') &
             call croak('unknown option ‘'//trim(case)//'’ (try ‘-h’)')

        ham1 = trim(case)//suf_ham//old; mom1 = trim(case)//suf_mommat//old
        ham2 = trim(case)//suf_ham//cur; mom2 = trim(case)//suf_mommat//cur
        ham3 = trim(case)//suf_ham//jnd; mom3 = trim(case)//suf_mommat//jnd

        inquire(FILE=mom1, EXIST=info); inquire(FILE=mom2, exist=mation)
        mommat = info .and. mation
     end select option
  case (2)
     call fetcharg(1,ham2); ham1 = trim(ham2)//old; ham3 = trim(ham2)//jnd
     call fetcharg(1,mom2); mom1 = trim(mom2)//old; mom3 = trim(mom2)//jnd
     mommat = .true.
  case (3)
     call fetcharg(1,ham1); call fetcharg(2,ham2); call fetcharg(3,ham3)
     mommat = .false.
  case (6)
     call fetcharg(1,ham1); call fetcharg(2,ham2); call fetcharg(3,ham3)
     call fetcharg(4,mom1); call fetcharg(5,mom2); call fetcharg(6,mom3)
     mommat = .true.
  case default
     call croak ('must have 1, 2, 3, or 6 arguments (see ‘-h’)')
  end select

  call maybin_open(unit_ham1, FILE=ham1, get=binham1, STATUS='old')
  call maybin_open(unit_ham2, FILE=ham2, get=binham2, STATUS='old')
  call maybin_open(unit_ham3, FILE=ham3, set=binham2, STATUS='replace')

  if (mommat) then
     open(unit_mom1, FILE=mom1, STATUS='old')
     open(unit_mom2, FILE=mom2, STATUS='old')
     open(unit_mom3, FILE=mom3, STATUS='replace')
  end if

!!! Splice ‘hk’
  call maybin_read(unit_ham1, Nkb1)
  call maybin_read(unit_ham2, Nkb2)

  if (Nkb1(2)/=Nkb2(2) .and. Nkb1(1)/=0) &
       call croak("Error: number of bands not consistent")

  if (binham2) allocate(vv(Nkb2(2)))

  call maybin_write(unit_ham3, (/ Nkb1(1)+Nkb2(1), Nkb2(2) /), fmt=fmt_vk_head)
  do i=1,Nkb1(1)
     call cpk_ham(unit_ham1, unit_ham3, Nkb2(2))
  enddo

  do i=1,Nkb2(1)
     call cpk_ham(unit_ham2, unit_ham3, Nkb2(2))
  enddo

  if (mommat) then
!!! Splice ‘mommat2’
     call copy(unit_mom1, unit_mom3, buf)
     do j=1,Nkb1(1)
        call cpk_mom(unit_mom1, unit_mom3)
     enddo

     read(unit_mom2,*)
     do j=1,Nkb2(1)
        call cpk_mom(unit_mom2, unit_mom3)
     enddo
  endif

contains
  subroutine copy(src, dst, buf)
    integer,          intent(in)  :: src, dst
    character(len=*), intent(out) :: buf

    read (src, '(A)') buf
    write(dst, '(A)') trim(buf)
  end subroutine copy

  subroutine cpk_ham(src, dst, Nb)
    integer, intent(in) :: src, dst, Nb
    integer :: j

    ! 21*2*#bands is large enough for a line from ‘hk’.
    character(len=42*Nkb2(2)) :: buf

    if (formatted(src)) then
       call copy(src, dst, buf)
       do j=1, Nb
          call copy(src, dst, buf)
       enddo
    else
       read(src) kw; write(dst) kw
       do j=1, Nb
          read(src) vv; write(dst) vv
       end do
    end if
  end subroutine cpk_ham

  subroutine cpk_mom(src, dst)
    integer, intent(in) :: src, dst

    integer :: nemin, nemax, i, j

    call copy(src, dst, buf)
    call copy(src, dst, buf); read(buf(28:38),"(2I5)") nemin, nemax
    call copy(src, dst, buf)

    do i = nemin, nemax
       do j = i, nemax
          call copy(src, dst, buf)
       enddo
    enddo
  end subroutine cpk_mom
end program combine_hamfiles


!! Time-stamp: <2015-11-10 16:16:44 assman@faepop36.tu-graz.ac.at>
