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
!!!           2014-2016 Elias Assmann
!!!

program combine_hamfiles
  use const,     only: BUFSZ, DPk
  use clio,      only: fetcharg, croak
  use util,      only: lowercase
  use woptic_io, only: fn_ham,    unit_hamcur => unit_ham,                 &
       &               fn_mommat, unit_momcur => unit_mommat, set_casename,&
       &               jnd => suf_jnd, old => suf_old
  use woptic,    only: fmt_vk_head
  use maybebin,  only: maybin_open, maybin_read, maybin_write, formatted

  implicit none

  character(*), parameter :: rev_str = "$version: v0.1.0-45-gf5b26c6$"
  character(*), parameter :: woptic_version = rev_str(11 : len (rev_str)-1)

  integer,          parameter :: unit_hamold=121, unit_momold=13
  integer,          parameter :: unit_hamnew=123, unit_momnew=16
  character(len=*), parameter :: cur=''

  integer          :: Nkb1(2), Nkb2(2), i,j
  logical          :: info, mation, mommat=.false.
  character(BUFSZ) :: hamold, hamcur, hamnew, momold, momcur, momnew, buf
  character(BUFSZ) :: file
  character(2)     :: updn=''

  ! The “HK” file may be binary (vvk!).  During the first woptic
  ! iteration, HK1 will be plain text even if HK2 is binary.
  logical :: binhamold, binhamcur

  ! These will be used to copy binary ‘vvk’ files.  Caveat: size(kw)=4
  ! assumes that only vvk, not hk or vk, files will ever be binary (3
  ! k-coordinates and ω).
  real(DPk)                 :: kw(4)
  complex(DPk), allocatable :: vv(:)

  args: select case (command_argument_count())
  case (1)                    ! joinham -h | joinham -v | joinham CASE
     call fetcharg(1, file)

     option: select case (file)
     case ('-h','-H','-help','--help')
        print '(A)', 'USAGE: joinham [--up|--dn] CASE'
        print '(A)', '   or  joinham HK MOMMAT'
        print '(A)', '   or  joinham HK1 HK2 HKOUT'
        print '(A)', '   or  joinham HK1 HK2 HKOUT MOM1 MOM2 MOMOUT'
        print '(A)', '   or  joinham { --help | --version }'
        print '(A)'
        print '(A)', 'Splice two ‘hk’ files together, and possibly two &
             &‘mommat2’ files.'
        print '(A)'
        print '(A)', 'OPTIONS:'
        print '(A)', '  --help, -h     output this message and exit'
        print '(A)', '  --version, -v  output version information and exit'
        print '(A)', '  --up|--dn      spin-polarized mode'
        call exit(0)

     case ('-v', '-V', '-version', '--version')
        print '("woptic_main ", A)', WOPTIC_VERSION
        call exit(0)

     case default
        if (file(1:1) == '-') &
             call croak('bad option ‘'//trim(file)//'’ (try ‘-h’)')

        call from_casename(file, updn)
     end select option

  case (2)                ! joinham HK MOMMAT | joinham [-up|-dn] CASE
     call fetcharg(1, hamcur)
     call fetcharg(2, momcur)

     spmode: select case(hamcur)
     case ('-up', '--up', '-dn', '--dn')
        updn = hamcur(len_trim(hamcur)-1 : )

        call from_casename(momcur, updn)

     case default
        hamold = trim(hamcur)//old; hamnew = trim(hamcur)//jnd
        momold = trim(momcur)//old; momnew = trim(momcur)//jnd
        mommat = .true.
     end select spmode

  case (3)                      ! joinham HK1 HK2 HK3
     call fetcharg(1,hamold);call fetcharg(2,hamcur);call fetcharg(3,hamnew)
     mommat = .false.

  case (6)                      ! joinham MOM1 MOM2 MOM3
     call fetcharg(1,hamold);call fetcharg(2,hamcur);call fetcharg(3,hamnew)
     call fetcharg(4,momold);call fetcharg(5,momcur);call fetcharg(6,momnew)
     mommat = .true.

  case default
     call croak ('must have 1, 2, 3, or 6 arguments (see ‘-h’)')
  end select args

  call maybin_open(unit_hamold, FILE=hamold, get=binhamold, STATUS='old')
  call maybin_open(unit_hamcur, FILE=hamcur, get=binhamcur, STATUS='old')
  call maybin_open(unit_hamnew, FILE=hamnew, set=binhamcur, STATUS='replace')

  if (mommat) then
     open(unit_momold, FILE=momold, STATUS='old')
     open(unit_momcur, FILE=momcur, STATUS='old')
     open(unit_momnew, FILE=momnew, STATUS='replace')
  end if

!!! Splice ‘hk’
  call maybin_read(unit_hamold, Nkb1)
  call maybin_read(unit_hamcur, Nkb2)

  if (Nkb1(2)/=Nkb2(2) .and. Nkb1(1)/=0) &
       call croak("Error: number of bands not consistent")

  if (binhamcur) allocate(vv(Nkb2(2)))

  call maybin_write(unit_hamnew, (/ Nkb1(1)+Nkb2(1), Nkb2(2) /), &
       &            fmt=fmt_vk_head)
  do i=1,Nkb1(1)
     call cpk_ham(unit_hamold, unit_hamnew, Nkb2(2))
  enddo

  do i=1,Nkb2(1)
     call cpk_ham(unit_hamcur, unit_hamnew, Nkb2(2))
  enddo

  if (mommat) then
!!! Splice ‘mommat2’
     call copy(unit_momold, unit_momnew, buf)
     do j=1,Nkb1(1)
        call cpk_mom(unit_momold, unit_momnew)
     enddo

     read(unit_momcur,*)
     do j=1,Nkb2(1)
        call cpk_mom(unit_momcur, unit_momnew)
     enddo
  endif

contains
  subroutine from_casename(file, updn)
    character(*) :: file
    character(2) :: updn

    call set_casename(file, UPDN=updn)

    hamold = trim(fn_ham)//old; momold = trim(fn_mommat)//old
    hamcur = trim(fn_ham)//cur; momcur = trim(fn_mommat)//cur
    hamnew = trim(fn_ham)//jnd; momnew = trim(fn_mommat)//jnd

    inquire(FILE=momold, EXIST=info); inquire(FILE=momcur, EXIST=mation)
    mommat = info .and. mation
  end subroutine from_casename

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


!! Time-stamp: <2016-02-15 20:16:55 assman@faepop36.tu-graz.ac.at>
