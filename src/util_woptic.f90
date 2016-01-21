!!! woptic/src/util_woptic.f90
!!!
!!!    Collection of routines for woptic
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2015 Elias Assmann
!!!

!!/=== Things that are exported here =============
!!
!! woptic: Things specific to woptic
!!
!!    get_mommat(), get_mommat_1k(), type(inwop_t), inwop_read(),
!!    MODE_*, fmt_*, woptic_version
!!
!!
!! maybebin: Read/write files that may be binary or plain text
!!
!!    maybin_open(), maybin_read(), maybin_write(),
!!    maybin_get(),  maybin_set()
!!
!!
!! woptic_io: Files and units for woptic
!!
!!    fn_*, unit_*, suf_*; print1or3(), set_casename()
!!
!!\===============================================


!---------------- Things specific to woptic            ----------------------
module woptic
  use const, only: DPk

  implicit none
  private
  public :: get_mommat, get_mommat_1k, inwop_t, inwop_read

  interface get_mommat
     module procedure get_mommat_unit, get_mommat_fname
  end interface get_mommat

  interface inwop_read
     module procedure inwop_read_fname, inwop_read_unit, inwop_read_argstr
  end interface inwop_read

!!! Matelmode: what to use for momentum matrix elements
!!!
!!! The tens place gives the “major mode”, the ones place is reserved
!!! for “sub-modes”.  Input given as single digit I will be equivalent
!!! to 10*I.
  integer, public, parameter :: &
       MODE_Peierls = 10, & ! Peierls approximation (band derivative)
       MODE_interp  = 20, & ! Wannier interpolated matrix elements
       MODE_optic   = 30, & ! matrix elements from ‘optic’
       MODE_F_OrigU = 31, & ! like MODE_optic, but take U(k) from ‘chk’
       MODE_Bloch   = 40, & ! LDA-only mode in Wannier basis
       MODE_LDA     = 50    ! LDA-only mode in KS basis, no interpolation

  integer, public, parameter :: MATELMODES(6) = & 
       (/ MODE_Peierls, MODE_interp, MODE_optic, MODE_F_OrigU, &
       &  MODE_Bloch, MODE_LDA /)

  integer, parameter :: LEN_MATEL=7
  character(LEN_MATEL), public, parameter :: MATELNAMES(6) = & 
       (/ "Peierls", "interp ", "optic  ", "OrigU  ", "Bloch  ", "LDA    " /)

  integer, parameter :: UNSUPPORTED_MATELMODES(3) = &
       (/ MODE_Peierls, MODE_Bloch, MODE_LDA /)

!!! Type representing ‘inwop’ files
  type inwop_t
!!!! Variables directly corresponding to entries in config file
     ! what to calculate (‘OPT’, ‘DOS’, or ‘JOINT’)
     character(4)          :: mode
     integer               :: matelmode
     character(LEN_MATEL)  :: matelname
     ! include intra-unit cell hopping (Peierls mode only)
     logical               :: intrahop
     ! maximum energy, energy spacing, broadening [eV]
     real(DPk)             :: Emax, dE, delterl
     integer               :: wint_dens   ! #internal freqs. per external freq.
     real(DPk)             :: wint_cutoff ! accuracy for wint range (→Nwx)
     ! outer and inner band windows
     integer               :: bmin_w2k,bmax_w2k, bmin,bmax
     real(DPk)             :: beta,chempot! inv. temp. [1/ev], μ [eV]
     real(DPk)             :: drudesep    ! “Drude” sumrule to this energy [eV]
     logical               :: orbresolv   ! compute orbitally resolved quant.?
     logical               :: selfE       ! use self-energy?
     integer,  allocatable :: iself(:)    ! indices of interacting orbitals
     logical               :: wfrot, shift! rotate WFs?, apply scissors shift?
     real(DPk)             :: Eshift      ! shift energy [eV]
     integer,  allocatable :: ishift(:)   ! indices of bands to be shifted
!!!! “Derived” variables
     real(DPk),allocatable :: w(:), E(:)  ! internal and external freq. grids
     real(DPk)             :: dw          ! int. freq. spacing
     integer               :: Nwx         ! number of extra int. freqs. on
                                          ! either side of [-Ω, 0]
     logical               :: optcond, dos! compute opt. cond.?, DOS?
     logical               :: joint       ! joint DOS mode
     logical               :: Peierls     ! Peierls approximation?
     logical               :: need_umatrix! copmute/read U(k)?
     logical               :: orig_umatrix! use U(k) from ‘case.chk’
     logical               :: mixed       ! include mixed transitions?
     logical               :: read_mommat ! need to read ‘case.mommat2’?
     logical               :: read_vk     ! need to read ‘case.vk?’ [interp.]?
     logical               :: read_vvk    ! need ‘.vvk??’ [mixed tr. interp.]?
     logical               :: read_energy ! need to read ‘case.energy’?
     logical               :: read_ham    ! need to read ‘case.hk’?
  end type inwop_t

  character(*), public, parameter :: &
       fmt_rweights    = "(15I5)",                                         &
       fmt_vr_head     = "('real space dipole elements for ', A)",         &
       fmt_vvr_head    = "('real space mixed dipole elements for ', A)",   &
       fmt_vr          = "(      5I5, 3(2X,2E14.6))",                      &
       fmt_vr_read     = "(3I5, 10X,  3(2X,2E14.6))",                      &
       fmt_vvr         = "(F8.3, 5I5, 6(2X,2E14.6))",                      &
       fmt_vvr_read    = "(      33X, 6(2X,2E14.6))",                      &
       fmt_vk_head     = "(2(I0, 1X), '	# Nk Nwann')",                     &
       fmt_vvk_head    = "(2(I0, 1X), '	# Nk*Nfreq, Nwann')",              &
       fmt_vk_kp       = "(3F12.8)",                                       &
       fmt_vvk_kp      = "(3F12.8, F20.8, '	# kx,ky,kz, ω')",          &
       fmt_vk          = "(200(E20.12))",                                  &
! Fortran 2008 “*()” does not work with older compilers, so we put a
! “ridiculously large” number for the number of bands …
!       fmt_vk          = "(*(E20.12))",                                    &
       fmt_vvk         = fmt_vk

  character(*),         parameter :: svn_rev_str  = '$Rev: 394 $'
  character(*), public, parameter :: woptic_version = &
       "0.98-r" // svn_rev_str(7:len(svn_rev_str)-2)

contains

  subroutine get_mommat_unit(unit, Vx, Vy, Vz, nk, nbmin, nbmax)
!!! Read header line off of ‘unit’ and call get_mommat_1k() [q.v.]
!!! ‘nk’ times.
    integer,      intent(in)                    :: unit, nk, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:,:) :: Vx, Vy, Vz

    integer :: jk

    read(unit,*)

    do jk=1,nk
       call get_mommat_1k(unit, &
            &             Vx(:,:,jk), Vy(:,:,jk), Vz(:,:,jk), &
            &             nbmin, nbmax)
    end do
  end subroutine get_mommat_unit

  subroutine get_mommat_1k(unit, Vkx, Vky, Vkz, nbmin, nbmax)
!!! Read selected momentum matrix elements corresponding to next
!!! k-point.
!!!
!!! unit            connected to a ‘mommat2’ file
!!!                 (read header off first!)
!!!
!!! Vk{x,y,z}       storage for matrix elements
!!!                 (only upper part is filled!)
!!!
!!! nbmin,nbmax     selected bands

    integer,      intent(in)                  :: unit, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:) :: Vkx, Vky, Vkz

    integer :: nbminloc, nbmaxloc, i, j, nb

    nb = nbmax - nbmin + 1

    read(unit,*)
    read(unit,"(T28, 2I5)") nbminloc, nbmaxloc
    read(unit,*)

    ! skip leading junk
    do i = 1, &               ! "little Gauss":
         & (nbmin - nbminloc)*(2*nbmaxloc-nbmin-nbminloc+3)/2
       read(unit,*)
    end do

    ! here comes the juicy part
    do i=nbmin,nbmax
       do j=i,nbmaxloc
          ! skip intermediate junk
          if ((j.lt.nbmin) .or. (j.gt.nbmax)) then
             read(unit, *)
             cycle
          end if

          ! read what we need
          read(unit, '(11X,6E13.6)') Vkx(i-nbmin+1, j-nbmin+1), &
               &                     Vky(i-nbmin+1, j-nbmin+1), &
               &                     Vkz(i-nbmin+1, j-nbmin+1)
       end do
    end do

    ! skip trailing junk up to next k-point
    do i = 1, &
         & (nbmaxloc-nbmax)*(nbmaxloc-nbmax+1)/2
       read(unit,*)
    end do
  end subroutine get_mommat_1k

!!! Wrapper for calling the above with a filename instead of an open
!!! unit
  subroutine get_mommat_fname(fname, Vx, Vy, Vz, nk, nbmin, nbmax)
    use util, only: newunit

    integer,      intent(in)                    :: nk, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:,:) :: Vx, Vy, Vz
    character(*), intent(in)                    :: fname

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_mommat_unit(lun, Vx, Vy, Vz, nk, nbmin, nbmax)
    close(lun)
  end subroutine get_mommat_fname


!!! Read and check ‘inwop’ file
  subroutine inwop_read_unit(lun, inw)
    use clio,  only: croak
    use util,  only: uppercase, string
    use const, only: BUFSZ

    integer,       intent(in)  :: lun
    type(inwop_t), intent(out) :: inw

    integer          :: Nshift, Nself, Nb, i, NE, wmin, wmax, ios
    character(BUFSZ) :: buf

    read(lun, '(A)') buf
    read(buf, *, IOSTAT=ios) inw%mode, inw%matelmode, inw%intrahop
    if (ios == 0) then          ! matelmode given as number
       if (inw%matelmode < 10) inw%matelmode = 10*inw%matelmode

       if (all(inw%matelmode /= MATELMODES)) &
            call croak("unknown matelmode: "//trim(string(inw%matelmode)))
    else                        ! matelmode given as string?
       read(buf, *) inw%mode, inw%matelname, inw%intrahop

       do i = 1, size(MATELNAMES)
          if (uppercase(inw%matelname) == uppercase(MATELNAMES(i))) then
             inw%matelmode = MATELMODES(i)
             goto 1001
          end if
       end do

       call croak("unknown matelmode: "//trim(inw%matelname))

1001   continue
    end if

    unsupported: if (any(inw%matelmode == UNSUPPORTED_MATELMODES)) then
       inw%matelname = "???"

       do i = 1, size(MATELMODES)
          if (inw%matelmode == MATELMODES(i)) then
             inw%matelname = MATELNAMES(i)
          end if
       end do

       call croak("FIXME: matelmode "//trim(string(inw%matelmode))// &
            &     " ("//trim(inw%matelname)//") is currently unsupported")
    end if unsupported

    read(lun,*) inw%Emax, inw%dE, inw%delterl, inw%wint_dens, inw%wint_cutoff

    read(lun, '(A)') buf
    read(buf, *, IOSTAT=ios) inw%bmin,inw%bmax, inw%bmin_w2k,inw%bmax_w2k
    if (ios /= 0) then
       read(buf, *) inw%bmin, inw%bmax
       inw%bmin_w2k = inw%bmin
       inw%bmax_w2k = inw%bmax
    end if

    if (inw%bmin_w2k == 0) inw%bmin_w2k = inw%bmin
    if (inw%bmax_w2k == 0) inw%bmax_w2k = inw%bmax
    if (inw%bmin     == 0) inw%bmin     = inw%bmin_w2k
    if (inw%bmax     == 0) inw%bmax     = inw%bmax_w2k
    read(lun,*) inw%beta, inw%chempot
    read(lun,*) inw%drudesep, inw%orbresolv

    read(lun,*) inw%selfE, Nself
    if (inw%selfE) then
       if (NSelf == 0) then
          NSelf = inw%bmax - inw%bmin + 1
          allocate(inw%iself(Nself))
          inw%iself = (/ (i, i=inw%bmin, inw%bmax) /)
          read(lun,*)
       else
          allocate(inw%iself(Nself))
          read(lun,*) inw%iself
       end if
    else
       Nself=0
       allocate(inw%iself(Nself))
       read(lun,*)
    end if

    read(lun,*) inw%wfrot, inw%shift, inw%Eshift, Nshift
    if (.not. inw%shift .or. inw%Eshift==0) Nshift=0
    allocate(inw%ishift(Nshift))
    if (Nshift>0) read(lun,*) inw%ishift

    !! construct frequency grids
    NE      = ceiling(inw%Emax/inw%dE)
    inw%dw  = inw%dE/inw%wint_dens
    inw%Nwx = find_nw_extra(inw%wint_cutoff, inw%beta, inw%dw)
    wmin    = -NE*inw%wint_dens-inw%Nwx
    wmax    = +NE*inw%wint_dens+inw%Nwx
    allocate(inw%E(   1:NE), &
         &   inw%w(wmin:wmax))
    inw%E = (/ (i*inw%dE, i=1,   NE)   /)
    inw%w = (/ (i*inw%dw, i=wmin,wmax) /)

    if (size(inw%w) == 0) then
       write(buf, "('internal frequency grid is empty' / 3(A,'=',G10.3))") &
            'Emax', inw%Emax, ', dE', inw%dE, ', dw', inw%dw
       call croak(trim(buf))
    end if

    !! consistency checks and derived options
    inw%optcond = .true.
    inw%dos     = .true.
    inw%joint   = .false.
    select case (uppercase(inw%mode))
    case('OPT') ;
    case('DOS') ;  inw%optcond = .false.
    case('JOINT'); inw%joint   = .true.
    case default; call croak("unknown MODE: "//inw%mode)
    end select

    inw%mixed = inw%bmin/=inw%bmin_w2k .or. inw%bmax/=inw%bmax_w2k

    inw%intrahop     = inw%matelmode == MODE_Peierls .and. inw%intrahop
    inw%read_energy  = inw%matelmode >= MODE_LDA      .or. inw%mixed 
    inw%read_mommat  = inw%matelmode >= MODE_optic    .or. &
         &            (inw%matelmode == MODE_interp  .and. inw%mixed)
    inw%read_vk      = inw%matelmode == MODE_interp
    inw%read_vvk     = inw%matelmode == MODE_interp  .and. inw%mixed
    inw%read_ham     = inw%matelmode >  MODE_Peierls .and. &
         &             inw%matelmode <  MODE_LDA
    inw%need_umatrix = inw%matelmode >  MODE_interp  .and. &
         &             inw%matelmode <  MODE_LDA
    inw%orig_umatrix = inw%matelmode == MODE_F_OrigU
    inw%Peierls      = inw%matelmode == MODE_Peierls

    if (Nshift==0) inw%Eshift = 0

    ! Now we discard the “sub-mode part” of matelmode because
    ! woptic_main does not know about sub-modes.
    inw%matelmode = (inw%matelmode / 10) * 10

    if (inw%Emax < 0) &
         call croak("Emax should be ≥0, not " // trim(string(inw%Emax)))
    if (inw%dE <= 0) &
         call croak("dE should be >0, not " // trim(string(inw%dE)))
    if (inw%beta < 0) &
         call croak("beta should be ≥0, not " // trim(string(inw%beta)))
    if (inw%delterl <= 0) &
         call croak("delterl should be >0, not" // trim(string(inw%delterl)))
    if ( inw%bmin_w2k<=0           .or. &
         inw%bmin<inw%bmin_w2k     .or. inw%bmax<inw%bmin .or. &
         inw%bmax_w2k<inw%bmin_w2k .or. inw%bmax>inw%bmax_w2k) &
         call croak("1 ≤ bmin_w2k ≤ bmin ≤ bmax ≤ bmax_w2k must hold")
    if (inw%drudesep < 0) &
         call croak("drudesep should be ≥0, not "//trim(string(inw%drudesep)))
    Nb = inw%bmax_w2k - inw%bmin_w2k + 1
    if (any(inw%ishift < inw%bmin_w2k .or. inw%ishift > inw%bmax_w2k)) &
         call croak("shift indices must be inside outer window")
    if (any(inw%iself  < inw%bmin     .or. inw%iself  > inw%bmax)) &
         call croak("interacting indices must be contained in inner window")
  end subroutine inwop_read_unit

  pure integer function find_nw_extra(cutoff, beta, dw)
    !! The integrals over internal frequency ω for a given external
    !! frequency Ω should span [-Ω-δ, +δ], where δ is determined by
    !! temperature broadening and the desired accuracy.
    !!
    !! This is encapsulated in the factor [f(ω) - f(ω+Ω)]/Ω for the
    !! dynamical, and β/4/cosh²(½βω) for the static quantities.  Here,
    !! we use a single constant δ which is determined by the static
    !! weight factor, since this converges slowest relative to its
    !! maximum value.  I.e., δ is determined by
    !!
    !!    β/4/cosh²(½β δ) = cutoff,
    !!
    !! and find_nw_extra returns the smallest value on the ω grid
    !! larger than δ.
    real(DPk), intent(in) :: cutoff, beta, dw
    real(DPk) :: w

    w = 2/beta * acosh(sqrt(beta/4/cutoff))

    find_nw_extra = ceiling(w/dw)
  end function find_nw_extra

!!! Wrappers for calling the above with a filename instead of an open
!!! unit
  subroutine inwop_read_fname(fname, inw)
    use util, only: newunit

    character(*),   intent(in)  :: fname
    type(inwop_t),  intent(out) :: inw

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call inwop_read_unit(lun, inw)
    close(lun)
  end subroutine inwop_read_fname

  subroutine inwop_read_argstr(arg, inw)
    use util, only: newunit
    use clio, only: argstr

    type(argstr),  intent(in)  :: arg
    type(inwop_t), intent(out) :: inw

    call inwop_read_fname(arg%s, inw)
  end subroutine inwop_read_argstr

!!! noninteracting spectral function
!!!
!!! To allow inlining, this should be contained in the user file.

!!!  pure elemental real(DPk) function spectral(omega, epsilon, delterl)
!!!    use const, only: PI
!!!
!!!    real(DPk), intent(in) :: omega, epsilon, delterl
!!!
!!!    spectral = delterl/PI / ((omega-epsilon)**2 + delterl**2)
!!!  end function spectral
end module woptic

!---------------- Read/write files that may be binary or plain text -----------
module maybebin
  use const, only: DPk
  implicit none
  private
  public :: maybin_open, maybin_read, maybin_write, maybin_get, maybin_set, &
       formatted

  logical         :: first=.true.
  logical, target :: binary

  interface maybin_read
     module procedure myread_, myread_i, myread_r, myread_c, myread_a
  end interface maybin_read

  interface maybin_write
     module procedure mywrite_, mywrite_i, mywrite_r, mywrite_c, mywrite_a
  end interface maybin_write

contains

  subroutine maybin_set(bin)
    logical, intent(in) :: bin

    binary = bin
    first  = .false.
  end subroutine maybin_set

  pure logical function maybin_get()
    maybin_get = binary
  end function maybin_get

  subroutine maybin_open(unit, file, get, set, status)
    use clio,  only: croak

    integer,      intent(in)                   :: unit
    character(*), intent(in)                   :: file
    logical,      intent(out),optional, target :: get
    logical,      intent(in), optional, target :: set
    character(*), intent(in), optional         :: status

    character(100)   :: buf
    character        :: c
    character(12)    :: f
    logical          :: detect
    logical, pointer :: b
    character(7)     :: stat

    if (present(get) .and. present(set)) &
         call croak("both `get_bin' and `set_bin' present in call to "// &
         &          "maybin_open()")

    if      (present(get)) then
       detect = .true.
       b      => get
    else if (present(set)) then
       detect = .false.
       b      => set
    else
       detect =  first
       b      => binary

       first = .false.
    end if

    if (detect) then
       b = .false.

       !! We consider the file plain text if the first N bytes are
       !! printable ASCII characters.
       !!
       !! FIXME: do something less hacky.
       open (unit=unit, file=file, status='old', access='stream')
       read (unit, END=1001) buf

       b = .not. printable(buf)
       goto 1002

1001   rewind(unit)
       do
          read(unit, END=1002) c
          if (.not. printable(c)) then
             b = .true.
             goto 1002
          end if
       end do

1002   close(unit)
    end if

    if (b) then
       f = 'unformatted'
    else
       f = 'formatted'
    end if

    if (present(status)) then
       stat=status
    else
       stat='old'
    end if

    open(unit=unit, file=file, status=stat, form=f)
  end subroutine maybin_open

  pure logical function printable(string)
    character(*), intent(in) :: string
    character                :: c
    integer                  :: i

    printable = .true.

    do i=1, len_trim(string)
       c = string(i:i)
       select case (iachar(c))
       case (0:8, 14:31, 126:)
          printable = .false.
          return
       end select
    end do
  end function printable

!!! Check if a file is opened in FORMATTED mode.
!!!
!!! NB1: Execution times reveal no difference between this check and a
!!! simple logical.
!!!
!!! NB2: If the File is not opened, this returns .false.  (Inquire()
!!! returns "UNDEFINED".)
  logical function formatted(unit)
    integer, intent(in) :: unit
    character(12)       :: f

    inquire(unit, FORM=f)

    formatted = f=='FORMATTED'
  end function formatted

  subroutine myread_(unit)
    integer,      intent(in)            :: unit

    if (.not. formatted(unit)) then
       read(unit)
    else
       read(unit,*)
    end if
  end subroutine myread_

  subroutine myread_i(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    integer,      intent(out)           :: data(:)

    if (.not. formatted(unit)) then
       read(unit) data
    else if (present(fmt)) then
       read(unit, fmt) data(:)
    else
       read(unit, *)   data(:)
    end if
  end subroutine myread_i

  subroutine myread_r(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    real(DPk),    intent(out)           :: data(:)

    if (.not. formatted(unit)) then
       read(unit) data
    else if (present(fmt)) then
       read(unit, fmt) data(:)
    else
       read(unit, *)   data(:)
    end if
  end subroutine myread_r

  subroutine myread_c(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    complex(DPk), intent(out)           :: data(:)

    if (.not. formatted(unit)) then
       read(unit) data
    else if (present(fmt)) then
       read(unit, fmt) data(:)
    else
       read(unit, *)   data(:)
    end if
  end subroutine myread_c

  subroutine myread_a(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    character(*), intent(out)           :: data

    if (.not. formatted(unit)) then
       read(unit) data
    else if (present(fmt)) then
       read(unit, fmt) data
    else
       read(unit, *)   data
    end if
  end subroutine myread_a

  subroutine mywrite_(unit)
    integer,      intent(in)            :: unit

    if (.not. formatted(unit)) then
       write(unit)
    else
       write(unit,*)
    end if
  end subroutine mywrite_

  subroutine mywrite_i(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    integer,      intent(in)            :: data(:)

    if (.not. formatted(unit)) then
       write(unit) data
    else if (present(fmt)) then
       write(unit, fmt) data(:)
    else
       write(unit, *)   data(:)
    end if
  end subroutine mywrite_i

  subroutine mywrite_r(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    real(DPk),    intent(in)            :: data(:)

    if (.not. formatted(unit)) then
       write(unit) data
    else if (present(fmt)) then
       write(unit, fmt) data(:)
    else
       write(unit, *)   data(:)
    end if
  end subroutine mywrite_r

  subroutine mywrite_c(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    complex(DPk), intent(in)            :: data(:)

    if (.not. formatted(unit)) then
       write(unit) data
    else if (present(fmt)) then
       write(unit, fmt) data(:)
    else
       write(unit, *)   data(:)
    end if
  end subroutine mywrite_c

  subroutine mywrite_a(unit, data, fmt)
    integer,      intent(in)            :: unit
    character(*), intent(in),  optional :: fmt
    character(*), intent(in)            :: data

    if (.not. formatted(unit)) then
       write(unit) data
    else if (present(fmt)) then
       write(unit, fmt) data
    else
       write(unit, *)   data
    end if
  end subroutine mywrite_a

end module maybebin

!--------------- Files and units for woptic          ------------------------
module woptic_io
  use const, only: BUFSZ
  implicit none
  public

  integer :: i

  integer,          parameter   :: unit_inwop=10
  character(*),     parameter   ::  suf_inwop='.inwop'
  character(BUFSZ)              ::   fn_inwop
  integer,          parameter   :: unit_fermi=11
  character(*),     parameter   ::  suf_fermi='.fermi'
  character(BUFSZ)              ::   fn_fermi
  integer,          parameter   :: unit_ham=12
  character(*),     parameter   ::  suf_ham='.hk'
  character(BUFSZ)              ::   fn_ham
  integer,          parameter   :: unit_intrahop=13
  character(*),     parameter   ::  suf_intrahop='.intrahop'
  character(BUFSZ)              ::   fn_intrahop
  integer,          parameter   :: unit_wfrot=14
  character(*),     parameter   ::  suf_wfrot='.wfrot'
  character(BUFSZ)              ::   fn_wfrot
  integer,          parameter   :: unit_hr=15
  character(*),     parameter   ::  suf_hr='_hr.dat'
  character(BUFSZ)              ::   fn_hr
  integer,          parameter   :: unit_chk=16
  character(*),     parameter   ::  suf_chk='.chk'
  character(BUFSZ)              ::   fn_chk
  integer,          parameter   :: unit_inwf=17
  character(*),     parameter   ::  suf_inwf='.inwf'
  character(BUFSZ)              ::   fn_inwf
  integer,          parameter   :: unit_struct=20
  character(*),     parameter   ::  suf_struct='.struct'
  character(BUFSZ)              ::   fn_struct
  integer,          parameter   :: unit_symop=21
  character(*),     parameter   ::  suf_symop='.symop'
  character(BUFSZ)              ::   fn_symop
  integer,          parameter   :: unit_energy=22
  character(*),     parameter   ::  suf_energy='.energy'
  character(BUFSZ)              ::   fn_energy
  integer,          parameter   :: unit_mommat=23
  character(*),     parameter   ::  suf_mommat='.mommat2'
  character(BUFSZ)              ::   fn_mommat
  integer,          parameter   :: unit_vr=27
  character(*),     parameter   ::  suf_vr='.vr'
  character(BUFSZ)              ::   fn_vr
  integer,          parameter   :: unit_vk(3)=(/ 271, 272, 273 /)
  character(*),     parameter   ::  suf_vk(3)=(/ '.vkx', '.vky', '.vkz' /)
  character(BUFSZ)              ::   fn_vk(3)
  integer,          parameter   :: unit_vvr=28
  character(*),     parameter   ::  suf_vvr='.vvr'
  character(BUFSZ)              ::   fn_vvr
  integer,          parameter   :: unit_vvk(6)=(/ (i, i=281,286) /)
  character(*),     parameter   ::  suf_vvk(6)= &
       (/ '.vvkxx', '.vvkxy', '.vvkxz', '.vvkyy', '.vvkyz', '.vvkzz' /)
  character(BUFSZ)              ::   fn_vvk(6)
  integer,          parameter   :: unit_klist=30
  character(*),     parameter   ::  suf_klist='.klist'
  character(BUFSZ)              ::   fn_klist
  integer,          parameter   :: unit_fklist=31
  character(*),     parameter   ::  suf_fklist='.klist_full'
  character(BUFSZ)              ::   fn_fklist
  integer,          parameter   :: unit_kadd=32
  character(*),     parameter   ::  suf_kadd='.klist_add'
  character(BUFSZ)              ::   fn_kadd
  integer,          parameter   :: unit_tet=40
  character(*),     parameter   ::  suf_tet='.tetra'
  character(BUFSZ)              ::   fn_tet
  integer,          parameter   :: unit_ftet=41
  character(*),     parameter   ::  suf_ftet='.tetra_full'
  character(BUFSZ)              ::   fn_ftet
  integer,          parameter   :: unit_map=42
  character(*),     parameter   ::  suf_map='.map'
  character(BUFSZ)              ::   fn_map
  integer,          parameter   :: unit_voe=43
  character(*),     parameter   ::  suf_voe='.voe'
  character(BUFSZ)              ::   fn_voe
  integer,          parameter   :: unit_optcond=51
  character(*),     parameter   ::  suf_optcond='.optcondw'
  character(BUFSZ)              ::   fn_optcond
  integer,          parameter   :: unit_optorb(6)=(/511,512,513,514,515,516/)
  character(*),     parameter   ::  suf_optorb(6)= &
       (/ suf_optcond//'_orbxx', suf_optcond//'_orbxy', suf_optcond//'_orbxz',&
       &  suf_optcond//'_orbyy', suf_optcond//'_orbyz', suf_optcond//'_orbzz'/)
  character(BUFSZ)              ::   fn_optorb(6)
  integer,          parameter   :: unit_contr=52
  character(*),     parameter   ::  suf_contr='.kcontribw'
  character(BUFSZ)              ::   fn_contr
  integer,          parameter   :: unit_wdos=53
  character(*),     parameter   ::  suf_wdos='.wdos'
  character(BUFSZ)              ::   fn_wdos
  integer,          parameter   :: unit_doscontr=54
  character(*),     parameter   ::  suf_doscontr='.wdoskcontribw'
  character(BUFSZ)              ::   fn_doscontr
  integer,          parameter   :: unit_K1=55
  character(*),     parameter   ::  suf_K1='.K1w'
  character(BUFSZ)              ::   fn_K1
  integer,          parameter   :: unit_selfE=56
  character(*),     parameter   ::  suf_selfE='.selfE'
  character(BUFSZ)              ::   fn_selfE
  integer,          parameter   :: unit_outwop=60
  character(*),     parameter   ::  suf_outwop='.outputwop'
  character(BUFSZ)              ::   fn_outwop
  integer,          parameter   :: unit_outref=61
  character(*),     parameter   ::  suf_outref='.outputref'
  character(BUFSZ)              ::   fn_outref
  integer,          parameter   :: unit_outvr=62
  character(*),     parameter   ::  suf_outvr='.outputvr'
  character(BUFSZ)              ::   fn_outvr
  integer,          parameter   :: unit_outvk=63
  character(*),     parameter   ::  suf_outvk='.outputvk'
  character(BUFSZ)              ::   fn_outvk
contains
  subroutine print1or3(lun, label, unit, tensor)
    use const, only: DPk
    implicit none

    integer,      intent(in) :: lun
    character(*), intent(in) :: label, unit
    real(DPk),    intent(in) :: tensor(3,3)

    character(*), parameter :: &
         fmt1 = '(A12, " [", A, "]", F12.3)', &
         fmt3 = '(A11, "(xx,yy,zz) [", A, "]", 3F12.3)'

    real(DPk) :: xx, yy, zz

    xx=tensor(1,1); yy=tensor(2,2); zz=tensor(3,3)

    if (max(abs((xx-yy)/xx), abs((xx-zz)/xx), abs((zz-yy)/xx)) > 1e-3) then
       write(lun, fmt3) label, unit, xx, yy, zz
    else
       write(lun, fmt1) label, unit, xx
    end if
  end subroutine print1or3

  subroutine set_casename(file, band)
    use clio,  only: argstr

    implicit none

    type(argstr), intent(in)           :: file
    logical,      intent(in), optional :: band

    character(5) :: suf
    integer      :: i
    suf=''
    if (present(band)) then
       if (band) suf='_band'
    end if

    fn_inwop    =trim(file%s)//trim(suf_inwop    )
    fn_fermi    =trim(file%s)//trim(suf_fermi    )
    fn_ham      =trim(file%s)//trim(suf_ham      )
    fn_intrahop =trim(file%s)//trim(suf_intrahop )
    fn_wfrot    =trim(file%s)//trim(suf_wfrot    )
    fn_hr       =trim(file%s)//trim(suf_hr       )
    fn_chk      =trim(file%s)//trim(suf_chk      )
    fn_inwf     =trim(file%s)//trim(suf_inwf     )
    fn_struct   =trim(file%s)//trim(suf_struct   )
    fn_symop    =trim(file%s)//trim(suf_symop    )
    fn_energy   =trim(file%s)//trim(suf_energy   )
    fn_mommat   =trim(file%s)//trim(suf_mommat   )
    fn_klist    =trim(file%s)//trim(suf_klist    )//suf
    fn_vr       =trim(file%s)//trim(suf_vr       )
    fn_vvr      =trim(file%s)//trim(suf_vvr      )
    fn_fklist   =trim(file%s)//trim(suf_fklist   )
    fn_kadd     =trim(file%s)//trim(suf_kadd     )
    fn_tet      =trim(file%s)//trim(suf_tet      )
    fn_ftet     =trim(file%s)//trim(suf_ftet     )
    fn_map      =trim(file%s)//trim(suf_map      )
    fn_voe      =trim(file%s)//trim(suf_voe      )
    fn_optcond  =trim(file%s)//trim(suf_optcond  )
    fn_contr    =trim(file%s)//trim(suf_contr    )//suf
    fn_wdos     =trim(file%s)//trim(suf_wdos     )//suf
    fn_doscontr =trim(file%s)//trim(suf_doscontr )
    fn_K1       =trim(file%s)//trim(suf_K1       )//suf
    fn_selfE    =trim(file%s)//trim(suf_selfE    )//suf
    fn_outwop   =trim(file%s)//trim(suf_outwop   )
    fn_outref   =trim(file%s)//trim(suf_outref   )
    fn_outvr    =trim(file%s)//trim(suf_outvr    )
    fn_outvk    =trim(file%s)//trim(suf_outvk    )
    do i=1,size(fn_optorb)
       fn_optorb(i)=trim(file%s)//trim(suf_optorb(i))
    end do
    do i=1,size(fn_vk)
       fn_vk(i)    =trim(file%s)//trim(suf_vk(i)    )
    end do
    do i=1,size(fn_vvk)
       fn_vvk(i)   =trim(file%s)//trim(suf_vvk(i)   )
    end do
  end subroutine set_casename
end module woptic_io


!! Time-stamp: <2015-06-02 12:09:35 assman@faepop23.tu-graz.ac.at>
