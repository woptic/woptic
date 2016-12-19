!!! woptic/src/compute_vr.f90
!!!
!!!    Computes V(r), the dipole matrix-elements in real space.
!!!
!!! Copyright 2012      Philipp Wissgott
!!!           2013-2016 Elias Assmann
!!!
!!! FIXME: This program follows the WET principle too much.
!!!

program compute_vr
  use util,      only: ptime, string
  use const,     only: DPk, Ryd_eV
  use clio,      only: croak, argstr, fetcharg
  use structmod, only: struct_t, struct_read
  use inwfmod,   only: inwf_t, inwf_read
  use maybebin,  only: maybin_open, maybin_read, maybin_write, &
       maybin_set, maybin_get
  use woptic,    only: get_mommat_1k, inwop_t, inwop_read, &
       fmt_rweights, fmt_vr_head, fmt_vvr_head, fmt_vr, fmt_vvr
  use Wannier90, only: chk_t, chk_read
  use woptic_io, only: set_casename, &
       ! input:
                    unit_mommat,          unit_hr, unit_wfrot,              &
        suf_inwop,   suf_mommat, suf_chk,  suf_hr,                          &
         fn_inwop,    fn_mommat,  fn_chk,   fn_hr,   fn_wfrot,              &
       ! input for mixed transitions / disentanglement:
       unit_energy, unit_fermi,                                             &
        suf_energy,  suf_fermi, suf_struct, suf_inwf,                       &
         fn_energy,   fn_fermi,  fn_struct,  fn_inwf,                       &
       ! output:
       unit_outvr,  unit_vr,    unit_vvr,                                   &
        suf_outvr,   suf_vr,     suf_vvr,                                   &
         fn_outvr,    fn_vr,      fn_vvr

  implicit none

  character(*), parameter :: rev_str = "$version: v0.1.0-70-g3e3f99f$"
  character(*), parameter :: woptic_version = rev_str(11 : len (rev_str)-1)

!!! Formats for log file
  character(*), parameter ::            &
       fmtII = "(A40, I6, I5)",         &
       fmtA  = "(A40, A6)"

!!! Formats for output
  character(*), parameter ::                                          &
       fmt_help   = '(A, T10, "case", A, T25, A, T50, A)',            &
       fmt_hlp2   = '(A, T10,         A, T25, A, T50, A)'

!!! Command line arguments and input files
  type(argstr)   :: arg, case
  type(inwop_t)  :: inw
  type(chk_t)    :: chk
  type(struct_t) :: stru        ! we need NNEQ for mixed transitions
  type(inwf_t)   :: inwf        ! we need bmin,bmax for disentanglement

!!! Options
  logical      :: no_umat=.false., have_so=.false.
  character(2) :: updn=''

!!! Loop variables
  integer :: i, v,w, k, R, f, ii, iarg

!!! Integer parameters
  ! Nb [outer window] ≥ Nwann [inner window]; Nb_dis: disentanglement window;
  ! Nf: internal frequencies; Nb_nonint = Nb - Nwann
  ! bmin,bmax: outer window in Wien2k band indices
  ! wmin,wmax: indices of the Wannier window inside the outer window
  integer :: Nb,bmin,bmax, Nwann,wmin,wmax, Nb_dis,vmin,vmax, Nf,fmin,fmax
  integer :: Nb_nonint, Nk, NR, Nbloc

!!! Temporaries
  real(DPk) :: rdotk, EFermi, mxr

  real(DPk),    allocatable, dimension(:)   :: tmprot
  real(DPk),    allocatable, dimension(:,:) :: rpt_cart, kpt_cart, bands
  complex(DPk), allocatable, dimension(:,:) :: rotmat, fac
  ! inonint: indices of the non-Wannier states inside the outer window
  integer,      allocatable, dimension(:)   :: rweights, inonint
  integer,      allocatable, dimension(:,:) :: rpt_frac

!!! Temporary matrices
  ! the matrix transfoming Bloch to Wannier states; U_iv(k) or
  ! V_in(k) * U_nv(k) [with disentanglement]
  complex(DPk), pointer,     dimension(:,:)     :: UVmat
  ! spectral function matrix A_ii(k, ω)
  real(DPk),    allocatable, dimension(:)       :: A

!!! Momentum matrix elements V
  ! V_ij(k) [as read]
  complex(DPk), allocatable, dimension(:,:)     :: Vx,    Vy,    Vz
  ! V_wv(k) = UVmat+ * V * UVmat [transformed into Wannier gauge]
  complex(DPk), allocatable, dimension(:,:,:)   :: UVUx,  UVUy,  UVUz
  ! W_wv(R) [in direct space]
  complex(DPk), allocatable, dimension(:,:,:)   :: VRx,   VRy,   VRz

!!! For mixed transitions:
  ! mixed momentum matrix elements V_iw
  complex(DPk), allocatable, dimension(:,:)     :: VWBx, VWBy, VWBz
  ! V_wi(k) V_iv(k)
  complex(DPk), allocatable, dimension(:,:,:)   :: &
       UVVUxx, UVVUxy, UVVUxz, UVVUyy, UVVUyz, UVVUzz
  ! W_wv(k, ω) = Sum_i V_wi(k) A_ii(k,ω) V_iv(k)
  complex(DPk), allocatable, dimension(:,:,:,:) :: &
       VAVxx,   VAVxy,  VAVxz,  VAVyy,  VAVyz,  VAVzz
  ! W_wv(R, ω) = F_k[W_wv(k,ω)]
  complex(DPk), allocatable, dimension(:,:,:,:) :: &
       VVRxx,   VVRxy,  VVRxz,  VVRyy,  VVRyz,  VVRzz

  call maybin_set(.true.)

!!! Argument parsing
  if (command_argument_count() < 1) &
       call croak('Usage: compute_vr [OPTIONS] case')

  do iarg=1,command_argument_count()
     call fetcharg(iarg, arg)
     select case(arg%s)
     case ('-h', '-H', '-help', '--help')
        write(*,'(A)') &
             "compute_vr: compute dipole matrix elements in direct space"
        write(*,*) 
        write(*,fmt_hlp2)"USAGE:", "compute_vr [--text] [--up|--dn] [--so]&
             & CASE"
        write(*,*) 
        write(*,fmt_help)"INPUT:", suf_inwop, "input file"
        write(*,fmt_help)" (W90)", suf_chk,   "Wannier90 checkpoint file"
        write(*,fmt_help)"",       suf_hr,    "Hamiltonian in direct space"
        write(*,fmt_help)" (W2k)", suf_mommat,"momentum matrix elements"
        write(*,fmt_help)"",       suf_struct,"Wien2k master input file",  &
             "[mixed transitions]"
        write(*,fmt_help)"",       suf_energy,"energies from lapw1",       &
             "[mixed transitions]"
        write(*,fmt_help)" (w2w)", suf_fermi, "Fermi energy",              &
             "[mixed transitions]"
        write(*,fmt_help)"",       suf_inwf,  "w2w input file",            &
             "[disentanglement]"
        write(*,*) 
        write(*,fmt_help)"OUTPUT:",suf_outvr,"log file"
        write(*,fmt_help)"",       suf_vr,   "U^+ V U"
        write(*,fmt_help)"",       suf_vvr,  "V A V", "[mixed transitions]"
        write(*,*) 
        write(*,fmt_hlp2)"OPTIONS:", "--help, -h"
        write(*,fmt_hlp2)"",         "--version, -v"
        write(*,fmt_hlp2)"",         "--text, -t", "plain-text case.vvr"
        write(*,fmt_hlp2)"",         "--up|--dn",  "spin-polarized mode"
        stop

     case ('-v', '-V', '-version', '--version')
        print '("compute_vr ", A)', WOPTIC_VERSION
        call exit(0)

     case('--noU', '-noU')
        stop "FIXME: --noU is currently broken."
        no_umat = .true.

     case('-t', '--text')
        call maybin_set(.false.)

     case ('-up', '--up')
        updn='up'
     case ('-dn', '--dn')
        updn='dn'
     case ('-so', '--so')
        have_so = .true.

     case default
        if (arg%s(1:1) == '-') &
             call croak('unknown option: '//trim(arg%s))

        call fetcharg(iarg, case)
     end select
  end do

  call set_casename(case, HAVE_SO=have_so, UPDN=updn)

  call ptime(UNIT=unit_outvr)


  open (unit_outvr, FILE=fn_outvr, STATUS='replace')
  write(unit_outvr, '(A)') &
       "------------- Dipole Matrix Elements in Direct Space -------------"

  call inwop_read(fn_inwop, inw)

  call chk_read(fn_chk, chk, READ_MMN=.false.)

  if (chk%have_disentangled .and. inw%mixed) &
       call croak('mixed transitions with disentanglement not implemented')

  Nk = chk%num_kpts

  allocate(kpt_cart(3, Nk))

  Nwann = inw%bmax - inw%bmin + 1

  if (chk%num_wann /= Nwann) &
       call croak("number of Wannier functions inconsistent (‘"&
       &          //trim(fn_chk)// "’ vs. ‘" //trim(fn_inwop)// "’)")

  if (chk%have_disentangled) then
     call inwf_read(fn_inwf, inwf)
     bmin   = inwf%bmin
     bmax   = inwf%bmax
     Nb_dis = bmax - bmin + 1

     ! indices for Wannier-Wannier part of V(k)
     vmin = 1
     vmax = Nb_dis
  else
     Nb_dis = 0
     bmin   = inw%bmin_w2k
     bmax   = inw%bmax_w2k

     ! These variables are only needed in the mixed case 
     wmin = inw%bmin - inw%bmin_w2k + 1
     wmax = wmin + Nwann - 1
     Nb_nonint = bmax-bmin+1 - Nwann

     vmin=wmin; vmax=wmax
  end if

  Nb = bmax - bmin + 1

  fmin = lbound(inw%w,1)
  fmax = merge(ubound(inw%w,1), fmin-1, inw%mixed)
  Nf   = fmax - fmin + 1

  if (inw%wfrot) then
     call croak('WF rotation not implemented')
     write(unit_outvr,*) "rotating...."
     open(unit_wfrot, FILE=fn_wfrot, STATUS='old')
     allocate(rotmat(Nwann,Nwann), tmprot(Nwann*2))
     do v=1,Nwann
        read(unit_wfrot,*) tmprot
        do w=1,Nwann
           rotmat(v,w) = cmplx(tmprot(2*w-1),tmprot(2*w), DPk)
        end do
     end do
     close(unit_wfrot)
  end if

  write(unit_outvr,fmtII) &
       "number of k-points/Wannier functions:", Nk, Nwann
  write(unit_outvr,fmtA ) "apply U matrices?",  merge(" no", "yes", no_umat)
  write(unit_outvr,fmtA ) "mixed transitions?", merge("yes", " no", inw%mixed)
  write(unit_outvr,fmtA ) "disentanglement?",   merge("yes", " no", &
       chk%have_disentangled)


!!! Things special to cases with outer window: read E_F, struct (for
!!! NNEQ), bands; allocate frequencies
  if (inw%mixed) then
     open (unit_fermi, FILE=fn_fermi, STATUS='old')
     read (unit_fermi,*) EFermi
     close(unit_fermi)

     call struct_read(fn_struct, stru)

     allocate(bands(Nb_nonint, Nk))
     open (unit_energy, FILE=fn_energy, STATUS='old')
     do ii = 1, stru%nneq
        read(unit_energy,*)
        read(unit_energy,*)
     enddo

     do k = 1, Nk
        read(unit_energy,"(73X,i6)") Nbloc
        if (Nbloc < bmax) &
             call croak('energy file does not span outer band window: '// &
             trim(string(Nbloc)) //' < '// trim(string(bmax)) &
             //" @ k-point #"// trim(string(k)))

        do i=1,bmin-1
           read(unit_energy,*)
        end do

        do i=1,wmin-1
           read(unit_energy,*) ii, mxr
           bands(i, k) = (mxr-EFermi)*Ryd_eV
        enddo

        do i=wmin,wmax
           read(unit_energy,*)
        end do

        do i=wmin, Nb-Nwann
           read(unit_energy,*) ii, mxr
           bands(i, k) = (mxr-EFermi)*Ryd_eV
        end do

        do i = bmax+1, Nbloc
           read(unit_energy,*)
        end do
     enddo
     close(unit_energy)

     write(unit_outvr,fmtII) "number of frequencies:", Nf
     write(unit_outvr,fmtA ) "binary output?", merge("yes"," no", maybin_get())
  end if
  write(unit_outvr,*)


!!! Read case_hr.dat for vectors and weights in direct space
  open(unit_hr, FILE=fn_hr, STATUS='old')
  read(unit_hr,*)
  read(unit_hr,*)
  read(unit_hr,*) NR
  allocate(rweights(NR), rpt_frac(NR,3), rpt_cart(NR,3))

  rweights = 0
  do R=1,NR/15
     read(unit_hr,*) rweights((R-1)*15+1:min(R*15,NR))
  end do
  if (mod(NR, 15) /= 0) then
     R = NR/15+1
     read(unit_hr,*) rweights((R-1)*15+1:min(R*15,NR))
  end if

  do R=1,NR
     do v=1,Nwann
        do w=1,Nwann
           read(unit_hr,*) rpt_frac(R,:)
        end do
     end do
  end do
  close(unit_hr)

  !! Transform fractional to cartesian lattice vectors.  We may as
  !! well work in Angstrom here, the factor falls out in k*R
  !! anyway.  NB: recip_lattice includes the factor of 2*pi
  kpt_cart = matmul(transpose(chk%recip_lattice), chk%kpt_latt)
  rpt_cart = matmul(rpt_frac, chk%real_lattice)


!!! Allocate dynamic storage
  allocate(fac(NR, Nk))

  !! k-space
  allocate(Vx(Nb, Nb), &
       &   Vy(Nb, Nb), &
       &   Vz(Nb, Nb))

  allocate(UVUx(Nwann, Nwann, Nk), &
       &   UVUy(Nwann, Nwann, Nk), &
       &   UVUz(Nwann, Nwann, Nk))

  if (inw%mixed) then
     allocate(UVVUxx(Nwann, Nb_nonint, Nwann), &
          &   UVVUxy(Nwann, Nb_nonint, Nwann), &
          &   UVVUxz(Nwann, Nb_nonint, Nwann), &
          &   UVVUyy(Nwann, Nb_nonint, Nwann), &
          &   UVVUyz(Nwann, Nb_nonint, Nwann), &
          &   UVVUzz(Nwann, Nb_nonint, Nwann))

     allocate(VAVxx(Nwann, Nwann, Nk, fmin:fmax))
     allocate(VAVxy(Nwann, Nwann, Nk, fmin:fmax))
     allocate(VAVxz(Nwann, Nwann, Nk, fmin:fmax))
     allocate(VAVyy(Nwann, Nwann, Nk, fmin:fmax))
     allocate(VAVyz(Nwann, Nwann, Nk, fmin:fmax))
     allocate(VAVzz(Nwann, Nwann, Nk, fmin:fmax))
  end if

  !! R-space
  allocate(VRx(Nwann,Nwann,NR)); VRx = 0
  allocate(VRy(Nwann,Nwann,NR)); VRy = 0
  allocate(VRz(Nwann,Nwann,NR)); VRz = 0

  if (inw%mixed) then
     allocate(VVRxx(Nwann, Nwann, NR, fmin:fmax)); VVRxx = 0
     allocate(VVRxy(Nwann, Nwann, NR, fmin:fmax)); VVRxy = 0
     allocate(VVRxz(Nwann, Nwann, NR, fmin:fmax)); VVRxz = 0
     allocate(VVRyy(Nwann, Nwann, NR, fmin:fmax)); VVRyy = 0
     allocate(VVRyz(Nwann, Nwann, NR, fmin:fmax)); VVRyz = 0
     allocate(VVRzz(Nwann, Nwann, NR, fmin:fmax)); VVRzz = 0
  end if

  if (chk%have_disentangled) allocate(UVmat(Nb_dis, Nwann))

  if (inw%mixed) then
     allocate(VWBx   (Nb_nonint, Nwann), &
          &   VWBy   (Nb_nonint, Nwann), &
          &   VWBz   (Nb_nonint, Nwann), &
          &   inonint(Nb_nonint),        &
          &   A      (Nb_nonint))
     inonint = (/ (i, i=1,wmin-1), (i, i=wmax+1,Nb) /)
  end if

  call ptime("reading files")


!!! Read case.mommat2, populating UVU[xyz]
  open(unit_mommat, FILE=fn_mommat, STATUS='old')
  read(unit_mommat,*)
  kpoint: do k=1,Nk
     call get_mommat_1k(unit_mommat, Vx, Vy, Vz, bmin, bmax)

     if (chk%have_disentangled) then
        UVmat = matmul(chk%u_matrix_opt(:,:,k), chk%u_matrix(:,:,k))
     else
        UVmat => chk%u_matrix(:,:,k)
     end if

     UVUx(:,:,k) = matmul(transpose(conjg(UVmat)), &
          &               hemul(Vx(vmin:vmax, vmin:vmax), UVmat))
     UVUy(:,:,k) = matmul(transpose(conjg(UVmat)), &
          &               hemul(Vy(vmin:vmax, vmin:vmax), UVmat))
     UVUz(:,:,k) = matmul(transpose(conjg(UVmat)), &
          &               hemul(Vz(vmin:vmax, vmin:vmax), UVmat))

     mixed: if (inw%mixed) then
        !! Complete V by hermiticity
        forall (i=1:Nb_nonint, inonint(i)<wmin)
           VWBx(i, 1:Nwann) = Vx(inonint(i), wmin:wmax)
           VWBy(i, 1:Nwann) = Vy(inonint(i), wmin:wmax)
           VWBz(i, 1:Nwann) = Vz(inonint(i), wmin:wmax)
        end forall

        forall (i=1:Nb_nonint, inonint(i)>wmax)
           VWBx(i, 1:Nwann) = conjg(Vx(wmin:wmax, inonint(i)))
           VWBy(i, 1:Nwann) = conjg(Vy(wmin:wmax, inonint(i)))
           VWBz(i, 1:Nwann) = conjg(Vz(wmin:wmax, inonint(i)))
        end forall

        VWBx = matmul(VWBx, UVmat)
        VWBy = matmul(VWBy, UVmat)
        VWBz = matmul(VWBz, UVmat)

        do v=1,Nwann
           do w=1,Nwann
              UVVUxx(v,:,w) = conjg(VWBx(:,v)) * VWBx(:,w)
              UVVUxy(v,:,w) = conjg(VWBx(:,v)) * VWBy(:,w)
              UVVUxz(v,:,w) = conjg(VWBx(:,v)) * VWBz(:,w)
              UVVUyy(v,:,w) = conjg(VWBy(:,v)) * VWBy(:,w)
              UVVUyz(v,:,w) = conjg(VWBy(:,v)) * VWBz(:,w)
              UVVUzz(v,:,w) = conjg(VWBz(:,v)) * VWBz(:,w)
           end do
        end do

        VAV_sum: do f = fmin,fmax
           A = spectral(inw%w(f), bands(:,k), inw%delterl)

           forall(v=1:Nwann, w=1:Nwann)
              VAVxx(v,w, k,f) = sum(UVVUxx(v, :, w) * A)
              VAVxy(v,w, k,f) = sum(UVVUxy(v, :, w) * A)
              VAVxz(v,w, k,f) = sum(UVVUxz(v, :, w) * A)
              VAVyy(v,w, k,f) = sum(UVVUyy(v, :, w) * A)
              VAVyz(v,w, k,f) = sum(UVVUyz(v, :, w) * A)
              VAVzz(v,w, k,f) = sum(UVVUzz(v, :, w) * A)
           end forall
        end do VAV_sum
     end if mixed
  end do kpoint
  close(unit_mommat)
  call ptime("reading mommat2 & applying U(k)")


!!! Fourier trafo to R
  furry_transformer: do R=1,NR
     do k=1,Nk
        rdotk    = dot_product(kpt_cart(:,k), rpt_cart(R,:))
        fac(R,k) = exp(-(0,1)*rdotk)/real(Nk,DPk)

        VRx(:,:,R) = VRx(:,:,R) + fac(R,k)*UVUx(:,:,k)
        VRy(:,:,R) = VRy(:,:,R) + fac(R,k)*UVUy(:,:,k)
        VRz(:,:,R) = VRz(:,:,R) + fac(R,k)*UVUz(:,:,k)
     end do
  end do furry_transformer
  call ptime('FT to V(R)')

  furrier_transformer: do f = fmin,fmax
     do R=1,NR
        do k=1,Nk
           VVRxx(:,:,R,f) = VVRxx(:,:,R,f) + fac(R,k)*VAVxx(:,:,k,f)
           VVRxy(:,:,R,f) = VVRxy(:,:,R,f) + fac(R,k)*VAVxy(:,:,k,f)
           VVRxz(:,:,R,f) = VVRxz(:,:,R,f) + fac(R,k)*VAVxz(:,:,k,f)
           VVRyy(:,:,R,f) = VVRyy(:,:,R,f) + fac(R,k)*VAVyy(:,:,k,f)
           VVRyz(:,:,R,f) = VVRyz(:,:,R,f) + fac(R,k)*VAVyz(:,:,k,f)
           VVRzz(:,:,R,f) = VVRzz(:,:,R,f) + fac(R,k)*VAVzz(:,:,k,f)
        end do
     end do
  end do furrier_transformer
  call ptime('FT to VAV(R, w)')


!!! Output
  open (unit_vr, file=fn_vr, status='replace')
  write(unit_vr, fmt_vr_head) trim(case%s)
  write(unit_vr, *) Nwann
  write(unit_vr, *) NR
  write(unit_vr, fmt_rweights) rweights
  do R=1,NR
     do v=1,Nwann
        do w=1,Nwann
           write(unit_vr, fmt_vr) &
                rpt_frac(R,:), v,w, VRx(v,w,R), VRy(v,w,R), VRz(v,w,R)
        end do
     end do
  end do
  close(unit_vr)

  if (inw%mixed) then
     call maybin_open (unit_vvr, file=fn_vvr, status='replace')
     call maybin_write(unit_vvr, trim(case%s),    fmt=fmt_vvr_head)
     call maybin_write(unit_vvr, (/ Nwann, Nf /), fmt='(2(1X, I0))')
     call maybin_write(unit_vvr, (/ NR /),        fmt='(I0)')
     call maybin_write(unit_vvr, rweights,        fmt=fmt_rweights)
     do f=fmin,fmax
        do R=1,NR
           do v=1,Nwann
              do w=1,Nwann
                 if (maybin_get()) then
                    write(unit_vvr) &
                         VVRxx(v,w,R,f), VVRxy(v,w,R,f), VVRxz(v,w,R,f), &
                         VVRyy(v,w,R,f), VVRyz(v,w,R,f), VVRzz(v,w,R,f)
                 else
                    write(unit_vvr, fmt_vvr) inw%w(f), rpt_frac(R,:), v,w, &
                         VVRxx(v,w,R,f), VVRxy(v,w,R,f), VVRxz(v,w,R,f), &
                         VVRyy(v,w,R,f), VVRyz(v,w,R,f), VVRzz(v,w,R,f)
                 end if
              end do
           end do
        end do
     end do
     close(unit_vvr)
  end if
  call ptime('writing V(R)')

contains
  pure elemental real(DPk) function spectral(omega, epsilon, delterl)
    use const, only: PI

    real(DPk), intent(in) :: omega, epsilon, delterl

    spectral = delterl/PI / ((omega-epsilon)**2 + delterl**2)
  end function spectral

  function hemul(he, ge)
    use util, only: string

    complex(DPk), intent(in), dimension(:,:) :: he
    complex(DPk), intent(in), dimension(:,:) :: ge

    complex(DPk), dimension(size(he,1), size(ge,2)) :: hemul

    integer :: M, N
    complex(DPk), parameter :: one=1, zero=0

    M = size(he,1)
    N = size(ge,2)

    if (any((/ size(he,2), size(ge,1) /) /= M)) &
         call croak('hemul: bad matrix shapes '//trim(string(shape(he)))// &
         ', '//trim(string(shape(ge))))

    call zhemm('L','U', M,N, one, he,M, ge,M, zero, hemul,M)
  end function hemul
end program compute_vr
