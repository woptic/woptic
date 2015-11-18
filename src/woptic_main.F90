!!! woptic/src/woptic_main.F90
!!!
!!!    Main program to compute the optical conductivity via adaptive
!!!    kmesh refinement
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2013-2015 Elias Assmann
!!!

program woptic_main
  use const,      only: DPk, PI, ERROR_UNIT
  use structmod,  only: struct_t, struct_read
  use kpoints,    only: get_kmesh_klist, count_kmesh_klist
  use util,       only: string, ptime, ptimer, ptick, ptock, ptot, inverse3x3
#define str(x) trim(string(x))
  use clio,       only: argstr, fetcharg, croak, carp
  use maybebin,   only: maybin_open, maybin_read
  use Wannier90,  only: chk_t, chk_read
  use woptic,     only: inwop_t, inwop_read, get_mommat, WOPTIC_VERSION, &
       MODE_Peierls, MODE_interp, MODE_optic, MODE_Bloch, MODE_LDA
  use woptic_io,  only: print1or3, set_casename, &
       & unit_intrahop, unit_outwop, unit_struct, unit_wfrot, unit_tet,     &
       & unit_doscontr, unit_mommat, unit_fklist, unit_klist, unit_map,     &
       &                unit_fermi,  unit_ham,   unit_K1,      &
       & unit_energy,   unit_contr,  unit_ftet,   unit_wdos,  unit_optcond, &
       & unit_optorb,   unit_selfE,  unit_vk,     unit_vvk,   unit_chk,     &
       &  suf_intrahop,  suf_outwop,  suf_inwop,   suf_wfrot,  suf_tet,     &
       &  suf_doscontr,  suf_mommat,  suf_fklist,  suf_klist,  suf_map,     &
       &  suf_struct,    suf_fermi,   suf_ham,    suf_K1,      &
       &  suf_energy,    suf_contr,   suf_ftet,    suf_wdos,   suf_optcond, &
       &  suf_optorb,    suf_selfE,   suf_vk,      suf_vvk,    suf_chk,     &
       &   fn_intrahop,   fn_outwop,   fn_inwop,    fn_wfrot,   fn_tet,     &
       &   fn_doscontr,   fn_mommat,   fn_fklist,   fn_klist,   fn_map,     &
       &   fn_struct,     fn_fermi,    fn_ham,     fn_K1,      &
       &   fn_energy,     fn_contr,    fn_ftet,     fn_wdos,    fn_optcond, &
       &   fn_optorb,     fn_selfE,    fn_vk,       fn_vvk,     fn_chk

  implicit none

  real(DPk), parameter :: KPT_TOL = 1e-10_DPk

!!! Naming conventions for the observables optcond, DOS, DCcond,
!!! Seebeck, K1:
!!!
!!!    The variable “O” has the k-resolved contributions to quantity
!!!    O; “O_tot” is k-summed; “_orb” means orbital-resolved
!!!    contributions (these are only computed if requested in
!!!    case.inwop).
!!!
!!!    The following exist:
!!!        optcond_tot, optcond, optcond_tot_orb, optcond_orb
!!!        DOS_tot,              DOS_tot_orb,     DOS_orb
!!!        DCcond_tot,  DCcond,  DCcond_tot_orb,  DCcond_orb
!!!        Seebeck_tot,          Seebeck_tot_orb
!!!        K1_tot,      K1,      K1_tot_orb,      K1_orb
!!!
!!!    In addition, DOS and optcond depend on external frequency E
!!!    (DOS has the full -maxw:maxw range, optcond only 1:maxw-1); and
!!!    optcond, DCcond, Seebeck, and K1 depend on dimensional indices
!!!    i,j.
!!!
!!!    Indices are odered as: i, j, bands, E, k.

  real   (DPk), allocatable, dimension(:)         :: fw, DOS_tot
  real   (DPk), allocatable, dimension(:,:)       :: DOS_tot_orb
  real   (DPk), allocatable, dimension(:,:,:)     :: &
       DCCond, K1, DCcond_tot_orb, K1_tot_orb,       &
       DOS_orb, optcond_tot, Seebeck_tot_orb
  real   (DPk), allocatable, dimension(:,:,:,:)   :: &
       DCcond_orb, K1_orb, optcond, optcond_tot_orb
  real   (DPk), allocatable, dimension(:,:,:,:,:) :: optcond_orb
  complex(DPk), allocatable, dimension(:,:)       :: SE, G
  complex(DPk), allocatable, dimension(:,:,:)     :: VA1, VA2, VAV, Awk

  ! mapping of dimensional indices i,j to single index
  integer,      parameter :: ij(6,2) = &
       reshape((/ 1,1,1, 2,2,3, 1,2,3, 2,3,3 /), shape(ij))
  integer,      dimension(4,8)        :: reftet
  real   (DPk), dimension(size(ij,1)) :: sumruledrude, sumrule
  real   (DPk), dimension(10)         :: wd
  real   (DPk), dimension(3,3)        :: DCcond_tot, K1_tot, Seebeck_tot
  real   (DPk), dimension(3,3)        :: S, St, D, K, O, xdc, xoc, xk1, &
       BR1, BR1inv

  integer      :: NWF            ! # inner-window (Wannier) states
  integer      :: WFmin, WFmax   ! inner window
  integer      :: NKS            ! # outer-window states
  integer      :: KSmin, KSmax   ! outer window
  integer      :: Nw             ! # internal frequencies
  integer      :: NE             ! # external frequencies
  integer      :: wmin, wmax     ! (maximal) w-window
  integer      :: Nk             ! # (new) k-points
  integer      :: Nt             ! # tetrahedra
  integer      :: Nsym           ! # symmetry operations
  
  integer      :: Nkold, Nk_echo
  integer      :: i,j, jj, jda, jw, jE, jwE
  integer      :: jt, T, jb,jb1,jb2, jd, jk,jk1,jk2, js, jarg, ii
  real   (DPk) :: mxr, sumtetraweights, df, convfac, x,xd,xk !, y, z
  complex(DPk) :: derivative(3)
  logical      :: have_file, band=.false. ! command line argument, --band option
  logical      :: vvk_binary              ! vvk files in binary (mixed trans.)

  character(*), parameter :: &
       fmt_subtime  = "('   of which ',A,T33, &
                      & '(sec):',F8.3,' wall;',F9.3,' CPU')", &
       fmt_K1DC     = "(I8,6E20.12)", &
       fmt_contr    = "(   6E20.12)", &
       fmt_doscontr = "( 100E20.12)" ! FIXME: this will fail for > 100
                                     ! bands :-)

  type(inwop_t)            :: inw
  type(struct_t)           :: stru
  type(chk_t), allocatable :: chk
  type(argstr)             :: arg, file
  type(ptimer)             :: &
       specfct_timer = ptimer(fmt=fmt_subtime), &
       optcond_timer = ptimer(fmt=fmt_subtime)

  integer,      allocatable, dimension(:)       :: inonint
  integer,      allocatable, dimension(:,:)     :: tetra
  real   (DPk), allocatable, dimension(:)       :: wtetra
  real   (DPk), allocatable, dimension(:,:)     :: bands, bands_full, kpts
  real   (DPk), allocatable, dimension(:,:,:)   :: symop
  complex(DPk), allocatable, dimension(:,:)     :: rotmat
#ifdef WANNIER_RANDOM_PHASE
  complex(DPk), allocatable, dimension(:,:)     :: phase
  real   (DPk), allocatable, dimension(:,:)     :: phi
#endif
  complex(DPk), allocatable, dimension(:,:,:)   :: u_matrix, Hk, Hkfull
  complex(DPk), allocatable, dimension(:,:,:,:) :: Hkd, Hkd_full

  real   (DPk), parameter :: echarge = 1.602176487e-19_DPk   ! [A·s]
  real   (DPk), parameter :: hbar    = 1.05457148e-34_DPk    ! [J·s]
  ! [hbar^3/emass^2*(m/ang)^5]
  real   (DPk), parameter :: factor1 = 1.413351709413265e+08_DPk
  character(*), parameter :: suf_new = '_new'
  character(*), parameter :: usage = '(                                       &
&"woptic_main: compute optical conductivity for a given k-mesh"              /&
&/"USAGE",                                                                    &
&T10,"woptic_main [--band] CASE",                                            /&
&/"FILES",                                                                    &
&T10,"(prefixed by CASE)",                                                   /&
&"*  '//suf_inwop    //'",T22,"woptic main input file"                       /&
&"*  '//suf_struct   //'",T22,"Wien2k master input file"                     /&
&"*  '//suf_klist    //'",T22,"symmetrized k-points"                         /&
&"*  '//suf_tet      //'",T22,"symmetrized tetrahedra"                       /&
&"   '//suf_energy   //'",T22,"energies from lapw1"                          /&
&"   '//suf_fermi    //'",T22,"Fermi energy"                                 /&
&"   '//suf_mommat   //'",T22,"matrix elements from optic"                  /&
&"   '//suf_chk      //'",T22,"Wannier90 checkpoint file"                    /&
&"   '//suf_vk(1)//','//suf_vk(2)//','//suf_vk(3)//' ",                       &
&                         T22,"Wannier-interpolated matrix elements"         /&
&"   '//suf_vvk(1)   //'",T22,"Wannier-interpolated mixed matrix elements"   /&
&"   '//suf_ham      //'",T22,"Wannier Hamiltonian H(R)"                     /&
&"   '//suf_selfe    //'",T22,"self-energy Σ(ω) (e.g. from DMFT)"            /&
&"   '//suf_wfrot    //'",T22,"Wannier function rotation matrix"             /&
&"P  '//suf_fklist   //'",T22,"unsymmetrized k-points"                       /&
&"P  '//suf_ftet     //'",T22,"unsymmetrized tetrahedra"                     /&
&"P  '//suf_map      //'",T22,"mapping of klist_full to klist"               /&
&"P  '//suf_intrahop //'",T22,"distance matrix for intra-u.c. hopping "                    /&
&"*U '//suf_contr    //'",T22,"function values for estimator"                /&
&"*U '//suf_K1       //'",T22,""                                             /&
&"*U '//suf_doscontr //'",T22,"stored DOS"                                   /&
&"*W '//suf_outwop   //'",T22,"diagnostic output"                            /&
&"*W '//suf_optcond  //'",T22,"optical conductivity"                         /&
&" W '//suf_optorb(1)//'",T22,"orbitally resolved optical conductivity"      /&
&"*W '//suf_wdos     //'",T22,"(joint) density of states"                    /&
&/"Files marked ‘W’ are written, and ‘U’, read and updated.  The updated"    /&
&"file ‘F’ is written to ‘F'//suf_new//'’.  Other files are read only.  "    /&
&"Precisely which of these files are used depends on OPTIONS and settings in"/&
&"‘CASE'//suf_inwop//'’ (*: always; P: Peierls mode)."                       /&
&/"OPTIONS",                                                                  &
&T10,"--help, -h"                                                            /&
&T10,"--version, -v"                                                         /&
&T10,"--band", T22, "use BZ path"                                             &
&)'

!!!------------- Argument handling              -----------------------------
  jarg=1
  arguments: do while (jarg <= command_argument_count())
     call fetcharg(jarg, arg)

     select case (arg%s)
     case ('--')
        call fetcharg(command_argument_count(), file)
        have_file=.true.
        exit arguments
     case ('-h','-H','-help','--help')
        print usage
        call exit(0)
     case ('-v', '-V', '-version', '--version')
        print '("woptic_main ", A)', WOPTIC_VERSION
        call exit(0)
     case ('-band', '--band')
        band = .true.
     case default
        if (arg%s(1:1) == '-') &
             call croak('unknown option ' // trim(arg%s))
        call fetcharg(jarg, file)
        have_file=.true.
     end select

     jarg = jarg+1
  end do arguments

  if (.not. have_file) &
       call croak('CASE argument must be given')

  call set_casename(file, band)

!!!------------- Open files for reading         -----------------------------
  call inwop_read(fn_inwop, inw)
                        open(unit_struct,   FILE=fn_struct,   STATUS='old')
                        open(unit_tet,      FILE=fn_tet,      STATUS='old')
                        open(unit_contr,    FILE=fn_contr,    STATUS='old')
                        open(unit_K1,       FILE=fn_K1,       STATUS='old')
                        open(unit_klist,    FILE=fn_klist,    STATUS='old')
  if (inw%DOS)          open(unit_doscontr, FILE=fn_doscontr, STATUS='old')
  if (inw%read_energy)  open(unit_energy,   FILE=fn_energy,   STATUS='old')
  if (inw%read_energy)  open(unit_fermi,    FILE=fn_fermi,    STATUS='old')
  if (inw%read_mommat)  open(unit_mommat,   FILE=fn_mommat,   STATUS='old')
  if (inw%read_ham)     open(unit_ham,      FILE=fn_ham,      STATUS='old')
  if (inw%selfE)        open(unit_selfE,    FILE=fn_selfE,    STATUS='old')
  if (inw%wfrot)        open(unit_wfrot,    FILE=fn_wfrot,    STATUS='old')
  if (inw%orig_umatrix) open(unit_chk,      FILE=fn_chk,      STATUS='old', &
       &                                    FORM='unformatted')
  if (inw%intrahop)     open(unit_intrahop, FILE=fn_intrahop, STATUS='old')
  if (inw%Peierls) then
                        open(unit_fklist,   FILE=fn_fklist,   STATUS='old')
                        open(unit_ftet,     FILE=fn_ftet,     STATUS='old')
                        open(unit_map,      FILE=fn_map,      STATUS='old')
  end if
  if (inw%read_vk) then
     do i=1,size(unit_vk)
                        open(unit_vk(i),    FILE=fn_vk(i),    STATUS='old')
     end do
  end if
  if (inw%read_vvk) then
     call        maybin_open(unit_vvk(1),   FILE=fn_vvk(1),   STATUS='old', &
          &                  GET=vvk_binary)
     do i=2,size(unit_vvk)
        call     maybin_open(unit_vvk(i),   FILE=fn_vvk(i),   STATUS='old', &
             &               SET=vvk_binary)
     end do
  end if


!!!------------- Preparations: read, allocate … -----------------------------
  call struct_read(unit_struct, stru); close(unit_struct)

  ! some shortcuts for oft-used limits
  KSmin = inw%bmin_w2k;    KSmax = inw%bmax_w2k;    NKS = KSmax - KSmin + 1
  WFmin = inw%bmin;        WFmax = inw%bmax;        NWF = WFmax - WFmin + 1
  wmin  = lbound(inw%w,1); wmax  = ubound(inw%w,1); Nw  =  wmax -  wmin + 1
  NE    = size  (inw%E);

  ! count / read k-points as needed
  get_Nk: if (inw%orig_umatrix) then
     allocate(chk)
     call chk_read(unit_chk, chk); close(unit_chk)
     call get_kmesh_klist(unit_klist, kpts, stru)
     Nk = size(kpts, 1)
  else
     Nk = count_kmesh_klist(unit_klist)
  end if get_Nk
  close(unit_klist)

  !allocation of some arrays
  Nsym = size(stru%rsym,3)

  allocate(inonint(WFmin-KSmin + KSmax-WFmax))
  allocate( Hk    (WFmin:WFmax, WFmin:WFmax,    Nk), &
       &    Hkd   (KSmin:KSmax, KSmin:KSmax, 3, Nk), &
       &    bands (KSmin:KSmax,                 Nk), &
#ifdef WANNIER_RANDOM_PHASE
       &    phi   (WFmin:WFmax,                 Nk), &
       &    phase (WFmin:WFmax,                 Nk), &
#endif
       &    symop (3, 3, Nsym)          &
  )
  inonint = (/ (i, i=KSmin,WFmin-1), (i, i=WFmax+1,KSmax) /)

#ifdef WANNIER_RANDOM_PHASE
  call init_random_seed()
  
  call random_number(phi)
  phase = exp((0,2)*PI*phi)
#endif

  BR1inv = stru%lat2car/2/PI
  call inverse3x3(BR1inv, BR1)
  do i=1, size(symop,3)
     symop(:,:,i) = matmul(matmul(BR1, stru%rsym(:,:,i)), &
          &                BR1inv)
  end do

  open (unit_outwop, FILE=fn_outwop, STATUS='replace')
  call ptime(UNIT=unit_outwop)
  call ptime(UNIT=unit_outwop, TIMER=specfct_timer)
  call ptime(UNIT=unit_outwop, TIMER=optcond_timer)
  write(unit_outwop,'(A)') &
       "----- OPTICAL CONDUCTIVITY WITH WANNIER FUNCTIONS -----"
  write(unit_outwop, '(A)') &
       "      This is woptic_main " // woptic_version
  write(unit_outwop,*)
  write(unit_outwop, '(I0," k points, ",I0," Wannier functions ")') Nk, NWF
  write(unit_outwop, '(I0," external frequencies up to ",F7.2)') NE, inw%E(NE)
  write(unit_outwop, '(I0," internal frequencies from ",SP,F7.3," to ",F7.3)')&
       Nw, inw%w(wmin), inw%w(wmax)
  write(unit_outwop, '("number of elements: ", I0)') stru%nneq
  write(unit_outwop, '("band windows: [",&
       &I0, "..", I0, "] ⊆ [", I0, "..", I0, "]")') WFmin,WFmax, KSmin,KSmax
  write(unit_outwop, "('unit cell volume: ',G10.6E1,' bohr³')") stru%vol

  call read_tetra(unit_tet, tetra, wtetra); close(unit_tet)

  write(unit_outwop,*)

!!!------------- Conditionally read files according to options --------------
  if (inw%wfrot) then
     allocate(rotmat(WFmin:WFmax, WFmin:WFmax))
     do jb=WFmin,WFmax
        read(unit_wfrot,*) rotmat(jb,:)
     end do
     close(unit_wfrot)
  end if

  if (inw%read_energy) then
     write(unit_outwop,'(A)') '>>> read energies from '//trim(fn_energy)
     if (inw%orig_umatrix) then
        allocate(bands_full(KSmin:KSmax, chk%num_kpts))
        call read_energy(unit_energy, bands_full, KSmin, KSmax, chk%num_kpts)
     else
        call read_energy(unit_energy, bands, KSmin, KSmax, Nk)
     end if
     write(unit_outwop,*)
     close(unit_energy)
  end if

  if (inw%read_ham) then
     call read_ham(unit_ham, Hk, WFmin, WFmax, Nk)
     close(unit_ham)
  end if
  
  if (inw%read_mommat) then
     write(unit_outwop,'(A)') &
          '>>> read Wien2k matrix elements from '//trim(fn_mommat)
     write(unit_outwop,*)

     if (inw%orig_umatrix) then
        allocate(Hkd_full(KSmin:KSmax, KSmin:KSmax, 3, chk%num_kpts))
        call get_mommat(unit_mommat,                                  &
             Hkd_full(:,:,1,:), Hkd_full(:,:,2,:), Hkd_full(:,:,3,:), &
             chk%num_kpts, KSmin, KSmax)
     else
        call get_mommat(unit_mommat, Hkd(:,:,1,:), Hkd(:,:,2,:), Hkd(:,:,3,:),&
             &          Nk, KSmin, KSmax)
     end if
  end if

  get_U: if (inw%orig_umatrix) then
     allocate(u_matrix(WFmin:WFmax, WFmin:WFmax, Nk))

     ! Search kpts's jk1-th k-point in chk%kpt_latt.  We assume the
     ! two are sorted in the same way, and that the k-point list for
     ! ‘mommat2’ and ‘energy’ is the same.
     jk1 = 1; jk2 = 1
     do
        if (all(abs(kpts(jk1,:) - chk%kpt_latt(:,jk2)) < KPT_TOL)) then
           u_matrix(:,:,jk1) = conjg(transpose(chk%u_matrix(:,:,jk2)))
           Hkd  (:,:,:, jk1) = Hkd_full (:,:,:, jk2)
           if (inw%read_energy) &
                bands(:,  jk1) = bands_full  (:,  jk2)

           do jb1 = KSmin, KSmax
              Hkd(jb1, jb1, :, jk1)    =       Hkd_full(jb1, jb1, :, jk2)
              do jb2 = jb1+1, KSmax
                 Hkd(jb1, jb2, :, jk1) =       Hkd_full(jb1, jb2, :, jk2)
                 Hkd(jb2, jb1, :, jk1) = conjg(Hkd_full(jb1, jb2, :, jk2))
              end do
           end do

           jk1 = jk1 + 1
           jk1 = jk1 + 1
           if (jk1 > Nk) exit
        end if

        jk2 = jk2 + 1
        if (jk2 > chk%num_kpts) call croak('k-point from ‘klist’ not found'//&
        'in ‘chk’: '//str(kpts(jk1,:)))
     end do

     deallocate(kpts, chk, Hkd_full)
     if (inw%read_energy) deallocate(bands_full)
  elseif (inw%need_umatrix) then
     allocate(u_matrix(WFmin:WFmax, WFmin:WFmax, Nk))
     do jk=1,Nk
        call get_umatrix(Hk(:,:,jk), NWF, u_matrix(:,:,jk))
     end do
  end if get_U

  complete_Hkd: if (inw%read_mommat) then
     forall (jb1=KSmin:KSmax, jb2=KSmin:KSmax, jb2>jb1)
        Hkd(jb2, jb1, :,:) = conjg(Hkd(jb1, jb2, :,:))
     end forall
  end if complete_Hkd


!!! Construct / read self-energy
  allocate(SE(WFmin:WFmax, wmin:wmax))
  SE = -(0,1)*inw%delterl
  if (inw%selfE) then
     write(unit_outwop,'(A)') &
          '>>> read self-energy from ‘'//trim(fn_selfE)//'’ …'
     write(unit_outwop,*)

     SE(inw%iself,:) = read_SE(unit_selfE, size(inw%iself), wmin, wmax)
     close(unit_selfE)

     write(unit_outwop,'(A)') &
          '    self-energy will be applied to orbitals '//str(inw%iself)
  end if

  if (inw%shift) then
     if(inw%read_vvk) call croak( &
          "FIXME: scissors shift with VAV interpolation is unimplemented")
     if (any(inw%ishift < WFmin .or. inw%ishift > WFmax)) call croak( &
          "FIXME: scissors shift for non-Wannier bands is unimplemented")

     write(unit_outwop,'(4X,A,F0.3,A)') 'applying scissors shift ', &
          inw%Eshift,' to orbitals '//str(inw%ishift)

     SE(inw%ishift(:), :) = SE(inw%ishift(:), :) - inw%Eshift
  end if

  if (inw%chempot /= 0 .and. inw%selfE) then
     write(unit_outwop,'(4X,A,F0.3,A)') 'applying chemical potential ', &
          inw%chempot,' to orbitals '//str(inw%iself)

     SE(inw%iself (:), :) = SE(inw%iself (:), :) - inw%chempot
  end if


!!!------------- Check matelmode and do appropriate things to V(k) ----------
  matelmode: select case (inw%matelmode)
  case (MODE_Peierls)
     write(unit_outwop,'(A)') &
          '>>> use matrix elements from Peierls approximation'
     write(unit_outwop,*)

     call get_hkderiv()

  case (MODE_interp)
     write(unit_outwop,'(A)') &
          '>>> use Wannier-interpolated matrix elements V(k)'
     write(unit_outwop,*)

     ! read_vk() will replace W-W part of Hkd
     call read_vk(unit_vk, Hkd(WFmin:WFmax, WFmin:WFmax,:,:))
     do jd=1,size(unit_vk)
        close(unit_vk(jd))
     end do

     if (inw%read_vvk) then
        write(unit_outwop,'(A)') &
             '>>> read Wannier-interpolated mixed matrix elements VAV(k)'
        write(unit_outwop,*)
        call read_vvk_heads(unit_vvk)
     end if

  case (MODE_optic)
     write(unit_outwop,'(A)') '>>> use rotated matrix elements'
     write(unit_outwop,'(A)') '    with Wannier90 Hamiltonian'
     write(unit_outwop,*)

     do jk=1,Nk
        do jd=1,3
           ! For the case including mixed transitions, we have to
           ! apply the U(k) also to (one side of) the mixed matrix
           ! elements.
           !
           ! This order of factors is different from Wannier90, but
           ! gives the same results as the old code
           Hkd(:, WFmin:WFmax, jd,jk) = &
                matmul(Hkd(:, WFmin:WFmax, jd,jk), &
                &      transpose(conjg(u_matrix(:,:,jk))))
           Hkd(WFmin:WFmax, :, jd,jk) = &
                matmul(u_matrix(:,:,jk), &
                &      Hkd(WFmin:WFmax, :, jd,jk))
        end do
     end do

  case (MODE_Bloch)          ! use diagonal Hamiltonian
     write(unit_outwop,'(A)') '>>> use original Wien2k matrix elements'
     write(unit_outwop,'(A)') '    with diagonal Hamiltonian'
     write(unit_outwop,*)
     write(unit_outwop,'(A)') '    “I''m so mean I make medicine sick!”'
     write(unit_outwop,*)

     !rotate Hamiltonian back to diagonal basis
     do jk=1,Nk
        Hk(:,:,jk) = matmul(transpose(conjg(u_matrix(:,:,jk))), &
             &              matmul(Hk(:,:,jk), u_matrix(:,:,jk)))
     end do

  case (MODE_LDA)            ! wien2k only mode
     write(unit_outwop,'(A)') '>>> use original Wien2k matrix elements'
     write(unit_outwop,'(A)') '    with diagonal Hamiltonian from Wien2k'
     write(unit_outwop,*)
     do jb=WFmin,WFmax
        Hk(jb,jb,:) = bands(jb,:)
     end do

  case default
     call croak("unknown matrix element mode: "// &
          str(inw%matelmode))
  end select matelmode

#ifdef WANNIER_RANDOM_PHASE
  write(unit_outwop,'(A)') '    with random phases manually added'
  write(unit_outwop,*)

  forall (jb=WFmin:WFmax, jk=1:Nk)
     Hkd(jb,:,1,jk) = Hkd(jb,:,1,jk) *       phase(jb,jk)
     Hkd(:,jb,1,jk) = Hkd(:,jb,1,jk) * conjg(phase(jb,jk))
  end forall
#endif

  if (allocated(u_matrix)) deallocate(u_matrix)

!!!------------- Let's go                       -----------------------------
  write(unit_outwop,'(A)') '>>> compute optical conductivity …'
  write(unit_outwop,*)

  if (inw%Peierls) then
     !Expl:    2*pi*elem.charge^2/hbar/unitcellvol[bohr^3] * ang/bohr*meter/ang*cm/meter (*sigma[bohr^2] from Ham.)
     convfac = 2*PI*echarge**2/hbar/stru%vol/0.529177_DPk*10E8_DPk!conv. to(Siemens/cm)
  else!wien2k matrix elements
     !Expl: echarge*hbar^3[Js]*(J/eV)^2/emass^2/bohr^5*(bohr^5/0.52^5 ang^5)*(10^50 ang^5/m^5)*(m/10^2 cm)
     !strange wien2k factor coming from joint.f: 64/pi???
     !spin is already included in matrix elements
     convfac = 64*factor1/stru%vol/(0.529177**5)/10E2_DPk!conv. to(Siemens/cm)
  end if

  write(unit_outwop,*)"mode: ",inw%matelname
  write(unit_outwop,*)"number of orbitals contributing to optics: ",NKS

  if (inw%joint) then
     if (inw%mixed) call croak("JOINT unimplemented for mixed transitions")
     ! FIXME: this is (probably) broken by VAV interpolation.  When
     ! reimplementing, it could be optimized based on matel1 == matel2
     ! always.
     write(unit_outwop,*)'JOINT modus: compute joint density of states'
     Hkd     = 1
     convfac = 1
  end if
  write(unit_outwop,*)"convfac: ",convfac

  allocate(Awk        (     WFmin:WFmax,WFmin:WFmax, wmin:wmax    ), &
       &   G          (     WFmin:WFmax,WFmin:WFmax               ), &
       &   VA1        (     WFmin:WFmax,WFmin:WFmax, wmin:wmax    ), &
       &   VA2        (     WFmin:WFmax,WFmin:WFmax, wmin:wmax    ), &
       &   optcond    (3,3,                          1:NE,      Nk), &
       &   optcond_tot(3,3,                          1:NE         ), &
       &   DCcond     (3,3,                                     Nk), &
       &   K1         (3,3,                                     Nk), &
       &   DOS_orb    (     KSmin:KSmax,             wmin:wmax, Nk), &
       &   DOS_tot_orb(     KSmin:KSmax,             wmin:wmax    ), &
       &   DOS_tot    (                              wmin:wmax    ), &
       &   fw         (                              wmin:wmax    ))
  if (inw%mixed) allocate( &
       &   VAV        (     WFmin:WFmax,WFmin:WFmax, wmin:wmax    ))

  DOS_orb = 0; optcond_tot = 0
  optcond = 0; DOS_tot_orb = 0

  fw = 1 / (exp(inw%beta * inw%w) + 1)

  write(unit_outwop,*)'parameters: [Emax dE delterl beta]'
  write(unit_outwop,"(2x,4F13.6)")inw%Emax,inw%dE,inw%delterl,inw%beta

  if (inw%orbresolv) then
     allocate(optcond_orb    (3,3, KSmin:KSmax, 1:NE,  Nk), &
          &   optcond_tot_orb(3,3, KSmin:KSmax, 1:NE     ), &
          &   DCcond_orb     (3,3, KSmin:KSmax,        Nk), &
          &   K1_orb         (3,3, KSmin:KSmax,        Nk), &
          &   DCcond_tot_orb (3,3, KSmin:KSmax           ), &
          &   Seebeck_tot_orb(3,3, KSmin:KSmax           ), &
          &   K1_tot_orb     (3,3, KSmin:KSmax           ))
     optcond_orb = 0;     DCcond_orb=0;       K1_orb = 0
     optcond_tot_orb = 0; DCcond_tot_orb = 0; K1_tot_orb = 0
  end if


!!!------------- Read old k-resolved contributions --------------------------
  Nkold = 0; DCcond = 0; K1 = 0
  if (.not. band) then
     read(unit_contr,*) Nkold
  end if
  do jk=1,Nkold
     read(unit_contr,fmt_K1DC) &
          i, ( DCcond(ij(jj,1), ij(jj,2), jk), jj=1,size(ij,1) )

     read(unit_K1,fmt_K1DC) &
          i, ( K1(ij(jj,1),ij(jj,2), jk), jj=1,size(ij,1) )

     do jE=1,NE
        read(unit_contr,fmt_contr) &
             ( optcond(ij(jj,1),ij(jj,2), jE, jk), jj=1,size(ij,1) )
     end do
     if (inw%DOS) then
        do jw=wmin,wmax
           read(unit_doscontr,fmt_doscontr) DOS_orb(:, jw, jk)
        end do
     end if
  end do
  close(unit_contr)
  close(unit_doscontr)
  close(unit_K1)

  call ptime("optcalc setup")


!!!------------- Compute new k-resolved contributions -----------------------
  write(unit_outwop,"(/, '+------------------------------+')")
  Nk_echo = min(100, (Nk-Nkold)/10+1)
  !$OMP parallel do
  k_points: do jk=Nkold+1,Nk
     call ptick(specfct_timer)
     compute_A_wk: do jw = wmin,wmax
        Awk(:,:, jw) = specmat(jw, jk)
        do jb=WFmin,WFmax
           DOS_orb(jb,jw,jk) = DOS_orb(jb,jw,jk) + real(Awk(jb,jb,jw))
        end do
        DOS_orb(inonint, jw, jk) = DOS_orb(inonint, jw, jk) &
             + specfun(inonint, jw, jk)
     end do compute_A_wk
     call ptock(specfct_timer)

     tensoridx: do jj=1,size(ij,1)
        i = ij(jj,1)
        j = ij(jj,2)

        compute_VA_wk: do jw = wmin,wmax
           VA1(:,:,jw) = matmul(Hkd(WFmin:WFmax,WFmin:WFmax,i,jk), Awk(:,:,jw))
           VA2(:,:,jw) = matmul(Hkd(WFmin:WFmax,WFmin:WFmax,j,jk), Awk(:,:,jw))
        end do compute_VA_wk

        get_VAV_wk: if (inw%read_vvk) then
           call read_vvk_kpt(unit_vvk(jj), VAV)
#ifdef WANNIER_RANDOM_PHASE
           forall (jb=WFmin:WFmax)
              VAV(jb,:,:) = VAV(jb,:,:) *       phase(jb,jk)
              VAV(:,jb,:) = VAV(:,jb,:) * conjg(phase(jb,jk))
           end forall
#endif
        elseif         (inw%mixed) then
           compute_VAV_wk: do jw = wmin,wmax
              VAV(:,:,jw) = 0
              do jb1=WFmin,WFmax
                 do jb2=WFmin,WFmax
                    VAV(jb1,jb2,jw) = VAV(jb1,jb2,jw) +      &
                         sum(Hkd(jb1, inonint,      i, jk) * &
                         &   specfun (inonint, jw,     jk) * &
                         &   Hkd     (inonint, jb2, j, jk))
                 end do
              end do
           end do compute_VAV_wk
        end if get_VAV_wk

        call ptick(optcond_timer)
        compute_optcond: if(inw%optcond) then
           xd=0; xk=0
           dc_cond: do jw = -inw%Nwx, +inw%Nwx
              df = -inw%beta * cosh(inw%beta*inw%w(jw)/2)**(-2)/4
              x=0

              bloch_dc: do ii=1, size(inonint)
                 jb = inonint(ii)
                 x = x + sum( &
                      &      specfun(jb,      jw,    jk) &
                      &*     specfun(inonint, jw,    jk) &
                      &*real(Hkd(jb, inonint,     i, jk) &
                      &*     Hkd(    inonint, jb, j, jk)))
              end do bloch_dc

              wannier_dc: do jb=WFmin,WFmax
                 x = x + real(sum(VA1(jb, :, jw) * VA2(:, jb, jw)))
              end do wannier_dc

              mixed_dc: if (inw%mixed) then
                 do jb = WFmin,WFmax
                    x = x + 2*real(sum(VAV(jb,:,jw) * Awk(:,jb, jw)))
                 end do
              end if mixed_dc

              xd = xd + inw%dw*df          *x
              xk = xk + inw%dw*df*inw%w(jw)*x
           end do dc_cond
           DCcond(i,j, jk) = DCcond(i,j, jk) - xd
           K1    (i,j, jk) = K1    (i,j, jk) + xk
           
           ext_freq: do jE = 1, NE
              jwE = inw%wint_dens*jE
              int_freq: do jw = -jwE-inw%Nwx, +inw%Nwx
                 x=0
                 bloch: do ii=1, size(inonint)
                    jb = inonint(ii)
                    x = x + sum(&
                         & specfun(inonint, jw,          jk)  &
                         * specfun(jb,      jw+jwE,      jk) &
                         * real(Hkd(jb,      inonint, i, jk)  &
                         *      Hkd(inonint, jb,      j, jk)))
                 end do bloch

                 wannier: do jb = WFmin,WFmax
                    x = x + real(sum( &
                         VA1(jb, :, jw) * VA2(:, jb, jw+jwE)))
                 end do wannier

                 mixed: if (inw%mixed) then
                    do jb = WFmin,WFmax
                       x = x &
                            + real(sum(VAV(jb,:, jw)       &
                            &         *Awk(:,jb, jw+jwE))) &
                            + real(sum(Awk(jb,:, jw)       &
                            &         *VAV(:,jb, jw+jwE)))
                    end do
                 end if mixed

                 optcond(i,j, jE, jk) = optcond(i,j, jE, jk) + &
                      inw%dw*(fw(jw) - fw(jw+jwE))/inw%E(jE) * x
              end do int_freq
           end do ext_freq
        end if compute_optcond
        call ptock(optcond_timer)

        orbresolv: if (inw%orbresolv) then
           call croak('FIXME: orbresolv still has to be adapted to VAV-related changes!')
           orb_optcond: if(inw%optcond) then
              orb_dc: do jw= wmin, wmax
                 df = -inw%beta * cosh(inw%beta*inw%w(jw)/2)**(-2)/4

                 do jb = 1,NKS
                    x = real(sum(VA1(jb, :, jw) * VA2(:, jb, jw)))
                    DCcond_orb(i,j, jb, jk) = DCcond_orb(i,j, jb, jk) &
                         - inw%dE*df*x
                    K1_orb(i,j, jb, jk) = K1_orb(i,j, jb, jk) &
                         + inw%dE*df*inw%w(jw)*x
                 end do
              end do orb_dc

              orb_ext: do jE = 1, NE
                 orb_int: do jw = -inw%wint_dens*jE-inw%Nwx, +inw%Nwx
                    do jb = 1,NKS
                       x = real(sum( &
                            VA1(jb, :, jw) * VA2(:, jb, jw+jE)))

                       optcond_orb(i,j, jb, jE, jk) = &
                            optcond_orb(i,j, jb, jE, jk) &
                            + inw%dE*(fw(jw) - fw(jw+jE))/inw%E(jE) * x
                    end do
                 end do orb_int
              end do orb_ext
           end if orb_optcond
        end if orbresolv
     end do tensoridx !loop over entries of the cond. matrix

     if (mod(jk, Nk_echo)==0) then 
        write(unit_outwop, &
             '("| k-point",I5," of ",I0," | ",F3.0,"%",T32,"|")') &
             jk, Nk, 100._DPk*jk/Nk
        flush(unit_outwop)
     end if
  end do k_points
  !$OMP end do
  write(unit_outwop,"('+------------------------------+', /)")

  deallocate(VA1, VA2, fw)
  do i=1,size(unit_vvk)
     close(unit_vvk(i))
  end do

  call ptime("optcalc loop")
  call ptot("Green function", specfct_timer)
  call ptot("optcond",        optcond_timer)


!!!------------- Open files for writing         -----------------------------
  open(unit_optcond,   FILE=fn_optcond,                  STATUS='replace')
  open(unit_contr,     FILE=trim(fn_contr    )//suf_new, STATUS='replace')
  open(unit_K1,        FILE=trim(fn_K1       )//suf_new, STATUS='replace')
  if (inw%DOS) then
  open(unit_wdos,      FILE=fn_wdos,                     STATUS='replace')
  open(unit_doscontr,  FILE=trim(fn_doscontr )//suf_new, STATUS='replace')
  end if
  if (inw%orbresolv) then; do i=1,size(unit_optorb)
  open(unit_optorb(i), FILE=fn_optorb(i),                STATUS='replace')
  end do; end if

  write(unit_contr,*) Nk, NE, size(ij,1), convfac

  do jk=1,Nk
     write(unit_contr, fmt_K1DC) &
          jk, ( DCcond(ij(jj,1), ij(jj,2), jk), jj=1,size(ij,1) )

     write(unit_K1, fmt_K1DC) &
          jk, ( K1(ij(jj,1),ij(jj,2), jk), jj=1,size(ij,1) )

     do jE=1,NE
        write(unit_contr, fmt_contr) &
             ( optcond(ij(jj,1),ij(jj,2), jE, jk), jj=1,size(ij,1))
     end do

     if (inw%DOS) then
        do jw=wmin,wmax
           write(unit_doscontr, fmt_doscontr) DOS_orb(:, jw, jk)
        end do
     end if
  end do
  close(unit_contr)
  close(unit_doscontr)
  close(unit_K1)

  call ptime("write contr")

  if (band) then
     stop
  end if
  
  do i=1,3
     do j=i+1,3
        DCcond (j,i,:)   = DCcond (i,j,:)
        K1     (j,i,:)   = K1     (i,j,:)
        optcond(j,i,:,:) = optcond(i,j,:,:)
        if (inw%orbresolv) then
           DCcond_orb (j,i,:,:)   = DCcond_orb (i,j,:,:)
           K1_orb     (j,i,:,:)   = K1_orb     (i,j,:,:)
           optcond_orb(j,i,:,:,:) = optcond_orb(i,j,:,:,:)
        end if
     end do
  end do


!!!------------- Symmetrized tetrahedral integration ------------------------
  write(unit_outwop,*)
  write(unit_outwop,'(A)') ">>> symmetrized tetrahedral integration"
  write(unit_outwop,*)
  sumtetraweights = sum(wtetra)

  DCcond_tot  = 0; K1_tot = 0
  tetrahedra: do jt=1,Nt
     reftet(:, 1) = tetra((/ 1, 5, 6,  7 /), jt)
     reftet(:, 2) = tetra((/ 2, 5, 8,  9 /), jt)
     reftet(:, 3) = tetra((/ 3, 6, 8, 10 /), jt)
     reftet(:, 4) = tetra((/ 4, 7, 9, 10 /), jt)
     reftet(:, 5) = tetra((/ 5, 6, 7, 10 /), jt)
     reftet(:, 6) = tetra((/ 5, 6, 8, 10 /), jt)
     reftet(:, 7) = tetra((/ 5, 7, 9, 10 /), jt)
     reftet(:, 8) = tetra((/ 5, 8, 9, 10 /), jt)

!!! Timings for SVO (no orbresolv), SE in t2g + outer window
!!!
!!! ifort + mkl:
!!!    Times for optcond               (sec):   0.013 wall;    0.013
!!!    Times for optcalc integration   (sec):   5.321 wall;    5.318
!!! gfortran44 -llapack:
!!!    Times for optcond               (sec):   0.070 wall;    0.070
!!!    Times for optcalc integration   (sec):  26.481 wall;   26.473
!!! orig: 
!!!    Times for optcond               (sec):   0.163 wall;    0.164
!!!    Times for optcalc integration   (sec):  59.650 wall;   59.641
     xdc = 0; xk1 = 0
     sym1: do js=1,Nsym
        S = symop(:,:, js)
        St = transpose(S)
        do jd=1,8
           do jda=1,4
              D = DCcond(:,:, reftet(jda,jd))
              K =     K1(:,:, reftet(jda,jd))

              xdc = xdc + matmul(St, matmul(D, S))
              xk1 = xk1 + matmul(St, matmul(K, S))
           end do
        end do
     end do sym1

     DCcond_tot  = DCcond_tot + xdc * wtetra(jt)/32/Nsym
     K1_tot      = K1_tot     + xk1 * wtetra(jt)/32/Nsym
     
     freq: do jw = 1, NE
        xoc = 0;
        do jd=1,8
           do jda=1,4
              O = optcond(:,:, jw, reftet(jda,jd))
              sym2: do js=1,Nsym
                 S = symop(:,:, js)
                 St = transpose(S)

                 xoc = xoc + matmul(St, matmul(O, S))
              end do sym2
           end do
        end do
        optcond_tot(:,:,jw) = optcond_tot(:,:,jw) + xoc * wtetra(jt)/32/Nsym
     end do freq
!     if (jt==1) call ptime("optcond")

     if (inw%DOS) then
        do jd=1,8
           DOS_tot_orb(:, :) = DOS_tot_orb(:, :) + &
                sum(DOS_orb(:, :, reftet(:, jd)), 3) * wtetra(jt)/32
        end do
     end if
!     if (jt==1) call ptime("DOS")
  end do tetrahedra

  call ptime("optcond integration")

  orb: if (inw%orbresolv) then
     bands_loop: do jb=KSmin,KSmax
        xdc = 0; xk1 = 0

        sym3: do js=1,Nsym
           S = symop(:,:, js)
           St = transpose(S)

           tetrahedra2: do jt=1,Nt
              wd(1: 4) = wtetra(jt)/20/Nsym
              wd(5:10) = wtetra(jt)/ 5/Nsym

              do jd=1,10
                 jk = tetra(jt,jd)

                 D = DCcond_orb(:,:, jb, jk)
                 K =     K1_orb(:,:, jb, jk)

                 xdc = xdc + matmul(St, matmul(D, S))*wd(jd)
                 xK1 = xK1 + matmul(St, matmul(K, S))*wd(jd)
              end do
           end do tetrahedra2
        end do sym3

        DCcond_tot_orb(:,:, jb)  = DCcond_tot_orb(:,:, jb) + xdc
        K1_tot_orb    (:,:, jb)  = K1_tot_orb    (:,:, jb) + xk1
     end do bands_loop

     freq2: do jw = 1, NE
        bands2: do jb=KSmin, KSmax
           xoc = 0

           sym4: do js=1,Nsym
              S = symop(:,:, js)
              St = transpose(S)

              tetrahedra3: do jt=1,Nt
                 wd(1: 4) = wtetra(jt)/20/Nsym
                 wd(5:10) = wtetra(jt)/ 5/Nsym

                 do jd=1,10
                    O = optcond_orb(:,:, jb, jw, tetra(jt,jd))

                    xoc = xoc + matmul(St, matmul(O, S))*wd(jd)
                 end do
              end do tetrahedra3
           end do sym4

           optcond_tot_orb(:,:, jb, jw)  = optcond_tot_orb(:,:, jb, jw) + &
                xoc
        end do bands2
     end do freq2
  end if orb
     
  call ptime("orbresolv integration")


!!!------------- Post-processing and output     -----------------------------
  !compute thermopower
  !explanation of factor: conv2muV*kboltz/echarge
  Seebeck_tot = 86.173427909006634_DPk * K1_tot / DCcond_tot * inw%beta
  if (inw%orbresolv) then
     Seebeck_tot_orb = &
          86.173427909006634_DPk * K1_tot_orb / DCcond_tot_orb * inw%beta
  end if

  !compute DOS
  if (inw%DOS) DOS_tot = sum(DOS_tot_orb, 1)

  !unit conversion to Siemens/cm
  DCcond_tot  = DCcond_tot  *convfac
  optcond_tot = optcond_tot *convfac
  if (inw%orbresolv) then
     DCcond_tot_orb  = DCcond_tot_orb  *convfac
     optcond_tot_orb = optcond_tot_orb *convfac
  end if

  call ptime("optcalc S/dc")

!!! General info
  call print1or3(unit_outwop, 'dc cond.', '1/Ω cm', DCcond_tot)
  call print1or3(unit_outwop, 'thermopower', 'μV/K', Seebeck_tot)

  write(unit_outwop,*) 'Step for intergral=', inw%dw
  write(unit_outwop,*) 'sum of tetrahedral weights=', sumtetraweights
  write(unit_outwop,*) "conversion factor:", convfac

!!! Sumrules
  ii = ceiling(inw%drudesep/inw%dE)
  do jj=1,size(ij,1)
     sumruledrude(jj) = sum(optcond_tot(ij(jj,1), ij(jj,2),    1:ii))*inw%dE
     sumrule     (jj) = sum(optcond_tot(ij(jj,1), ij(jj,2), ii+1:  ))*inw%dE
  end do

  write(unit_outwop,*) "sumrules [Drude,interband]:", &
       sumruledrude(1), sumrule(1)

  do jj=1,size(ij,1)
     sumruledrude(jj) = sumruledrude(jj) + &
          DCcond_tot(ij(jj,1), ij(jj,2)) / sumtetraweights*inw%dE
  end do
  write(unit_outwop,*) "sumrule +dc [Drude]       :", sumruledrude(1)

!!! Optical conductivity
  write(unit_optcond,*)"# Optical conductivity in Ω-1 cm-1"
  write(unit_optcond,*)"# ω σ_xx σ_xy σ_xz σ_yy σ_yz σ_zz"
  write(unit_optcond,'(1E20.12,100E20.12)') &
       0._DPk, ( DCcond_tot(ij(jj,1),ij(jj,2)), jj=1,size(ij,1) )
  do jE=1,NE
     write(unit_optcond,'(1E20.12,100E20.12)') &
          inw%E(jE), ( optcond_tot(ij(jj,1), ij(jj,2), jE), jj=1,size(ij,1) )
  end do
  close(unit_optcond)

!!! Density of states
  if (inw%DOS) then
     write(unit_wdos,*)"# spectral function in eV^-1"
     write(unit_wdos,*)"# ω total orbital1 orbital2 …"
     do jw=wmin,wmax
        write(unit_wdos,'(200F13.8)') inw%w(jw), DOS_tot(jw), DOS_tot_orb(:,jw)
     end do
  end if
  close(unit_wdos)

!!! Orbitally resolved quantities
  if (inw%orbresolv) then
     do jj=1,size(ij,1)
        write(unit_optorb(jj),'(A2,100E20.12)') &
             "# ", Seebeck_tot_orb(ij(jj,1), ij(jj,2), :)

        write(unit_optorb(jj),'(1E20.12,100E20.12)') &
             0_DPk, DCcond_tot_orb(ij(jj,1), ij(jj,2), :)

        do jE=1,NE
           write(unit_optorb(jj),'(1E20.12,100E20.12)') &
                inw%E(jE), optcond_tot_orb(ij(jj,1), ij(jj,2), :, jE)
        end do
        close(unit_optorb(jj))
     end do
  end if

  call ptime("optcalc output")

contains

subroutine get_hkderiv()
  real(DPk) :: tmpr_Nb(2*NWF)
  integer   :: Nkfull, Nk2, Nb, Ntfull

  complex(DPk), allocatable, dimension(:,:,:,:) :: Hkdfull
  real   (DPk), allocatable, dimension(:,:,:)   :: distmatrix
  real   (DPk), allocatable, dimension(:,:)     :: kfull
  integer,      allocatable, dimension(:,:)     :: tetrafull, map
  integer,      allocatable, dimension(:)       :: mapcount, patches

  !we need the full kmesh only in the Peierls case
  call get_kmesh_klist(unit_fklist, kfull, stru)
  close(unit_fklist)
  Nkfull = size(kfull, 1)

  allocate( Hkdfull (KSmin:KSmax, KSmin:KSmax, 3, Nkfull), &
       &    Hkfull  (KSmin:KSmax, KSmin:KSmax, Nkfull), &
       &    mapcount(Nk),      &
       &    map     (2, Nkfull) )

  Hkdfull = 0; mapcount = 0

  read(unit_ftet,*) Ntfull
  allocate(tetrafull(Ntfull,10))
  do jt=1,Ntfull
     read(unit_ftet,*) tetrafull(jt,:),mxr
  end do
  close(unit_ftet)

  do jk=1,Nkfull
     read(unit_map,*) map(:,jk)
  end do
  close(unit_map)

  !in case of Peierls the full Hamilonian without symmetries is read in
  read(unit_ham,*)  Nk2,Nb
  if (Nkfull /= Nk2) &
       call croak("number of k-points inconsistent, &
       & check ‘.klist_full’ and ‘.hk’ files")
  if (Nb /= NWF) &
       call croak("number of bands inconsistent, &
       & check ‘.inwop’ and ‘.hk’ files")

  allocate(patches(Nkfull))
  !     read-in Hamiltonian
  do jk=1,Nkfull
     read(unit_ham,*)
     do jb = WFmin,WFmax
        read(unit_ham,*) tmpr_Nb; Hkfull(jb, WFmin:WFmax, jk) = &
             & cmplx(tmpr_Nb(1::2), tmpr_Nb(2::2), DPk)
     end do
  end do
  close(unit_ham)

  !now fill up the Hamiltonian with band energies if required
  do jk=1,Nkfull
     do jb=KSmin,WFmin-1
        if (bands(jb, map(1,jk)) == 0) &
             call carp("Warning: lower index of bands is smaller &
             & than overall wien2k index")

        Hkfull(jb,jb,jk) = bands(jb, map(1,jk))
     end do
     do jb=WFmax+1, KSmax
        if (bands(jb, map(1,jk)) == 0) &
             call carp("Warning: upper index of bands is larger &
             & than overall wien2k index")

        Hkfull(jb, jb, jk) = bands(jb, map(1,jk))
     end do
  end do
  write(unit_outwop,'(A)') '>>> finished reading Hamiltonian'
  write(unit_outwop,'(A)') '    compute group velocity via derivative of &
       & Hamiltonian …'
  write(unit_outwop,*)
  !compute the derivative on the full kmesh
  patches = 0
  tetrahedra: do jt=1,Ntfull
     do jd=1,10
        T = tetrafull(jt,jd)
        do jb1=KSmin,KSmax
           do jb2=KSmin,KSmax
              call compute_tetraderiv(tetrafull(jt,:), kfull,          &
                   &                  Hkfull(jb1,jb2,tetrafull(jt,:)), &
                   &                  jd, derivative)

              Hkdfull(jb1,jb2,:,T) = Hkdfull(jb1,jb2,:,T) + derivative(:)
           end do
        end do
        patches(T) = patches(T) + 1
     end do
  end do tetrahedra

  normalize_Hfull: do jk=1,Nk
     Hkdfull(:,:,:,jk) = Hkdfull(:,:,:,jk) / patches(jk)
  end do normalize_Hfull

  !now map to symmetrized mesh
  do jk=1,Nkfull
     if (mapcount(map(1,jk)) /= 0) cycle

     do jb1=KSmin,KSmax
        do jb2=KSmin,KSmax
           derivative(:) = Hkdfull(jb1,jb2,:,jk)
           do jd=1,3
              Hkd(jb1,jb2, jd, map(1,jk)) = &
                   sum(symop(jd, :,map(2,jk)) * derivative(:))
           end do
        end do
     end do

     Hk(:,:,  map(1,jk)) = Hkfull(:,:,jk)
     mapcount(map(1,jk)) = mapcount(map(1,jk)) + 1
  end do

  normalize_H: do jk=1,Nk
     Hkd(:,:,:,jk) = Hkd(:,:,:,jk)/mapcount(jk)
  end do normalize_H

  !internal hopping
  if (inw%intrahop) then
     write(unit_outwop,'(A)') '>>> add internal unitcell hopping …'
     write(unit_outwop,*)
     !read-in distance matrix
     allocate(distmatrix(WFmin:WFmax, WFmin:WFmax, 3))

     do jd=1,3
        do jb=WFmin,WFmax
           read(unit_intrahop,*) distmatrix(jb, :, jd)
        end do
     end do
     close(unit_intrahop)

     do jk=1,Nk
        do jd=1,3
           Hkd(WFmin:WFmax,WFmin:WFmax, jd,jk) =     &
                Hkd(WFmin:WFmax,WFmin:WFmax, jd,jk) &
                - (0,1)*distmatrix(:,:,jd)*Hk(:,:,jk)
        end do
     end do

     deallocate(distmatrix)
  end if
end subroutine get_hkderiv

subroutine read_ham(lun, ham, bmin, bmax, Nk)
  use const, only: DPk

  integer,      intent(in)  :: lun, bmin, bmax, Nk
  complex(DPk), intent(out) :: ham(bmin:bmax, bmin:bmax, 1:Nk)

  real(DPk) :: tmpr_Nb(2*(bmax-bmin+1))
  integer   :: jk, jb, Nb, Nk2

  read(lun,*)  Nk2,Nb
  if (Nk /= Nk2) &
       call croak("number of k-points inconsistent, &
       & check ‘.klist’ and ‘.hk’ files")
  if (Nb /= bmax-bmin+1) &
       call croak("number of bands inconsistent, &
       & check ‘.inwop’ and ‘.hk’ files")

  !     read-in Hamiltonian
  do jk=1,Nk
     read(lun,*)
     do jb=bmin,bmax
        read(lun,*) tmpr_Nb;
        ham(jb, bmin:bmax, jk) = cmplx(tmpr_Nb(1::2), tmpr_Nb(2::2), DPk)
     end do
  end do

  if (inw%wfrot) then
     write(unit_outwop,'(A)') ">>> rotating …"
     write(unit_outwop,*)
     do jk=1,Nk
        ham(:,:, jk) = &
             matmul(conjg(transpose(rotmat)), matmul(ham(:,:, jk), rotmat))
     end do
     deallocate(rotmat)
  end if
end subroutine read_ham

subroutine read_energy(lun, ene, bmin, bmax, Nk)
  use const, only: Ryd_eV, DPk

  integer,   intent(in)  :: lun, bmin, bmax, Nk
  real(DPk), intent(out) :: ene(bmin:bmax, 1:Nk)

  integer   :: ii, jk, jb, Nbloc
  real(DPk) :: E, EFermi

  read (unit_fermi,*) EFermi
  close(unit_fermi)
  write(unit_outwop, "('Wien2k Fermi energy: ',G10.6E1,' Ry')") EFermi

  do ii=1,stru%nneq
     read(lun,*)
     read(lun,*)
  end do
  do jk=1,size(ene,2)
     read(lun,"(73X,i6)") Nbloc
     if (Nbloc < bmax) &
          call croak('energy file does not span outer band window')

     do jb=1,bmin-1
        read(lun,*)
     end do

     do jb=bmin,bmax
        read(lun,*) ii, E
        ene(jb, jk) = (E-EFermi)*Ryd_eV
     end do

     do jb = bmax+1, Nbloc
        read(lun,*)
     end do
  end do
end subroutine read_energy

function read_SE(lun, Nb, wmin, wmax) result(SE)
  use const, only: DPk

  integer,      intent(in)  :: lun, Nb, wmin, wmax
  complex(DPk)              :: SE(1:Nb, wmin:wmax)

  complex(DPk), allocatable :: SEin(:,:)
  real   (DPk), allocatable :: wSE(:)
  real   (DPk)              :: tmpr_Ninter(2*Nb)
  real   (DPk)              :: weight
  integer                   :: jw, jw1, jw2, jwhi, jwlo, Nw, Nb2, Nw_in

  Nw = size(SE,2)
  read(lun,*) Nw_in, Nb2
  if (Nb2 /= Nb) &
       call croak("#bands in ‘.selfE’ is inconsistent with ‘.inwop’")

  allocate(SEin(Nb, Nw_in), wSE(Nw_in))
  do jw=1,Nw_in
     read(lun,*) wSE(jw), tmpr_Ninter
     SEin(:,jw) = cmplx(tmpr_Ninter(1::2), tmpr_Ninter(2::2), DPk)
  end do

  jw1 = 2

  !! interpolate self-energy at our frequencies linearly from input
  !! self-energy (at arbitrary frequencies)
  jwlo = max(ceiling(wSE(1)    /inw%dw), wmin)
  jwhi = min(floor  (wSE(Nw_in)/inw%dw), wmax)
  do jw2 = jwlo,jwhi
     !! find input frequencies that enclose the target frequency
     do
        if (wSE(jw1) >= inw%w(jw2)) exit
        jw1 = jw1 + 1
     end do

     weight = (wSE(jw1) - inw%w(jw2)) / (wSE(jw1) - wSE(jw1-1))

     SE(:, jw2) = SEin(:,jw1-1) * weight + &
          &       SEin(:,jw1  ) * (1-weight) 
  end do

  if (inw%w(wmin) < wSE(1) .or. inw%w(wmax) > wSE(Nw_in)) then
     call carp("input self-energy does not cover necessary frequencies [" &
          //str(inw%w(wmin))//", "&
          //str(inw%w(wmax))//"].")
     write(ERROR_UNIT, '(A)') &
          'Will use boundary values for outer frequencies'

     do jw = wmin, jwlo-1
        SE(:, jw) = SE(:, jwlo)
     end do
     do jw = jwhi+1, wmax
        SE(:, jw) = SE(:, jwhi)
     end do
  end if

  deallocate(SEin, wSE)

  where(abs(aimag(SE)) < inw%delterl) &
       SE = cmplx(real(SE), -inw%delterl, DPk)
end function read_SE

subroutine read_tetra(lun, tet, weigh)
  use const, only: DPk

  integer,   intent(in)  :: lun
  integer,   intent(out), allocatable :: tet(:,:)
  real(DPk), intent(out), allocatable :: weigh(:)

  integer :: jt, kdim(3)

  read (lun,*) Nt, kdim
  write(unit_outwop, "('reciprocal increments are [',3G10.4E1,']')") &
       1.0_DPk/kdim

  allocate(tet(10,Nt), weigh(Nt))
  do jt=1,Nt
     read(lun,*) tet(:,jt), weigh(jt)
  end do
end subroutine read_tetra

subroutine read_vk(lun, Vk)
  integer,      intent(in)  :: lun(:)
  complex(DPk), intent(out) :: Vk(:,:,:,:)

  integer   :: jk, jd, jb, Nk
  real(DPk) :: tmpr_Nb(2*size(Vk,1))

  do jd=1,size(lun)
     read(lun(jd),*) Nk
     if (Nk /= size(Vk, 4)) call croak("number of k-points inconsistent, &
          & check ‘.klist’ and ‘.vk’ files")
  end do

  do jk=1,Nk
     do jd=1,size(lun)
        read(lun(jd),*)
     end do

     do jd=1,size(lun)
        do jb=1,size(Vk,1)
           read(lun(jd),*) tmpr_Nb
           Vk(jb,:,jd,jk) = cmplx(tmpr_Nb(1::2), tmpr_Nb(2::2), DPk)
        end do
     end do
  end do
end subroutine read_vk

subroutine read_vvk_heads(lun)
  integer, intent(in) :: lun(:)

  integer :: tmp(2)

  do i=1,size(lun)
     call maybin_read(lun(i), tmp)
     if (tmp(1) /= Nk*Nw) call croak("#(k-points) × #(internal frequencies) &
          &inconsistent, check ‘.klist’ and ‘.vvk’ files")
     if (tmp(2) /= NWF) call croak("number of bands inconsistent, &
          &check ‘.inwop’ and ‘.vvk’ files")
  end do
end subroutine read_vvk_heads

subroutine read_vvk_kpt(lun, VAV)
  use const, only: DPk

  integer,      intent(in)  :: lun
  complex(DPk), intent(out) :: VAV(:,:,:)

  integer   :: jb, jw
  real(DPk) :: tmpr_Nb(2*size(VAV,1))

  do jw = lbound(VAV,3), ubound(VAV,3)
     call maybin_read(lun)

     do jb = lbound(VAV,1), ubound(VAV,1)
        call maybin_read(lun, tmpr_Nb)
        VAV(jb, :, jw) = cmplx(tmpr_Nb(1::2), tmpr_Nb(2::2), DPk)
     end do
  end do
end subroutine read_vvk_kpt

pure elemental real(DPk) function specfun(jb, jw, jk)
  use const, only: PI

  integer, intent(in) :: jb, jw, jk

  specfun = inw%delterl/PI / ((inw%w(jw) - bands(jb,jk))**2 + inw%delterl**2)
end function specfun

function specmat(jw, jk)
  use const, only: DPk

  integer, intent(in) :: jw, jk
  complex(DPk)        :: specmat(WFmin:WFmax, WFmin:WFmax)

  complex(DPk) :: G(WFmin:WFmax, WFmin:WFmax)
  integer      :: jb

  !! prepare Green function inverse
  G = -Hk(:,:,jk)
  do jb=WFmin,WFmax
     G(jb,jb) = G(jb,jb) + inw%w(jw) - SE(jb,jw)
  end do

  call invert_matrix(G)
  specmat = (G - transpose(conjg(G))) * (0,.5_DPk)/PI
end function specmat

subroutine invert_matrix(A)
  use const, only: DPk
  use clio,  only: croak
  use util,  only: string

  complex(DPk), intent(inout) :: A(:,:)

  complex(DPk) :: work(64*size(A,1))
  integer      :: ipiv(size(A,1)), n, info

  n = size(A,1)

  call zgetrf(n,n, A, n, ipiv, info)
  if (info /= 0) call croak("error in zgetrf(), info="//str(info))
  call zgetri(n,   A, n, ipiv, work, size(work), info)
  if (info /= 0) call croak("error in zgetri(), info="//str(info))
end subroutine invert_matrix

subroutine get_umatrix(H,nw,U)
  use iso_fortran_env, only: ERROR_UNIT
  use util,            only: string
  use clio,            only: carp

  implicit none

  integer nw,info,j,i
  complex(DPk) :: H(nw,nw),work(2*nw),u(nw,nw)
  real   (DPk) :: eigval(nw),rwork(7*nw)

  do i=1,nw
     do j=1,nw
        U(i,j) = H(i,j)
     end do
  end do
  call ZHEEV( 'V', 'U', nw, U, nw, eigval, work, 2*nw, rwork,info )

  if (info.ne.0) then
     call carp("Error: in zheev lapack routine, info:" // str(info))
     write(ERROR_UNIT, '("Hamiltonian:")')
     do i=1,nw
        do j=1,nw
           write(ERROR_UNIT,*)i,j,H(i,j)
        end do
     end do
     call croak()
  end if
end subroutine get_umatrix


subroutine compute_tetraderiv(tetra, k, Hloctet, node, derivative)
  !computes the derivative within a tetrahedron along cartesian coordinate dimidx

  use const, only: DPk
  use util,  only: string
  use clio,  only: croak

  implicit none

  integer,      intent(in)  :: tetra(10), node
  real   (DPk), intent(in)  :: k(:,:)
  complex(DPk), intent(in)  :: Hloctet(10)
  complex(DPk), intent(out) :: derivative(3)

  integer      :: dNdt(10,4,10), ipiv(4), j1, j2, info
  complex(DPk) :: RHS(4,3), func(4), A(4,4)
  real   (DPk) :: kk(10,3)

  kk(:,:) = k(tetra(:), :)

  !assemble dNdt
  dNdt = 0
  !at node1        !at node2        !at node3         !at node4       
  dNdt(1,1,1)= 3;  dNdt(1,1,2)=-1;  dNdt( 1,1,3)=-1;  dNdt( 1,1,4)=-1;
  dNdt(2,2,1)=-1;  dNdt(5,1,2)= 4;  dNdt( 7,1,3)= 4;  dNdt( 8,1,4)= 4;
  dNdt(5,2,1)= 4;  dNdt(2,2,2)= 3;  dNdt( 2,2,3)=-1;  dNdt( 2,2,4)=-1;
  dNdt(3,3,1)=-1;  dNdt(3,3,2)=-1;  dNdt( 6,2,3)= 4;  dNdt( 9,2,4)= 4;
  dNdt(7,3,1)= 4;  dNdt(6,3,2)= 4;  dNdt( 3,3,3)= 3;  dNdt( 3,3,4)=-1;
  dNdt(4,4,1)=-1;  dNdt(4,4,2)=-1;  dNdt( 4,4,3)=-1;  dNdt(10,3,4)= 4;
  dNdt(8,4,1)= 4;  dNdt(9,4,2)= 4;  dNdt(10,4,3)= 4;  dNdt( 4,4,4)= 3;

  !at node5        !at node6         !at node7         !at node8       
  dNdt(5,1,5)= 2;  dNdt( 1,1,6)=-1;  dNdt( 1,1,7)= 1;  dNdt( 1,1,8)= 1;
  dNdt(1,1,5)= 1;  dNdt( 5,1,6)= 2;  dNdt( 7,1,7)= 2;  dNdt( 8,1,8)= 2;
  dNdt(2,2,5)= 1;  dNdt( 7,1,6)= 2;  dNdt( 2,2,7)=-1;  dNdt( 2,2,8)=-1;
  dNdt(5,2,5)= 2;  dNdt( 2,2,6)= 1;  dNdt( 5,2,7)= 2;  dNdt( 5,2,8)= 2;
  dNdt(3,3,5)=-1;  dNdt( 6,2,6)= 2;  dNdt( 6,2,7)= 2;  dNdt( 9,2,8)= 2;
  dNdt(6,3,5)= 2;  dNdt( 3,3,6)= 1;  dNdt( 3,3,7)= 1;  dNdt( 3,3,8)=-1;
  dNdt(7,3,5)= 2;  dNdt( 6,3,6)= 2;  dNdt( 7,3,7)= 2;  dNdt( 7,3,8)= 2;
  dNdt(4,4,5)=-1;  dNdt( 4,4,6)=-1;  dNdt( 4,4,7)=-1;  dNdt(10,3,8)= 2;
  dNdt(8,4,5)= 2;  dNdt( 9,4,6)= 2;  dNdt( 8,4,7)= 2;  dNdt( 4,4,8)= 1;
  dNdt(9,4,5)= 2;  dNdt(10,4,6)= 2;  dNdt(10,4,7)= 2;  dNdt( 8,4,8)= 2;

  !at node9         !at node10        
  dNdt( 1,1,9)=-1;  dNdt( 1, 1,10)=-1;
  dNdt( 5,1,9)= 2;  dNdt( 7, 1,10)= 2;
  dNdt( 8,1,9)= 2;  dNdt( 8, 1,10)= 2;
  dNdt( 2,2,9)= 1;  dNdt( 2, 2,10)=-1;
  dNdt( 9,2,9)= 2;  dNdt( 6, 2,10)= 2;
  dNdt( 3,3,9)=-1;  dNdt( 9, 2,10)= 2;
  dNdt( 6,3,9)= 2;  dNdt( 3, 3,10)= 1;
  dNdt(10,3,9)= 2;  dNdt(10, 3,10)= 2;
  dNdt( 4,4,9)= 1;  dNdt(4 , 4,10)= 1;
  dNdt( 9,4,9)= 2;  dNdt(10, 4,10)= 2;

  do j1=1,4
     func(j1) = sum(dNdt(:,j1,node) * Hloctet(:))
  end do

  A = 0
  A(1,:) = 1
  do j1=1,3
     do j2=1,4
        A(1+j1,j2) = sum(dNdt(:, j2, node) * kk(:,j1))
     end do
  end do

  RHS = 0
  do j1=1,3
     RHS(j1+1,j1) = 1
  end do

  call zgetrf(4,4,A,4,ipiv,info)
  if (info.ne.0) &
       call croak("ERROR in lapack zgetrf, INFO="//str(info))

  call zgetrs('N',4,3,A,4,ipiv,RHS,4,info)
  if (info.ne.0) &
       call croak("ERROR in lapack zgetrs, INFO="//str(info))

  do j2=1,3
     derivative(j2) = 1 + sum(RHS(:,j2)*func(:))
  end do
end subroutine compute_tetraderiv

end program

#ifdef WANNIER_RANDOM_PHASE
! from the GCC documentation,
! <https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html>,
! adapted for old ifort
subroutine init_random_seed()
  ! use iso_fortran_env, only: int64
  use ifport, only: getpid
  use util, only: newunit
  implicit none
  integer, parameter :: int64=8
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit(un), file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
#endif

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-11-10 15:12:52 assman@faepop36.tu-graz.ac.at>
