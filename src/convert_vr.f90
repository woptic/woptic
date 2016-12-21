!!! woptic/src/convert_vr.f90
!!!
!!!    Fourier-transforms dipole matrix-elements V(r) back to k-space,
!!!    V(k)
!!!
!!! Copyright 2012      Philipp Wissgott
!!!           2013-2016 Elias Assmann
!!!

PROGRAM convert_vr
  use util,      only: ptime, string
  use const,     only: DPk, PI
  use clio,      only: croak, argstr, fetcharg
  use structmod, only: struct_t, struct_read
  use maybebin,  only: maybin_open, maybin_read, maybin_write
  use woptic,    only: inwop_t, inwop_read, &
       fmt_rweights, fmt_vr_read, fmt_vvr_read, fmt_vk_head, fmt_vk_kp, &
       fmt_vk, fmt_vvk_head, fmt_vvk_kp, fmt_vvk
  use kpoints,   only: get_kmesh_klist
  use woptic_io, only: set_casename, &
       ! input:
       unit_vr,    unit_vvr,                                       &
        suf_vr,     suf_vvr,  suf_struct,  suf_inwop,  suf_klist,  &
         fn_vr,      fn_vvr,   fn_struct,   fn_inwop,   fn_klist,  &
       ! output:
       unit_outvk, unit_vk,  unit_vvk,                             &
        suf_outvk,                                                 &
         fn_outvk,   fn_vk,    fn_vvk

  implicit none

  character(*), parameter :: rev_str = "$version: v0.1.0-81-gc103383$"
  character(*), parameter :: woptic_version = rev_str(11 : len (rev_str)-1)

  integer        :: v, w, k, r, o, i, id, omin, omax
  integer        :: Nb, Nr, Nk, Nw, iarg, check(3)
  real(DPk)      :: imre(12), rdotk
  logical        :: binout=.true., binin
  character(2)   :: updn=''

  type(struct_t) :: stru
  type(argstr)   :: arg, case
  type(inwop_t)  :: inw

  complex(DPk), allocatable, dimension(:,:,:,:)   :: vr,  vk
  complex(DPk), allocatable, dimension(:,:,:,:,:) :: vvr, vvk
  complex(DPk), allocatable :: fac(:,:)
  real(DPk),    allocatable :: kpts(:,:)
  integer,      allocatable :: rweights(:), rvec(:,:)

  !! formats for log file
  character(*), parameter ::            &
       fmtII = "(A40, I6, I5)",         &
       fmtA  = "(A40, A6)"

  !! formats for output
  character(*), parameter ::                                      &
       fmt_help = '(A, T10, "case", A, T25, A, T50, A)',          &
       fmt_hlp2 = '(A, T10,         A, T25, A, T50, A)'

  !! just a string ...
  character(*), parameter :: A = "(A)"


!!! Argument parsing
  if (command_argument_count() < 1) &
       call croak('Usage: convert_vr [OPTIONS] case')

  do iarg=1,command_argument_count()
     call fetcharg(iarg, arg)
     select case (arg%s)
     case ('-h','-H','-help','--help')
        write(*,A) "convert_vr: convert direct-space dipole &
             &matrix elements back to k-space"
        write(*,*) 
        write(*,fmt_hlp2) "USAGE:", "convert_vr [--text] [--up|--dn] CASE"
        write(*,*) 
        write(*,fmt_help) "INPUT:", suf_struct, "Wien2k master input file"
        write(*,fmt_help) "", suf_klist, "target k-points"
        write(*,fmt_help) "", suf_inwop, 'woptic input file'
        write(*,fmt_help) "", suf_vr, "direct-space matrix elements V(R)"
        write(*,fmt_help) "", suf_vvr, "direct-space matrix elements VAV", &
             '[mixed transitions]'
        write(*,*) 
        write(*,fmt_help) "OUTPUT:", suf_outvk, "log file"
        write(*,fmt_help) "",  '.vk{x,y,z}', "k-space matrix elements V(k)"
        write(*,fmt_help) "", '.vvkIJ', "k-space matrix elements VAV"
        write(*,*) 
        write(*,fmt_hlp2) "OPTIONS:", "-h, --help"
        write(*,fmt_hlp2) "", "-v, --version"
        write(*,fmt_hlp2) "", "-t, --text", "read/write VAV in plain text"
        write(*,fmt_hlp2) "", "--up|--dn", "spin-polarized mode"
        call exit(0)

     case ('-v', '-V', '-version', '--version')
        print '("convert_vr ", A)', WOPTIC_VERSION
        call exit(0)

     case('-t', '--text')
        binout = .false.
     case ('-up', '--up')
        updn='up'
     case ('-dn', '--dn')
        updn='dn'

     case default
        if (arg%s(1:1) == '-') &
             call croak('unknown option: '//trim(arg%s))

        call fetcharg(iarg, case)
     end select
  end do

  call set_casename(case, UPDN=updn)

!!! Read inwop, struct
  call struct_read(fn_struct, stru)
  call  inwop_read(fn_inwop,  inw)
  omin = lbound(inw%w,1)
  omax = merge(ubound(inw%w,1), omin-1, inw%mixed)
  Nw   = omax-omin+1

  if (size(inw%ishift) /= 0 .and. inw%mixed) call croak( &
       "FIXME: scissors shift with VAV interpolation is unimplemented")

!!! Read klist
  call get_kmesh_klist(fn_klist, kpts, stru)
  nk = size(kpts, 1)

!!! Open output file
  open (unit_outvk, FILE=fn_outvk, STATUS='replace')
  write(unit_outvk, A) &
       "------------- k-Interpolated Dipole Matrix Elements -------------"
  write(unit_outvk,fmtA ) "mixed transitions?", &
       merge("yes", " no", inw%mixed)
  if (inw%mixed) &
       write(unit_outvk,fmtII) "number of frequencies:", Nw
  write(unit_outvk,*)
  call ptime(UNIT=unit_outvk)

  flush(unit_outvk)

!!! Read V(R), [VAV(R)]
  open(unit_vr, file=fn_vr, status='old')
  if (inw%mixed) then
     call maybin_open(unit_vvr, file=fn_vvr, status='old', get=binin)
     binout = binout .and. binin
     write(unit_outvk,fmtA) "binary output?", merge("yes", " no", binout)
     flush(unit_outvk)
  end if
  read(unit_vr, *)   ; if (inw%mixed) call maybin_read(unit_vvr)! label
  read(unit_vr, *) Nb; if (inw%mixed) call maybin_read(unit_vvr, check(1:2))
  read(unit_vr, *) Nr; if (inw%mixed) call maybin_read(unit_vvr, check(3:3))

  if (Nb/=(inw%bmax-inw%bmin+1)) &
       call croak("`inwop' and `vr' do not agree on #bands: "   &
       &          //trim(string(inw%bmax-inw%bmin+1))//' /= ' &
       &          //trim(string(Nb)))

  if (inw%mixed .and. any(check /= (/ Nb, Nw, Nr /))) &
       call croak("`vvr' file inconsistent -- (Nb,Nw,Nr) = (" &
       &          //trim(string(check(1)))//','               &
       &          //trim(string(check(2)))//','               &
       &          //trim(string(check(3)))//') vs ('          &
       &          //trim(string(Nb))//','                     &
       &          //trim(string(Nw))//','                     &
       &          //trim(string(Nr))//')')

  write(unit_outvk,fmtII) &
       "number of k-points/Wannier functions:", Nk, Nb
  flush(unit_outvk)

  allocate(rweights(nr),             &
       &   fac     (nr, nk),         &
       &   rvec    (3, nr),          &
       &   vr      (3, nb,nb, nr),   &
       &   vk      (3, nb,nb, nk))

  allocate(vvr(6, nb,nb, nr, omin:omax), &
       &   vvk(6, nb,nb, nk, omin:omax))

  ! For binary ‘vvr’, rweights is *one* record.
  if (inw%mixed .and. binin) read(unit_vvr) rweights
  do r=1,nr/15
     read(unit_vr, fmt_rweights) rweights((r-1)*15+1 : min(r*15,nr))
     if (inw%mixed .and. .not. binin) read(unit_vvr,*)
  enddo
  if (mod(nr, 15) /= 0) then
     r = nr/15+1
     read(unit_vr, fmt_rweights) rweights((r-1)*15+1 : min(r*15,nr))
     if (inw%mixed .and. .not. binin) read(unit_vvr,*)
  endif

  call ptime('read preparation')

  read_vr: do r=1,Nr
     do v=1,Nb
        do w=1,Nb
           read(unit_vr, fmt_vr_read) rvec(:,r), imre(1:2*size(vr,1))

           do id=1,size(vr,1)
              vr(id, v,w, r) = cmplx(imre(2*id-1), imre(2*id), DPk)
           end do
        enddo
     enddo
  enddo read_vr
  call ptime('reading V(R)')

  read_vvr: do o=omin,omax
     do r=1,Nr
        do v=1,Nb
           do w=1,Nb
              call maybin_read(unit_vvr, imre, fmt=fmt_vvr_read)

              do id=1,size(vvr,1)
                 vvr(id, v,w, r, o) = cmplx(imre(2*id-1), imre(2*id), DPk)
              end do
           end do
        end do
     end do
  end do read_vvr
  call ptime('reading VAV(R, w)')

  vk=0
  furry_transformer: do k=1,Nk
     do r=1,Nr
        rdotk    = 2*PI* dot_product(kpts(k,:), &
             &                       rvec(:,r))
        fac(r,k) = exp((0,1)*rdotk) / rweights(r)

        vk(:,:,:,k) = vk(:,:,:,k) + fac(r,k)*vr(:,:,:,r)
     end do
  end do furry_transformer
  call ptime('FT to V(k)')

  vvk=0
  furrier_transformer: do o=omin,omax
     do k=1,Nk
        do r=1,Nr
           vvk(:,:,:,k,o) = vvk(:,:,:,k,o) + fac(r,k)*vvr(:,:,:,r,o)
        end do
     end do
  end do furrier_transformer
  call ptime('FT to VAV(k, w)')

  write_vk: do i=1,size(vk,1)
     open (unit_vk(i), file=fn_vk(i), status='replace')
     write(unit_vk(i), fmt_vk_head) Nk, Nb
     
     do k=1,Nk
        write(unit_vk(i), fmt_vk_kp) kpts(k,:)
        do w=1,Nb
           write(unit_vk(i), fmt_vk) vk(i, w, :, k)
        end do
     end do

     close(unit_vk(i))
  end do write_vk

  write_vvk: do i=1,merge(size(vvk,1), 0, inw%mixed)
     call maybin_open (unit_vvk(i), fn_vvk(i), status='replace', set=binout)
     call maybin_write(unit_vvk(i), (/ Nk*Nw, Nb /), fmt=fmt_vvk_head)

     do k=1,Nk
        do o=omin,omax
           call maybin_write(unit_vvk(i), (/ kpts(k,:),inw%w(o) /), fmt_vvk_kp)
           do w=1,Nb
              call maybin_write(unit_vvk(i), vvk(i, w, :, k, o), fmt_vvk)
           end do
        end do
     end do
     close(unit_vvk(i))
  end do write_vvk

  call ptime('writing V(k)')
end program convert_vr
