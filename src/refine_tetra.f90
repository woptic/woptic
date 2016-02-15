!!! woptic/src/refine_tetra.f90
!!!
!!!    Adaptive refinement of tetrahedral k-mesh 
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2016 Elias Assmann
!!!

PROGRAM refine_tmesh
  use util,      only: inverse3x3, string
  use const,     only: DPk
  use structmod, only: struct_t, struct_read
  use kpoints,   only: count_kmesh_klist
  use clio,      only: fetcharg, argstr, croak
  use woptic_io, only: set_casename, &
       unit_klist, unit_fklist, unit_tet,   unit_ftet, unit_contr, unit_voe,&
       unit_struct,unit_map,    unit_inwop, unit_kadd, unit_outref,         &
       suf_klist,  suf_fklist,  suf_tet,    suf_ftet,  suf_contr,  suf_voe, &
       suf_struct, suf_map,     suf_inwop,  suf_kadd,  suf_outref,          &
       fn_klist,   fn_fklist,   fn_tet,     fn_ftet,   fn_contr,   fn_voe,  &
       fn_struct,  fn_map,      fn_inwop,   fn_kadd,   fn_outref,           &
       suf_rfd

  implicit none

  character(*), parameter :: rev_str = "$version: v0.1.0-45-gf5b26c6$"
  character(*), parameter :: woptic_version = rev_str(11 : len (rev_str)-1)

  integer   :: Nk, Nkfull, Nt, Nev, Nvoe, Nnewt, NE, Nsig, Nnewk, Nnewsk, Nkp
  integer   :: iarg, iv, iE, ik, jk, ir, it, jd1, tcount, kdiv, jsig
  integer   :: ndim(3), kdim(3)
  real(DPk) :: dE, dwei, emin, emax, rvec(3)
  
  real(DPk) :: rdum
  integer   :: idum

  real(DPk), allocatable, dimension(:) :: &
       nabk(:), tmp(:), newwtetra(:), tvariance(:), wtetra(:)
  integer,   allocatable, dimension(:,:) :: &
       k, kfull, kfulltmp, newk, tetra, newtetra, map, newmap, VOEidx
  real(DPk), allocatable :: kfull_int(:,:), kcontribw(:,:,:)
  integer,   allocatable :: tclass(:), newtclass(:), VOE(:,:,:)
  logical,   allocatable :: trefine(:)

  real(DPk), dimension(3,3) :: RR, GG

  logical        :: have_file=.false.
  logical        :: inter=.false., init=.false.
  character(2)   :: updn=''
  real(DPk)      :: theta=0.5_DPk
  type(argstr)   :: arg, file
  integer        :: init_steps
  type(struct_t) :: stru

  character(*), parameter :: usage = '(                                      &
&"refine_tetra: adaptive refinement of tetrahedral k-mesh"                 /&
&/"USAGE",                                                                   &
&T10,"refine_tetra [OPTIONS] CASE",                                        /&
&/"FILES",                                                                  &
&T10,"(prefixed by CASE)",                                                 /&
&"   '//suf_inwop //'",T20,"woptic main input file (for dE)"               /&
&"   '//suf_struct//'",T20,"Wien2k master input file"                      /&
&"   '//suf_contr //'",T20,"function values for estimator (w.r.t klist)"   /&
&"   '//suf_outref//'",T20,"log file"                                      /&
&" * '//suf_klist //'",T20,"symmetrized k-points"                          /&
&" * '//suf_fklist//'",T20,"unsymmetrized k-points"                        /&
&" * '//suf_tet   //'",T20,"symmetrized tetrahedra"                        /&
&" * '//suf_ftet  //'",T20,"unsymmetrized tetrahedra"                      /&
&" * '//suf_voe   //'",T20,"list of k-points on tetrahedral edges"         /&
&" * '//suf_map   //'",T20,"internal mapping of klist_full to klist"       /&
&/"Files marked ‘*’ are read and updated (except with --init).  The"       /&
&"updated file ‘F’ is written to ‘F'//suf_rfd//'’.  The list of added"     /&
&"k-points is written to ‘CASE'//suf_kadd//'’."                            /&
&/"OPTIONS",                                                                &
&T10,"--theta Θ", T21,"(0≤Θ≤1) defines the ‘harshness’ of refinement"      /&
&                 T20,"Θ=0: uniform;   Θ=1: most adaptive refinement"      /&
&T10,"--init N",                                                            &
&T20,"initial refinement with N steps (in general, 2 … 4)"                 /&
&T10,"--inter",                                                             &
&T20,"give larger weight to higher-energy contributions"                   /&
&T10,"--up|--dn",                                                           &
&T20,"spin-polarized calculation"                                          /&
&T10,"--help, --version"                                                    &
&)'


!!!------------- Argument handling              -----------------------------
  iarg=1
  arguments: do while (iarg <= command_argument_count())
     call fetcharg(iarg, arg)

     select case (arg%s)
     case ('--')
        call fetcharg(command_argument_count(), file)
        have_file=.true.
        exit arguments
     case ('-h','-H','-help','--help')
        print usage
        call exit(0)
     case ('-inter', '--inter')
        inter = .true.
     case ('-th', '--th', '-theta', '--theta')
        iarg=iarg+1
        call fetcharg(iarg, theta)
        if (theta<0 .or. theta>1) &
             call croak('0 ≤ Θ ≤ 1 must hold')
     case ('-init', '--init')
        iarg=iarg+1
        init=.true.
        call fetcharg(iarg, init_steps)
     case ('-v', '-V', '-version', '--version')
        print '("refine_tetra ", A)', WOPTIC_VERSION
        call exit(0)
     case ('-up', '--up')
        updn='up'
     case ('-dn', '--dn')
        updn='dn'
     case default
        if (arg%s(1:1) == '-') &
             call croak('unknown option ' // trim(arg%s))
        call fetcharg(iarg, file)
        have_file=.true.
     end select

     iarg = iarg+1
  end do arguments

  if (.not. have_file) &
       call croak('CASE argument must be given')

  call set_casename(file, UPDN=updn)


!!!------------- Open files for reading         -----------------------------
  if (.not. init) then
  open(unit_klist , FILE=fn_klist , STATUS='old')
  open(unit_fklist, FILE=fn_fklist, STATUS='old')
  open(unit_ftet  , FILE=fn_ftet  , STATUS='old')
  open(unit_voe   , FILE=fn_voe   , STATUS='old')
  open(unit_contr , FILE=fn_contr , STATUS='old')
  open(unit_map   , FILE=fn_map   , STATUS='old')
  end if
  open(unit_struct, FILE=fn_struct, STATUS='old')
  open(unit_inwop , FILE=fn_inwop , STATUS='old')

!!!------------- Open log file                  -----------------------------
  open(unit_outref, FILE=fn_outref, STATUS='replace')


!!!------------- Read input files               -----------------------------
  read (unit_inwop,*)
  read (unit_inwop,*) rdum, dE
  close(unit_inwop)

  call struct_read(unit_struct, stru)
  close(unit_struct)

  if (stru%ortho) then
     RR(1,:) = stru%brlat(1,:) / stru%a
     RR(2,:) = stru%brlat(2,:) / stru%a
     RR(3,:) = stru%brlat(3,:) / stru%a

     call inverse3x3(RR, GG)
     write(unit_outref,*) &
          "transformation matrix between wien2k and internal representation"
     write(unit_outref,"(3F12.5)") GG(1,:)
     write(unit_outref,"(3F12.5)") GG(2,:)
     write(unit_outref,"(3F12.5)") GG(3,:)
  end if

  nvoe = 0

  kdiv=0

  init_: if (init) then
     !get initial mesh
     nt = 48*8**(init_steps-1) 
     allocate(tvariance(nt),tclass(nt),newtclass(8*nt))
     nkfull = (3*2**init_steps)**3
     allocate(tetra(nt,10),wtetra(nt),newk(nkfull,3),newtetra(48*8**(init_steps-1),10),&
          newwtetra(48*8**(init_steps-1)),trefine(48*8**(init_steps-1)))
     allocate(VOE((3*2**init_steps)**3,300,2),VOEidx(6*48*8**(init_steps),3))
     VOE = 0
     VOEidx = 0
     call get_initialmesh(init_steps,tcount,newtetra,newwtetra,newtclass,nkfull,newk,nvoe,VOE,VOEidx)
     ndim = 2**(init_steps+1)
     nnewk = 0
     kdim = 2**init_steps
     nk = 0
  else
     !standard run: load k-mesh
     read(unit_ftet,*)nt,ndim
     allocate(tetra(nt,10),wtetra(nt),newk(48*nt,3),newtetra(8*nt,10),newwtetra(8*nt),trefine(nt))
     allocate(tvariance(nt),tclass(nt),newtclass(8*nt))
     do it=1,nt
        read(unit_ftet,*)tetra(it,:),wtetra(it),tclass(it)
     enddo
     close(unit_ftet)

     nk = count_kmesh_klist(unit_klist)
     allocate(k(nk,3),nabk(nk))
     do jk=1,nk
        if(jk.eq.1) then 
           read(unit_klist,1523) idum,(k(jk,ik),ik=1,3),kdiv, &
                dwei,emin,emax,nkp,kdim
        else
           read(unit_klist,1520) idum, (k(jk,ik),ik=1,3),kdiv,dwei
        endif
     enddo
     close(unit_klist)

     nkfull = count_kmesh_klist(unit_fklist)
     allocate(kfull(nkfull,3),kfulltmp(nkfull,3),map(nkfull,2),kfull_int(nkfull,3))
     do jk=1,nkfull
        if(jk.eq.1) then 
           read(unit_fklist,1523) idum,(kfull(jk,ik),ik=1,3),kdiv, &
                dwei,emin,emax,nkp,kdim
        else
           read(unit_fklist,1520) idum, (kfull(jk,ik),ik=1,3),kdiv,dwei
        endif
     enddo
     close(unit_fklist)
1523 FORMAT(I10,4I10,3f5.1,4x,i6,10x,3i3,1x) 
1520 FORMAT(I10,4I10,f5.1)     

     !map first to internal coordinates(real numerators) and than
     !back to integer values
     do jk=1,nkfull
        do jd1=1,3
           rvec(jd1) = real(kfull(jk,jd1), DPk)/real(ndim(jd1), DPk)
        enddo
        if (stru%ortho) rvec = matmul(RR, rvec)

        do jd1=1,3
           kfull(jk,jd1) = int(rvec(jd1)*ndim(jd1))
        enddo
     enddo

     !load mesh information
     read(unit_voe,*)nvoe
     allocate(VOE(48*nt,300,2),VOEidx(nvoe + 48*nt,3))
     VOE = 0
     VOEidx = 0
     do iv=1,nvoe
        read(unit_voe,*)VOEidx(iv,:)
        nev = 1
        do while (VOE(VOEidx(iv,1),nev,1).ne.0)
           nev = nev + 1
        enddo
        VOE(VOEidx(iv,1),nev,1) = VOEidx(iv,2)
        VOE(VOEidx(iv,1),nev,2) = VOEidx(iv,3)
     enddo
     do jk=1,nkfull
        nev = 0
        do while (VOE(jk,nev+2,1).ne.0)
           nev = nev + 1
        enddo
        VOE(jk,1,1) = nev
     enddo
     close(unit_voe)
  endif init_

  if (init) then
     NE = 0
     nsig = 1
     allocate(kcontribw(nk,0:NE,nsig),tmp(nsig))
     kcontribw = 1d0
  else
     write(unit_outref,*)"standard mode"
     read(unit_contr,*)idum,NE,nsig
     if (idum.ne.nk) call croak("number of k-points inconsistent, &
          &check .kcontribw and .klist")

     allocate(kcontribw(nk,0:NE,nsig),tmp(nsig))
     do jk=1,nk
        read(unit_contr,3000)idum, tmp
        do jsig=1,nsig
           kcontribw(jk,0,jsig) = tmp(jsig)  
        enddo
        do iE=1,NE
           read(unit_contr,3001)tmp(1:nsig)
           do jsig=1,nsig
              kcontribw(jk,iE,jsig) = tmp(jsig)  
           enddo
        enddo
     enddo
     close(unit_contr)
3000 FORMAT(I8,6E20.12)
3001 FORMAT(6E20.12)
  endif

  if (.not. init) then
     write(unit_outref,*)">>> refinement of tetrahedral mesh"
     write(unit_outref,*)"  theta=",theta
     write(unit_outref,*)"  number of initial tetrahedra:",nt
     write(unit_outref,*)"  number of initial full mesh k-points:",nkfull
     write(unit_outref,*)"  number of initial symmetrized k-points:",nk
     write(unit_outref,*)"  sum of tetrahedron weights before refinement:", &
          sum(wtetra)

     !get symmetry mesh information  
     kfulltmp = kfull

     do jk=1,nkfull
        read(unit_map,*)map(jk,:)
     enddo
     close(unit_map)

     !compute error estimators and mark elements
     call nullestimator(nt,nkfull,nk,tetra,tclass,wtetra,map,NE,nsig,&
          kcontribw,theta,trefine,inter,dE)

     newk(1:nkfull,:) = kfull !remember old k points
     !refine marked elements
     call refine_tetra(nt,tetra,tclass,wtetra,nkfull,48*nt,newk,&
          &            nvoe,VOE,VOEidx,trefine,8*nt,&
          &            newtetra,newwtetra,newtclass,nnewk,tcount) 
     ndim = ndim*2 !doubled divisor for new mesh
  endif

!!!------------- Open files for writing         -----------------------------
  open(unit_klist , FILE=trim(fn_klist )//suf_rfd, STATUS='replace')
  open(unit_fklist, FILE=trim(fn_fklist)//suf_rfd, STATUS='replace')
  open(unit_tet   , FILE=trim(fn_tet   )//suf_rfd, STATUS='replace')
  open(unit_ftet  , FILE=trim(fn_ftet  )//suf_rfd, STATUS='replace')
  open(unit_voe   , FILE=trim(fn_voe   )//suf_rfd, STATUS='replace')
  open(unit_map   , FILE=trim(fn_map   )//suf_rfd, STATUS='replace')
  open(unit_kadd  , FILE=trim(fn_kadd  ),          STATUS='replace')

!!!------------- Write updated files            -----------------------------
  write(unit_voe,*) nvoe
  do iv=1,nvoe
     write(unit_voe,*)VOEidx(iv,:)
  enddo
  close(unit_voe)

  !compute the shapes of the tetrahedra
  call compute_shapeparameters(nkfull+nnewk, tcount, newk(1:nkfull+nnewk,:), &
       &                       newtetra(1:tcount,:), ndim)

  do jk=1,nkfull+nnewk
     !transform to conventional reciprocal lattice vectors
     if (stru%ortho) then
        rvec = matmul(GG, real(newk(jk,:), DPk))
     else
        rvec = newk(jk,:)
     end if

     if (jk.eq.1) then
        write(unit_fklist,1533) jk, nint(rvec),ndim(1), &
             rdum,-7.,1.5,nkfull+nnewk,kdim
     else
        write(unit_fklist,1530) jk, nint(rvec),ndim(1),rdum
     endif
  enddo
  write(unit_fklist,1521)
  close(unit_fklist)

  write(unit_ftet,*) tcount,ndim 
  do it=1,tcount
     write(unit_ftet,"(10I10,E20.12,1I10)")(newtetra(it,ir),ir=1,10),newwtetra(it),newtclass(it)
  enddo
  close(unit_ftet)

  !symmetrize kmesh
  allocate(newmap(nkfull+nnewk,2))
  newmap = 0
  if (.not. init) then
     newmap(1:nkfull,:) = map !remember the old map
     k = 2*k !going to new integer mesh with doubled divisor
  endif
  write(unit_outref,*) "  symmetrize k-point mesh "
  call symmetrize_kmesh(nk, nkfull+nnewk, newk(1:nkfull+nnewk,:), k, &
       &                size(stru%ksym,3), stru%ksym, ndim, nnewsk, newmap)
  write(unit_outref,*) &
       "  number of full mesh k-points after refinement: ",nkfull+nnewk
  write(unit_outref,*) &
       "  number of symmetrized k-points after refinement: ",nnewsk
  write(unit_outref,*) "  symmetrize tetra mesh "

  !symmetrize tetrahedra
  call symmetrize_tetra(nkfull+nnewk,newmap,tcount,newtetra(1:tcount,:),newwtetra(1:tcount),nnewt)
  write(unit_outref,*) &
       "  number of full mesh tetrahedra after refinement: ",tcount
  write(unit_outref,*) &
       "  number of symmetrized tetrahedra after refinement: ",nnewt
  write(unit_outref,*) &
       "  sum of tetrahedron weights after refinement", sum(newwtetra)

  !write-out symmetrized mesh
  do jk=1,nnewsk
     !transform to conventional reciprocal vectors,
     !but only the new ones, since the old ones have already been transformed
     if (jk.gt.nk) then
        if (stru%ortho) then
           rvec = matmul(GG, real(newk(jk,:), DPk))
        else
           rvec = newk(jk,:)
        end if
     else
        rvec = k(jk,:)
     endif

     if (jk.eq.1) then
        write(unit_klist,1533) jk, nint(rvec),ndim(1), &
             rdum,-7.,1.5,nnewsk,kdim
     else
        write(unit_klist,1530) jk, nint(rvec),ndim(1),rdum
     endif
  enddo
  write(unit_klist,1521)
  close(unit_klist)

  !write-out map
  do jk=1,nkfull+nnewk
     write(unit_map,*)newmap(jk,:)
  enddo
  close(unit_map)

  do jk=nk+1,nnewsk
     !transform to conventional reciprocal lattice vectors
     if (stru%ortho) then
        rvec = matmul(GG, real(newk(jk,:), DPk))
     else
        rvec = newk(jk,:)
     end if

     if (jk.eq.nk+1) then
        write(unit_kadd,1533) jk,nint(rvec),2*kdiv, &
             rdum,-7.,1.5,nnewsk,kdim
     else
        write(unit_kadd,1530) jk, nint(rvec),2*kdiv,rdum
     endif
  enddo
  write(unit_kadd,1521)
  close(unit_kadd)

  write(unit_tet,*) nnewt,ndim
  do it=1,nnewt
     write(unit_tet,"(10I10,E20.12)")(newtetra(it,ir),ir=1,10),newwtetra(it) 
  enddo
  close(unit_tet)

1530 FORMAT(I10,4I10,f5.1)                                              
1521 format('END',/)
1533 FORMAT(I10,4I10,3f5.1,4x,i6,' k, div: (',3i3,')')    

CONTAINS
  FUNCTION kmean(kp1,kp2,ndim) result(outm)

    implicit none

    integer kp1(3),kp2(3),ndim(3),outm(3),jd

    outm = 1
    do jd=1,3
       if (kp1(jd).gt.ndim(jd)) then
          kp1(jd) = 2*ndim(jd)-kp1(jd)
       endif
       if (kp2(jd).gt.ndim(jd)) then
          kp2(jd) = 2*ndim(jd)-kp2(jd)
       endif
       outm(jd) = (kp1(jd)+kp2(jd))/2
    enddo
  END FUNCTION kmean



END PROGRAM refine_tmesh

SUBROUTINE desym(nk,k,tetra,ndim,outk) 

  implicit none

  integer nk,k(nk,3),tetra(4),ndim(3),outk(4,3)
  integer jd,idx1(4),idx2(4),jt

  outk = 0
  !   tmpk = 0
  do jt=1,4
     outk(jt,:) = k(tetra(jt),:)
  enddo
  do jd=1,3
     idx1 = 0
     idx2 = 0
     do jt=1,4
        if (outk(jt,jd).gt.ndim(jd)/2) then
           idx1(jt) = 1
        elseif (outk(jt,jd).eq.0) then
           idx2(jt) = 1
        endif
     enddo
     if ((any(idx1.eq.1)).and.(any(idx2.eq.1))) then
        do jt=1,4
           if (idx2(jt).eq.1) then
              outk(jt,jd) = ndim(jd)
           endif
        enddo
     endif
  enddo

END SUBROUTINE desym


!     *****************************************************************
FUNCTION find_kpoint2(nk,k,kvec) result(outp)
  implicit none

  integer outp,nk,k(max(nk,1),3),kvec(3),j

  if (nk.eq.0) then
     outp = -1
     return
  endif
  outp = 0
  do j=1,nk
     if ((k(j,1).eq.kvec(1)).and.(k(j,2).eq.kvec(2)).and.(k(j,3).eq.kvec(3))) then
        outp = j
        return
     endif
  enddo

END FUNCTION find_kpoint2


SUBROUTINE symmetrize_kmesh(nskold,nk,k,skold,nop,symop,Ndim,nsk,map)
  !on input k: unsymmetrized klist
  !on output k: symmetrized klist with rest zeros
  !the plan is to use the symmetry operations on the full kmesh
  !than order the kpoints and then take the first k point
  ! as representative for the symmetrized list
  ! It is thus necessary to remember the original index during the ordering
  use util,      only: string
  use clio,      only: croak
  use woptic_io, only: unit_outref

  implicit none

  integer :: nk,jk,js,jd1,jd2,nop
  integer :: k(nk,3),tmpk(nk*nop,5),tmpk2(nk*nop,5),nsk
  integer :: tmpk3(nk*nop,5),lock,nskold,mapidx(nk*nop),mapcount
  integer :: symop(3,3,nop)
  integer :: rotvec(3),ndim(3),skold(nskold,3)
  integer :: map(nk,2),counter,jk2
  integer :: idx(nk*nop),flux(ndim(1)+2),flux2(ndim(1)+2),fidx,fidx2,tmp(nk*nop),nf,nf2
  logical :: isold(nk)

  write(unit_outref,*)"  number of symmetry operations:",nop
  tmpk = 0
  nsk = 0
  isold = .false.

  counter =1
  tmpk = 0
  !act all symmetry operations on all unsymmetrized k-points
  do js=1,nop
     do jk=1,nk
        rotvec = 0   
        do jd1=1,3
           do jd2=1,3
              rotvec(jd1) = rotvec(jd1) + 2*int(symop(jd1,jd2,js))*k(jk,jd2)
           enddo
           rotvec(jd1) = mod(rotvec(jd1),2*Ndim(jd1))
           rotvec(jd1) = rotvec(jd1) + (1-isign(1,rotvec(jd1)))*Ndim(jd1)
        enddo

        tmpk(counter,1:3) = rotvec(1:3)/2
        tmpk(counter,4) = jk
        tmpk(counter,5) = js
        counter = counter + 1
     enddo
  enddo

  !order according to first dimensional entry
  call indexx(tmpk,idx,nop*nk)
  tmpk2 = 0
  do jk=1,nop*nk
     tmpk2(jk,:) = tmpk(idx(jk),:)
  enddo

  !now go for the other entries
  !the klist now has intervals of equal 1st entries
  !which have to be ordered according to the second entry
  flux(1) = 0
  fidx = 1
  do jk=2,nop*nk
     !find the number of intervals and their boundaries
     if (tmpk2(jk-1,1).ne.tmpk2(jk,1)) then
        fidx = fidx + 1
        flux(fidx) = jk - 1
     endif
  enddo
  flux(fidx+1) = nop*nk
  tmpk = 0
  tmpk3 = 0
  do jk=2,fidx+1
     nf = flux(jk) - flux(jk-1)
     tmp(1:nf) = tmpk2(flux(jk-1)+1:flux(jk),2)
     call indexx(tmp,idx,nf)
     tmpk(flux(jk-1)+1:flux(jk),:) = tmpk2(flux(jk-1)+idx(1:nf),:)

     !now the same game for the third entry
     fidx2 = 1
     flux2 = 0
     flux2(1) = flux(jk-1) 
     do jk2=flux(jk-1)+2,flux(jk)
        if (tmpk(jk2-1,2).ne.tmpk(jk2,2)) then
           fidx2 = fidx2 + 1
           flux2(fidx2) = jk2 - 1
        endif
     enddo
     flux2(fidx2+1) = flux(jk)
     do jk2=2,fidx2 + 1
        nf2 = flux2(jk2) - flux2(jk2-1)
        tmp(1:nf2) = tmpk(flux2(jk2-1)+1:flux2(jk2),3)
        call indexx(tmp,idx,nf2)
        tmpk3(flux2(jk2-1)+1:flux2(jk2),:) = tmpk(flux2(jk2-1)+idx(1:nf2),:)
     enddo
  enddo

  nsk = nskold

  k(1:nsk,:) = skold
  tmpk = 0
  tmpk(1,:) = tmpk3(1,:)
  lock = 1

  mapcount = 1
  mapidx = 0 
  mapidx(1) = map(tmpk3(1,4),1)!map(tmpk3(1,4),1)

  !now take a first look at the ordered list
  !and check whether the old symmetrized klist 
  !already covers a k point
  do jk = 2,nop*nk
     if ((tmpk3(jk,1).ne.tmpk3(jk-1,1)).or.&
          (tmpk3(jk,2).ne.tmpk3(jk-1,2)).or.&
          (tmpk3(jk,3).ne.tmpk3(jk-1,3))) then
        mapcount = mapcount + 1
     endif
     if (mapidx(mapcount).eq.0) then
        mapidx(mapcount) = map(tmpk3(jk,4),1)
     endif
  enddo


  !here comes the real mapping
  if (mapidx(1).eq.0) then
     nsk = nsk + 1
     map(tmpk3(1,4),1) = nsk
     mapidx(1) = nsk
  else
     map(tmpk3(1,4),1) = mapidx(1)
  endif
  map(tmpk3(1,4),2) = tmpk3(1,5)
  mapcount = 1
  do jk = 2,nop*nk
     if ((tmpk3(jk,1).ne.tmpk3(jk-1,1)).or.&
          (tmpk3(jk,2).ne.tmpk3(jk-1,2)).or.&
          (tmpk3(jk,3).ne.tmpk3(jk-1,3))) then
        mapcount = mapcount + 1
     endif
     if (map(tmpk3(jk,4),1).eq.0) then
        !k-point has no mapping
        if (mapidx(mapcount).eq.0) then
           !no other k-point in this symmetry group 
           !has a mapping -> it has to be a new one
           nsk = nsk + 1
           k(nsk,:) = tmpk3(jk,1:3)
           mapidx(mapcount) = nsk
           map(tmpk3(jk,4),1) = nsk
           map(tmpk3(jk,4),2) = tmpk3(jk,5)
        else
    !another k-point in this symmetry group 
    !has a mapping -> map to a kpoint from the old symmetrized klist 
           map(tmpk3(jk,4),1) = mapidx(mapcount)
           map(tmpk3(jk,4),2) = tmpk3(jk,5)
        endif
     endif
  enddo

  do jk=1,nk
     if (map(jk,1).eq.0) &
          call croak("Error in symmetrize k-mesh "//trim(string(jk)))
  enddo

END SUBROUTINE symmetrize_kmesh

SUBROUTINE symmetrize_tetra(maxk,map,nt,tetra,wtetra,countsymtetra)
  use const, only: DPk

  implicit none

  integer maxk,map(maxk),nt,tetra(nt,10),compvec,jc
  real(DPk) wtetra(nt),tmpw(nt),r1
  integer jt,jd,jd1,jd2,val1,val2,idx(nt),tmpt(nt,10)
  integer :: ichang,loctetra(10),it(10)
  integer :: countsymtetra

  !first map the tetrahedra vertices to symmetrized kmesh
  do jt=1,nt
     do jd=1,10       
        tetra(jt,jd) = map(tetra(jt,jd))
     enddo
  enddo

  !now order vertices for each tetrahedron from lowest to highest value
  do jt=1,nt
     do jd1=1,3                                                       
        do jd2=jd1+1,4                                                     
           val1=tetra(jt,jd1)                                                  
           val2=tetra(jt,jd2)                                                  
           tetra(jt,jd1)= min(val1,val2)                                     
           tetra(jt,jd2)= max(val1,val2)
        enddo
     enddo
  enddo
  do jt=1,nt
     do jd1=5,9                                                       
        do jd2=jd1+1,10                                                     
           val1=tetra(jt,jd1)                                                  
           val2=tetra(jt,jd2)                                                  
           tetra(jt,jd1)= min(val1,val2)                                     
           tetra(jt,jd2)= max(val1,val2)
        enddo
     enddo
  enddo

  !reorder tetra list according to the first vertex
  !heapsort from wien2k
  call indexx(tetra,idx,nt)
  do jt=1,nt
     tmpt(jt,:) = tetra(idx(jt),:)
     tmpw(jt) = wtetra(idx(jt))
  enddo


  do jd=2,10
     ichang = 1
     do while (ichang.eq.1)
        ichang = 0
        do jt=2,nt
           compvec = 0
           do jc=1,jd-1
              if (tmpt(jt,jc).eq.tmpt(jt-1,jc)) then
                 compvec = compvec + 1
              endif
           enddo
           if ((compvec.eq.jd-1).and.(tmpt(jt,jd).lt.tmpt(jt-1,jd))) then
              do jc=jd,10
                 it(jc)=tmpt(jt-1,jc)   
                 tmpt(jt-1,jc)=tmpt(jt,jc)
                 tmpt(jt,jc)=it(jc)
              enddo
              r1=tmpw(jt-1)
              tmpw(jt-1)=tmpw(jt)
              tmpw(jt)=r1
              ichang=1
           end if
        enddo
     enddo
  enddo

  !count the tetrahedra and sum up the weight
  countsymtetra = 1
  jt=1
  loctetra = tmpt(1,:)
  tetra = 0
  wtetra = 0d0
  do while (jt.le.nt)
     if ((tmpt(jt,1).eq.loctetra(1)).and.(tmpt(jt,2).eq.loctetra(2)).and.&
          (tmpt(jt,3).eq.loctetra(3)).and.(tmpt(jt,4).eq.loctetra(4)).and.& 
          (tmpt(jt,4).eq.loctetra(4)).and.(tmpt(jt,5).eq.loctetra(5)).and.&
          (tmpt(jt,6).eq.loctetra(6)).and.(tmpt(jt,7).eq.loctetra(7)).and.&
          (tmpt(jt,8).eq.loctetra(8)).and.(tmpt(jt,9).eq.loctetra(9)).and.&
          (tmpt(jt,10).eq.loctetra(10))) then
        wtetra(countsymtetra) = wtetra(countsymtetra) + tmpw(jt)
     else
        tetra(countsymtetra,:) = loctetra
        countsymtetra = countsymtetra + 1
        loctetra = tmpt(jt,:)
        wtetra(countsymtetra) = tmpw(jt)

     endif
     jt=jt+1
  enddo
  tetra(countsymtetra,:) = loctetra


END SUBROUTINE symmetrize_tetra

FUNCTION test1(k) result(outk)

  implicit none

  integer k,outk

  if (k.gt.16) then
     outk = 32 - k
  else
     outk = k
  endif
END FUNCTION test1

SUBROUTINE nullestimator(nt,nkfull,nk,tetra,tclass,wtetra,map,NE,nsig,kcontribw,theta,trefine,inter,dE)
  use const,           only: DPk
  use iso_fortran_env, only: ERROR_UNIT
  use woptic_io,       only: unit_outref

  implicit none

  integer   :: nt,nkfull,nk,jt,jd,tclass(nt)
  integer   :: tetra(nt,10),map(nkfull),NE,nsig,jsig,iE,reftet(8,4)
  integer   :: jref
  real(DPk) :: kcontribw(nk,0:NE,nsig),estimator(nt),maxest,theta,sumest
  real(DPk) :: wtetra(nt),est1,est2,est3,contrw,dE!,mmatrix(nkfull,nkfull),!refmat(4,4)!,RHS(nkfull)
  logical   :: trefine(nt),inter
  !   integer, allocatable ::  idx(:),ipiv(:)
  !   real(DPk), allocatable :: mmatrixred(:,:),RHSred(:)


  write(unit_outref,*)" error estimators for frequencies [w errorestimator]"
  estimator = 0d0
  do iE=0,NE
     contrw = 0d0
     do jsig=1,nsig
        do jt=1,nt
           reftet(1,1:4) = tetra(jt, (/ 1, 5, 6, 7  /))
           reftet(2,1:4) = tetra(jt, (/ 2, 5, 8, 9  /))
           reftet(3,1:4) = tetra(jt, (/ 3, 6, 8, 10 /))
           reftet(4,1:4) = tetra(jt, (/ 4, 7, 9, 10 /))
           if (tclass(jt).eq.1) then !central octahedron
              !corresp. to class 1 in Endres paper
              reftet(5,1:4) = tetra(jt, (/ 5, 6, 7, 10 /))
              reftet(6,1:4) = tetra(jt, (/ 5, 6, 8, 10 /))
              reftet(7,1:4) = tetra(jt, (/ 5, 7, 9, 10 /))
              reftet(8,1:4) = tetra(jt, (/ 5, 8, 9, 10 /))
           else
              !corresp. to class 2 in Endres paper
              reftet(5,1:4) = tetra(jt, (/ 5, 6, 7, 9  /))
              reftet(6,1:4) = tetra(jt, (/ 5, 6, 8, 9  /))
              reftet(7,1:4) = tetra(jt, (/ 6, 7, 9, 10 /))
              reftet(8,1:4) = tetra(jt, (/ 6, 8, 9, 10 /))
           endif
           est1 = 0d0
           est2 = 0d0
           est3 = 0d0
           do jd=1,4
              est1 = est1 + wtetra(jt)*kcontribw(map(tetra(jt,jd)),iE,jsig)/4d0
              est2 = est2 - wtetra(jt)*kcontribw(map(tetra(jt,jd)),iE,jsig)/20d0
           enddo

           do jd=5,10
              est2 = est2 +&
                   wtetra(jt)*kcontribw(map(tetra(jt,jd)),iE,jsig)/5d0
           enddo

           do jref=1,8
              do jd=1,4
                 est3 = est3 + wtetra(jt)*kcontribw(map(reftet(jref,jd)),iE,jsig)/32d0
              enddo
           enddo

           if (inter) then

              estimator(jt) = estimator(jt) + dabs(est1-est3)/dble(nsig)*iE*dE!**2/dabs(est2-est3)/dble(nsig)
              contrw = contrw + dabs(est1-est3)/dble(nsig)*iE*dE!**2/dabs(est2-est3)/dble(nsig)
           else
              estimator(jt) = estimator(jt) + dabs(est1-est3)/dble(nsig)!**2/dabs(est2-est3)/dble(nsig)
              contrw = contrw + dabs(est1-est3)/dble(nsig)!**2/dabs(est2-est3)/dble(nsig)
           endif
        enddo
     enddo
     write(unit_outref,"(2F12.5)")iE*dE,contrw
  enddo

  !find maximal entry 
  maxest = 0d0
  do jt=1,nt
     if (estimator(jt).gt.maxest) then
        maxest = estimator(jt)
     endif
  enddo

  write(unit_outref,*) "  maximal error estimator: ",maxest
  !mark elements to refine
  trefine = .false.
  do jt=1,nt
     if (estimator(jt) >= theta*maxest) then
        trefine(jt) = .true.
     endif
  enddo

  sumest = 0
  write(unit_outref,*) " error estimators for elements &
       &[elementidx errorestimator marked weight vertices]"
  do jt=1,nt
     sumest = sumest + estimator(jt)
     write(unit_outref,"(I10,F13.8,L3,F13.8,4I10)") &
          jt,estimator(jt),trefine(jt),wtetra(jt),tetra(jt,1:4)
  enddo
  write(unit_outref,*) "estimator = ",sumest

END SUBROUTINE nullestimator

SUBROUTINE refine_tetra(nt,tetra,tclass,wtetra,nkfull,maxk,newk,nvoe,VOE,VOEidx,trefine,maxt, &
     newtetra,newwtetra,newtclass,nnewk,tcount) 
  use const,           only: DPk
  use util,            only: string
  use clio,            only: croak, carp
  use iso_fortran_env, only: ERROR_UNIT

  implicit none

  integer maxt,nt,maxk
  integer nkfull,tetra(nt,10),tclass(nt),newtetra(maxt,10),newtclass(maxt),tcount
  integer jt,npidx(12,2),VOE(maxk,300,2),nvoe,VOEidx(nvoe + 48*nt,3)
  integer jd,jk,kidx(6),locv,nev,jv,nnewk,newk(maxk,3)
  integer eidx1,eidx2,eidx11,eidx22
  real(DPk) :: wtetra(nt),newwtetra(maxt)
  logical   :: trefine(maxt), tchange

  tcount = nt
  newtetra  = 0
  newwtetra = 0
  newtclass = 0
  newtetra(1:nt,:) = tetra
  newwtetra(1:nt) = wtetra
  newtclass(1:nt) = tclass
  nnewk = 0
  !update the original k points
  do jk=1,nkfull
     newk(jk,:) = 2*newk(jk,:)
  enddo

  npidx(1,:) = (/ 1, 2 /)
  npidx(2,:) = (/ 1, 3 /)
  npidx(3,:) = (/ 1, 4 /)
  npidx(4,:) = (/ 2, 3 /)
  npidx(5,:) = (/ 2, 4 /)
  npidx(6,:) = (/ 3, 4 /)

  tchange = .true.
  tchange_loop: do while (tchange)
     tchange = .false.
     do jt=1,nt
        if (trefine(jt)) then
           if ( tetra(jt,1)==tetra(jt,2) .or. &
                tetra(jt,2)==tetra(jt,3) .or. &
                tetra(jt,3)==tetra(jt,4)) then
              call carp("Error: multiple vertices "//trim(string(jt)))
              write(ERROR_UNIT, *) tetra(jt,:)
              call croak
           else
              !built 8 new tetrahedron and weights
              newtetra(jt,1:4)       = tetra(jt, (/ 1, 5,  6,  7 /))
              newtetra(tcount+1,1:4) = tetra(jt, (/ 5, 2,  8,  9 /))
              newtetra(tcount+2,1:4) = tetra(jt, (/ 6, 8,  3, 10 /))
              newtetra(tcount+3,1:4) = tetra(jt, (/ 7, 9, 10,  4 /))
              newtclass(jt)       = tclass(jt)
              newtclass(tcount+1) = tclass(jt)
              newtclass(tcount+2) = tclass(jt)
              newtclass(tcount+3) = tclass(jt)
              if (tclass(jt).eq.1) then !central octahedron
                 !corresp. to class 1 in Endres paper
                 newtetra(tcount+4,1:4) = tetra(jt, (/  9, 7, 10,  5 /))
                 newtetra(tcount+5,1:4) = tetra(jt, (/ 10, 6,  5,  7 /))
                 newtetra(tcount+6,1:4) = tetra(jt, (/  9, 8,  5, 10 /))
                 newtetra(tcount+7,1:4) = tetra(jt, (/  6, 5,  8, 10 /))
                 newtclass(tcount+4) = -1
                 newtclass(tcount+5) = -1
                 newtclass(tcount+6) = +1
                 newtclass(tcount+7) = +1
              else
                 !corresp. to class 2 in Endres paper
                 newtetra(tcount+4,1:4) = tetra(jt, (/ 10, 6, 8, 9 /))
                 newtetra(tcount+5,1:4) = tetra(jt, (/ 10, 9, 7, 6 /))
                 newtetra(tcount+6,1:4) = tetra(jt, (/  5, 8, 6, 9 /))
                 newtetra(tcount+7,1:4) = tetra(jt, (/  5, 7, 9, 6 /))
                 newtclass(tcount+4) = -1
                 newtclass(tcount+5) = -1
                 newtclass(tcount+6) = +1
                 newtclass(tcount+7) = +1
              endif
              do jd=0,7
                 if (jd.eq.0) then
                    newwtetra(jt) = wtetra(jt)/8
                    npidx(1,:) = newtetra(jt, (/ 1, 2 /))
                    npidx(2,:) = newtetra(jt, (/ 1, 3 /))
                    npidx(3,:) = newtetra(jt, (/ 1, 4 /))
                    npidx(4,:) = newtetra(jt, (/ 2, 3 /))
                    npidx(5,:) = newtetra(jt, (/ 2, 4 /))
                    npidx(6,:) = newtetra(jt, (/ 3, 4 /))
                 else
                    newwtetra(tcount+jd) = wtetra(jt)/8
                    npidx(1,:) = newtetra(tcount+jd, (/ 1, 2 /))
                    npidx(2,:) = newtetra(tcount+jd, (/ 1, 3 /))
                    npidx(3,:) = newtetra(tcount+jd, (/ 1, 4 /))
                    npidx(4,:) = newtetra(tcount+jd, (/ 2, 3 /))
                    npidx(5,:) = newtetra(tcount+jd, (/ 2, 4 /))
                    npidx(6,:) = newtetra(tcount+jd, (/ 3, 4 /))
                 endif

                 do jk=1,6
                    kidx(jk) = 0
                    eidx1 = min(npidx(jk,1),npidx(jk,2))
                    eidx2 = max(npidx(jk,1),npidx(jk,2))
                    nev = VOE(eidx1,1,1)
                    do jv=1,nev
                       if (VOE(eidx1,jv+1,1).eq.eidx2) then
                          kidx(jk) = VOE(eidx1,jv+1,2)
                          exit
                       endif
                    enddo
                    if (kidx(jk).eq.0) then
                       nnewk = nnewk + 1
                       kidx(jk) = nkfull+nnewk
                       newk(kidx(jk),:) = (newk(npidx(jk,1),:)+newk(npidx(jk,2),:))/2
                       VOE(eidx1,1,1) = nev + 1
                       VOE(eidx1,nev+2,1) = eidx2
                       VOE(eidx1,nev+2,2) = kidx(jk)

                       nvoe = nvoe + 1

                       VOEidx(nvoe,:) = (/ min(npidx(jk,1),npidx(jk,2)) ,max(npidx(jk,1),npidx(jk,2)), kidx(jk)/)
                    endif
                    if (jd.eq.0) then
                       newtetra(jt,4+jk) = kidx(jk)
                    else
                       newtetra(tcount+jd,4+jk) = kidx(jk)
                    endif
                 enddo
              enddo
              tcount = tcount + 7
           endif
        endif
     enddo

     !homogenize mesh
     do jt=1,nt
        if (.not.trefine(jt)) then
           npidx(1,:) = (/ newtetra(jt,1), newtetra(jt,5)/)
           npidx(2,:) = (/ newtetra(jt,1), newtetra(jt,6)/)
           npidx(3,:) = (/ newtetra(jt,1), newtetra(jt,7)/)
           npidx(4,:) = (/ newtetra(jt,2), newtetra(jt,8)/)
           npidx(5,:) = (/ newtetra(jt,2), newtetra(jt,9)/)
           npidx(6,:) = (/ newtetra(jt,3), newtetra(jt,10)/)
           npidx(7,:) = (/ newtetra(jt,2), newtetra(jt,5)/)
           npidx(8,:) = (/ newtetra(jt,3), newtetra(jt,6)/)
           npidx(9,:) = (/ newtetra(jt,4), newtetra(jt,7)/)
           npidx(10,:) = (/ newtetra(jt,3), newtetra(jt,8)/)
           npidx(11,:) = (/ newtetra(jt,4), newtetra(jt,9)/)
           npidx(12,:) = (/ newtetra(jt,4), newtetra(jt,10)/)
           do jk=1,12
              locv = 0
              eidx1 = min(npidx(jk,1),npidx(jk,2))
              eidx2 = max(npidx(jk,1),npidx(jk,2))
              nev = VOE(eidx1,1,1)
              do jv=1,nev
                 if (VOE(eidx1,jv+1,1).eq.eidx2) then
                    locv = VOE(eidx1,jv+1,2)
                    exit
                 endif
              enddo
              if (locv.ne.0) then
                 eidx11 = min(npidx(jk,1),locv)
                 eidx22 = max(npidx(jk,1),locv)
                 nev = VOE(eidx11,1,1)
                 do jv=1,nev
                    if (VOE(eidx11,jv+1,1).eq.eidx22) then
                       trefine(jt) = .true.
                       tchange = .true.
                       exit
                    endif
                 enddo
                 eidx11 = min(npidx(jk,2),locv)
                 eidx22 = max(npidx(jk,2),locv)
                 nev = VOE(eidx11,1,1)
                 do jv=1,nev
                    if (VOE(eidx11,jv+1,1).eq.eidx22) then
                       trefine(jt) = .true.
                       tchange = .true.
                       exit
                    endif
                 enddo
              endif
           enddo
        else
           trefine(jt) = .false.
        endif
     enddo
  enddo tchange_loop
END SUBROUTINE refine_tetra


SUBROUTINE get_initialmesh(init_steps,nt,tetra,wtetra,tclass,nkfull,kfull,nvoe,VOE,VOEidx)
  use const,           only: DPk
  use iso_fortran_env, only: ERROR_UNIT
  use woptic_io,       only: unit_outref
  use util,            only: string
  use clio,            only: croak, carp

  implicit none

  integer init_steps,counter,j1,j2,j3,nkfull,kmax,maxt
  integer kfull((3*2**init_steps)**3,3),Ti(4),krot(4,3)
  integer nt,tetra(48*8**(init_steps-1),10),tclass(48*8**(init_steps-1)),newpoints(6,3),kidx(6)
  integer newtetra(48*8**(init_steps-1),10),newtclass(48*8**(init_steps-1))
  integer VOE((3*2**init_steps)**3,300,2),VOEidx(6*48*8**(init_steps),3)
  real(DPk) :: S(3,3,48),wtetra(48*8**(init_steps-1)),newwtetra(48*8**(init_steps-1))
  integer eidx1,eidx2,jd,jk,jr,js,jt,nnewk,nvoe,tcount,find_kpoint2,npidx(6,2)
  integer :: jsign, jperm
  logical :: trefine(48*8**(init_steps-1))

  integer, parameter :: axes(3,3) = reshape( (/ &
       & 0, 0, 1, &
       & 0, 1, 0, &
       & 1, 0, 0  &
       /), (/3,3/) )
  integer :: perm(0:2) = (/ 1, 2, 3 /)


  kmax=(3*2**init_steps)**3
  maxt = 48*8**(init_steps-1)
  kfull = 0

  counter = 1 
  do j1=1,3
     do j2=1,3
        do j3=1,3
           kfull(counter,:) = (/j1-1,j2-1,j3-1/)  
           counter =  counter + 1
        enddo
     enddo
  enddo
  nkfull = counter - 1

  Ti = (/1,2,5,14/)

  do jperm = 0,5
     do jsign = 1,8
        do jd = 0,2
           S(:, jd+1, jsign + jperm*8) = (-1)**(jsign/2**jd) * axes(:,perm(jd))
        end do
     end do

     perm(mod((/ jperm, jperm+1 /), 3)) = perm(mod((/ jperm+1, jperm /), 3))
  end do

  tetra = 0
  do js=1,48
     krot = 0
     do jt=1,4
        do j1=1,3
           do j2=1,3
              krot(jt,j1) = krot(jt,j1) +int(S(j1,j2,js))*kfull(Ti(jt),j2)-int(S(j1,j2,js))
           enddo
           krot(jt,j1) = krot(jt,j1) + 1
        enddo
     enddo
     do jt=1,4
        do jk=1,nkfull
           if ((krot(jt,1).eq.kfull(jk,1)).and.(krot(jt,2).eq.kfull(jk,2))&
                .and.(krot(jt,3).eq.kfull(jk,3))) then
              tetra(js,jt) = jk
           endif
        enddo
     enddo
     tclass(js) = int(S(1,1,js)*S(2,2,js)*S(3,3,js) + S(1,2,js)*S(2,3,js)*S(3,1,js) &
          + S(1,3,js)*S(2,1,js)*S(3,2,js) - S(1,1,js)*S(2,3,js)*S(3,2,js) &
          - S(1,2,js)*S(2,1,js)*S(3,3,js) - S(1,3,js)*S(2,2,js)*S(3,1,js))
  enddo



  wtetra = 0.125d0/6d0
  nnewk = 0
  VOE = 0
  VOEidx = 0
  nvoe = 0
  npidx(1,:) = (/ 1, 2 /)
  npidx(2,:) = (/ 1, 3 /)
  npidx(3,:) = (/ 1, 4 /)
  npidx(4,:) = (/ 2, 3 /)
  npidx(5,:) = (/ 2, 4 /)
  npidx(6,:) = (/ 3, 4 /)
  nt = 48
  kfull = 2*kfull

  do jt=1,nt
     if ( tetra(jt,1)==tetra(jt,2) .or. &
          tetra(jt,2)==tetra(jt,3) .or. &
          tetra(jt,3)==tetra(jt,4)) then
        call carp("multiple vertices "//trim(string(jt)))
        write(ERROR_UNIT, *) tetra(jt,:)
        call croak()
     else
        newpoints(1,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,2),:))/2
        newpoints(2,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,3),:))/2
        newpoints(3,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,4),:))/2
        newpoints(4,:) =(kfull(tetra(jt,2),:)+kfull(tetra(jt,3),:))/2
        newpoints(5,:) =(kfull(tetra(jt,2),:)+kfull(tetra(jt,4),:))/2
        newpoints(6,:) =(kfull(tetra(jt,3),:)+kfull(tetra(jt,4),:))/2
     endif
     do jd=1,6
        kidx(jd)= find_kpoint2(nkfull+nnewk,kfull(1:nkfull+nnewk,:),newpoints(jd,:))
        if (kidx(jd).le.0) then
           nnewk = nnewk + 1
           kidx(jd) = nkfull+nnewk
           kfull(kidx(jd),:) = newpoints(jd,:)
           eidx1 = min(tetra(jt,npidx(jd,1)),tetra(jt,npidx(jd,2)))
           eidx2 = max(tetra(jt,npidx(jd,1)),tetra(jt,npidx(jd,2)))
           VOE(eidx1,1,1) = VOE(eidx1,1,1) + 1
           VOE(eidx1,VOE(eidx1,1,1)+2,1) = eidx2
           VOE(eidx1,VOE(eidx1,1,1)+2,2) = kidx(jd)
           nvoe = nvoe + 1
           VOEidx(nvoe,:) = (/eidx1  ,eidx2 , kidx(jd) /)
        endif
        tetra(jt,4+jd) = kidx(jd)
     enddo
  enddo



  nkfull = nkfull + nnewk
  nnewk = 0

  do jr=1,init_steps-1 
     trefine = .true.
     call refine_tetra(nt,tetra(1:nt,:),tclass(1:nt),wtetra(1:nt),nkfull,kmax,kfull,&
          nvoe,VOE,VOEidx(1:nvoe + 48*nt,:),trefine,maxt,&
          newtetra,newwtetra,newtclass,nnewk,tcount) 
     nkfull = nkfull + nnewk
     nt = tcount
     tetra = newtetra
     wtetra = newwtetra
     tclass = newtclass 
  enddo

  write(unit_outref,*)">>> initial tetrahedral mesh"
  write(unit_outref,*)"  number of tetrahedra:",nt
  write(unit_outref,*)"  number of full mesh k points:",nkfull
END SUBROUTINE get_initialmesh

SUBROUTINE compute_shapeparameters(nk,nt,k,tetra,ndim)
  use const,     only: DPk
  use woptic_io, only: unit_outref

  implicit none

  integer   :: nk,nt,jt,ndim(3),counter
  integer   :: tetra(nt,10),k(nk,3)
  real(DPk) :: vol,crad,qual(nt),a(3),b(3),c(3),quals(nt)
  real(DPk) :: crossab(3),crossca(3),crossbc(3),tmp(3)

  do jt=1,nt
     a = dble(k(tetra(jt,2),:) - k(tetra(jt,1),:))/ndim(1)
     b = dble(k(tetra(jt,3),:) - k(tetra(jt,1),:))/ndim(1)
     c = dble(k(tetra(jt,4),:) - k(tetra(jt,1),:))/ndim(1)

     crossab(1) = a(2)*b(3) - a(3)*b(2)
     crossab(2) = a(3)*b(1) - a(1)*b(3)
     crossab(3) = a(1)*b(2) - a(2)*b(1)
     crossca(1) = c(2)*a(3) - c(3)*a(2)
     crossca(2) = c(3)*a(1) - c(1)*a(3)
     crossca(3) = c(1)*a(2) - c(2)*a(1)
     crossbc(1) = b(2)*c(3) - b(3)*c(2)
     crossbc(2) = b(3)*c(1) - b(1)*c(3)
     crossbc(3) = b(1)*c(2) - b(2)*c(1)
     vol = dabs(dot_product(a,crossbc))/6d0
     tmp = dot_product(a,a)*crossbc+dot_product(b,b)*crossca &
          +dot_product(c,c)*crossab
     crad = dsqrt(tmp(1)**2+tmp(2)**2+tmp(3)**2)/12d0/vol
     qual(jt) = crad**3/vol
  enddo

  quals = 0d0
  counter = 0
  do jt=1,nt
     if (.not.any(quals(1:counter).eq.qual(jt))) then
        counter = counter + 1
        quals(counter) = qual(jt)
     endif
  enddo

  write(unit_outref,*)"  there are ",counter," different shapes"
  write(unit_outref,*)"  qualities ",quals(1:counter)
END SUBROUTINE compute_shapeparameters


!! Time-stamp: <2016-02-15 18:35:46 assman@faepop36.tu-graz.ac.at>
