!!! woptic/src/inwopcheck.f90
!!!
!!!    Parses ‘inwop‘ file, checks for consistency, and outputs
!!!    variables needed by ‘woptic’
!!!
!!! Copyright 2014-2015 Elias Assmann
!!!

program check_inwop
  use woptic,    only: inwop_t, inwop_read
  use woptic_io, only: set_casename, scratch,                             &
       fn_klist, fn_fklist, fn_voe, fn_map, fn_contr, fn_doscontr, fn_K1, &
       fn_energy, fn_mommat, fn_ham, fn_vk, fn_vvk, fn_optcond, fn_wdos,  &
       fn_outref, fn_hr, fn_vr, fn_vvr, fn_ksym, fn_inop, fn_fermi,       &
       fn_outwop, fn_kadd, fn_hist, fn_klist, fn_tet, fn_ftet,            &
       suf_new, suf_sym, suf_old, suf_rfd, suf_jnd
  use clio,      only: argstr, fetcharg, croak
  use util,      only: basename

  type(inwop_t) :: inw
  type(argstr)  :: fname, arg
  integer       :: i, jarg
  logical       :: runlapw1, runoptic, convert_vr, mixed_vr, peierls, need_ham
  character(2)  :: updn = ''
  logical       :: have_band=.false., have_so=.false., have_file=.false.

  call fetcharg(1, fname)

  call inwop_read(fname, inw)

  print '("matelmode=", I0)', inw%matelmode

  runlapw1   = inw%read_mommat .or. inw%read_energy
  runoptic   = inw%read_mommat
  convert_vr = inw%read_vk
  mixed_vr   = inw%read_vvk
  peierls    = inw%peierls
  need_ham   = inw%read_ham

  if (inw%orig_umatrix) then
     runlapw1 = .false.
     runoptic = .false.
  end if

#define prlog(name) print '(A, "=", A)', #name, merge("yes", "   ", ##name)
  prlog(runlapw1  )
  prlog(runoptic  )
  prlog(convert_vr)
  prlog(mixed_vr  )
  prlog(peierls   )
  prlog(need_ham  )

  print_filenames: if (command_argument_count() > 1) then
     jarg=2
     arguments_rest: do while (jarg <= command_argument_count())
        call fetcharg(jarg, arg)

        select case (arg%s)
        case ('-band', '--band')
           have_band = .true.
        case ('-up', '--up')
           updn='up'
        case ('-dn', '--dn')
           updn='dn'
        case ('-so', '--so')
           have_so = .true.
        case default
           if (arg%s(1:1) == '-') &
                call croak('unknown option ' // trim(arg%s))
           call fetcharg(jarg, fname)
           have_file=.true.
        end select

        jarg = jarg+1
     end do arguments_rest

     if (.not. have_file) &
          call croak('casename not supplied')

     call set_casename(fname, UPDN=updn, HAVE_SO=have_so, HAVE_BAND=have_band)

#define prvar(name)    print '(A, "=", A)',  'fn_'//#name, trim( fn_##name)
#define prsuf(name)    print '(A, "=", A)', 'suf_'//#name, trim(suf_##name)
#define prvas(name, n) print '(A, "=", ''"'', '//#n//'(A, " "), ''"'')',  \
     'fn_'//#name, ( trim( fn_##name(i)), i=1,##n )
#define prbas(name)    print '(A)',   trim(basename(fn_##name))
#define prbss(name, n) print '(A)', ( trim(basename(fn_##name(i))), i=1,##n )

     prvar(hist  );   prvar(contr   )
     prvar(klist );   prvar(doscontr)
     prvar(tet   );   prvar(K1      )
     prvar(ftet  );   prvar(optcond )
     prvar(klist );   prvar(wdos    )
     prvar(kadd  );   prvar(outref  )
     prvar(fklist);   prvar(fermi   )
     prvar(voe   );   prvar(ksym    )
     prvar(map   );   prvar(outwop  )

     if (runlapw1  )  prvar(energy  ); prsuf(new)
     if (runoptic  )  prvar(mommat  ); prsuf(sym)
     if (runoptic  )  prvar(inop    ); prsuf(old)
     if (need_ham  )  prvar(ham     ); prsuf(jnd)
     if (need_ham  )  prvar(hr      ); prsuf(rfd)
     if (convert_vr)  prvas(vk,  3  )
     if (convert_vr)  prvar(vr      )
     if (mixed_vr  )  prvas(vvk, 6  )
     if (mixed_vr  )  prvar(vvr     )

     print '("SCRATCH=", A)', trim(scratch%s)

     print '(A, ''="'')', 'scratchfiles'
     prbas(ftet)
     prbas(tet)
     prbas(map)
     prbas(voe)
     prbas(doscontr)
     prbas(contr)
     prbas(K1)
     print '(''"'')'
  end if print_filenames

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

!! Time-stamp: <2016-02-23 19:15:51 assman@faepop36.tu-graz.ac.at>
