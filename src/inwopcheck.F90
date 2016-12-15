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

  character(*), parameter :: fmt_var = '(A, "=", A)'

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

  print '(A, "=", A)', 'runlapw1'  , merge("yes", "   ", runlapw1  )
  print '(A, "=", A)', 'runoptic'  , merge("yes", "   ", runoptic  )
  print '(A, "=", A)', 'convert_vr', merge("yes", "   ", convert_vr)
  print '(A, "=", A)', 'mixed_vr'  , merge("yes", "   ", mixed_vr  )
  print '(A, "=", A)', 'peierls'   , merge("yes", "   ", peierls   )
  print '(A, "=", A)', 'need_ham'  , merge("yes", "   ", need_ham  )

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

     if (runlapw1  ) print fmt_var, 'fn_energy',   trim(fn_energy  )
     if (runoptic  ) print fmt_var, 'fn_mommat',   trim(fn_mommat  )
     if (runoptic  ) print fmt_var, 'fn_inop'  ,   trim(fn_inop    )
     if (need_ham  ) print fmt_var, 'fn_ham'   ,   trim(fn_ham     )
     if (need_ham  ) print fmt_var, 'fn_hr'    ,   trim(fn_hr      )
     if (convert_vr) print fmt_mul(3), 'fn_vk' , ( trim(fn_vk (i)  ), i=1,3 )
     if (convert_vr) print fmt_var, 'fn_vr'    ,   trim(fn_vr      )
     if (mixed_vr  ) print fmt_mul(6), 'fn_vvk', ( trim(fn_vvk(i)  ), i=1,6 )
     if (mixed_vr  ) print fmt_var, 'fn_vvr'   ,   trim(fn_vvr     )

     print                 fmt_var, 'fn_hist'    , trim(fn_hist    )
     print                 fmt_var, 'fn_klist'   , trim(fn_klist   )
     print                 fmt_var, 'fn_tet'     , trim(fn_tet     )
     print                 fmt_var, 'fn_ftet'    , trim(fn_ftet    )
     print                 fmt_var, 'fn_klist'   , trim(fn_klist   )
     print                 fmt_var, 'fn_kadd'    , trim(fn_kadd    )
     print                 fmt_var, 'fn_fklist'  , trim(fn_fklist  )
     print                 fmt_var, 'fn_voe'     , trim(fn_voe     )
     print                 fmt_var, 'fn_map'     , trim(fn_map     )
     print                 fmt_var, 'fn_contr'   , trim(fn_contr   )
     print                 fmt_var, 'fn_doscontr', trim(fn_doscontr)
     print                 fmt_var, 'fn_K1'      , trim(fn_K1      )
     print                 fmt_var, 'fn_optcond' , trim(fn_optcond )
     print                 fmt_var, 'fn_wdos'    , trim(fn_wdos    )
     print                 fmt_var, 'fn_outref'  , trim(fn_outref  )
     print                 fmt_var, 'fn_fermi'   , trim(fn_fermi   )
     print                 fmt_var, 'fn_ksym'    , trim(fn_ksym    )
     print                 fmt_var, 'fn_outwop'  , trim(fn_outwop  )

     print                 fmt_var, 'suf_new'    , trim(suf_new    )
     print                 fmt_var, 'suf_sym'    , trim(suf_sym    )
     print                 fmt_var, 'suf_old'    , trim(suf_old    )
     print                 fmt_var, 'suf_jnd'    , trim(suf_jnd    )
     print                 fmt_var, 'suf_rfd'    , trim(suf_rfd    )

     print                 fmt_var, 'SCRATCH'    , trim(scratch%s  )

     print fmt_mul(7), 'scratchfiles', trim(basename(fn_K1 )), &
          trim(basename(fn_ftet    )), trim(basename(fn_tet)), &
          trim(basename(fn_contr   )), trim(basename(fn_voe)), &
          trim(basename(fn_doscontr)), trim(basename(fn_map))
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

contains
  pure character(len=30) function fmt_mul(N)
    use util, only: string
    implicit none

    integer, intent(in) :: N

    fmt_mul = '(A, "=''", '//trim(string(N))//'(/A), /"''")'
  end function fmt_mul
end program check_inwop

!! Time-stamp: <2016-09-14 14:33:19 assman@faepop71.tu-graz.ac.at>
