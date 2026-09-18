MODULE mosart_lake_snowice_emulator_mod
!-------------------------------------------------------------------------------------------
! Online snow/ice emulator for MOSART-Lake (one-way).
!
! Runs the frozen causal LSTM (T23: single layer, unidirectional, hidden 96, input 15,
! output 4) inside MOSART, driven by an EXTERNAL monthly forcing file whose path is a
! namelist setting. Nothing about the emulator depends on where that forcing came from,
! so a user can supply their own as long as it follows the documented schema.
!
! The forcing file is read in object space (49,708 objects x 10 features per month,
! about 4 MB) on every task through PIO and indexed locally: an object is exactly a
! (routing cell, main-or-tributary) pair, which is already how MOSART is decomposed, so
! no object-space domain decomposition is needed. Each object's LSTM is independent of
! every other one, so results are invariant to the MPI layout by construction.
!
! STATE AND RESTART.  The LSTM hidden/cell state is NOT written to the restart file.
! Instead it is rebuilt by replaying `lake_snowice_spinup_months` months of forcing from
! zero before the first month the model needs. The packaged cold-start diagnostic shows
! the state is effectively converged after 12-24 months (mean |df_ice| 0.043 at 0 months,
! 0.0022 at 12, 0.0013 at 24), so a 24-month replay costs a few seconds at init and buys
! freedom from a 76 MB restart field that RtmIO could not write anyway (it has no
! decomposed 2-D path). The price is that a split run is not bit-for-bit with a
! continuous one; the difference is bounded by the numbers above and is reported in the log.
!
! TIME SEMANTICS.  Monthly values are mid-month anchors interpolated linearly in time,
! identical to the prescribed-driver path, so that running the same case with
! lake_snowice_source='driver' and ='emulator' differs only in where the four numbers
! came from. That requires holding the bracketing pair of months, i.e. the LSTM is
! stepped one month ahead of the model calendar. That is legitimate here: the forcing is
! prescribed and fully known in advance, exactly as it is for the offline driver file.
!-------------------------------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use RtmFileUtils  , only : getfil
  use RtmIO         , only : ncd_pio_openfile, ncd_pio_closefile, ncd_io, file_desc_t, var_desc_t
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog, lake_snowice_model_file, lake_snowice_forcing_file, &
                             lake_snowice_spinup_months
  use RunoffMod     , only : rtmCTL, TLake_r, TLake_t, TUnit_lake_r, TUnit_lake_t
  use pio
  use mosart_lake_snowice_kernel_mod, only : NDYN, NSTA, NIN, NHID, NTGT, NGAT, &
                                             snowice_lstm_step, snowice_head_post

  implicit none
  private

  public :: lake_snowice_emulator_init
  public :: lake_snowice_emulator_bracket   ! fill the four fields for a bracketing month pair
  public :: lake_snowice_emulator_range     ! first/last month the forcing file covers

  ! --- frozen model, replicated on every task (about 0.35 MB) ---
  real(r8), allocatable, save :: w_ih(:,:)      ! (NIN ,NGAT)  transposed by NetCDF ordering
  real(r8), allocatable, save :: w_hh(:,:)      ! (NHID,NGAT)
  real(r8), allocatable, save :: b_lstm(:)      ! (NGAT)   already bias_ih + bias_hh
  real(r8), allocatable, save :: w_head(:,:)    ! (NHID,NTGT)
  real(r8), allocatable, save :: b_head(:)      ! (NTGT)
  real(r8), allocatable, save :: dyn_mean(:), dyn_std(:)   ! (NDYN)

  ! --- object bookkeeping ---
  integer , save              :: nobj_glob = 0
  integer , allocatable, save :: gobj_main(:)   ! (numr) global object index, 0 = none
  integer , allocatable, save :: gobj_trib(:)
  real(r8), allocatable, save :: stat_r(:,:)    ! (NSTA, begr:endr) normalized statics
  real(r8), allocatable, save :: stat_t(:,:)

  ! --- LSTM state, carried on the local decomposition; never restarted ---
  real(r8), allocatable, save :: h_r(:,:), c_r(:,:)   ! (NHID, begr:endr)
  real(r8), allocatable, save :: h_t(:,:), c_t(:,:)

  ! --- cached monthly predictions for the bracketing pair ---
  real(r8), allocatable, public, save :: p_r0(:,:), p_r1(:,:)  ! (NTGT, begr:endr)
  real(r8), allocatable, public, save :: p_t0(:,:), p_t1(:,:)
  integer , save :: ym_p0 = -999999, ym_p1 = -999999
  integer , save :: ym_state = -999999          ! month the LSTM state has been stepped through

  ! --- forcing file ---
  type(file_desc_t), save :: fncid
  type(var_desc_t) , save :: fvdesc
  integer, allocatable, save :: forc_yyyymm(:)
  integer, save :: nt_forc = 0
  logical, save :: initialized = .false.

contains

!-------------------------------------------------------------------------------------------
  subroutine lake_snowice_emulator_init()

    character(len=256) :: locfn
    integer :: ier, dimid, n, g, ncell, nloc_r, nloc_t, spin
    integer :: nobj_f, ndyn_f, nhid_f, nin_f, ntgt_f
    logical :: found
    type(file_desc_t) :: mncid
    integer , allocatable :: obj_cell0(:), obj_is_main(:), tmp_cell0(:), tmp_is_main(:)
    real(r8), allocatable :: stat_all(:,:)
    character(len=*), parameter :: subname = '(lake_snowice_emulator_init)'

    if (initialized) return

    if (len_trim(lake_snowice_model_file) == 0) &
         call shr_sys_abort(subname//' lake_snowice_model_file is empty')
    if (len_trim(lake_snowice_forcing_file) == 0) &
         call shr_sys_abort(subname//' lake_snowice_forcing_file is empty')

    !------------------------------------------------------------------ model package
    call getfil(lake_snowice_model_file, locfn, 0)
    if (masterproc) then
       write(iulog,*) trim(subname),' reading emulator model package: ',trim(locfn)
       call shr_sys_flush(iulog)
    end if
    call ncd_pio_openfile(mncid, trim(locfn), 0)

    call get_dimlen(mncid, 'object', nobj_glob, subname)
    call get_dimlen(mncid, 'input' , nin_f    , subname)
    call get_dimlen(mncid, 'hidden', nhid_f   , subname)
    call get_dimlen(mncid, 'target', ntgt_f   , subname)
    call get_dimlen(mncid, 'dynamic', ndyn_f  , subname)
    if (nin_f /= NIN .or. nhid_f /= NHID .or. ntgt_f /= NTGT .or. ndyn_f /= NDYN) then
       write(iulog,*) trim(subname),' model shape mismatch: input/hidden/target/dynamic = ', &
                      nin_f, nhid_f, ntgt_f, ndyn_f,' expected ',NIN,NHID,NTGT,NDYN
       call shr_sys_abort(subname//' emulator model package has an unexpected shape')
    end if

    allocate(w_ih(NIN,NGAT), w_hh(NHID,NGAT), b_lstm(NGAT), &
             w_head(NHID,NTGT), b_head(NTGT), dyn_mean(NDYN), dyn_std(NDYN))
    call rd_r2(mncid,'lstm_weight_ih', w_ih  , subname)
    call rd_r2(mncid,'lstm_weight_hh', w_hh  , subname)
    call rd_r1(mncid,'lstm_bias'     , b_lstm, subname)
    call rd_r2(mncid,'head_weight'   , w_head, subname)
    call rd_r1(mncid,'head_bias'     , b_head, subname)
    call rd_r1(mncid,'dynamic_mean'  , dyn_mean, subname)
    call rd_r1(mncid,'dynamic_std'   , dyn_std , subname)
    if (any(dyn_std <= 0._r8)) call shr_sys_abort(subname//' non-positive dynamic_std')
    if (.not. all(is_finite2(w_ih)) .or. .not. all(is_finite2(w_hh)) .or. &
        .not. all(is_finite1(b_lstm))) &
         call shr_sys_abort(subname//' non-finite value in the LSTM weights')

    allocate(obj_cell0(nobj_glob), obj_is_main(nobj_glob))
    call rd_i1(mncid,'obj_cell0'  , obj_cell0  , subname)
    call rd_i1(mncid,'obj_is_main', obj_is_main, subname)

    allocate(stat_all(NSTA,nobj_glob))
    call rd_r2(mncid,'static_input_normalized', stat_all, subname)

    call ncd_pio_closefile(mncid)

    !------------------------------------------------------------------ forcing file
    call getfil(lake_snowice_forcing_file, locfn, 0)
    if (masterproc) then
       write(iulog,*) trim(subname),' reading emulator forcing: ',trim(locfn)
       call shr_sys_flush(iulog)
    end if
    call ncd_pio_openfile(fncid, trim(locfn), 0)

    call get_dimlen(fncid, 'time'  , nt_forc, subname)
    call get_dimlen(fncid, 'object', nobj_f , subname)
    call get_dimlen(fncid, 'feature', ndyn_f, subname)
    if (nobj_f /= nobj_glob .or. ndyn_f /= NDYN) then
       write(iulog,*) trim(subname),' forcing has ',nobj_f,' objects and ',ndyn_f, &
                      ' features; model package has ',nobj_glob,' and ',NDYN
       call shr_sys_abort(subname//' forcing and model package disagree on shape')
    end if

    allocate(forc_yyyymm(nt_forc))
    call rd_i1(fncid,'yyyymm', forc_yyyymm, subname)
    do n = 2, nt_forc
       if (month_add(forc_yyyymm(n-1),1) /= forc_yyyymm(n)) then
          write(iulog,*) trim(subname),' forcing yyyymm is not contiguous at index ',n, &
                         ': ',forc_yyyymm(n-1),' -> ',forc_yyyymm(n)
          call shr_sys_abort(subname//' emulator forcing has a gap in its month axis')
       end if
    end do

    ! The object axis must agree item for item, otherwise every lake gets another lake's
    ! weather. This is the one cross-file check that can never be relaxed.
    allocate(tmp_cell0(nobj_glob), tmp_is_main(nobj_glob))
    call rd_i1(fncid,'obj_cell0'  , tmp_cell0  , subname)
    call rd_i1(fncid,'obj_is_main', tmp_is_main, subname)
    if (any(tmp_cell0 /= obj_cell0) .or. any(tmp_is_main /= obj_is_main)) &
         call shr_sys_abort(subname//' forcing and model package have different object mappings')
    deallocate(tmp_cell0, tmp_is_main)

    ier = pio_inq_varid(fncid, 'forcing', fvdesc)
    if (ier /= PIO_NOERR) call shr_sys_abort(subname//' forcing file has no "forcing" variable')

    !------------------------------------------------------------------ cell -> object lookup
    allocate(gobj_main(rtmCTL%numr), gobj_trib(rtmCTL%numr))
    gobj_main = 0 ; gobj_trib = 0
    do g = 1, nobj_glob
       ncell = obj_cell0(g) + 1                 ! file is 0-based
       if (ncell < 1 .or. ncell > rtmCTL%numr) then
          write(iulog,*) trim(subname),' object ',g,' has cell ',ncell,' outside 1..',rtmCTL%numr
          call shr_sys_abort(subname//' object cell index is outside the routing grid')
       end if
       if (obj_is_main(g) /= 0) then
          if (gobj_main(ncell) /= 0) call shr_sys_abort(subname//' duplicate main object on a cell')
          gobj_main(ncell) = g
       else
          if (gobj_trib(ncell) /= 0) call shr_sys_abort(subname//' duplicate tributary object on a cell')
          gobj_trib(ncell) = g
       end if
    end do

    !------------------------------------------------------------------ local statics and state
    allocate(stat_r(NSTA,rtmCTL%begr:rtmCTL%endr), stat_t(NSTA,rtmCTL%begr:rtmCTL%endr))
    allocate(h_r(NHID,rtmCTL%begr:rtmCTL%endr), c_r(NHID,rtmCTL%begr:rtmCTL%endr))
    allocate(h_t(NHID,rtmCTL%begr:rtmCTL%endr), c_t(NHID,rtmCTL%begr:rtmCTL%endr))
    allocate(p_r0(NTGT,rtmCTL%begr:rtmCTL%endr), p_r1(NTGT,rtmCTL%begr:rtmCTL%endr))
    allocate(p_t0(NTGT,rtmCTL%begr:rtmCTL%endr), p_t1(NTGT,rtmCTL%begr:rtmCTL%endr))
    stat_r = 0._r8 ; stat_t = 0._r8
    p_r0 = 0._r8 ; p_r1 = 0._r8 ; p_t0 = 0._r8 ; p_t1 = 0._r8

    nloc_r = 0 ; nloc_t = 0
    do n = rtmCTL%begr, rtmCTL%endr
       g = gobj_main(rtmCTL%gindex(n))
       if (g > 0) then
          stat_r(:,n) = stat_all(:,g) ; nloc_r = nloc_r + 1
       end if
       g = gobj_trib(rtmCTL%gindex(n))
       if (g > 0) then
          stat_t(:,n) = stat_all(:,g) ; nloc_t = nloc_t + 1
       end if
    end do
    deallocate(stat_all, obj_cell0, obj_is_main)

    spin = max(0, lake_snowice_spinup_months)
    initialized = .true.

    if (masterproc) then
       write(iulog,*) trim(subname),' emulator ready: ',nobj_glob,' global objects, forcing ', &
                      forc_yyyymm(1),' .. ',forc_yyyymm(nt_forc),' (',nt_forc,' months)'
       write(iulog,*) trim(subname),' spin-up replay = ',spin,' months; LSTM state is NOT ', &
                      'carried in the restart file, so a split run is not bit-for-bit with ', &
                      'a continuous one (bounded by the packaged cold-start diagnostic)'
       call shr_sys_flush(iulog)
    end if
    write(iulog,*) trim(subname),' local objects on this task: main ',nloc_r,' tributary ',nloc_t

  end subroutine lake_snowice_emulator_init

!-------------------------------------------------------------------------------------------
  subroutine lake_snowice_emulator_bracket(ym0, ym1)
    ! Make predictions for months ym0 and ym1 (consecutive) available in p_*0 / p_*1 and
    ! copy them into TLake_%. The caller does the time interpolation.

    integer, intent(in) :: ym0, ym1

    integer :: spin, ymstart, ym
    character(len=*), parameter :: subname = '(lake_snowice_emulator_bracket)'

    if (.not. initialized) call lake_snowice_emulator_init()

    if (ym0 == ym_p0 .and. ym1 == ym_p1) return

    if (ym0 < forc_yyyymm(1) .or. ym1 > forc_yyyymm(nt_forc)) then
       write(iulog,*) trim(subname),' need months ',ym0,'..',ym1,' but forcing covers ', &
                      forc_yyyymm(1),'..',forc_yyyymm(nt_forc)
       call shr_sys_abort(subname//' emulator forcing does not cover the requested months')
    end if

    ! Advance by one month if that is all that is needed, otherwise rebuild from a replay.
    if (ym_p1 == ym0 .and. ym_state == ym0) then
       p_r0 = p_r1 ; p_t0 = p_t1
       call step_month(ym1)
       call harvest(p_r1, p_t1)
    else
       spin = max(0, lake_snowice_spinup_months)
       ymstart = month_add(ym0, -spin)
       if (ymstart < forc_yyyymm(1)) ymstart = forc_yyyymm(1)
       if (masterproc) then
          write(iulog,*) trim(subname),' replaying LSTM from ',ymstart,' to ',ym1, &
                         ' to reach bracket ',ym0,' -> ',ym1
          call shr_sys_flush(iulog)
       end if
       h_r = 0._r8 ; c_r = 0._r8 ; h_t = 0._r8 ; c_t = 0._r8
       ym = ymstart
       do
          call step_month(ym)
          if (ym == ym0) call harvest(p_r0, p_t0)
          if (ym == ym1) then
             call harvest(p_r1, p_t1)
             exit
          end if
          ym = month_add(ym, 1)
       end do
    end if

    ym_p0 = ym0 ; ym_p1 = ym1

  end subroutine lake_snowice_emulator_bracket

!-------------------------------------------------------------------------------------------
  subroutine lake_snowice_emulator_range(ym_first, ym_last)
    integer, intent(out) :: ym_first, ym_last
    if (.not. initialized) call lake_snowice_emulator_init()
    ym_first = forc_yyyymm(1)
    ym_last  = forc_yyyymm(nt_forc)
  end subroutine lake_snowice_emulator_range

!-------------------------------------------------------------------------------------------
  subroutine step_month(yyyymm)
    ! One LSTM step for every local object, using the forcing of calendar month yyyymm.

    integer, intent(in) :: yyyymm

    real(r8), allocatable :: fbuf(:,:)        ! (NDYN, nobj_glob) one month, all objects
    integer  :: start(3), count(3), ier, tidx, n, g
    real(r8) :: x(NIN)
    character(len=*), parameter :: subname = '(step_month)'

    tidx = time_index(yyyymm)
    if (tidx < 1) then
       write(iulog,*) trim(subname),' no forcing record for ',yyyymm
       call shr_sys_abort(subname//' missing month in the emulator forcing file')
    end if

    ! NetCDF declares forcing(time,object,feature); Fortran sees (feature,object,time).
    allocate(fbuf(NDYN,nobj_glob))
    start = (/ 1, 1, tidx /)
    count = (/ NDYN, nobj_glob, 1 /)
    ier = pio_get_var(fncid, fvdesc, start, count, fbuf)
    if (ier /= PIO_NOERR) call shr_sys_abort(subname//' failed to read the forcing slice')

    do n = rtmCTL%begr, rtmCTL%endr
       g = gobj_main(rtmCTL%gindex(n))
       if (g > 0) then
          x(1:NDYN)     = (fbuf(:,g) - dyn_mean) / dyn_std
          x(NDYN+1:NIN) = stat_r(:,n)
          call snowice_lstm_step(w_ih, w_hh, b_lstm, x, h_r(:,n), c_r(:,n))
       end if
       g = gobj_trib(rtmCTL%gindex(n))
       if (g > 0) then
          x(1:NDYN)     = (fbuf(:,g) - dyn_mean) / dyn_std
          x(NDYN+1:NIN) = stat_t(:,n)
          call snowice_lstm_step(w_ih, w_hh, b_lstm, x, h_t(:,n), c_t(:,n))
       end if
    end do

    deallocate(fbuf)
    ym_state = yyyymm

  end subroutine step_month

!-------------------------------------------------------------------------------------------
  subroutine harvest(pr, pt)
    ! Read the four targets out of the current hidden state for every local object.

    real(r8), intent(inout) :: pr(:,rtmCTL%begr:), pt(:,rtmCTL%begr:)

    integer :: n

    do n = rtmCTL%begr, rtmCTL%endr
       if (gobj_main(rtmCTL%gindex(n)) > 0) then
          call snowice_head_post(w_head, b_head, h_r(:,n), pr(:,n))
       else
          pr(:,n) = 0._r8
       end if
       if (gobj_trib(rtmCTL%gindex(n)) > 0) then
          call snowice_head_post(w_head, b_head, h_t(:,n), pt(:,n))
       else
          pt(:,n) = 0._r8
       end if
    end do

  end subroutine harvest

!-------------------------------------------------------------------------------------------
  integer function time_index(yyyymm) result(tidx)
    integer, intent(in) :: yyyymm
    integer :: i
    tidx = -1
    do i = 1, nt_forc
       if (forc_yyyymm(i) == yyyymm) then ; tidx = i ; return ; end if
    end do
  end function time_index

  integer function month_add(yyyymm, k) result(ym)
    integer, intent(in) :: yyyymm, k
    integer :: y, m, tot
    y = yyyymm/100 ; m = mod(yyyymm,100)
    tot = y*12 + (m-1) + k
    ym  = (tot/12)*100 + mod(tot,12) + 1
  end function month_add

!-------------------------------------------------------------------------------------------
  subroutine get_dimlen(ncid, name, n, caller)
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: name, caller
    integer          , intent(out)   :: n
    integer :: ier, dimid
    ier = pio_inq_dimid(ncid, trim(name), dimid)
    if (ier /= PIO_NOERR) call shr_sys_abort(trim(caller)//' missing dimension '//trim(name))
    ier = pio_inq_dimlen(ncid, dimid, n)
    if (ier /= PIO_NOERR .or. n <= 0) &
         call shr_sys_abort(trim(caller)//' bad dimension '//trim(name))
  end subroutine get_dimlen

  subroutine rd_r1(ncid, name, a, caller)
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: name, caller
    real(r8)         , intent(inout) :: a(:)
    logical :: found
    call ncd_io(ncid=ncid, varname=trim(name), flag='read', data=a, readvar=found)
    if (.not. found) call shr_sys_abort(trim(caller)//' missing variable '//trim(name))
  end subroutine rd_r1

  subroutine rd_r2(ncid, name, a, caller)
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: name, caller
    real(r8)         , intent(inout) :: a(:,:)
    logical :: found
    call ncd_io(ncid=ncid, varname=trim(name), flag='read', data=a, readvar=found)
    if (.not. found) call shr_sys_abort(trim(caller)//' missing variable '//trim(name))
  end subroutine rd_r2

  subroutine rd_i1(ncid, name, a, caller)
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: name, caller
    integer          , intent(inout) :: a(:)
    logical :: found
    call ncd_io(ncid=ncid, varname=trim(name), flag='read', data=a, readvar=found)
    if (.not. found) call shr_sys_abort(trim(caller)//' missing variable '//trim(name))
  end subroutine rd_i1

  pure function is_finite1(a) result(m)
    real(r8), intent(in) :: a(:)
    logical :: m(size(a))
    m = (a == a) .and. (abs(a) < huge(1._r8))
  end function is_finite1

  pure function is_finite2(a) result(m)
    real(r8), intent(in) :: a(:,:)
    logical :: m(size(a,1),size(a,2))
    m = (a == a) .and. (abs(a) < huge(1._r8))
  end function is_finite2

END MODULE mosart_lake_snowice_emulator_mod
