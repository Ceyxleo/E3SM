MODULE mosart_lake_snowice_mod
!-------------------------------------------------------------------------------------------
! Prescribed lake snow/ice for MOSART-Lake (one-way).
!
! The driver file holds monthly fields produced offline by the snow/ice emulator and
! repacked onto the routing grid (see perlmutter_script/repack_snowice_driver.py).
! MOSART never feeds anything back: these are boundary conditions, not prognostic state.
!
! Monthly values are treated as MID-MONTH anchors and interpolated linearly in time.
! The emulator targets are monthly means, so anchoring them at the end of a month would
! shift the whole seasonal cycle by about 15 days - larger than the emulator's own
! freeze-up error. Interpolating (rather than holding each month constant) also gives a
! continuous dh_ice/dt, which a later prognostic-latent-heat version needs.
!
! Nothing here is carried in the restart file: the bracketing months are recovered from
! the absolute calendar date on every startup.
!-------------------------------------------------------------------------------------------

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort, shr_sys_flush
  use RtmFileUtils  , only : getfil
  use RtmIO         , only : ncd_pio_openfile, ncd_io, file_desc_t
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog, lake_snowice_driver_file, lake_snowice_source
  use mosart_lake_snowice_emulator_mod, only : lake_snowice_emulator_init,    &
                                              lake_snowice_emulator_bracket, &
                                              lake_snowice_emulator_range,   &
                                              p_r0, p_r1, p_t0, p_t1
  use RunoffMod     , only : rtmCTL, TLake_r, TLake_t
  use RtmTimeManager, only : get_curr_date, get_calendar
  use pio           , only : pio_inq_dimid, pio_inq_dimlen, PIO_NOERR

  implicit none
  private

  public :: lake_snowice_init
  public :: lake_snowice_update
  public :: lake_ice_snow_bc

  ! Optical and thermal properties of the ice/snow cover. Compile-time constants: they are
  ! shared by the main-channel and tributary sides, so changing one changes both.
  real(r8), parameter :: ice_albedo    = 0.35_r8
  real(r8), parameter :: snow_albedo   = 0.75_r8
  real(r8), parameter :: ice_sw_trans  = 0.20_r8   ! shortwave transmitted through bare ice
  real(r8), parameter :: snow_sw_trans = 0.05_r8   ! shortwave transmitted through snow-covered ice
  real(r8), parameter :: k_ice         = 2.2_r8    ! W/m/K
  real(r8), parameter :: k_snow        = 0.3_r8    ! W/m/K
  real(r8), parameter :: emiss_ice     = 0.97_r8
  real(r8), parameter :: stefan_b      = 5.67e-8_r8
  real(r8), parameter :: min_resist    = 1.0e-3_r8 ! numerical floor only; the physical
                                                   ! atmospheric resistance is computed below
  real(r8), parameter :: snowice_tol   = 1.0e-10_r8
  real(r8), parameter :: TINYVALUE     = 1.0e-14_r8

  ! One calendar month of the four prescribed fields, on the local decomposition.
  type snowice_buf
     real(r8), pointer :: fice(:) => null()   ! ice areal fraction (-)
     real(r8), pointer :: hice(:) => null()   ! ice thickness (m)
     real(r8), pointer :: fsnw(:) => null()   ! snow-covered fraction (-)
     real(r8), pointer :: dsnw(:) => null()   ! snow depth (m)
  end type snowice_buf

  type(file_desc_t), save :: ncid
  integer, allocatable, save :: driver_yyyymm(:)
  integer, save :: ntime_driver = 0
  type(snowice_buf), save :: r0, r1, t0, t1     ! bracketing months, main / tributary
  integer, save :: cur_ym0 = -999999            ! yyyymm currently held in r0/t0
  integer, save :: cur_ym1 = -999999            ! yyyymm currently held in r1/t1
  logical, save :: initialized = .false.

contains

!-------------------------------------------------------------------------------------------
  subroutine lake_snowice_init()

    character(len=256) :: locfn
    integer :: dimid, ier
    logical :: found
    character(len=*), parameter :: subname = '(lake_snowice_init)'

    if (initialized) return

    ! The online emulator owns its own files; the prescribed driver is not opened at all.
    if (trim(lake_snowice_source) == 'emulator') then
       call lake_snowice_emulator_init()
       initialized = .true.
       return
    end if

    if (len_trim(lake_snowice_driver_file) == 0) then
       call shr_sys_abort(subname//' lake_snowice_driver_file is empty')
    end if

    call getfil(lake_snowice_driver_file, locfn, 0)
    if (masterproc) then
       write(iulog,*) trim(subname),' reading lake snow/ice driver: ',trim(locfn)
       call shr_sys_flush(iulog)
    end if
    call ncd_pio_openfile(ncid, trim(locfn), 0)

    ier = pio_inq_dimid(ncid, 'time', dimid)
    if (ier /= PIO_NOERR) call shr_sys_abort(subname//' driver has no time dimension')
    ier = pio_inq_dimlen(ncid, dimid, ntime_driver)
    if (ier /= PIO_NOERR .or. ntime_driver <= 0) &
         call shr_sys_abort(subname//' driver has a bad time dimension')

    allocate(driver_yyyymm(ntime_driver))
    call ncd_io(ncid=ncid, varname='yyyymm', flag='read', data=driver_yyyymm, readvar=found)
    if (.not. found) call shr_sys_abort(subname//' driver has no yyyymm variable')

    call alloc_buf(r0); call alloc_buf(r1)
    call alloc_buf(t0); call alloc_buf(t1)

    initialized = .true.
    if (masterproc) then
       write(iulog,*) trim(subname),' driver covers ',driver_yyyymm(1),' .. ', &
                      driver_yyyymm(ntime_driver),' (',ntime_driver,' months)'
       call shr_sys_flush(iulog)
    end if

  end subroutine lake_snowice_init

!-------------------------------------------------------------------------------------------
  subroutine lake_snowice_update()
    ! Refresh the interpolated snow/ice state for the current model date.
    ! Both sources deliver monthly values on the same mid-month anchors, so the time
    ! interpolation below is shared and a 'driver' run and an 'emulator' run differ only
    ! in where the four numbers came from.

    integer  :: yr, mon, day, tod
    integer  :: ym0, ym1, yr0, mo0, yr1, mo1, i0, i1, n
    integer  :: ym_first, ym_last
    real(r8) :: pos, half, span, w
    real(r8) :: dt_span
    logical  :: emul
    character(len=*), parameter :: subname = '(lake_snowice_update)'

    if (.not. initialized) call lake_snowice_init()

    emul = (trim(lake_snowice_source) == 'emulator')

    call get_curr_date(yr, mon, day, tod)

    ! Position inside the current month, in days, and the mid-month anchor.
    pos  = real(day-1,r8) + real(tod,r8)/86400._r8
    half = 0.5_r8*real(days_in_month(yr,mon),r8)

    if (pos < half) then
       call prev_month(yr, mon, yr0, mo0)
       yr1 = yr ; mo1 = mon
       span = 0.5_r8*real(days_in_month(yr0,mo0),r8) + half
       w    = (0.5_r8*real(days_in_month(yr0,mo0),r8) + pos) / span
    else
       yr0 = yr ; mo0 = mon
       call next_month(yr, mon, yr1, mo1)
       span = half + 0.5_r8*real(days_in_month(yr1,mo1),r8)
       w    = (pos - half) / span
    end if
    w = max(0._r8, min(1._r8, w))

    ym0 = yr0*100 + mo0
    ym1 = yr1*100 + mo1

    if (emul) then
       call lake_snowice_emulator_range(ym_first, ym_last)
    else
       ym_first = driver_yyyymm(1)
       ym_last  = driver_yyyymm(ntime_driver)
    end if

    ! Clamp to the end points: outside the covered range the nearest month is held.
    ! Half a month of zero-order hold at each end is acceptable; a case that runs wholly
    ! outside the covered period is not, and is caught here.
    if (ym1 < ym_first .or. ym0 > ym_last) then
       write(iulog,*) trim(subname),' model date ',yr*10000+mon*100+day, &
                      ' is outside the snow/ice range ',ym_first,' .. ',ym_last
       call shr_sys_abort(subname//' case period is not covered by the snow/ice source')
    end if
    if (ym0 < ym_first) then ; ym0 = ym_first ; w = 1._r8 ; end if
    if (ym1 > ym_last ) then ; ym1 = ym_last  ; w = 0._r8 ; end if

    if (ym0 /= cur_ym0 .or. ym1 /= cur_ym1) then
       if (emul) then
          call lake_snowice_emulator_bracket(ym0, ym1)
       else
          i0 = time_index(ym0)
          i1 = time_index(ym1)
          if (i0 < 1 .or. i1 < 1) then
             write(iulog,*) trim(subname),' missing driver record for ',ym0,' or ',ym1
             call shr_sys_abort(subname//' missing monthly record in the snow/ice driver')
          end if
          call read_month('r', r0, i0) ; call read_month('r', r1, i1)
          call read_month('t', t0, i0) ; call read_month('t', t1, i1)
       end if
       cur_ym0 = ym0 ; cur_ym1 = ym1
       if (masterproc) then
          write(iulog,*) trim(subname),' snow/ice bracket now ',ym0,' -> ',ym1, &
                         ' (source ',trim(lake_snowice_source),')'
          call shr_sys_flush(iulog)
       end if
    end if

    ! dh_ice/dt is constant across the bracket. It is unused in this version but is the
    ! one quantity a prognostic latent-heat term would need, and deriving it here keeps it
    ! consistent with the h_ice actually being applied.
    dt_span = span*86400._r8

    if (emul) then
       do n = rtmCTL%begr, rtmCTL%endr
          TLake_r%lake_ice_fraction(n) = (1._r8-w)*p_r0(1,n) + w*p_r1(1,n)
          TLake_r%lake_ice_depth(n)    = (1._r8-w)*p_r0(2,n) + w*p_r1(2,n)
          TLake_r%lake_snow_fraction(n)= (1._r8-w)*p_r0(3,n) + w*p_r1(3,n)
          TLake_r%lake_snow_depth(n)   = (1._r8-w)*p_r0(4,n) + w*p_r1(4,n)
          TLake_r%lake_dhice_dt(n)     = (p_r1(2,n) - p_r0(2,n)) / dt_span

          TLake_t%lake_ice_fraction(n) = (1._r8-w)*p_t0(1,n) + w*p_t1(1,n)
          TLake_t%lake_ice_depth(n)    = (1._r8-w)*p_t0(2,n) + w*p_t1(2,n)
          TLake_t%lake_snow_fraction(n)= (1._r8-w)*p_t0(3,n) + w*p_t1(3,n)
          TLake_t%lake_snow_depth(n)   = (1._r8-w)*p_t0(4,n) + w*p_t1(4,n)
          TLake_t%lake_dhice_dt(n)     = (p_t1(2,n) - p_t0(2,n)) / dt_span
       end do
    else
       do n = rtmCTL%begr, rtmCTL%endr
          TLake_r%lake_ice_fraction(n) = (1._r8-w)*r0%fice(n) + w*r1%fice(n)
          TLake_r%lake_ice_depth(n)    = (1._r8-w)*r0%hice(n) + w*r1%hice(n)
          TLake_r%lake_snow_fraction(n)= (1._r8-w)*r0%fsnw(n) + w*r1%fsnw(n)
          TLake_r%lake_snow_depth(n)   = (1._r8-w)*r0%dsnw(n) + w*r1%dsnw(n)
          TLake_r%lake_dhice_dt(n)     = (r1%hice(n) - r0%hice(n)) / dt_span

          TLake_t%lake_ice_fraction(n) = (1._r8-w)*t0%fice(n) + w*t1%fice(n)
          TLake_t%lake_ice_depth(n)    = (1._r8-w)*t0%hice(n) + w*t1%hice(n)
          TLake_t%lake_snow_fraction(n)= (1._r8-w)*t0%fsnw(n) + w*t1%fsnw(n)
          TLake_t%lake_snow_depth(n)   = (1._r8-w)*t0%dsnw(n) + w*t1%dsnw(n)
          TLake_t%lake_dhice_dt(n)     = (t1%hice(n) - t0%hice(n)) / dt_span
       end do
    end if

  end subroutine lake_snowice_update

!-------------------------------------------------------------------------------------------
  subroutine lake_ice_snow_bc(ice_fraction, ice_depth, snow_fraction, snow_depth, &
       sw_down, water_albedo, q_open_non_sw, water_temp, air_temp, wind,          &
       fopen, sh_net_to_water, non_sw_flux_to_water, wind_mix_factor, t_ice_surf)
    ! Partition the lake surface into open water / bare ice / snow-covered ice and return
    ! the shortwave and non-shortwave fluxes reaching the water, per unit TOTAL lake area.
    !
    ! Ice surface temperature is not carried as a state. The ice is represented as a
    ! thermal resistance in series with the atmosphere, so the skin temperature follows in
    ! closed form and is returned as a diagnostic. Fluxes over ice are therefore a bulk
    ! conductance to air, not a resolved surface energy balance: no separate longwave,
    ! sublimation or over-ice saturation vapour pressure.

    real(r8), intent(in)  :: ice_fraction, ice_depth, snow_fraction, snow_depth
    real(r8), intent(in)  :: sw_down          ! downward shortwave (W/m2)
    real(r8), intent(in)  :: water_albedo     ! open-water albedo (-)
    real(r8), intent(in)  :: q_open_non_sw    ! non-shortwave flux per unit OPEN WATER area (W/m2)
    real(r8), intent(in)  :: water_temp       ! lake surface water temperature (K)
    real(r8), intent(in)  :: air_temp         ! air temperature (K)
    real(r8), intent(in)  :: wind             ! wind speed (m/s)
    real(r8), intent(out) :: fopen, sh_net_to_water, non_sw_flux_to_water
    real(r8), intent(out) :: wind_mix_factor, t_ice_surf

    real(r8) :: fbareice, fsnowice, hice, hsnow
    real(r8) :: lambda_atm, r_atm, r_bareice, r_snowice, r_eff, ficetot

    fopen    = max(0._r8, min(1._r8, 1._r8 - ice_fraction))
    fsnowice = max(0._r8, min(1._r8, snow_fraction))
    fbareice = max(0._r8, min(1._r8, ice_fraction - snow_fraction))
    hice     = max(0._r8, ice_depth)
    hsnow    = max(0._r8, snow_depth)

    sh_net_to_water = max(sw_down, 0._r8) * (                        &
         fopen    * (1._r8 - water_albedo)                        +  &
         fbareice * (1._r8 - ice_albedo ) * ice_sw_trans          +  &
         fsnowice * (1._r8 - snow_albedo) * snow_sw_trans )

    ! Atmospheric resistance from the model's own exchange coefficients: the linearised
    ! longwave term plus the bulk sensible coefficient used in the open-water budget. A
    ! fixed floor here would give an unphysical flux as the ice thins (a 0.01 m2K/W floor
    ! implies a 100 W/m2/K exchange coefficient, several times the realistic value).
    lambda_atm = 4._r8*emiss_ice*stefan_b*air_temp**3 + 1.5701_r8*wind
    r_atm      = 1._r8/max(lambda_atm, TINYVALUE)
    r_bareice  = hice/k_ice + r_atm
    r_snowice  = hice/k_ice + hsnow/k_snow + r_atm

    non_sw_flux_to_water = fopen*q_open_non_sw                          &
         + fbareice*(air_temp - water_temp)/max(r_bareice, min_resist)  &
         + fsnowice*(air_temp - water_temp)/max(r_snowice, min_resist)

    ! Wind-driven mixing acts only on the open fraction.
    wind_mix_factor = fopen

    ficetot = fbareice + fsnowice
    if (ficetot > snowice_tol) then
       r_eff = (fbareice*r_bareice + fsnowice*r_snowice)/ficetot
       t_ice_surf = water_temp + (air_temp - water_temp)*(r_eff - r_atm)/max(r_eff, TINYVALUE)
    else
       t_ice_surf = water_temp
    end if

  end subroutine lake_ice_snow_bc

!-------------------------------------------------------------------------------------------
! Helpers
!-------------------------------------------------------------------------------------------
  subroutine alloc_buf(b)
    type(snowice_buf), intent(inout) :: b
    integer :: begr, endr
    begr = rtmCTL%begr ; endr = rtmCTL%endr
    allocate(b%fice(begr:endr), b%hice(begr:endr), b%fsnw(begr:endr), b%dsnw(begr:endr))
    b%fice = 0._r8 ; b%hice = 0._r8 ; b%fsnw = 0._r8 ; b%dsnw = 0._r8
  end subroutine alloc_buf

  subroutine read_month(side, b, tidx)
    character(len=1),  intent(in)    :: side      ! 'r' (main channel) or 't' (tributary)
    type(snowice_buf), intent(inout) :: b
    integer,           intent(in)    :: tidx
    call read_field(side//'_lake_ice_fraction',       b%fice, tidx)
    call read_field(side//'_lake_ice_thickness',      b%hice, tidx)
    call read_field(side//'_lake_snow_cover_fraction',b%fsnw, tidx)
    call read_field(side//'_lake_snow_depth',         b%dsnw, tidx)
  end subroutine read_month

  subroutine read_field(varname, field, tidx)
    character(len=*),  intent(in)    :: varname
    real(r8), pointer, intent(inout) :: field(:)
    integer,           intent(in)    :: tidx
    logical :: found
    character(len=*), parameter :: subname = '(read_field)'
    call ncd_io(ncid=ncid, varname=trim(varname), data=field, dim1name='allrof', &
                flag='read', nt=tidx, readvar=found)
    if (.not. found) call shr_sys_abort(subname//' missing driver variable '//trim(varname))
  end subroutine read_field

  integer function time_index(yyyymm)
    integer, intent(in) :: yyyymm
    integer :: k
    time_index = -1
    do k = 1, ntime_driver
       if (driver_yyyymm(k) == yyyymm) then
          time_index = k
          return
       end if
    end do
  end function time_index

  integer function days_in_month(yr, mon)
    integer, intent(in) :: yr, mon
    integer, parameter :: dm(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    days_in_month = dm(mon)
    if (mon == 2 .and. trim(get_calendar()) == 'GREGORIAN') then
       if ((mod(yr,4) == 0 .and. mod(yr,100) /= 0) .or. mod(yr,400) == 0) days_in_month = 29
    end if
  end function days_in_month

  subroutine prev_month(yr, mon, pyr, pmon)
    integer, intent(in)  :: yr, mon
    integer, intent(out) :: pyr, pmon
    if (mon == 1) then
       pyr = yr - 1 ; pmon = 12
    else
       pyr = yr ; pmon = mon - 1
    end if
  end subroutine prev_month

  subroutine next_month(yr, mon, nyr, nmon)
    integer, intent(in)  :: yr, mon
    integer, intent(out) :: nyr, nmon
    if (mon == 12) then
       nyr = yr + 1 ; nmon = 1
    else
       nyr = yr ; nmon = mon + 1
    end if
  end subroutine next_month

end MODULE mosart_lake_snowice_mod
