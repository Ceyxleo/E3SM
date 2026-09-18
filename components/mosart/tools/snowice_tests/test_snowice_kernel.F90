program test_snowice_kernel
!-------------------------------------------------------------------------------------------
! Golden-vector check for mosart_lake_snowice_kernel_mod.
!
! Reads the frozen model package with plain NetCDF (no PIO, no MOSART), replays the eight
! packaged golden objects over their 24 months from h=c=0, and compares against the
! packaged hidden state, cell state, transformed output and post-processed output.
!
!   ftn -ffree-line-length-none -o test_snowice_kernel \
!       shr_kind_mod.F90 MOSART_lake_snowice_kernel_mod.F90 test_snowice_kernel.F90 -lnetcdff
!-------------------------------------------------------------------------------------------
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mosart_lake_snowice_kernel_mod
  use netcdf
  implicit none

  character(len=*), parameter :: MODELF = &
    '/pscratch/sd/l/lil325/mosart_lake_ml_emulator/outputs/deployment/causal_lstm_snowice_v2.nc'

  integer :: ncid, vid, dimid, nobj_g, ntime_g, o, t, ier
  real(r8), allocatable :: w_ih(:,:), w_hh(:,:), b(:), w_head(:,:), b_head(:)
  real(r8), allocatable :: gin(:,:,:), gout(:,:,:), gpp(:,:,:), gh(:,:,:), gc(:,:,:)
  real(r8) :: h(NHID), c(NHID), y(NTGT), yraw(NTGT)
  real(r8) :: eh, ec, eo, ep, mh, mc, mo, mp
  logical  :: bad

  call chk(nf90_open(MODELF, NF90_NOWRITE, ncid), 'open')
  call chk(nf90_inq_dimid(ncid,'golden_object',dimid),'d1')
  call chk(nf90_inquire_dimension(ncid,dimid,len=nobj_g),'d1len')
  call chk(nf90_inq_dimid(ncid,'golden_time',dimid),'d2')
  call chk(nf90_inquire_dimension(ncid,dimid,len=ntime_g),'d2len')
  write(*,'(a,i0,a,i0,a)') ' golden vectors: ',nobj_g,' objects x ',ntime_g,' months'

  allocate(w_ih(NIN,NGAT), w_hh(NHID,NGAT), b(NGAT), w_head(NHID,NTGT), b_head(NTGT))
  allocate(gin(NIN,ntime_g,nobj_g), gout(NTGT,ntime_g,nobj_g), gpp(NTGT,ntime_g,nobj_g))
  allocate(gh(NHID,ntime_g,nobj_g), gc(NHID,ntime_g,nobj_g))

  call rd('lstm_weight_ih', w_ih);  call rd('lstm_weight_hh', w_hh)
  call rd1('lstm_bias', b);         call rd('head_weight', w_head)
  call rd1('head_bias', b_head)
  call rd3('golden_model_input', gin)
  call rd3('golden_output_r8', gout)
  call rd3('golden_output_postprocessed_r8', gpp)
  call rd3('golden_hidden_r8', gh)
  call rd3('golden_cell_r8', gc)
  call chk(nf90_close(ncid),'close')

  mh = 0._r8 ; mc = 0._r8 ; mo = 0._r8 ; mp = 0._r8
  do o = 1, nobj_g
     h = 0._r8 ; c = 0._r8
     do t = 1, ntime_g
        call snowice_lstm_step(w_ih, w_hh, b, gin(:,t,o), h, c)
        call snowice_head_post(w_head, b_head, h, y, yraw)
        eh = maxval(abs(h - gh(:,t,o)))      ; mh = max(mh, eh)
        ec = maxval(abs(c - gc(:,t,o)))      ; mc = max(mc, ec)
        eo = maxval(abs(yraw - gout(:,t,o))) ; mo = max(mo, eo)
        ep = maxval(abs(y    - gpp(:,t,o)))  ; mp = max(mp, ep)
     end do
  end do

  write(*,'(a,es12.4)') ' max |hidden      - golden| = ', mh
  write(*,'(a,es12.4)') ' max |cell        - golden| = ', mc
  write(*,'(a,es12.4)') ' max |output      - golden| = ', mo
  write(*,'(a,es12.4)') ' max |postprocess - golden| = ', mp
  bad = (mh > 1.e-10_r8) .or. (mc > 1.e-10_r8) .or. (mo > 1.e-10_r8) .or. (mp > 1.e-10_r8)
  if (bad) then
     write(*,*) ' FAIL: tolerance 1e-10 exceeded'
     stop 1
  else
     write(*,*) ' PASS: Fortran kernel reproduces the packaged golden vectors to 1e-10'
  end if

contains
  subroutine chk(s, w)
    integer, intent(in) :: s ; character(len=*), intent(in) :: w
    if (s /= NF90_NOERR) then
       write(*,*) 'netcdf error at ',trim(w),': ',trim(nf90_strerror(s)) ; stop 2
    end if
  end subroutine chk
  subroutine rd(n, a)
    character(len=*), intent(in) :: n ; real(r8), intent(out) :: a(:,:)
    call chk(nf90_inq_varid(ncid,n,vid),n) ; call chk(nf90_get_var(ncid,vid,a),n)
  end subroutine rd
  subroutine rd1(n, a)
    character(len=*), intent(in) :: n ; real(r8), intent(out) :: a(:)
    call chk(nf90_inq_varid(ncid,n,vid),n) ; call chk(nf90_get_var(ncid,vid,a),n)
  end subroutine rd1
  subroutine rd3(n, a)
    character(len=*), intent(in) :: n ; real(r8), intent(out) :: a(:,:,:)
    call chk(nf90_inq_varid(ncid,n,vid),n) ; call chk(nf90_get_var(ncid,vid,a),n)
  end subroutine rd3
end program test_snowice_kernel
