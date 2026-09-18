MODULE mosart_lake_snowice_kernel_mod
!-------------------------------------------------------------------------------------------
! Numerical kernel of the frozen snow/ice LSTM (T23), with no MOSART dependencies.
!
! Kept separate from mosart_lake_snowice_emulator_mod on purpose: the emulator module owns
! file I/O, the object mapping and the MOSART state, none of which a correctness test for
! the arithmetic should have to link against. The golden vectors packaged with the model
! are checked against THIS module by tools/test_snowice_kernel.F90, so the tested code and
! the running code are the same code.
!
! Contract (verified against causal_lstm_snowice_v2.nc golden vectors to 1.4e-15):
!   gates = W_ih x + W_hh h + b        gate order i,f,g,o; b is PyTorch bias_ih + bias_hh
!   c     = sigmoid(f)*c + sigmoid(i)*tanh(g)
!   h     = sigmoid(o)*tanh(c)
!   raw   = W_head h + b_head
!   f_ice = sigmoid(raw1) ; h_ice = softplus(raw2)
!   f_snow= f_ice*sigmoid(raw3) ; d_snow = softplus(raw4)
!   then the export post-processing cascade with a 1e-4 zero threshold.
!
! Array shapes follow what NetCDF hands Fortran, i.e. transposed relative to the Python
! declarations: lstm_weight_ih is (384,15) in the file and (15,384) here.
!-------------------------------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: NDYN, NSTA, NIN, NHID, NTGT, NGAT
  public :: snowice_lstm_step
  public :: snowice_head_post
  public :: snowice_sigmoid
  public :: snowice_softplus

  integer, parameter :: NDYN = 10          ! hurs,huss,sfcwind,tasmax,tasmin,rlds,rsds,prsn,ps,pr
  integer, parameter :: NSTA = 5           ! packaged statics, already normalized
  integer, parameter :: NIN  = NDYN + NSTA ! 15
  integer, parameter :: NHID = 96
  integer, parameter :: NTGT = 4           ! f_ice, h_ice, f_snow, d_snow
  integer, parameter :: NGAT = 4*NHID      ! 384

  real(r8), parameter :: ZERO_THRESH = 1.0e-4_r8  ! export post-processing threshold
  real(r8), parameter :: SP_THRESH   = 20.0_r8    ! torch softplus(beta=1, threshold=20)

contains

!-------------------------------------------------------------------------------------------
  subroutine snowice_lstm_step(w_ih, w_hh, b, x, h, c)

    real(r8), intent(in)    :: w_ih(NIN,NGAT)
    real(r8), intent(in)    :: w_hh(NHID,NGAT)
    real(r8), intent(in)    :: b(NGAT)
    real(r8), intent(in)    :: x(NIN)
    real(r8), intent(inout) :: h(NHID), c(NHID)

    real(r8) :: g(NGAT)
    integer  :: k

    g = matmul(x, w_ih) + matmul(h, w_hh) + b
    do k = 1, NHID
       c(k) = snowice_sigmoid(g(NHID+k))*c(k) + snowice_sigmoid(g(k))*tanh(g(2*NHID+k))
       h(k) = snowice_sigmoid(g(3*NHID+k))*tanh(c(k))
    end do

  end subroutine snowice_lstm_step

!-------------------------------------------------------------------------------------------
  subroutine snowice_head_post(w_head, b_head, h, y, y_raw)
    ! y     = post-processed targets (what MOSART uses)
    ! y_raw = transformed but un-thresholded targets, which is what the packaged
    !         golden_output_r8 holds. Optional so the run path can skip it.

    real(r8), intent(in)  :: w_head(NHID,NTGT)
    real(r8), intent(in)  :: b_head(NTGT)
    real(r8), intent(in)  :: h(NHID)
    real(r8), intent(out) :: y(NTGT)
    real(r8), intent(out), optional :: y_raw(NTGT)

    real(r8) :: raw(NTGT), fice, hice, fsnw, dsnw

    raw = matmul(h, w_head) + b_head

    fice = snowice_sigmoid(raw(1))
    hice = snowice_softplus(raw(2))
    fsnw = fice * snowice_sigmoid(raw(3))
    dsnw = snowice_softplus(raw(4))

    if (present(y_raw)) y_raw = (/ fice, hice, fsnw, dsnw /)

    if (fice < ZERO_THRESH) fice = 0._r8
    if (fice <= 0._r8) then
       hice = 0._r8
    else if (hice < ZERO_THRESH) then
       hice = 0._r8
    end if
    fsnw = min(fsnw, fice)
    if (fsnw < ZERO_THRESH) fsnw = 0._r8
    if (fsnw <= 0._r8) then
       dsnw = 0._r8
    else if (dsnw < ZERO_THRESH) then
       dsnw = 0._r8
    end if

    y(1) = min(max(fice,0._r8),1._r8)
    y(2) = max(hice,0._r8)
    y(3) = min(max(fsnw,0._r8), y(1))
    y(4) = max(dsnw,0._r8)

  end subroutine snowice_head_post

!-------------------------------------------------------------------------------------------
  pure real(r8) function snowice_sigmoid(z) result(s)
    real(r8), intent(in) :: z
    if (z >= 0._r8) then
       s = 1._r8/(1._r8 + exp(-z))
    else
       s = exp(z)/(1._r8 + exp(z))
    end if
  end function snowice_sigmoid

  pure real(r8) function snowice_softplus(z) result(s)
    real(r8), intent(in) :: z
    if (z > SP_THRESH) then
       s = z
    else
       s = log(1._r8 + exp(-abs(z))) + max(z,0._r8)
    end if
  end function snowice_softplus

END MODULE mosart_lake_snowice_kernel_mod
