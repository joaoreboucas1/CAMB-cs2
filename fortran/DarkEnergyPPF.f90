    module DarkEnergyPPF
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyPPF
        real(dl) :: c_Gamma_ppf = 0.4_dl
    contains
    procedure :: ReadParams => TDarkEnergyPPF_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyPPF_PythonClass
    procedure :: Init => TDarkEnergyPPF_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyPPF_PerturbedStressEnergy
    procedure :: diff_rhopi_Add_Term => TDarkEnergyPPF_diff_rhopi_Add_Term
    procedure, nopass :: SelfPointer => TDarkEnergyPPF_SelfPointer
    procedure, private :: setcgammappf
    end type TDarkEnergyPPF

    public TDarkEnergyPPF
    contains

    subroutine TDarkEnergyPPF_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyPPF) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
    this%cs2_0 = Ini%Read_Double('cs2_0', 1.d0)
    if (this%cs2_0 /= 1._dl) error stop 'cs2_0 supported by PPF model'
    call this%setcgammappf

    end subroutine TDarkEnergyPPF_ReadParams

    function TDarkEnergyPPF_PythonClass()
    character(LEN=:), allocatable :: TDarkEnergyPPF_PythonClass

    TDarkEnergyPPF_PythonClass = 'DarkEnergyPPF'
    end function TDarkEnergyPPF_PythonClass


    subroutine TDarkEnergyPPF_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyPPF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyPPF_SelfPointer

    subroutine TDarkEnergyPPF_Init(this, State)
    use classes
    use results
    use config
    class(TDarkEnergyPPF), intent(inout) :: this
    class(TCAMBdata), intent(inout), target :: State
    ! JVR MOD BEGIN: adding variables for integrating alpha_B
    real(dl), parameter :: a_ini = 1e-5, alpha_B_ini = 0d0
    real(dl) :: a, dalpha_B, dlog_a, rho_de, rho_m, rho_gamma, rho_tot, w_tot, last_term, w_de
    integer :: i
    ! JVR MOD END

    call this%TDarkEnergyEqnOfState%Init(State)
    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        this%num_perturb_equations = 1
    end if
    ! JVR MOD BEGIN: disabled this error
    ! if (this%cs2_0 /= 1._dl) &
    !     call GlobalError('DarkEnergyPPF does not support varying sound speed',error_unsupported_params)
    ! JVR MOD END

    ! JVR MOD BEGIN: populating array of alpha_B and log_a
    select type(State)
    class is (CAMBData)
        if (State%CP%use_cs2) then
            State%CP%alpha_B(1) = alpha_B_ini
            State%CP%log_a(1) = log(a_ini)
            dlog_a = -State%CP%log_a(1)/(alpha_B_len-1)
            do i = 1, alpha_B_len-1
                a = exp(State%CP%log_a(i))
                rho_de = State%Omega_de * a**(-3.0_dl*(1.0 + this%w_lam + this%wa)) * exp(-3.0_dl*this%wa*(1.0_dl - a))
                rho_m  = (State%grhoc + State%grhob)/State%grhocrit * a**(-3.0_dl)
                rho_gamma  = (State%grhog + State%grhornomass)/State%grhocrit * a**(-4.0_dl)
                rho_tot = rho_gamma + rho_m + rho_de
                w_de = this%w_lam + this%wa*(1.0_dl - a)
                w_tot = (rho_gamma/3.0_dl + w_de*rho_de)/rho_tot
                last_term = (4.0*rho_gamma/3.0 + rho_m)/rho_tot
                dalpha_B = this%cs2_0*(State%CP%alpha_K + 1.5_dl*State%CP%alpha_B(i)**(2.0_dl)) \
                           + (State%CP%alpha_B(i) - 2.0_dl)*(1.5_dl*(1.0_dl + w_tot) + 0.5_dl*State%CP%alpha_B(i)) + 3.0_dl*last_term
                State%CP%log_a(i+1) = State%CP%log_a(i) + dlog_a
                State%CP%alpha_B(i+1) = State%CP%alpha_B(i) + dalpha_B*dlog_a
            end do
        end if
    end select
    ! JVR MOD END

    end subroutine TDarkEnergyPPF_Init

    subroutine setcgammappf(this)
    class(TDarkEnergyPPF) :: this

    this%c_Gamma_ppf = 0.4_dl * sqrt(this%cs2_0)

    end subroutine setcgammappf


    function TDarkEnergyPPF_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    !Get derivative of anisotropic stress
    class(TDarkEnergyPPF), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, w, grhok, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh

    if (this%is_cosmological_constant .or. this%no_perturbations) then
        ppiedot = 0
    else
        hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
        ppiedot = 3._dl * dgrhoe + dgqe * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t * (1 + w) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
        ppiedot = ppiedot * adotoa / Kf1
    end if

    end function TDarkEnergyPPF_diff_rhopi_Add_Term


    subroutine TDarkEnergyPPF_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyPPF), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: Gamma, S_Gamma, ckH, Gammadot, Fa, sigma
    real(dl) :: vT, grhoT, k2

    if (this%no_perturbations) then
        dgrhoe=0
        dgqe=0
        return
    end if

    k2=k**2
    !ppf
    grhoT = grho - grhov_t
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t * (1 + w) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = this%c_Gamma_ppf * k / adotoa

    if (ckH * ckH > 1000) then
        ! Was ckH^2 > 30 originally, but this is better behaved (closer to fluid)
        ! for some extreme models (thanks Yanhui Yang, Simeon Bird 2024)
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe = S_Gamma - Gammadot / adotoa - Gamma
    dgqe = -dgqe / Fa * 2._dl * k * adotoa + vT * grhov_t * (1 + w)
    dgrhoe = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe

    end subroutine TDarkEnergyPPF_PerturbedStressEnergy



    end module DarkEnergyPPF
