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
    real(dl) :: a, dalpha_B, dlog_a, w_tot, last_term, w_de, alpha_K
    real(dl) :: grho_de, grho_no_de_t, grho_tot, gpres_no_de, grho_nu, gpres_nu
    integer :: i, j, nu_i
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
            alpha_K = State%CP%alpha_K_0
            if (State%CP%alpha_K_parametrization .eq. 1) then
                call this%BackgroundDensityAndPressure(state%grhov, a_ini, grho_de)
                grho_tot = State%grho_no_de(a_ini)/(a_ini**2) + grho_de
                alpha_K = alpha_K*grho_de/grho_tot/State%Omega_DE
            end if
            State%CP%mu(1) = 1.0_dl + State%CP%alpha_B(1)**(2) \
                             / (2.0_dl*this%cs2_0*(alpha_K + 1.5_dl*State%CP%alpha_B(1)**(2)))
            dlog_a = -State%CP%log_a(1)/(alpha_B_len-1)
            do i = 1, alpha_B_len-1
                if (State%CP%alpha_B(i)**2 > 1e6*alpha_K) then
                    ! JVR NOTE: for many cases, \alpha_B just diverges (i.e. becomes too big and positive)
                    ! This is not a problem since \mu has a well-defined limit when \alpha_B -> \inf
                    ! In practice, I enforce this with the threshold defined above in the `if` statement
                    ! And then I just fill the rest of the arrays with the last values and break out of the integration loop
                    do j = i, alpha_B_len
                        State%CP%log_a(j) = State%CP%log_a(1) + (j-1)*dlog_a
                        State%CP%alpha_B(j) = State%CP%alpha_B(i)
                        State%CP%mu(j) = State%CP%mu(i)
                    end do
                    exit
                end if
                a = exp(State%CP%log_a(i))
                ! NOTE: in CAMB convention, grho = 8*pi*G*a^2*rho
                call this%BackgroundDensityAndPressure(State%grhov, a, grho_de)
                w_de = this%w_lam + this%wa*(1.0_dl - a)
                grho_no_de_t = State%grho_no_de(a)/a/a ! NOTE: grho_no_de returns 8*pi*G*a^4*rho
                grho_tot = grho_no_de_t + grho_de

                gpres_no_de = 0.0 ! NOTE: counts 8*pi*G*a^2*P, the pressure of photons, massless and massive nu
                if (State%CP%Num_Nu_Massive > 0) then
                    do nu_i = 1, State%CP%nu_mass_eigenstates
                        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i), grho_nu, gpres_nu)
                        gpres_no_de = gpres_no_de + gpres_nu
                    end do
                end if

                gpres_no_de = gpres_no_de + (State%grhog + State%grhornomass)/3.0_dl/a**2

                w_tot = (gpres_no_de + w_de*grho_de)/grho_tot

                last_term = (grho_no_de_t + gpres_no_de)/grho_tot
                alpha_K = State%CP%alpha_K_0
                if (State%CP%alpha_K_parametrization .eq. 1) then
                    alpha_K = alpha_K*grho_de/grho_tot/State%Omega_DE
                end if

                dalpha_B = this%cs2_0*(alpha_K + 1.5_dl*State%CP%alpha_B(i)**(2)) \
                           + (State%CP%alpha_B(i) - 2.0_dl)*(1.5_dl*(1.0_dl + w_tot) + 0.5_dl*State%CP%alpha_B(i))\
                           + 3.0_dl*last_term
                State%CP%log_a(i+1) = State%CP%log_a(i) + dlog_a
                State%CP%alpha_B(i+1) = State%CP%alpha_B(i) + dalpha_B*dlog_a
                State%CP%mu(i+1) = 1.0_dl + State%CP%alpha_B(i+1)**(2) \
                             / (2.0_dl*this%cs2_0*(alpha_K + 1.5_dl*State%CP%alpha_B(i+1)**(2)))
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
