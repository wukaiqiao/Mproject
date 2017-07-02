!  Module name: URS2	                                               !
!                                                                      !
!  Purpose: User Defined functions to implement corrections Sundaresan !
!                                                                      !
!  This routine is called from the calc_coeff and transport_prop       !
!                                                                      !
!  Author: M. van Dijk                                Date: 05-JUNE-17 !
!                                                                      !
!  Comments:                                                           !
!  
!  GRANULAR_ENERGY = .TRUE.                                         C
!        FRICTION = .TRUE.                                             C
!           EP_s(IJK,M) > EPS_f_min  -->  friction + viscous(pde)      C
!           EP_s(IJK,M) < EP_f_min   -->  viscous (pde)                C
!        FRICTION = .FALSE.                                            C
!           EP_g < EP_star  -->  friction_schaeffer + viscous(pde)     C
!           EP_g >= EP_star -->  viscous (pde)   
!
! Please disble blending   
! Can we do visc_s then muplitply by a term?
! Otherwise the stress has to be added as a source term.in solve v starred
!
!IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
!                Mu_star = eta0
!Suggest to call after calc_mu_default, then it changes epmu_s naturally                                                                  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

! Changes in source code:
! Created new RDF, defined as USR_RDF
! Used correlation of Torquato for USR_RDF
! Deleted standard call for USR2
! Added calls for USR2 in calc_coeff.f and transport_prop.f

SUBROUTINE USR2(INDEX2)

! initialize variables from other routines
!---------------------------------------------------------------------//

  USE rdf
  USE compar
  USE coeff
  USE constant
  use discretelement, only: des_mmax
  USE fldvar
  USE functions
  USE geometry
  USE indices
  USE param
  USE param1
  USE physprop
  USE visc_s !get solid viscosity variables
  USE kintheory
  USE run
  USE toleranc
  USE usr
  !USE tau_s, only : TAU_U_S, TAU_V_S, TAU_W_S
  USE Trace !Wu
  IMPLICIT NONE

! index to check if tau or dissipation will be adjusted
  INTEGER, INTENT(IN) :: INDEX2

! Indices
  INTEGER :: IJK, M

! variables for GD model
  DOUBLE PRECISION :: press_star, c_star, zeta0_star, &
                          nu_gamma_star, &
                          lambda_num, cd_num, zeta1
  DOUBLE PRECISION :: D_PM, EP_SM
  DOUBLE PRECISION :: nu0, Chi

! variables for Sunderesan implementation
  DOUBLE PRECISION :: eps_c, par_fric, eta_k, eta_c, eta_b, eta_star
  DOUBLE PRECISION :: Psi, fun_mu, e_eff, Mu_ys
  DOUBLE PRECISION :: H_ce, J_ce, K_ce, I_I, K_K, M_CE
  DOUBLE PRECISION :: Chi_I, beta, temp_var 
  DOUBLE PRECISION :: Delta_tau, Dleta_gamma,DeltaX
  DOUBLE PRECISION :: gamma_s, I_0, alpha_0, alpha_1, alpha_2  
  DOUBLE PRECISION :: RateStrain(3,3), VelGrad(3,3)
  DOUBLE PRECISION :: TAU_SHEAR, ETA_SHEAR, P_gd

!---------------------------------------------------------------------//

! Initiate loop to calculate the adjustments according to Sundaresan
  DO IJK = ijkstart3, ijkend3
    IF ( FLUID_AT(IJK) ) THEN

   ! local aliases
      M=MMAX                !only one solid phase
      Chi = G_0(IJK,M,M)    !Assgin the radial distribution function to Chi
      EP_SM = EP_s(IJK,M)   !Assgin the solid fraction to EP_SM
      D_PM = D_P(IJK,M)     !Assgin the solid diameter to D_PM 

    !You should add ROP_SM = ROP_S(IJK,M) to call density
   ! local parameters from Sundaresan for par_fric = 0.1 (function! will be better)
      par_fric = 0.1
      ETA_SHEAR = 0.268 !yield stress ratio

   ! **************u-independent parameters
      I_0 = 0.2 ! Ineartia number 0
      alpha_0 = 0.36 
      alpha_1 = 0.06 
      alpha_2 = 0.58
   ! **************u-independent parameters

!Wu: 
! calculation of Sundaresan correction factors
     ! inertia contributions to viscosity (I think all of this should be taken from GD?)
      eta_k = (1.d0-(2.d0/5.d0)*(1.d0+C_e)*(1.d0-3.d0*C_e)* &
        EP_SM*Chi)/((1.d0-(1.d0/4.d0)*(1.d0-C_e)**2.d0-(5.d0/24.d0)* & 
        (1.d0-C_e**2.d0))*Chi)

      eta_c = (4.d0/5.d0)*(1.d0+C_e)*EP_SM*Chi*eta_k !calculated in pde_gd

      eta_b = (384.d0/(25.d0*pi))*(1.d0+C_e)*EP_SM**2.d0*Chi !calculated in pde_gd

      eta_star = eta_k + eta_c + eta_b !calculated in pde_gd
      
      psi = 1.d0+(3.d0/10.d0)*(1.d0-C_e**2.d0)**(-2.d0/3.d0)*  &
       (1.d0-exp(-8.d0*par_fric)) 

    ! function for the effective restitution coefficient (!this can be a constant used!!!!!)
      fun_mu = (3.d0/2.d0)*par_fric*exp(-3.d0*par_fric)

    ! effective restitution coefficient (!this can be a constant used!!!!!!!)
      e_eff = C_e - fun_mu  

    ! dimensionless functions
      H_ce = EP_SM*(1.d0+2.d0*(1.d0+C_e)*EP_SM*Chi)
 
      J_ce = (5.d0*DSQRT(pi)/96.d0)*eta_star

      K_ce = (12.d0/DSQRT(pi))*EP_SM**2.d0*Chi*(1.d0-C_e**2.d0)

     ! adjusted energy dissipation functions with effective restitution
      K_K = K_ce*((1.d0-e_eff**2.d0)/(1.d0-C_e**2.d0)) !(K_prime)
    
      M_CE = max(J_ce/K_K, (alpha_1/DSQRT(EP_star_array(IJK)-EP_SM)))

   !WU:*************************************************
 
 !**************Calcaute shear rate modulus**************
 !Directly use trD_s2 in trace module which is Ds_m:Ds_m (calculated in calc_mu_s)
 !shear rate = 2* strain rate modulus, defined for 3D domain in Sun 2012 paper
 !strain rate modulus = sqrt(0.5*Ds_m:Ds_m)

      gamma_s =2.d0*dsqrt(0.5d0*trD_s2(IJK,M)) 

 !missing a solid fraction term compared to the Sun equation 6, here added
 !as it required to calcute the inetial number, but not sure if correct!!!!
 !need to check other paper!!!
 ! P_gd = P_s_c(IJK,M) 
      P_gd = P_s_c(IJK,M)*EP_SM 
 
 ! inertial number as Sun paper
      I_I = (gamma_s*D_PM)/(EP_SM*dsqrt(P_gd/ROP_S(IJK,M)))

 ! correction coefficient for inertial stress 
      DeltaX = 1.d0/((I_0/I_I)**1.5d0+1) !(It calls X?) 

 ! DeltaX (Check the constant zero definded in which module)
 ! Mu_ys is the yield stress contribution on the Mu_ys
 ! Comment: in this way, the equation become Tau_ys = mu_ys* D, rather than
 ! the deviator of D (which is S). Otherwise, it should modify the diffusion coefficient in the
 ! by adding - 1/3 trace_s_c * mu_ys
 
      IF (((DeltaX>ZERO).AND.(DeltaX<ONE)).AND.(trD_s2(IJK,M)>ZERO)) THEN !(Prevent diverging)
  !      Mu_ys = (1-DeltaX) * ETA_SHEAR * P_gd / dsqrt(0.5*trD_s2(IJK,M)) 
         Mu_ys = (1-DeltaX) * ETA_SHEAR * P_s_c(IJK,M) / dsqrt(0.5*trD_s2(IJK,M)) 
      ELSE
        Mu_ys = 0.d0
      ENDIF 
  
  !*****adjust the solid viscosity/stress*****

  ! transition function
      beta = EP_SM*psi*J_ce*dsqrt(K_ce/(K_K*H_ce))
   
    ! correction factor for the inertia shear stress (Delat_tau)
      Delta_tau = alpha_0/beta+(1.d0-alpha_0/beta)*Chi_I

    !****correct the mu_s_c(KT_mu) and mu_s (solid mu)***
    
    !In Sun's paper, it mentions only shear stress. If they treate solid as 
    !incompressiable, then bluk viscosity goes to zero, and D=S. They do not 
    !mention any lamda in theri paper, which is debatable. Therefore, the shear
    !stress is then mu_s_c * (2D-2/3trd(D)) and change this part only. Here we 
    !change the entire viscous stress. 2*mu_s_c *D+lamda*trd(D)
    !Otherwise, we set lambda to be zero then.....
    !In this way, the KT stress is: Delta_tau*(2*mu_s_c*D+lamda*trd(D))
      IF (INDEX2==1) THEN
  !Subroutine called at end of CALC_DEFAULT_MUS(M)
    
          mu_s_c(IJK,M) = mu_s_c(IJK,M)*Delta_tau
     !    lambda_s_c(IJK,M) = lambda_s_c(IJK,M)*Delta_tau
          lambda_s_c(IJK,M) = ZERO  ! test to see if the lambda zero runs
          ! IF (IJK==1000) write(*,*) 'mu_s_c and lambda_s_c changed successfully', mu_s_c(IJK,M),lambda_s_c(IJK,M)
      ENDIF
  
      IF (INDEX2==2) THEN
  !Subroutine called at end of CALC_DEFAULT_MUS(M)
            IF (mu_ys/=0) THEN
            ITERMU = mu_ys
            u_muprime = mu_ys - ITERMU
            usr_ratio = mu_ys / mu_s(IJK,M)
            u_muprime= 1 - exp(-50000.d0* dsqrt(0.5*trD_s2(IJK,M)) / ((1-DeltaX) * ETA_SHEAR * P_s_c(IJK,M)))

              IF (usr_ratio > 1000)  THEN
               ITERMU = u_muprime * mu_ys
              ENDIF
            mu_s(IJK,M) = mu_s_c(IJK,M) + ITERMU 
!           IF (IJK==1000) write(*,*) 'mu_s_v changed successfully',mu_ys,ITERMU

            USRrdfg(IJK,M) = mu_ys 
            ENDIF
           
           
      ENDIF
       
       
!****************************************************** 
  ! correction factor for the dissipation energy (Delat_tau)
  ! To be insert at calc_trasport_prop

      Dleta_gamma = beta*dsqrt(H_CE)/(K_K*M_ce)*Delta_tau
    
      IF (INDEX2 == 3) THEN
   ! If USR2 is called after dissipation calc., correct diss. terms
   !*******************missing power of 2???, what is specie?????**************     
   ! old diss. term = old diss. term * effective restitution * correction
   ! EDT_s_ip = energy diss. with difference in species gran. T. transfer
          EDT_s_ip(IJK,M,M) = EDT_s_ip(IJK,M,M) * &
          (1.d0-e_eff**2.d0)/(1.d0-C_e**2.d0)*Dleta_gamma
 
      ! EDvel_sM_ip = eta(1,1) writen in GD paper, but it neglected in Sun paper
          EDvel_sM_ip(IJK,M,M) = EDvel_sM_ip(IJK,M,M) * &
          (1.d0-e_eff**2.d0)/(1.d0-C_e**2.d0)*Dleta_gamma
      ENDIF
    !  ENDIF !end if (INDEX2==2)

!******************************************************
 ! Ajust the solid viscosity 
 ! You don't call right afte calc_default, epmu_s has to be changed as well.
    ENDIF
  ENDDO
  
  return 

END SUBROUTINE USR2


