C**************************************************************
C*********** MAT. TANG. CALCS FOR EQ. PLASTIC STRAIN **********  
C**************************************************************   
      subroutine get_stress( 
     1  pl_props, hist_n, total_strain, eps_eq_pl_ip,
     2  lc, hist_new, pl_stress, dStress_dStran, eps_eq_drive, KINC,
     3  DTIME, dStress_dEpsPl, dStress_dEpsBar_V, dEpsPl_dEps_V, 
     4  dEpsPl_dEpsPl_bar, PNEWDT, eps_eq_pl_ip_neg, 
     5  eps_eq_drive_neg, dStress_dEpsBar_neg_V, dEpsPl_neg_dEps_V,
     6  dEpsPl_dEpsPl_neg_bar, dEpsPl_neg_dEpsPl_bar, 
     7  dEpsPl_neg_dEpsPl_neg_bar, G_TENS, G_COMP, dg_dEps_V, dgNeg_dEps_V,
     8  delta_eps_eq_neg, delta_eps_eq_pos)    
! C************************** INPUT *********************************
! C* pl_props: Material properties (10 variables. Not all of them are used)
! C   E = pl_props(1), nu = pl_props(2), lc = pl_props(3)
! C   aa = pl_props(4), bb = pl_props(5), sigma_t_inf = pl_props(6)
! C   nu_p = pl_props(7), e1p(1) = pl_props(8), e1p(2) = pl_props(9)
! C   e1p(3) = pl_props(10)
! C   sigma_t_inf and e1p are not used.
! C* hist_n: History variables from previous step
! C   hist_n(1-6) = Total plastic strain from previous step
! C   hist_n(7) = local equivalent plastic strain
! C   hist_n(8) = load flag (Utilizing this tension and compression loads are distingushed)
! C   hist_n(9) = non-local equivalent plastic strain
! C   hist_n(10) = damage
! C* total_strain: Current strain
! C* eps_eq_pl_ip: Current non-local equivalent-plastic-strain (solution of non-local part)    
! C* lc: Internal length scale     
! C************************** OUTPUT ********************************
! C* hist_new: Updated history variables (14 variables)
! C* pl_stress: Updated stress tensor
! C* dStress_dStran: Material tangent
! C* eps_eq_drive: (LOCAL equivalent plastic strain) driving parameter for non-local solution
! C* dStress_dEpsPl: Derivative of stress with respect to Plastic Strain (? bak bir)
! C* dStress_dEpsBar_V: Derivative of stress with respect to non-local solution variable
! C* dEpsPl_dEps_V: Derivative of (LOCAL) equivalent-plastic-strain w.r.t. total strain
! C* dEpsPl_dEpsPl: Derivative of (LOCAL) equivalent-plastic-strain w.r.t. (LOCAL) equivalent-plastic-strain
! C************************** MODIFIED ******************************
! C* DTIME: Time step (Abaqus variable)
! C******* new terms form three field approach
! C eps_pl_eq_ip_neg negative eq.pl.strain
! C****************************** NOTES ******************************  
! C* PNEWDT: If this value is updated, then analysis continues with the last
! C* converged increment with the new DTIME value. The new DTIME value is
! C* computed as DTIME = DTIME * PNEWDT
      implicit none
      INCLUDE 'KARC.blc'   		  
! C** INPUTS
      real*8, intent(in) :: pl_props(10), hist_n(14), total_strain(6),
     1  eps_eq_pl_ip, lc, eps_eq_pl_ip_neg 
      integer, intent(in) :: KINC         
! C** OUTPUTS
      real*8, intent(out) :: hist_new(14), pl_stress(6), 
     1 dStress_dStran(6,6), eps_eq_drive, dStress_dEpsPl(6,6),
     2 dStress_dEpsBar_V(6), dEpsPl_dEps_V(6), dEpsPl_dEpsPl_bar,
     3 eps_eq_drive_neg, dStress_dEpsBar_neg_V(6),
     4 dEpsPl_neg_dEps_V(6), dEpsPl_dEpsPl_neg_bar,
     5 dEpsPl_neg_dEpsPl_bar, dEpsPl_neg_dEpsPl_neg_bar, G_TENS, G_COMP,
     6 dg_dEps_V(6), dgNeg_dEps_V(6), delta_eps_eq_neg, delta_eps_eq_pos	 
! C** Stress Linearization terms
! C** dStress_dStran, C1
! C** dStress_dEpsBar_V, C2  
! C** dStress_dEpsBar_neg_V, C2_neg
! C**   
! C** Local Positive Equivalent plastic strain linearization  
! C** dEpsPl_dEps_V, C3, D3   
! C** dEpsPl_dEpsPl_bar, C4, D4  
! C** dEpsPl_dEpsPl_neg_bar, C5, D5  
! C**
! C** Local Positive Equivalent plastic strain linearization  
! C** dEpsPl_neg_dEps_V, C3_neg, D3_neg 
! C** dEpsPl_neg_dEpsPl_bar, C4_neg, D4_neg
! C** dEpsPl_neg_dEpsPl_neg_bar, C5_neg, D5_neg           
! C MODIFIED ABAQUS VARIABLES (BOTH INPUT AND OUTPUT)
      real*8, intent(inout) :: DTIME, PNEWDT  
!C** LOCAL VARIABLES
      real*8  stress_eff(6), gamma_cur, delta_eps_p(3,3)
      real*8 a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_el2D(3,3), 
     1  eps_pl_n2D(3,3), eps_pl_n(6), eps_eq_n, stress_eff_F(3,3), 
     2  I_1(3,3), II(3,3,3,3), PP(3,3,3,3), E, nu, sigma_c_0,
     3  sigma_t_0, sigma_c_inf, sigma_t_inf, nu_p, lame, mu, kappa,
     4  alpha, kk, e1p(3), total_strain2D(3,3), stress_dev_tr(3,3), 
     5  nr_stres_dev_tr, pres_tr, stress_tr(3,3), sigma_c, sigma_t,
     6  phi_tr, delta_eps_eq, eps_eq_new,
     7  phi, phi_pr, stress_dev_new(3,3), pres_new, eps_pl_new2D(3,3), 
     8  stress_new(3,3), d_flag, eps_el(6), stress_pres_eff, F_m, 
     9  kappa_str, stress_dev_eff(3,3), I1_eff, J2_eff, d_m, d_m_neg,
     1  eps_pl_fin(6),  eps_pl_new(6), kappa_0, aa, bb, kappa_old, 
     1  elastic_strain(6), Cel_4D(3,3,3,3), d_mOLD, kappa_old_neg,
     1  kappa_str_neg, d_mOLD_neg, eps_eq_n_neg, eps_eq_new_neg,
     1  c_pos, c_neg
      integer i,j,k,l,m,n
!C TEST VAR.
      real*8 dD_dKappa, dKappa_dEpsPlBar, K4, dStress_dEpsBar(3,3)
!C principal strains
      real*8 AN(3,3), fact_neg_stress     
      real*8 comp_damage, tens_damage 
!C** history variables for non-local plastic strains	
	  real*8 eps_bar_pos_n, eps_bar_neg_n
!C** Old non-local equivalent plastic strains
	  real*8 eps_eq_pl_ip_n, eps_eq_pl_ip_neg_n
!C** Principals Stresses and Hyd. stress
      real*8 PS(3), S_H, GTENS_A, GTENS_B, norm_stress	  
! C** temprorary: initiation of damages kappa_T and kappa_C
	  real*8 kappa_C, kappa_T	
	  real*8 dLim
	  real*8 dm_new, dm_neg_new
	  integer t_flag 
	  real*8 DT_1, DT_2, DT_3, DT_4, DT_5
C** 
	  dLim = 0.2d0
	  kappa_T = 3.00e-2
	  kappa_C = 9.98e-2
	  G_TENS = 1.d0 ! prevent 1/NaN erros.
      G_COMP = 1.d0
! C
! C INITIALIZE OUTPUT, prevents NaN errors
! C
      dEpsPl_dEpsPl_bar = 0.0d0
      dEpsPl_dEpsPl_neg_bar = 0.d0 
      pl_stress = 0.d0
      stress_eff = 0.d0  
      dStress_dEpsBar_V = 0.d0
      dEpsPl_dEps_V = 0.d0	 
      dStress_dStran = 0.00d0
      dStress_dEpsPl = 0.d0	  
      hist_new = 0.d0

      eps_eq_drive = 0.d0
      eps_eq_drive_neg = 0.d0
	  
	  delta_eps_p = 0.d0
C
	  GTENS_A = 1.d0 ! A*(sh/norm(s))**B
	  GTENS_B = 1.d0 ! A*(sh/norm(s))**B
C      
      gamma_cur = 0.d0
	  
	  dgNeg_dEps_V = 0.d0
	  dg_dEps_V = 0.d0
! C The amount of reduction in negative princpal stresses (NOT USED ANYMORE)     
      fact_neg_stress = 0.10d0
      fact_neg_stress = 1.0d0
! C#### Parameters for the fit of sigma_t and sigma_c
! C## Be aware that, this curve is fit for MPa. If the units are changed
! C## then this parameters must be changed as well.
      CALL stressCurve(a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c)
! C **
! C ** LOAD HISTORY VARIABLES FROM PREVIOUS STEP
! C **      
! C Load plastic strain
      do i=1, 6
        eps_pl_n(i) = hist_n(i)
      enddo
      CALL convert2FullStrain(eps_pl_n, eps_pl_n2D)
C
      eps_eq_n = hist_n(7) ! LOCAL (driving term for 2nd PDE) POSITIVE equivalent plastic strain
      eps_eq_n_neg = hist_n(8) ! LOCAL (driving term for 3rd PDE) NEGATIVE equivalent plastic strain
      kappa_old = hist_n(9) ! POSITIVE NON-LOCAL (result of the 2nd PDE) from previous step
      kappa_str = eps_eq_pl_ip ! CURRENT POSITIVE non-local kappa from second PDE
      kappa_old_neg = hist_n(10) ! NEGATIVE NON-LOCAL (result of the 3rd PDE) from previous step
      kappa_str_neg = eps_eq_pl_ip_neg ! CURRENT NEGATIVE non-local kappa from third PDE      
      d_mOLD =  hist_n(11) ! POSITIVE Damage value from previous step
      d_mOLD_neg =  hist_n(12) ! NEGATIVE Damage value from previous step  
	  eps_eq_pl_ip_n = hist_n(13) ! Positive non-local eqivalent plastic strain from previous step
	  eps_eq_pl_ip_neg_n = hist_n(14) ! Negative non-local eqivalent plastic strain from previous step  
! C **     
! C ** MATERIAL PROPERTIES
! C **
      E = pl_props(1)
      nu = pl_props(2)
      aa = pl_props(4)
      bb = pl_props(5)
      nu_p = pl_props(7)
! C
      alpha = (9.d0/2.d0)*(1.d0-2.d0*nu_p)/(1.d0+nu_p)
      kk = 1.d0/(1.d0+2.d0*(nu_p**2.d0))
! C
      lame = E*nu / ((1.d0+nu)*(1.d0-2.d0*nu))
      mu = E / (2.d0*(1.d0+nu))
      kappa = E/(3.d0*(1.d0-2.d0*nu))
! C
      CALL getIandPterms(I_1, II, PP) ! Identitty (2nd, 4th) and Deviatoric Projection
      CALL convert2FullStrain(total_strain, total_strain2D)
!C **
!C ** CHECK FOR DAMAGE EVOLUTION
!C **
!C Check for POSITIVE damage evolution
      if(kappa_str.GT.kappa_old) then
        CALL exp_damage(kappa_str, d_m) 
      else
        d_m =  hist_n(11) 
      endif  
!C Check for NEGATIVE damage evolution
      if(kappa_str_neg.GT.kappa_old_neg) then
        CALL exp_damage_comp(kappa_str_neg, d_m_neg) 
      else
        d_m_neg =  hist_n(12) 
      endif  
      
	  DMAX = max(DMAX, max(d_m, d_mOLD))
      ! CALL updateDTIME(DTIME, DMAX)
   
      kappa_str = MAX(kappa_str, kappa_old)
      kappa_str_neg = MAX(kappa_str_neg, kappa_old_neg)	
!C**
!C** Assuming no variation in plastic strain compute trial quantities.
!C** Be aware that damage is included in TRIAL part. (ELASTIC PREDICTION)
!C**
      CALL getTrial(total_strain2D, eps_pl_n2D, PP, I_1, mu, kappa,
     1  eps_eq_n, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, 
     2  stress_dev_tr, nr_stres_dev_tr, pres_tr, stress_tr, 
     3  sigma_c, sigma_t, phi_tr, d_m, d_m_neg, eps_eq_n_neg)    
C
      if(phi_tr.LT.0.d0)then ! PREDICTION IS CORRECT      
! C** simpy means that, there is no evolution in plastic strain. Hence, no modification is required in trial value        
       CALL get_D_abaqus2(E, nu, mu, dStress_dStran)	! Compute elastic stiffness tensor 
! C
! C** Update driving terms for 2nd and 3rd equations (no evolution in those terms)
       eps_eq_drive = eps_eq_n ! previous eps_pl_pos, positive local eq. pl. strain
       eps_eq_drive_neg = eps_eq_n_neg ! previous eps_pl_neg, negative local eq. pl. strain
! C          
! C** Update stress
       pl_stress(1) = stress_tr(1,1)
       pl_stress(2) = stress_tr(2,2)
       pl_stress(3) = stress_tr(3,3)
       pl_stress(4) = stress_tr(1,2)
       pl_stress(5) = stress_tr(1,3)
       pl_stress(6) = stress_tr(2,3)
! C **
	   ! DMAX = max(DMAX, max(d_m, d_mOLD))
          ! CALL updateDTIME(DTIME, DMAX)	  		 
! C** 	  
! C** Update history variables
       hist_new(1)=eps_pl_n2D(1,1)
       hist_new(2)=eps_pl_n2D(2,2)
       hist_new(3)=eps_pl_n2D(3,3)
       hist_new(4)=2.d0*eps_pl_n2D(1,2)
       hist_new(5)=2.d0*eps_pl_n2D(1,3)
       hist_new(6)=2.d0*eps_pl_n2D(2,3)
       hist_new(7) = eps_eq_n
       hist_new(8) = eps_eq_n_neg
       hist_new(9) = kappa_str
       hist_new(10) = kappa_str_neg
       hist_new(11) = d_m ! old idi 
       hist_new(12) = d_m_neg
! C hist_new(13) and hist_new(14) will be previous step non-local strains in the next increment.	  
       hist_new(13) = eps_eq_pl_ip ! Positive eqivalent plastic strain from this step 
       hist_new(14) = eps_eq_pl_ip_neg ! Negative eqivalent plastic strain from this step		
       do i=1,12
        hist_new(i) = hist_n(i)
       enddo		   
      elseif(phi_tr.GE.0.d0)then
! C **
! C ** PLASTIC CORRECTOR
! C ** SIMPLY MEANS, ELASTIC PREDICTION IS NOT CORRECT.
! C ** THERE IS EVOLUTION IN PLASTIC STRAIN.
! C ** NEW PLASTIC STRAIN AND RELATED TERMS (STRESS, TANGENT etc.)
! C ** MUST BE COMPUTED. ALL THE TRIAL VALUES MUST BE MODIFIED.
! C **  
       CALL plasticCorrector(nr_stres_dev_tr, mu, pres_tr, alpha,
     1  kappa, kk, stress_dev_tr, I_1, a_t, b_t, c_t, d_t,a_c,b_c,c_c,
     2  d_c, eps_eq_n, eps_pl_n2D, phi_tr, delta_eps_p, delta_eps_eq,
     3  eps_eq_new, sigma_c, sigma_t, phi_pr, gamma_cur, eps_pl_new2D,
     4  stress_new, stress_dev_new, pres_new, d_m, PNEWDT,
     5  fact_neg_stress, eps_eq_n_neg, eps_eq_new_neg, c_pos, c_neg, 
     6  d_m_neg, KINC,
     7  delta_eps_eq_neg, delta_eps_eq_pos, DTIME)  
!	   dg_dEps_V = 0.d0 ! those are possibly wrong..
!	   dgNeg_dEps_V = 0.d0
! C Compute driving forces for 2nd and 3rd PDEs
	   eps_eq_drive = eps_eq_new
	   eps_eq_drive_neg = eps_eq_new_neg 
! C                  
! C Plastic strain tensor in Voigt notation
       eps_pl_new(1) = eps_pl_new2D(1,1)
       eps_pl_new(2) = eps_pl_new2D(2,2)
       eps_pl_new(3) = eps_pl_new2D(3,3)
       eps_pl_new(4) = 2.d0*eps_pl_new2D(1,2)
       eps_pl_new(5) = 2.d0*eps_pl_new2D(1,3)
       eps_pl_new(6) = 2.d0*eps_pl_new2D(2,3)
C          
	   elastic_strain = total_strain-eps_pl_new
C          
       CALL get_stress_eff(E, nu, mu, elastic_strain, stress_eff)
C          
       pl_stress(1)=stress_new(1,1)
       pl_stress(2)=stress_new(2,2)
       pl_stress(3)=stress_new(3,3)
       pl_stress(4)=stress_new(1,2)
       pl_stress(5)=stress_new(1,3)
       pl_stress(6)=stress_new(2,3) 

!C** UPDATE ELEMENT STIFFNESS MATRIX	
C** Compute terms required in the linearization part (Element subroutines)    
C** Stress Linearization terms
C** dStress_dStran, C1
C** dStress_dEpsBar_V, C2  
C** dStress_dEpsBar_neg_V, C2_neg
C**   
C** Local Positive Equivalent plastic strain linearization  
C** dEpsPl_dEps_V, C3, D3   
C** dEpsPl_dEpsPl_bar, C4, D4  
C** dEpsPl_dEpsPl_neg_bar, C5, D5  
C**
C** Local Positive Equivalent plastic strain linearization  
C** dEpsPl_neg_dEps_V, C3_neg, D3_neg 
C** dEpsPl_neg_dEpsPl_bar, C4_neg, D4_neg
C** dEpsPl_neg_dEpsPl_neg_bar, C5_neg, D5_neg                   
C    
       CALL get_mat_tang(d_m, stress_dev_new, kappa, mu, 
     1  phi_pr, stress_eff, sigma_c, sigma_t, kappa_str, 
     2  gamma_cur, alpha, kk, stress_new, pres_new, dStress_dStran,
     3  dStress_dEpsBar_V, total_strain2D, eps_pl_new2D, 
     4  delta_eps_p,dEpsPl_dEps_V, dEpsPl_dEpsPl_bar, eps_eq_new,
     5  eps_eq_n, aa, bb, fact_neg_stress,
     6  d_m_neg, c_neg, c_pos, dStress_dEpsBar_neg_V, 
     7  dEpsPl_neg_dEps_V, dEpsPl_dEpsPl_neg_bar, eps_eq_new_neg,
     8  dEpsPl_neg_dEpsPl_bar, dEpsPl_neg_dEpsPl_neg_bar, eps_eq_n_neg, kappa_str_neg,
     9  pl_stress, d_mOLD, d_mOLD_neg, kappa_old, kappa_old_neg, total_strain) ! Update tangent
!C** END OF THE UPDATE ELEMENT STIFFNESS MATRIX	
	   ! DMAX = max(DMAX, max(d_m, d_mOLD))
          ! CALL updateDTIME(DTIME, DMAX)		  
! Store new damage values as history variables	 
       hist_new(11) = d_m
       hist_new(12) = d_m_neg  			  

!C** Finally update history variables.
       hist_new(1)=eps_pl_new2D(1,1)
       hist_new(2)=eps_pl_new2D(2,2)
       hist_new(3)=eps_pl_new2D(3,3)
       hist_new(4)=2.d0*eps_pl_new2D(1,2)
       hist_new(5)=2.d0*eps_pl_new2D(1,3)
       hist_new(6)=2.d0*eps_pl_new2D(2,3)
       hist_new(7) = eps_eq_new
       hist_new(8) = eps_eq_new_neg
       hist_new(9) = kappa_str
       hist_new(10) = kappa_str_neg   
! C hist_new(13) and hist_new(14) will be previous step non-local strains in the next increment.	  
       hist_new(13) = eps_eq_pl_ip ! Positive eqivalent plastic strain from this step 
       hist_new(14) = eps_eq_pl_ip_neg ! Negative eqivalent plastic strain from this step		
      else ! ERROR STATE, covers NaN errors
C** If the increment is large, then some NaN errors are obtained.
C** With this block, in such cases step size is reduced and
C** analysis step back to last converged step. However this 
C** time size size is DT*PNEWDT.
C** DT-->Previous step size          
C** PNEWDT-->Current step size = Previous step size * PNEWDT         
       print*, 'KINC', KINC
       print*, 'd_flag', d_flag
       print*, 'total_strain', total_strain
       print*, 'eps_eq_pl_ip', eps_eq_pl_ip
       print*, 'd_m', d_m
       print*, 'NaN error.'
       print*, 'eps_eq_new', eps_eq_new
       print*, 'eps_eq_new_neg', eps_eq_new_neg
	   print*, 'ERR 1'
C         analysis stopped, increment is reduced bu 0.001
C         with this new increment continue from the LAST CONVERGED step
C         check log file for this condition
		 CALL XIT()	
       endif ! ELASTIC PREDICTOR && PLASTIC CORRECTOR 
!       CALL updateDTIME(DTIME, DMAX)	   
      end subroutine get_stress
C
C
C
      subroutine getQandSigma(eps_eq_n, eps_eq_n_neg, a_t, b_t, c_t, d_t,
     1   a_c, b_c, c_c, d_c, d_m, d_m_neg, sigma_c, sigma_t )
C**
C** SUBROUTINE UPDATES YIELDS STRESSES IN TENSION AND COMPRESSION     
C**
        implicit none
C INPUT        
        real*8, intent(in) :: eps_eq_n, a_t, b_t, c_t, d_t, a_c,
     1 b_c, c_c, d_c, d_m, d_m_neg, eps_eq_n_neg
C OUTPUT     
        real*8, intent(out) ::sigma_c, sigma_t
C Following the work of mediavilla, yield stresses are reduced by (1-d)
C tension --> (1-dt), compression --> (1-dc)
        sigma_t = (1.d0-d_m)*
     1  (a_t*EXP(b_t*eps_eq_n)+c_t*EXP(d_t*eps_eq_n))
        sigma_c = (1.d0-d_m_neg)*
     1  (a_c*EXP(b_c*eps_eq_n_neg)+c_c*EXP(d_c*eps_eq_n_neg))         
      endsubroutine getQandSigma
C
C
C
      subroutine dSigma_dD(eps_eq_n, a_t, b_t, c_t, d_t,
     1   a_c, b_c, c_c, d_c, d_m, dSt_dD, dSc_dD, eps_eq_n_neg )
!CCCC SUBROUTINE RETURNS derivatives of yield stresses w.r.t. Damage   
C tension --> d()/d(Dt), compressin --> d()/d(Dc)  
        implicit none
C INPUT        
        real*8, intent(in) :: eps_eq_n, a_t, b_t, c_t, d_t, a_c,
     1 b_c, c_c, d_c, d_m, eps_eq_n_neg
C OUTPUT     
        real*8, intent(out) ::dSt_dD, dSc_dD
        dSt_dD = -1.0d0*
     1  (a_t*EXP(b_t*eps_eq_n)+c_t*EXP(d_t*eps_eq_n))
        dSc_dD = -1.0d0*
     1  (a_c*EXP(b_c*eps_eq_n_neg)+c_c*EXP(d_c*eps_eq_n_neg)  )         
      endsubroutine dSigma_dD      
C
C
C
      subroutine getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr, d_m, c_pos, c_neg,
     3 eps_eq_new_neg, d_m_neg)
C**
C** SUBROUTINE RETURNS PHI_PR VALUE (dPHI_dGamma)
C**
        implicit none
C INPUT        
        real*8, intent(in) :: nr_stres_dev_tr, mu, gamma_cur, 
     1  pres_tr, alpha, kappa, kk, delta_eps_p(3,3), delta_eps_eq,
     2  stress_dev_tr(3,3), I_1(3,3), a_t, b_t, c_t, d_t, a_c,
     3  b_c, c_c, d_c, eps_eq_new, sigma_c, sigma_t, d_m, c_pos,
     4  c_neg, eps_eq_new_neg, d_m_neg
        integer, intent(in) :: iter_num
C OUTPUT        
        real*8, intent(out) :: phi_pr        
C LOCAL      
        real*8 dStressDev_dDeltaGamma(3,3), dStressPress_dDeltaGamma, 
     @  dSt_dEpsEq_pos, dSc_dEpsEq_neg, const1,dEpsEq_dDeltaEps(3,3), 
     @  stress_dev_new(3,3), pres_new, dDeltaEps_dGamma(3,3),
     @  dSt_dGamma, dSc_dGamma, dEpsEq_pos_dEpsEq, dEpsEq_neg_dEpsEq
C     
        integer i,j,k,l
C
	   dStressDev_dDeltaGamma = -(1.d0/6.d0)*mu*
     @  stress_dev_tr / ((mu*gamma_cur + (1.d0/6.d0))**2.d0) 
C
        dStressPress_dDeltaGamma = 
     @    -0.5d0*kappa*alpha*pres_tr / 
     @    (kappa*alpha*gamma_cur + 0.5d0)**2.d0
C
        dSt_dEpsEq_pos = (1.d0-d_m)* ! d(sigma_t)/(dEpsEq_pos)
     @    (a_t*b_t*EXP(b_t*eps_eq_new) + 
     @    c_t*d_t*EXP(d_t*eps_eq_new))
C
        dSc_dEpsEq_neg = (1.d0-d_m_neg)*
     @    (a_c*b_c*EXP(b_c*eps_eq_new_neg) +
     @    c_c*d_c*EXP(d_c*eps_eq_new_neg))
C
        const1 = 0.d0
        do i=1,3
          do j=1,3
            const1 =  const1 + kk*delta_eps_p(i,j)*delta_eps_p(i,j)
          enddo
        enddo
        const1 = SQRT(const1) ! const1 is __WHOLE__ eps_eq_new    
C
C** Prevents division by zero errors
		dEpsEq_dDeltaEps = 0.d0
        if(iter_num.EQ.1.OR.abs(const1).LT.1e-12)then
		  dEpsEq_dDeltaEps = 0.d0 ! 3x3 luk zero matrisi yapar.
!		  https://stackoverflow.com/questions/43433004/is-there-an-intrinsic-function-for-initializing-arrays-to-zero-in-fortran
        else
		  dEpsEq_dDeltaEps = kk*delta_eps_p/const1
        endif          
C
		stress_dev_new = stress_dev_tr/(1.d0+6.d0*mu*gamma_cur)
C
        pres_new = pres_tr/(1.d0+2.d0*gamma_cur*alpha*kappa)
C
		dDeltaEps_dGamma = 
     @    3.d0*stress_dev_new + (2.d0/3.d0)*alpha*pres_new*I_1	 

!        if(ABS(gamma_cur).GT.1e-6)then
C** if gamma_cur is > 0, some contribution must come from this term
C** otherwise ignore the effect
          dDeltaEps_dGamma = dDeltaEps_dGamma +
     @  gamma_cur * (3.d0*dStressDev_dDeltaGamma + 
     @  (2.d0/3.d0) * alpha * dStressPress_dDeltaGamma * I_1)	
	 
!        endif
C
        dSt_dGamma = 0.d0
        dSc_dGamma = 0.d0
C
        dEpsEq_pos_dEpsEq = c_pos
        dEpsEq_neg_dEpsEq = c_neg
        do i=1,3
          do j=1,3
            dSt_dGamma = dSt_dGamma + 
     @  dSt_dEpsEq_pos * dEpsEq_pos_dEpsEq *
     @  dEpsEq_dDeltaEps(i,j) * dDeltaEps_dGamma(i,j)
            dSc_dGamma = dSc_dGamma +
     @  dSc_dEpsEq_neg * dEpsEq_neg_dEpsEq * 
     @  dEpsEq_dDeltaEps(i,j) * dDeltaEps_dGamma(i,j)     
          enddo
        enddo
C
C** Finally sum everything up to calculate d(phi)/d(Gamma)
        phi_pr = 0.d0
        do i=1,3
          do j=1,3
            phi_pr =  phi_pr + 
     @  6.d0*dStressDev_dDeltaGamma(i,j)*stress_dev_new(i,j)
          enddo
        enddo
        phi_pr = phi_pr + 
     @   6.d0*dStressPress_dDeltaGamma*(sigma_c-sigma_t) +
     @   6.d0*pres_new*(dSc_dGamma - dSt_dGamma) -
     @   2.d0*dSc_dGamma*sigma_t -
     @   2.d0*sigma_c*dSt_dGamma
      end subroutine getPhiPr
C
C
C      
      subroutine getTrial(stran_tot, eps_pl_n2D, PP, I_1, mu, kappa,
     1  eps_eq_n, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, 
     2  stress_dev_tr, nr_stres_dev_tr, pres_tr, stress_tr, 
     3  sigma_c, sigma_t, phi_tr, d_m, d_m_neg, eps_eq_n_neg)
C**
C** CACLUALTES TRIAL STRESS and PHI_TRIAL
C** PERFORMS CALCULATIONS FOR ELASTIC PREDICTOR
C**     
        implicit none
C** INPUT
        real*8, intent(in) ::  stran_tot(3,3), eps_pl_n2D(3,3), 
     1   PP(3,3,3,3), I_1(3,3), mu, kappa, eps_eq_n, a_t,
     2   b_t, c_t, d_t, a_c, b_c, c_c, d_c, d_m, d_m_neg,
     3   eps_eq_n_neg
C** OUTPUT
        real*8, intent(out) ::  stress_dev_tr(3,3), nr_stres_dev_tr, 
     1  pres_tr, stress_tr(3,3), sigma_c, sigma_t, phi_tr
C** LOCAL VARIABLES
        real*8  eps_el_tr(3,3), eps_el_vol_tr, eps_el_dev_tr(3,3)
        integer i,j,k,l,m,n
C** CODE

		eps_el_tr = stran_tot - eps_pl_n2D
C		
        eps_el_vol_tr = eps_el_tr(1,1)+eps_el_tr(2,2)+eps_el_tr(3,3)
C
		eps_el_dev_tr = 0.d0 ! 3x3
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            eps_el_dev_tr(i,j) = eps_el_dev_tr(i,j) +
     1       PP(i,j,k,l) * eps_el_tr(k,l)
           enddo
          enddo
         enddo
        enddo
C
		stress_dev_tr =2.d0*mu*eps_el_dev_tr
C
        CALL normOf(stress_dev_tr, nr_stres_dev_tr)
C
        pres_tr = kappa*eps_el_vol_tr
C
		stress_tr = stress_dev_tr + I_1*pres_tr
C
        CALL getQandSigma(eps_eq_n, eps_eq_n_neg, a_t, b_t, c_t, d_t,
     1    a_c, b_c, c_c, d_c, d_m, d_m_neg, sigma_c, sigma_t )
C
        CALL getPhi(nr_stres_dev_tr,pres_tr,sigma_c,sigma_t,phi_tr)
      end subroutine getTrial
C
C
C
      subroutine getPhi(nr_stres_dev_tr,pres_tr,sigma_c,sigma_t,phi_tr)
        implicit none
C**
C** CALCULATES YIELD CRITERION (PHI)
C**        
        real*8, intent(in)::nr_stres_dev_tr,pres_tr,sigma_c,sigma_t
        real*8, intent(out)::phi_tr
C
        phi_tr = 3.0d0*(nr_stres_dev_tr**2.0d0) + 6.0d0*pres_tr*
     @ (sigma_c - sigma_t) - 2.0d0*sigma_c*sigma_t  
      end subroutine getPhi
C
C
C
      subroutine plasticCorrector(nr_stres_dev_tr, mu, pres_tr, 
     1  alpha, kappa, kk, stress_dev_tr, I_1, a_t, b_t, c_t,
     2  d_t, a_c, b_c, c_c, d_c, eps_eq_n, eps_pl_n2D, phi_tr,
     3  delta_eps_p, delta_eps_eq, eps_eq_new, sigma_c, sigma_t,
     4  phi_pr, gamma_cur, eps_pl_new2D, stress_new, stress_dev_new,
     5  pres_new, d_m, PNEWDT, fact_neg_stress, 
     6  eps_eq_n_neg, eps_eq_new_neg, c_pos, c_neg, d_m_neg,
     7  KINC,
     8  delta_eps_eq_neg, delta_eps_eq_pos, DTIME)
        implicit none
! C**
! C** USING RETURN MAPPING ALGORITHM, CALCULATES PLASTIC STRAIN
! C** AND UPDATES STRESSES
! C**        
! CCCC INPUT
        real*8, intent(in) ::  nr_stres_dev_tr, mu, pres_tr, alpha,
     1    kappa, kk, stress_dev_tr(3,3), I_1(3,3), a_t, b_t,
     2    c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_n, eps_pl_n2D(3,3),
     3    phi_tr, d_m, fact_neg_stress, eps_eq_n_neg, d_m_neg
! CCCC OUTPUT
        real*8, intent(out) ::  delta_eps_p(3,3), delta_eps_eq, 
     1   eps_eq_new, sigma_c, sigma_t, phi_pr, gamma_cur,
     2   eps_pl_new2D(3,3), stress_new(3,3), stress_dev_new(3,3),
     3   pres_new, eps_eq_new_neg, c_pos, c_neg,
     4   delta_eps_eq_neg, delta_eps_eq_pos
! C     
        integer, intent(in) :: KINC
! C        
        real*8, intent(inout) :: PNEWDT, DTIME
! CCCC LOCAL VAR.
        real*8 diff, TOL, phi, nr_stres_dev_new, eps_pl_new_v(6)
! CCCC Principal strain calculation values
        real*8 PS(3), AN(3,3), delta_eps_p_v(6), PSpos(3), PSneg(3) 
!        real*8 delta_eps_eq_neg, delta_eps_eq_pos ! may be changed to output later  
        real*8 delta_eps_p_pos(3,3), delta_eps_p_neg(3,3)   
        integer iter_num, iter_max, i, j, k, l, m, n
! CCCC Modifications in delta_eps_eq_neg and delta_eps_eq_neg  
		real*8 S1, S2, S3, S_H, norm_stress
		real*8 sig_pr(3), sig_cos(3,3), s1_pr, s2_pr, s3_pr,
     1   dg_dS1, dg_dS2, dg_dS3, n1(3), n2(3), n3(3),
     2   dS1pr_dStress(3,3), dS3pr_dStress(3,3), dS2pr_dStress(3,3),
     3   Cel_4D(3,3,3,3), dg_dEps(3,3)	 
	    real*8 dgNeg_dEps(3,3), counter_1
		integer s1_ind, s2_ind, s3_ind
! CCCC CODE
        diff = 1.0d0
        TOL = 1e-6
		TOL = 1e-3
        ! TOL = 1e-9
        iter_num = 1    
        iter_max = 50000     
        gamma_cur = 0.d0
        c_pos = 1.0d0
        c_neg = 1.0d0
! with plasticity DTIME is 1e-3		
!		DTIME = 1e-3
! C        
        delta_eps_eq = 0.d0
        eps_eq_new = eps_eq_n
        eps_eq_new_neg = eps_eq_n_neg
        phi = phi_tr
! C
! C Define sigma_c and sigma_t values
        CALL getQandSigma(eps_eq_new, eps_eq_new_neg, a_t, b_t, c_t, d_t,
     1   a_c, b_c, c_c, d_c, d_m, d_m_neg, sigma_c, sigma_t )        
! C
	    stress_dev_new = 0.0d0! new deviatoric stress
		delta_eps_p = 0.d0 ! increment in plastic strain (Delta_eps_p)
! C
		counter_1 = 1.d0
        do while(diff.GT.TOL.AND.iter_num.LT.iter_max) ! local N-R
		 counter_1 = counter_1 + 1.d0
         CALL getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr, d_m, c_pos, c_neg, 
     4 eps_eq_new_neg, d_m_neg)!COMPUTES D(PHI)/D(GAMMA)
! SABAH TEST ET..
	    ! if(ABS(phi_pr)<1e-6)then ! increases numerical stability. 1e-3 is too big ??
		 ! phi_pr = 1e-6
		! endif 	
		! phi_pr = phi_pr + 1e-3
! C
! C** Typcial N-R, x = x - _ALPHA_ * f(x)/f'(x)
         gamma_cur = gamma_cur - phi/phi_pr
		 if(gamma_cur.LT.0.d0)then
		  gamma_cur = 1e-2
		 endif 
! TODO: Use secant or bi-section to find the root
! secant first
		 
! C**
! C** Update all stresses and plastic strains 
! C**
         stress_dev_new = stress_dev_tr / (1.d0+6.d0*mu*gamma_cur)		  
! C
         pres_new = pres_tr / (1.d0+2.d0*gamma_cur*alpha*kappa)
! C** Update stress
         stress_new = stress_dev_new + I_1*pres_new				
! C		
         delta_eps_p= gamma_cur *
     1 (3.d0*stress_dev_new+(2.d0/3.d0)*alpha*pres_new*I_1)		  
! C* -------------------------------------------------------         
! C** To distinguish tension and compression forces
! C** principal strain values of delta_eps_p must be computed
! C** To use built-in function (SPRIND), delta_eps_p must 
! C** be put into Voigt form
! C** For the next version, distinguishing damage will be based on the hydraustatic stress
! C** IF s_h < 0 --> Compression case, else --> Tension case
         delta_eps_p_v(1)=delta_eps_p(1,1)
         delta_eps_p_v(2)=delta_eps_p(2,2)
         delta_eps_p_v(3)=delta_eps_p(3,3)
         delta_eps_p_v(4)=2.d0*delta_eps_p(1,2)
         delta_eps_p_v(5)=2.d0*delta_eps_p(1,3)
         delta_eps_p_v(6)=2.d0*delta_eps_p(2,3)
! C* SPRIND finds Principal quantities (PS) and related directions (AN)     
         CALL SPRIND(delta_eps_p_v, PS, AN, 2, 3, 3)
! C
         PSpos = 0.0d0 ! Positive part of delta_eps_p_v
         PSneg = 0.0d0 ! Negative part of delta_eps_p_v
! C
         do i=1,3
          if(PS(i).GE.0.d0)then
           PSpos(i) = PS(i)
          else
           PSneg(i) = PS(i)
          endif
         enddo
! C
! delta_eps_eq is incremental equivalent plastic strain regarding whole delta_eps_p_v
         delta_eps_eq = SQRT(kk*(PS(1)**2.d0 + PS(2)**2.d0 + PS(3)**2.d0))
! C
! C** delta_eps_eq_pos : POSITIVE incremental equivalent plastic strain
! C** delta_eps_eq_neg : NEGATIVE incremental equivalent plastic strain
! C** delta_eps_eq_pos and delta_eps_eq_neg are defined first
! C** then c_pos and c_neg defined, app 0
         delta_eps_eq_pos = SQRT(kk*(PSpos(1)**2.d0 + PSpos(2)**2.d0 + PSpos(3)**2.d0))
         delta_eps_eq_neg = SQRT(kk*(PSneg(1)**2.d0 + PSneg(2)**2.d0 + PSneg(3)**2.d0))	  
! C
! C 
! MOD: G_TENS and G_COMP eklendi
         eps_eq_new = eps_eq_n + delta_eps_eq_pos  ! positive driving force
		 eps_eq_new_neg = eps_eq_n_neg + delta_eps_eq_neg ! negative diriving force
! **
! ** c_pos and c_neg will be used in the calculation of the material tangent matrix (DDSDDE)
         c_pos = delta_eps_eq_pos / delta_eps_eq
         c_neg = delta_eps_eq_neg / delta_eps_eq		  
! C
! C** Update sigma_t and sigma_c with new incremental equivalent plastic strain values
         CALL getQandSigma(eps_eq_new, eps_eq_new_neg, a_t, b_t, c_t, d_t,
     1  a_c, b_c, c_c, d_c, d_m, d_m_neg, sigma_c, sigma_t)
! C
! C** Compute Phi
         CALL normOf(stress_dev_new, nr_stres_dev_new)
         CALL getPhi(nr_stres_dev_new,pres_new,sigma_c,sigma_t,phi)
! C
		 diff = ABS(phi)
         iter_num = iter_num + 1
        enddo !end of local N-R    
! compute phi_pr with new gamma_cur        
		CALL getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr, d_m, c_pos, c_neg, 
     4 eps_eq_new_neg, d_m_neg)!COMPUTES D(PHI)/D(GAMMA)		
! C**
! C** Check for errors, Basically try-catch block     
! C**
        if(ABS(diff).GT.TOL)then ! FIX: The case where DIFF>TOL is caught
		 print*, 'Diff is large.'
		 ! print*, 'Diff:', diff
		 ! print*, 'Increment:', KINC
         ! print*, 'eps_eq_new', eps_eq_new
         ! print*, 'eps_eq_new_neg', eps_eq_new_neg
         ! print*, 'dm:', d_m
		 ! print*, 'G_TENS:', G_TENS
		 ! print*, 'G_COMP:', G_COMP	
         ! print*, 'iter_num', iter_num		
		 ! print*, 'phi_pr', phi_pr		 
		 ! print*, 'gamma_cur', gamma_cur
		 ! print*, 'Stopping...'		 
		 ! PNEWDT = 1e-2
		 ! CALL XIT
		 PNEWDT = 1e-2
		endif 
		
        if(ABS(phi_pr).LT.1e-18)then
          print*, 'phi_pr is almost zero...'
		  print*, 'Increment:', KINC		  
          print*, 'phi_pr:', phi_pr
          PNEWDT = 1e-3
		  DTIME = DTIME / 10.d0
! phi_pr = phi_pr + 1e-5 yapsan olur ??		  
          CALL XIT()
        endif
! C
        ! if(gamma_cur.LT.0.d0.AND.KINC.GT.1)then
          ! print*, 'gamma_cur is negative', gamma_cur
          ! PNEWDT = 0.05d0
          ! CALL XIT()
        ! endif
! C
        if(iter_num.EQ.iter_max)then
          print*, '*** EROOR LOCAL NR 1***'
		  ! print*, 'Increment:', KINC		  
          ! print*, 'This is diff:', diff
          ! print*, 'iter_num:', iter_num
          ! print*,'Too much iteration..'
          ! print*, 'eps_eq_new', eps_eq_new
          ! print*, 'eps_eq_new_neg', eps_eq_new_neg
          ! print*, 'dm:', d_m
		  ! print*, 'G_TENS:', G_TENS
		  ! print*, 'G_COMP:', G_COMP		 
          PNEWDT = 1e-2
! if we come up with such error, then increment (DTIME) is reduced by 1e-3 (PNEWDT). 
! Afterwads, analysis starts with this new increment from the last converged increment.		  
          ! CALL XIT()
        endif
! C**        
! C** Regarding gamma value compute new stresses and delta_eps_p
! C** Update Plastic strain
! C**
        eps_pl_new2D = eps_pl_n2D + delta_eps_p	
! C** Update plastic strain
        eps_pl_new_v(1)=eps_pl_new2D(1,1)
        eps_pl_new_v(2)=eps_pl_new2D(2,2)
        eps_pl_new_v(3)=eps_pl_new2D(3,3)
        eps_pl_new_v(4)=2.d0*eps_pl_new2D(1,2)
        eps_pl_new_v(5)=2.d0*eps_pl_new2D(1,3)
        eps_pl_new_v(6)=2.d0*eps_pl_new2D(2,3)   
! compute d_Gplus_d_sigma
! s1_pr, s2_pr, s3_pr, dg_dS1, dg_dS2, dg_dS3, n1(3), n2(3), n3(3),
! dS1pr_dStress(3,3), dS3pr_dStress(3,3), dS2pr_dStress(3,3),		
      end subroutine plasticCorrector
! C
! C
! C
      subroutine convert2FullStrain(strainVoigt, strainFull)
          implicit none
C INPUT
          real*8, intent(in) ::  strainVoigt(6)
C OUTPUT
          real*8, intent(out) ::  strainFull(3,3)
          integer i, j
C CODE
		strainFull = 0.d0
C
        strainFull(1,1)=strainVoigt(1)
        strainFull(2,2)=strainVoigt(2)
        strainFull(3,3)=strainVoigt(3)
        strainFull(1,2)=0.5d0*strainVoigt(4)
        strainFull(1,3)=0.5d0*strainVoigt(5)
        strainFull(2,3)=0.5d0*strainVoigt(6)
        strainFull(2,1)=strainFull(1,2)
        strainFull(3,1)=strainFull(1,3)
        strainFull(3,2)=strainFull(2,3)
      end subroutine convert2FullStrain
C
C
C
      subroutine fullStress2Voigt(stress_full, stress_voigt)
          implicit none
C INPUT
          real*8, intent(in) ::  stress_full(3,3)
C OUTPUT
          real*8, intent(out) ::  stress_voigt(6)
C CODE
          stress_voigt(1) = stress_full(1,1)
          stress_voigt(2) = stress_full(2,2)
          stress_voigt(3) = stress_full(3,3)
          stress_voigt(4) = stress_full(1,2)
          stress_voigt(5) = stress_full(1,3)
          stress_voigt(6) = stress_full(2,3)
      end subroutine fullStress2Voigt   
C
C
C      
      subroutine get_D_abaqus(E, nu, mu, d_m, D, d_m_neg)
C**
C** Material tanget matrix as defined in Melro et al. 
C**  ABAQUS documentation ile aynı sonucu veriyor   
C**     
        implicit none
! INPUT
        real*8, intent(in) :: E, nu, mu, d_m, d_m_neg  
! OUTPUT
        real*8, intent(out) ::  D(6,6)
! LOCAL
        real*8 t1, t2, t3, Gd, Ld
        integer i, j
C
        t1 = E*(1.d0-d_m)*(1.0d0-d_m_neg)*
     1    (1.d0-nu*(1.d0-d_m)*(1.0d0-d_m_neg))
        t2 = (1.d0+nu*(1.d0-d_m)*(1.0d0-d_m_neg))*
     1    (1.d0-2.d0*nu*(1.d0-d_m)*(1.0d0-d_m_neg))
        t3 = E*nu*(((1.d0-d_m)*(1.0d0-d_m_neg))**2.d0)
C
        Gd = t1/t2
        Ld = t3/t2
C        
		D = 0.d0
C        
        D(1,1) = Gd
        D(2,2) = Gd
        D(3,3) = Gd
        D(4,4) = mu*(1.d0-d_m)*(1.0d0-d_m_neg)
        D(5,5) = mu*(1.d0-d_m)*(1.0d0-d_m_neg)
        D(6,6) = mu*(1.d0-d_m)*(1.0d0-d_m_neg)
        D(1,2)= Ld
        D(1,3)= Ld
        D(2,3)= Ld
        D(2,1)=D(1,2)
        D(3,1)=D(1,3)
        D(3,2)=D(2,3)                
      end subroutine get_D_abaqus  
! C Above is not used.	  
      subroutine get_D_abaqus2(E, nu, mu, D)
C**
C** ABAQUS documentation ile aynı sonucu veriyor  
C** If nu12, nu13 and nu23 will be used then this will be ready 
C**     
        implicit none
! INPUT
        real*8, intent(in) :: E, nu, mu 
! OUTPUT
        real*8, intent(out) ::  D(6,6)
! LOCAL
        real*8 TT
		integer i, j
			

		D = 0.d0
		
		TT = 1.d0 / (1.d0 - 3.d0 * nu**2.d0 - 2.d0 * nu**3.d0)	
		D(1,1) = E * (1.d0 - nu**2.d0)*TT	
		D(2,2) = E * (1.d0 - nu**2.d0)*TT	
		D(3,3) = E * (1.d0 - nu**2.d0)*TT	
		D(1,2) = E * (nu + nu**2.d0)*TT
		D(1,3) = E * (nu + nu**2.d0)*TT
		D(2,3) = E * (nu + nu**2.d0)*TT		
		D(4,4) = mu
		D(5,5) = mu
		D(6,6) = mu		
		
		do i=1,6
		  do j=1,6
		    D(j,i) = D(i,j)
		  enddo
		enddo	
      end subroutine get_D_abaqus2   	  
C
C
C      
      subroutine get_stress_eff(E, nu, mu, elastic_strain, stress_eff)    
C**
C** RETURN EFFECTIVE STRESS (STRESS WITHOUT DAMAGE TERMS)   
C**     
        implicit none
! INPUT
        real*8, intent(in) :: E, nu, mu, elastic_strain(6) 
! OUTPUT
        real*8, intent(out) ::  stress_eff(6)
! LOCAL
        real*8 D_el(6,6)
        integer i, j
C
        ! CALL get_D_abaqus(E, nu, mu, 0.d0, D_el, 0.0d0) ! no damage--> Elastic tangent tensor
		CALL get_D_abaqus2(E, nu, mu, D_el)	
C     
		stress_eff = 0.d0 
C             
		stress_eff = MATMUL(D_el, elastic_strain)
      end subroutine get_stress_eff       
C
C
C      
      subroutine get_mat_tang(d_m, stress_dev, kappa, mu, phi_pr,
     1    stress_eff, sigma_c, sigma_t, kappa_str, gamma_cur,
     2    alpha, kk, stress_tot, stress_pres, dStress_dStran,
     3    dStress_dEpsBar_V, total_strain2D, eps_pl_new2D, delta_eps_p,
     4    dEpsPl_dEps_V, dEpsPl_dEpsPl_bar, eps_eq_new, eps_eq_n,
     5    aa_dam, bb_dam, fact_neg_stress, d_m_neg, c_neg, c_pos,
     6    dStress_dEpsBar_neg_V, dEpsPl_neg_dEps_V, dEpsPl_dEpsPl_neg_bar,
     7    eps_eq_new_neg, dEpsPl_neg_dEpsPl_bar, dEpsPl_neg_dEpsPl_neg_bar, eps_eq_n_neg,
     7    kappa_str_neg, pl_stress, d_mOLD, d_mOLD_neg, kappa_old, kappa_old_neg,
     8    total_strain)
! C** 
! C** CALCULATES MATERIAL TANGENT AND THE PARTIAL DERIVATIVE
! C** TERMS WHICH ARE REQUIRED IN ELEMENT STFIFFNESS TENSOR     
! C** CALCULATIONS ARE IN FULL FORM      
! C** MAIN SUBROUTINE FOR ELEMENT STIFFNESS TENSORS
! C** TODO: 2.0d0 carpanı gelebilecek yerelre bak.
! C**
! C** modifications in OCTOBER 2021
! C** Stress Linearization terms
! C** dStress_dStran, C1 **
! C** dStress_dEpsBar_V, C2  
! C** dStress_dEpsBar_neg_V, C2_neg
! C**   
! C** Local Positive Equivalent plastic strain linearization  
! C** dEpsPl_dEps_V, C3, D3   
! C** dEpsPl_dEpsPl_bar, C4, D4 ** 
! C** dEpsPl_dEpsPl_neg_bar, C5, D5  
! C**
! C** Local Positive Equivalent plastic strain linearization  
! C** dEpsPl_neg_dEps_V, C3_neg, D3_neg 
! C** dEpsPl_neg_dEpsPl_bar, C4_neg, D4_neg
! C** dEpsPl_neg_dEpsPl_neg_bar, C5_neg, D5_neg **
! C**
! C** 
! C** With the new formulation modification C1, C3 and C4_neg parameters must be provided
! C** Other terms may be neglected.
        implicit none
! INPUT
        real*8, intent(in) :: d_m, stress_dev(3,3), kappa, mu, phi_pr,
     1  stress_eff(3,3), sigma_c, sigma_t, kappa_str, gamma_cur, alpha, 
     2  stress_tot(3,3), stress_pres, total_strain2D(3,3),
     3  eps_pl_new2D(3,3), kk, delta_eps_p(3,3), eps_eq_new, eps_eq_n,
     4  aa_dam, bb_dam, fact_neg_stress, d_m_neg, c_neg, c_pos, 
     5  eps_eq_new_neg, eps_eq_n_neg, kappa_str_neg, total_strain(6),
     6  pl_stress(6), d_mOLD, d_mOLD_neg, kappa_old, kappa_old_neg
! OUTPUT
        real*8,intent(out)::dStress_dStran(6,6),dStress_dEpsBar_V(6), 
     1   dEpsPl_dEps_V(6), dEpsPl_dEpsPl_bar, dStress_dEpsBar_neg_V(6), 
     2   dEpsPl_neg_dEps_V(6), dEpsPl_dEpsPl_neg_bar, 
     3   dEpsPl_neg_dEpsPl_bar, dEpsPl_neg_dEpsPl_neg_bar
! LOCAL
        integer i, j, k, l, m, n, p, t ! counters
        real*8 K1(3,3), K2, K3, nr_stress_dev,
     1    I_1(3,3), II(3,3,3,3), PP(3,3,3,3), Cel_4D(3,3,3,3)      
        real*8 p1(3,3), dStressDev_dStress(3,3,3,3), K4,
     1    dStress_dEps(3,3,3,3), press_eff, dev_stress_eff(3,3),
     2    nr_stres_dev, dD_dKappa, dKappa_dEpsPlBar, dPhi_dD,
     3    I1, J2, N_star(3,3), t1(3,3,3,3), NN(3,3,3,3),
     4    dStress_dStran_4D(3,3,3,3), dPress_dStress(3,3),
     5    dStress_dD(3,3), dI1_dEps(3,3)
        real*8 dStress_dEpsBar(3,3), dStressDev_dEps(3,3,3,3)
        real*8 const1, dDEpsEqPl_dDeltaEpsPl(3,3), 
     1    dDeltaEpsPl_dStressDev(3,3,3,3), dStress_dEpsPl(3,3,3,3),
     2    dDeltaEpsPl_dStress(3,3,3,3), AA(3,3), K1_p1(3,3), 
     3    K1_p2(3,3), dDeltaEpsPl_dStressPres(3,3), 
     4    dStressPres_dStress(3,3), dEpsPl_dEps(3,3), 
     5    dPress_dEps(3,3), dPress_dD, dStressDev_dD(3,3),
     6    stress_tot2(3,3), invMat(3,3), dPhi_dStressDev(3,3),
     7    dPhi_dstressPress, dNstar_dStressDev(3,3,3,3), 
     8    dNstar_dStressPress(3,3)
        integer ind_1(6), ind_2(6)
C terms for new K3        
        real*8 a_t, b_t, c_t, d_t,
     1   a_c, b_c, c_c, d_c, dSt_dD, dSc_dD,
     2   dStressPres_dD, const2, const3, const4
C new terms for K1
        real*8 dSigmaT_dEpsEq_pos, dSigmaC_dEpsEq_neg, dEpsEq_dDeltaEpsEqPl, 
     1    dSigmaT_dEps(3,3), dSigmaC_dEps(3,3),K1_p3(3,3), K1_p4(3,3), 
     2    dEpsEq_pos_dEpsEq, dEpsEq_neg_dEpsEq
        integer  p2, t2
C principal stress terms
        real*8 delta_eps_p_v(6), PS(3), AN(3,3)     
C** new tangent terms        
        real*8 dStress_dD_neg(3,3), dStressDev_dD_neg(3,3) , 
     1   dStressPres_dD_neg, dPhi_dD_neg, dKappa_dEpsPlBar_neg,
     2   K3_neg, dD_neg_dKappa, K4_neg, dStress_dEpsBar_neg(3,3),
     3   AA_neg(3,3), dEpsPl_neg_dEps(3,3)
C**        
C**  INITIATION OF MATRICES
C**
        dEpsPl_dEpsPl_bar = 0.d0
        dEpsPl_dEpsPl_neg_bar = 0.d0
        dPress_dEps=0.d0
        K1_p1 = 0.d0
        K1_p2 = 0.d0
        K1_p3 = 0.d0
        K1_p4 = 0.d0
        dSigmaT_dEps=0.d0
        dSigmaC_dEps=0.d0
        K1=0.d0
        dStress_dEpsBar = 0.d0
		dStress_dEpsBar_neg = 0.d0
        dStress_dD=0.d0
        dStress_dD_neg = 0.d0
        dI1_dEps=0.d0
        AA = 0.d0
        AA_neg = 0.d0
        dEpsPl_dEps = 0.d0
        dStressDev_dD = 0.d0
        dStressDev_dD_neg = 0.0d0
        stress_tot2 = 0.d0
		dPress_dStress= 0.0d0
		dDEpsEqPl_dDeltaEpsPl = 0.0d0
		dDeltaEpsPl_dStressPres = 0.0d0
		N_star = 0.0d0
		dEpsPl_neg_dEps=0.d0	
		
       dEpsPl_dEps_V = 0.d0
       dStress_dEpsBar_V = 0.d0
       dStress_dStran = 0.d0		
	   delta_eps_p_v=0.0d0	
	   
       NN = 0.d0
       dStressDev_dEps = 0.d0
       dStressDev_dStress = 0.d0
	   dStress_dEps = 0.0d0
	   dDeltaEpsPl_dStressDev = 0.0d0
	   dStress_dStran_4D = 0.d0
	   dStress_dEpsPl = 0.0d0	   	
! C**        
! C** START FORUMATLION      
! C** CALL UTILITY FUNCTIONS  
        CALL getIandPterms(I_1, II, PP) 
! C** Elastic stiffness tensor in 4D       
        CALL getCel4D(kappa, mu, Cel_4D)    
! C**         
! C** DERIVATIONS FOR K1          
! C**           
		dStressDev_dStress= PP
 
		dStress_dEps= Cel_4D
! C        
        do i=1,3
         do j=1,3 
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3 ! C derivative of deviatoric stress w.r.t. total strain               
			  dStressDev_dEps(i,j,m,n) = dStressDev_dEps(i,j,m,n) +			
     2  dStressDev_dStress(i,j,k,l) * dStress_dEps(k,l,m,n)
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
! C 
! C         
		dPress_dStress = (1.d0/3.d0)*I_1
! C 
! C         
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3 ! derivative of stress w.r.t. total strain
            dPress_dEps(k,l) = dPress_dEps(k,l) +
     1  dPress_dStress(i,j)*dStress_dEps(i,j,k,l)
           enddo
          enddo
         enddo
        enddo
C
CCC Calculate dSigmaC_dEps and dSigmaT_dEps
        CALL stressCurve(a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c)
! C 1. Derivative of yield stress w.r.t. equivalent plastic strain
! C dSigmaC_dEpsEq & dSigmaT_dEpsEq
! C due to damage evolution step size is reduced. For the gauss points
! C where no/very-little damage is formed this derivative may be a problem
! C with this if block such problem is prevented
! C dSigmaC_dEpsEq_neg : d(sigma_neg)/d(EpsPl_neg)
! C dSigmaT_dEpsEq_pos : d(sigma_pos)/d(EpsPl_pos)
! C Since it is stress related D_n values are used. (Constant damage step)
        if(ABS(eps_eq_new-eps_eq_n).GT.1e-12)then ! evolution in positive eq strain
         dSigmaT_dEpsEq_pos = (1.d0-d_m)*
     1    ( a_t*b_t*EXP(b_t*eps_eq_new) + c_t*d_t*EXP(d_t*eps_eq_new) )   
        else ! no evolution in plastic eq strain
         dSigmaT_dEpsEq_pos = 0.0d0
        endif
! C**		
        if(ABS(eps_eq_new_neg-eps_eq_n_neg).GT.1e-12)then ! evolution in negative eq strain
         dSigmaC_dEpsEq_neg = (1.0d0-d_m_neg)*
     1    ( a_c*b_c*EXP(b_c*eps_eq_new_neg) + c_c*d_c*EXP(d_c*eps_eq_new_neg) )
        else ! no evolution in plastic eq strain
         dSigmaC_dEpsEq_neg = 0.0d0
        endif		
! C 2. dEpsEq_dDeltaEpsEqPl
! C derivative of equivalent plastic strain w.r.t.
! C increment in equivalent plastic strain
        dEpsEq_dDeltaEpsEqPl = 1.d0
! C 3. dDEpsEqPl_dDeltaEpsPl (3,3)
! C derivative of increment in equivalent plastic strain w.r.t.
! C increment in plastic strain 
! C Here principal strains are used. Hence derivative must include
! C principal strains. 
! C        
        delta_eps_p_v(1)=delta_eps_p(1,1) ! convert to Voigt form (for principal strain calculations)
        delta_eps_p_v(2)=delta_eps_p(2,2)
        delta_eps_p_v(3)=delta_eps_p(3,3)
        delta_eps_p_v(4)=2.d0*delta_eps_p(1,2)
        delta_eps_p_v(5)=2.d0*delta_eps_p(1,3)
        delta_eps_p_v(6)=2.d0*delta_eps_p(2,3)
        CALL SPRIND(delta_eps_p_v, PS, AN, 2, 3, 3) ! principal values --> PS
! C
        const1 = SQRT(kk*(PS(1)**2.d0 + PS(2)**2.d0 + PS(3)**2.d0)) ! const1: delta_eps_pl   
! C check if const1 is 0.0d0 If it is zero, division by zero case may be seen    
! C        
        if(ABS(const1).LT.1e-12)then
		 dDEpsEqPl_dDeltaEpsPl = 0.d0
        else
		 dDEpsEqPl_dDeltaEpsPl = kk*delta_eps_p/const1	
        endif
! C       
! C 4. dDeltaEpsPl_dStressDev (3,3,3,3)
! C derivative of increment in plastic strain w.r.t
! C deviatoric stress      
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dDeltaEpsPl_dStressDev(i,j,k,l) = 
     1  3.d0*gamma_cur*I_1(i,k)*I_1(j,l)
           enddo
          enddo
         enddo
        enddo
! C        
! C 7. dDeltaEpsPl_dStressPres (3,3)
! C derivative of increment in plastic strain w.r.t pressure 
! C      
	    dDeltaEpsPl_dStressPres = (2.d0/3.d0)*gamma_cur*alpha*I_1
! C 9. dSigmaT_dEps (3,3)
        dEpsEq_pos_dEpsEq = c_pos ! since dEps_EQ_pl_pos = dEps_EQ_pl * c_pos
        dEpsEq_neg_dEpsEq = c_neg ! since dEps_EQ_pl_neg = dEps_EQ_pl * c_neg
! C      
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              do p=1,3
               do t=1,3
                dSigmaT_dEps(p,t) = dSigmaT_dEps(p,t) + dSigmaT_dEpsEq_pos *
     1           dEpsEq_pos_dEpsEq * dEpsEq_dDeltaEpsEqPl * 
     2           dDEpsEqPl_dDeltaEpsPl(i,j) * dDeltaEpsPl_dStressDev(i,j,k,l) *
     3           dStressDev_dStress(k,l,m,n) * dStress_dEps(m,n,p,t)     
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo 
         enddo
        enddo 
! C 
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              dSigmaT_dEps(m,n) = dSigmaT_dEps(m,n) + dSigmaT_dEpsEq_pos *
     1         dEpsEq_pos_dEpsEq * dEpsEq_dDeltaEpsEqPl * 
     2         dDEpsEqPl_dDeltaEpsPl(i,j) * dDeltaEpsPl_dStressPres(i,j)*
     3         dPress_dStress(k,l) * dStress_dEps(k,l,m,n)
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
! C       
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              do p=1,3
               do t=1,3
                dSigmaC_dEps(p,t) = dSigmaC_dEps(p,t) + dSigmaC_dEpsEq_neg *
     1           dEpsEq_neg_dEpsEq * dEpsEq_dDeltaEpsEqPl * 
     2           dDEpsEqPl_dDeltaEpsPl(i,j) * dDeltaEpsPl_dStressDev(i,j,k,l) *
     3           dStressDev_dStress(k,l,m,n) * dStress_dEps(m,n,p,t)       
                enddo
               enddo
              enddo
             enddo
            enddo
           enddo 
         enddo
        enddo 
! C 
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              dSigmaC_dEps(m,n) = dSigmaC_dEps(m,n) + dSigmaC_dEpsEq_neg * 
     1         dEpsEq_neg_dEpsEq * dEpsEq_dDeltaEpsEqPl * 
     2         dDEpsEqPl_dDeltaEpsPl(i,j) * dDeltaEpsPl_dStressPres(i,j) *
     3         dPress_dStress(k,l) * dStress_dEps(k,l,m,n)
             enddo
            enddo
           enddo
          enddo
         enddo 
        enddo       
! C  
! C 
! C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            K1_p1(k,l) = K1_p1(k,l) + 
     1       6.d0*dStressDev_dEps(i,j,k,l)*stress_dev(i,j)
           enddo
          enddo
         enddo
        enddo
! C         
		K1_p2 =6.d0*dPress_dEps*(sigma_c-sigma_t)
! C        
        K1_p3 = 6.d0 * stress_pres * (dSigmaC_dEps - dSigmaT_dEps)		
! C         
        K1_p4 = 2.d0*dSigmaC_dEps*sigma_t + 2.d0*dSigmaT_dEps*sigma_c  		
! C        
        K1 = K1_p1 + K1_p2 + K1_p3 - K1_p4	
! C K2 is already avalible   
! C  K2 dPhi_dGammaCur     
        K2 = phi_pr
		if(ABS(phi_pr).LT.1e-6)then
		 print*, 'phi_pr is almost zero'
		 CALL XIT
		endif 
! C 
! C  CALCULATIONS FOR N_star
! C 
! C         
		N_star = 3.d0*stress_dev + (2.d0/3.d0)*alpha*stress_pres*I_1
        
! C NN is dNstar_dEps  
! dStressDev_dEps
! dPress_dEps   
! C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            NN(i,j,k,l) = 3.d0 * dStressDev_dEps(i,j,k,l) +
     1  (2.d0/3.d0)*alpha*dPress_dEps(k,l)*I_1(i,j)
           enddo
          enddo
         enddo
        enddo
! C Construct dStress_dEps       
! C        
		dStress_dStran_4D = Cel_4D
! C
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              dStress_dStran_4D(i,j,m,n) = dStress_dStran_4D(i,j,m,n) + 
     1         (K2**(-1.d0)) * Cel_4D(i,j,k,l) * (N_star(k,l)*K1(m,n))
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
! C
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              dStress_dStran_4D(i,j,m,n) = dStress_dStran_4D(i,j,m,n) +
     1         gamma_cur*Cel_4D(i,j,k,l)*NN(k,l,m,n) !FIX: + idi - oldu 
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo         
! C
        ind_1 = reshape( (/1, 2, 3, 1, 1, 2/), (/ 6/) )
        ind_2 = reshape( (/1, 2, 3, 2, 3, 3/), (/ 6/) )
C
        do i=1,6 !CONVERT DDSDDE TO VOIGT FORM (3,3,3,3)-->(6,6)
         do j=1,6
          dStress_dStran(i,j) = dStress_dStran_4D(ind_1(i),ind_2(i),ind_1(j),ind_2(j))
         enddo
        enddo   ! C1, D1
C        
C CALCULATIONS FOR K3   
C** K3_pos --> K3
C** K3_neg --> K3_neg   
C
        CALL stressCurve(a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c)
C        
C dSt_dD --> d(sigma_pos)/d(d_m)       
C dSc_dD --> d(sigma_neg)/d(d_m_neg)       
        dSt_dD = -(a_t*EXP(b_t*eps_eq_new)+c_t*EXP(d_t*eps_eq_new))
        dSc_dD = -(a_c*EXP(b_c*eps_eq_new_neg)+c_c*EXP(d_c*eps_eq_new_neg))  		
C
		dPhi_dD = 6.d0*stress_pres*(-dSt_dD) - 2.d0*dSt_dD*sigma_c  
		dPhi_dD_neg =  6.d0*stress_pres * dSc_dD - 2.d0*dSc_dD*sigma_t ! FIX: 6*stress_pres*(-dSc_dD) idi
C
        dKappa_dEpsPlBar = 1.d0 
        dKappa_dEpsPlBar_neg = 1.0d0
C
	    CALL getDd_Dr_exp_damage(d_m, kappa_str, dD_dKappa)		
		CALL getDd_Dr_exp_damage_comp(d_m_neg, kappa_str_neg, dD_neg_dKappa)
C
        K3 = dPhi_dD*dD_dKappa*dKappa_dEpsPlBar
        K3_neg = dPhi_dD_neg*dD_neg_dKappa*dKappa_dEpsPlBar_neg
C
C        
C** dStress_dEpsBar calculations
C** abbr. C2, D2    
C** Since it is staggered solution those terms are zero 
        dStress_dEpsBar = 0.d0 ! 3by3
C
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dStress_dEpsBar(i,j) = dStress_dEpsBar(i,j) +
     1       (K2**(-1.d0))*K3 *Cel_4D(i,j,k,l)*N_star(k,l)
           enddo
          enddo
         enddo
        enddo
		
        dStress_dEpsBar_V(1) = dStress_dEpsBar(1,1)
        dStress_dEpsBar_V(2) = dStress_dEpsBar(2,2)
        dStress_dEpsBar_V(3) = dStress_dEpsBar(3,3)
        dStress_dEpsBar_V(4) = dStress_dEpsBar(1,2)
        dStress_dEpsBar_V(5) = dStress_dEpsBar(1,3)
        dStress_dEpsBar_V(6) = dStress_dEpsBar(2,3) ! C2, D2		

	    dStress_dEpsBar_neg = 0.d0		
C           
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dStress_dEpsBar_neg(i,j) = dStress_dEpsBar_neg(i,j) +
     1       (K2**(-1.d0))*K3_neg *Cel_4D(i,j,k,l)*N_star(k,l)
           enddo
          enddo
         enddo
        enddo
C
C Convert to voigt form       
        dStress_dEpsBar_neg_V(1) = dStress_dEpsBar_neg(1,1)
        dStress_dEpsBar_neg_V(2) = dStress_dEpsBar_neg(2,2)
        dStress_dEpsBar_neg_V(3) = dStress_dEpsBar_neg(3,3)
        dStress_dEpsBar_neg_V(4) = dStress_dEpsBar_neg(1,2)
        dStress_dEpsBar_neg_V(5) = dStress_dEpsBar_neg(1,3)
        dStress_dEpsBar_neg_V(6) = dStress_dEpsBar_neg(2,3) ! C2_neg    
C
C** Linearization of local eq.pl.strain
C** C3, D3
C
! C AA := derivative of local positive equivalent plastic strain w.r.t. new TOTAL plastic strain       
! C AA := d(eps_pl_pos)/d(eps_pl_n+1)
! C Similarly, AA_neg is negative counter part.
! C since derivative of total local equivalent positive plastic strain w.r.t. 
! C incremental local equivalent positive plastic strain is 1.d0, it is not included in formulations
! C** dEpsEq_pos_dEpsEq = derivative of incremental local equivalent positive plastic strain w.r.t. 
! C incremental total equivalent plastic strain, c_pos
! C** dDEpsEqPl_dDeltaEpsPl = derivative of incremental total equivalent plastic strain incremental total  plastic strain
! C** dDeltaEpsPl_dStressDev = derivative of incremental total  plastic strain w.r.t. deviatoric stress
! C** dStressDev_dStress = derivative of deviatoric stress w.r.t. total stress 
! C** dPress_dStress = derivative of pressure w.r.t. total stress 
! C** dStress_dEpsPl = derivative of total stress w.r.t. total plastic strain tensor
C
! dEpsEq_pos_dEpsEq
! dEpsEq_neg_dEpsEq 
		dStress_dEpsPl = 0.d0
        
	    dStress_dEpsPl = Cel_4D

		AA = 0.0d0
C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              do p=1,3
               do t=1,3
                AA(p,t) = AA(p,t) + dEpsEq_pos_dEpsEq * 
     1  dDEpsEqPl_dDeltaEpsPl(i,j)*dDeltaEpsPl_dStressDev(i,j,k,l)*
     2  dStressDev_dStress(k,l,m,n)*dStress_dEpsPl(m,n,p,t)  
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
C
! dEpsEq_pos_dEpsEq
! dEpsEq_neg_dEpsEq
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              AA(m,n) = AA(m,n) + dEpsEq_pos_dEpsEq * 
     1  dDEpsEqPl_dDeltaEpsPl(i,j)*dDeltaEpsPl_dStressPres(i,j)*
     2  dPress_dStress(k,l)*dStress_dEpsPl(k,l,m,n)  
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo   
! C dEpsPl_dEps
		dEpsPl_dEps = 0.d0
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dEpsPl_dEps(k,l) = dEpsPl_dEps(k,l) -
     1  (K2**-1.d0)*AA(i,j) * (N_star(i,j)*K1(k,l)) 
           enddo
          enddo
         enddo
        enddo
! C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dEpsPl_dEps(k,l) = dEpsPl_dEps(k,l) + 
     1  gamma_cur*AA(i,j)*NN(i,j,k,l)
           enddo
          enddo
         enddo
        enddo
!C Voigt Form     
!C C3 is d(epsPlEq)/d(eps)   
        dEpsPl_dEps_V(1) = dEpsPl_dEps(1,1)
        dEpsPl_dEps_V(2) = dEpsPl_dEps(2,2)
        dEpsPl_dEps_V(3) = dEpsPl_dEps(3,3)
        dEpsPl_dEps_V(4) = dEpsPl_dEps(1,2)
        dEpsPl_dEps_V(5) = dEpsPl_dEps(1,3)
        dEpsPl_dEps_V(6) = dEpsPl_dEps(2,3) ! C3
C**
C** C3neg
C**
C
		AA_neg = 0.0d0
C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              do p=1,3
               do t=1,3
                AA_neg(p,t) = AA_neg(p,t) + dEpsEq_neg_dEpsEq * 
     1  dDEpsEqPl_dDeltaEpsPl(i,j)*dDeltaEpsPl_dStressDev(i,j,k,l)*
     2  dStressDev_dStress(k,l,m,n)*dStress_dEpsPl(m,n,p,t)  
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
C
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            do m=1,3
             do n=1,3
              AA_neg(m,n) = AA_neg(m,n) + dEpsEq_neg_dEpsEq * 
     1  dDEpsEqPl_dDeltaEpsPl(i,j)*dDeltaEpsPl_dStressPres(i,j) *
     2  dPress_dStress(k,l)*dStress_dEpsPl(k,l,m,n)  
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo   
C dEpsPl_dEps
		dEpsPl_neg_dEps=0.d0
C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dEpsPl_neg_dEps(k,l) = dEpsPl_neg_dEps(k,l) -
     1  (K2**-1.d0)*AA_neg(i,j) * (N_star(i,j)*K1(k,l)) 
           enddo
          enddo
         enddo
        enddo
C        
        do i=1,3
         do j=1,3
          do k=1,3
           do l=1,3
            dEpsPl_neg_dEps(k,l) = dEpsPl_neg_dEps(k,l) + 
     1  gamma_cur*AA_neg(i,j)*NN(i,j,k,l)
           enddo
          enddo
         enddo
        enddo
C Voigt Form     
C C3 is d(epsPlEq)/d(eps)   
        dEpsPl_neg_dEps_V(1) = dEpsPl_neg_dEps(1,1)
        dEpsPl_neg_dEps_V(2) = dEpsPl_neg_dEps(2,2)
        dEpsPl_neg_dEps_V(3) = dEpsPl_neg_dEps(3,3)
        dEpsPl_neg_dEps_V(4) = dEpsPl_neg_dEps(1,2)
        dEpsPl_neg_dEps_V(5) = dEpsPl_neg_dEps(1,3)
        dEpsPl_neg_dEps_V(6) = dEpsPl_neg_dEps(2,3) ! C3_neg
C**
C**
C**
C dEpsPl_dEpsPl, C4, D4
C dEpsPl_dEpsPl_neg, C4_neg, D4_neg
        dEpsPl_dEpsPl_bar = 0.d0
        dEpsPl_dEpsPl_neg_bar = 0.d0
        do i=1,3
         do j=1,3
          dEpsPl_dEpsPl_bar = dEpsPl_dEpsPl_bar -
     1     (K2**-1.d0)*K3*AA(i,j)*N_star(i,j) ! C4
         enddo
        enddo

        do i=1,3
         do j=1,3
          dEpsPl_dEpsPl_neg_bar = dEpsPl_dEpsPl_neg_bar -
     1  (K2**-1.d0)*K3_neg*AA(i,j)*N_star(i,j) ! C5
         enddo
        enddo  
C**
C**
C**
C dEpsPl_dEpsPl_neg, C5, D5
C dEpsPl_neg_dEpsPl, C5_neg, D5_neg  
        dEpsPl_neg_dEpsPl_bar = 0.d0
        dEpsPl_neg_dEpsPl_neg_bar = 0.d0
        do i=1,3
         do j=1,3
          dEpsPl_neg_dEpsPl_bar = dEpsPl_neg_dEpsPl_bar -
     1  (K2**-1.d0)*K3*AA_neg(i,j)*N_star(i,j) ! C4_neg
         enddo
        enddo

        do i=1,3
         do j=1,3
          dEpsPl_neg_dEpsPl_neg_bar = dEpsPl_neg_dEpsPl_neg_bar -
     1  (K2**-1.d0)*K3_neg*AA_neg(i,j)*N_star(i,j) ! C5_neg
         enddo
        enddo 	
! TODO: derivaties of g+ and g- 
! TODO g+ and g- must be output   
      end subroutine get_mat_tang
C
C
C
      subroutine getCel4D(kappa, mu, Cel_4D)
C**
C** RETURN ELASTIC STIFFNESS TENSOR IN 4D        
        implicit none
! INPUT
        real*8, intent(in) :: kappa, mu
! OUTPUT
        real*8, intent(out) :: Cel_4D(3,3,3,3) 
! LOCAL
        integer i, j, k, l, m, n
        real*8 I_1(3,3), II(3,3,3,3), PP(3,3,3,3)

        CALL getIandPterms(I_1, II, PP)

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                Cel_4D(i,j,k,l)=kappa*I_1(i,j)*I_1(k,l) + 
     1      2.d0*mu*(II(i,j,k,l)-(1.d0/3.d0)*I_1(i,j)*I_1(k,l))
              enddo
            enddo
          enddo
        enddo
      end subroutine getCel4D
C
C
C
      subroutine stressCurve(a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c)
C**        
C** Returns constants for strain-stress fit   
C** t, c : Tension, compression 
C**
        implicit none
C
        real*8, intent(out)::a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c      
! C
! C OLD VALUE of Tension
! C          
        a_t=86.29d0
        b_t=3.775d0
        c_t=-61.29d0
        d_t=-241.d0  
! C
! C NEW VALUE of Compression
! C        
        a_c=121.0d0
        b_c=0.3869d0
        c_c=-53.57
        d_c=-353.7       

! PYTHON SCIPY FIT
! C revision in tension data 14.12.2021 (v4 bununla çözüldü)
        ! a_t= -6.16968110e+01
        ! b_t = -1.80849865e+02
        ! c_t = 9.40949615e+01
        ! d_t = 5.67112134e-02 
! MATLAB FIT slope 2.00	
       ! a =       93.35  (92.83, 93.87)
       ! b =      0.3544  (0.3466, 0.3623)
       ! c =      -61.07  (-62.08, -60.06)
       ! d =      -184.5  (-192.1, -176.8)	
        ! a_t= 93.35
        ! b_t = 0.3544
        ! c_t = -61.07
        ! d_t = -184.5	

! ! fit with slope of 1.01
        ! a_t= 94.19d0
        ! b_t = 0.01817d0
        ! c_t = -61.78d0
        ! d_t = -180.4d0	   
       
! C
! C NEW VALUE of Compression
! C        
        ! a_c = 1.231993e+02
        ! b_c = 0.d0
        ! c_c = -5.076497e+01
        ! d_c = -2.701449e+02		
       ! a =       122.1  (121.5, 122.7)
       ! b =      0.3574  (0.3486, 0.3662)
       ! c =      -49.71  (-51.33, -48.09)
       ! ! d =      -278.5  (-300.4, -256.6)
        ! a_c = 122.1d0
        ! b_c = 0.3574d0
        ! c_c = -49.71d0
        ! d_c = -278.5d0 		   
		
      end subroutine stressCurve 

