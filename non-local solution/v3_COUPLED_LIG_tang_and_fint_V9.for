      subroutine get_tang_fint(x, y, z, hist_n, hist_new, pl_props, 
     1  U_v, K_el, fint_tot, lc, KINC, DTIME, PNEWDT)
        implicit none
        INCLUDE 'KARC.blc'        
C INPUTS        
        real*8, intent(in) :: x(8), y(8), z(8),
     1  hist_n(112), pl_props(10), U_v(40), lc   
        integer, intent(in) :: KINC    
C OUTPUTS
        real*8, intent(out) :: K_el(40, 40), hist_new(112),
     1  fint_tot(40)   
C MODIFIED ABAQUS VARIABLE (BOTH INPUT AND OUTPUT)
        real*8, intent(inout) :: DTIME, PNEWDT
C LOCAL VARIABLES
        real*8 func_uu(24, 24), func_up(24,8), func_pu(8,24), 
     2   func_pp(8,8), dStress_dStran(6,6), pl_stress(6),
     3   xi_v(8), eta_v(8), tau_v(8), w_v(8), dd_dKappa, detJ, 
     4   Jac(3,3), B_ip(6,24), Ni(8), eps_eq_drive, 
     5   ip_strain(6), xi_v_NL(8), eta_v_NL(8), tau_v_NL(8), 
     6   w_v_NL(8), Jac_NL(3,3), B_ip_NL(3,8), detJ_NL,
     7   stress_eff(6), E, mu, nu, kappa, I_1(3,3),
     8   II(3,3,3,3), PP(3,3,3,3), sigma_c_0, sigma_t_0, 
     9   sigma_c_inf, sigma_t_inf, stress_pres_eff, D_el(6,6),
     1   dStress_dEpsPl(6,6)
C
        real*8 U_disp(24), d_m, eps_pl_eq_ip, U_phi(8), kappa_str, 
     1    g_val, dg_dd, gamma_cur, delta_eps_p(3,3), nu_p, alpha, kk,
     2    dStress_dEpsBar_V(6), dEpsPl_dEps_V(6), dEpsPl_dEpsPl_bar,
     3    d_mOLD, aa, bb, kappa_str_neg, d_m_neg, dg_dd_neg, g_val_neg,
     4    dd_dKappa_neg,d_m_negOLD  
C
        integer i_gp, i, j, k, l, m, n, start_, end_,
     1    n_hist, ind_1(40)

! C* NEW TERMS for three-field approach
        real*8 U_phi_neg(8), fint_phi_neg(8), f1_neg(8), f2_neg(8), fd_neg(8)
        real*8 func_upNeg(24,8), func_ppNeg(8,8), func_pNegu(8,24), 
     1   func_pNegp(8,8), func_pNegpNeg(8,8)
        real*8 eps_pl_eq_ip_neg, eps_eq_drive_neg, dStress_dEpsBar_neg_V(6),
     1   dEpsPl_neg_dEps_V(6), dEpsPl_dEpsPl_neg_bar, dEpsPl_neg_dEpsPl_bar, 
     2   dEpsPl_neg_dEpsPl_neg_bar, eps_eq_new, eps_eq_new_neg  
C
! C* vars for internal force column
        real*8 func(24), fint(24), fint_phi(8), f1(8), f2(8), fd(8)
        real*8 V_TOT, G_TENS, G_COMP, dg_dEps_V(6), dgNeg_dEps_V(6)
        real*8 delta_eps_eq_neg, delta_eps_eq_pos		
C
        V_TOT = 0.d0
C Initialize output variables
C
		fint_tot = 0.0d0
		hist_new = 0.d0
		K_el = 0.0d0
C                
        n_hist = 14     
C        
        aa = pl_props(4)
        bb = pl_props(5)
C
!CCCC create vector related displacement and non-local variable phi
        j = 1
        k = 1
        do i=1,36,5
          U_disp(j) = U_v(i)
          U_disp(j+1) = U_v(i+1)
          U_disp(j+2) = U_v(i+2)
          U_phi(k) = U_v(i+3) ! positive
          U_phi_neg(k) = U_v(i+4)
C          
          k = k+1
          j = j+3
        enddo
C
        CALL getIandPterms(I_1, II, PP)
C **
C   Placeholders for Element Internal Force Column
C **
! https://stackoverflow.com/questions/43433004/is-there-an-intrinsic-function-for-initializing-arrays-to-zero-in-fortran  
C
	    fint = 0.d0
        fint_phi = 0.d0
        fint_phi_neg = 0.d0
C
        f1=0.d0
        f2=0.d0
        fd=0.d0
C
        f1_neg=0.d0
        f2_neg=0.d0
        fd_neg=0.d0  		
C **   
C Placeholer for Element Stiffness matrix  
C **   
C 1st Row KUU, KUP, KUPneg
        func_uu = 0.d0
        func_up = 0.d0
        func_upNeg = 0.d0		
C
C ** 2nd row KPU, KPP, KPPneg
		func_pu = 0.d0
		func_pp = 0.d0
		func_ppNeg = 0.d0 
C
C ** 3rd row KPnegU, KPnegP, KPnegPneg   
		func_pNegu = 0.d0
		func_pNegp = 0.d0
		func_pNegpNeg = 0.d0
C       
C
        CALL getXi(xi_v, eta_v, tau_v, w_v)
        CALL getXi_NL(xi_v_NL, eta_v_NL, tau_v_NL, w_v_NL)

C
        do i_gp=1,8 ! Loop over integration points
!C
          CALL getJacobian(x, y, z, xi_v(i_gp),
     1         eta_v(i_gp), tau_v(i_gp), Jac, B_ip, detJ)
          CALL getJacobian_NL(x, y, z, xi_v_NL(i_gp), eta_v_NL(i_gp),
     1        tau_v_NL(i_gp), Jac_NL, B_ip_NL, detJ_NL, Ni)
!C
!C ** Strain at integration poinrt	  
		  ip_strain = MATMUL(B_ip, U_disp)
		  
!C ** non-local variable (equivalent_plastic_strain) at integration point
		  eps_pl_eq_ip = DOT_PRODUCT(Ni, U_phi)
		  eps_pl_eq_ip_neg = DOT_PRODUCT(Ni, U_phi_neg)
		  
! C ** To keep track of history var.s start_ and end_ are used
! C each i.p. has related start and end numbers
! C i.p. 1: hist(1, 14)
! ...
! C i.p. 7: hist(85, 98)
! C i.p. 8: hist(99, 112)
          start_ = n_hist*(i_gp-1)+1
          end_ = start_+n_hist-1
!C          
!C ** Basically UMAT subroutine
          CALL get_stress( pl_props, hist_n(start_:end_),
     1     ip_strain, eps_pl_eq_ip, lc, hist_new(start_:end_),
     2     pl_stress, dStress_dStran, eps_eq_drive, KINC,
     3     DTIME, dStress_dEpsPl, dStress_dEpsBar_V, dEpsPl_dEps_V,
     4     dEpsPl_dEpsPl_bar, PNEWDT, eps_pl_eq_ip_neg, 
     5     eps_eq_drive_neg, dStress_dEpsBar_neg_V, dEpsPl_neg_dEps_V,
     6     dEpsPl_dEpsPl_neg_bar, dEpsPl_neg_dEpsPl_bar, 
     7     dEpsPl_neg_dEpsPl_neg_bar, G_TENS, G_COMP, dg_dEps_V, dgNeg_dEps_V,
     8     delta_eps_eq_neg, delta_eps_eq_pos)
!C		  
          d_mOLD = hist_n(start_+10) ! damage variable (based on non-local solution)
          d_m_negOLD = hist_n(start_+11) ! damage variable (based on non-local solution)
		  
          kappa_str = hist_new(start_+8) ! non-local kappa from second PDE
          kappa_str_neg = hist_new(start_+9) ! non-local kappa from third PDE
          d_m = hist_new(start_+10) ! damage variable (based on non-local solution)
          d_m_neg = hist_new(start_+11) ! damage variable (based on non-local solution)
!C **
!C ** Store integration point data for dummy elements
!C **
          do i=1,6
            STRESS_k(noel_k, i_gp, i) = pl_stress(i)
            STRAIN_k(noel_k, i_gp, i) = ip_strain(i)
          enddo
!C**
          EPS_PL_K(noel_k, i_gp,1)=hist_new(start_) ! components of plastic strain tensor
          EPS_PL_K(noel_k, i_gp,2)=hist_new(start_+1)
          EPS_PL_K(noel_k, i_gp,3)=hist_new(start_+2)
          EPS_PL_K(noel_k, i_gp,4)=hist_new(start_+3)
          EPS_PL_K(noel_k, i_gp,5)=hist_new(start_+4)
          EPS_PL_K(noel_k, i_gp,6)=hist_new(start_+5)
!C
          EPS_EQ_K(noel_k, i_gp)=hist_new(start_+6)! LOCAL equivalent plasstic strain (+)
          EPS_EQ_N_K(noel_k, i_gp)=hist_new(start_+7)! LOCAL equivalent plasstic strain (-)
!C
          DM_K(noel_k, i_gp)=hist_new(start_+10)!  Dt
          DM_N_K(noel_k, i_gp)=hist_new(start_+11)! Dc 
!C
!C **
!C ** COMPUTE Element Internal force column
!C **                      
		  fint = fint - w_v(i_gp) * MATMUL(TRANSPOSE(B_ip), pl_stress)*detJ

!C  Positive non-local eq.pl.strain       
!C  fint_phi = f1 + f2 - fd
!C	  
		  f1 = f1 + w_v_NL(i_gp) * Ni * DOT_PRODUCT(Ni,U_phi) * detJ
!C
          CALL g_func(d_m, g_val) ! localizing gradient, interaction function
		  f2 = f2 + w_v_NL(i_gp)*g_val*(lc**2.d0) *
     1     MATMUL(TRANSPOSE(B_ip_NL), MATMUL(B_ip_NL, U_phi)) * detJ
!C  
		  fd = fd + w_v_NL(i_gp) * Ni * eps_eq_drive * detJ
!C **
!C  Negative non-local eq.pl.strain         
!C  fint_phi_neg = f1_neg + f2_neg - fd_neg
!C
		  f1_neg = f1_neg + w_v_NL(i_gp) * Ni * DOT_PRODUCT(Ni,U_phi_neg) * detJ
          CALL g_func(d_m_neg, g_val_neg) ! localizing gradient, interaction function
!C
		  f2_neg = f2_neg + w_v_NL(i_gp) * g_val_neg * (lc**2.d0) *
     1     MATMUL(TRANSPOSE(B_ip_NL), MATMUL(B_ip_NL, U_phi_neg)) * detJ
!C
		  fd_neg = fd_neg + w_v_NL(i_gp) * Ni * eps_eq_drive_neg * detJ
C **
C Compute Element Tangent
! C **
! C ** 1st row
! C **
C KU,U
C use d(stress)/d(strain)
          CALL get_func_uu(w_v(i_gp), B_ip, dStress_dStran, detJ, func_uu) ! K_UU, C1
C 
C KU,E    
C use d(stress)/d(EpsBar_pos)
C derivative of stress w.r.t. Positive local equivalent plastic strain             
        CALL get_func_up(w_v(i_gp), B_ip, Ni, detJ, dStress_dEpsBar_V, func_up) ! K_UE, C2
C
C KU,Eneg
C use d(stress)/d(EpsBar_neg)
C derivative of stress w.r.t. Negative local equivalent plastic strain  
        CALL get_func_upNeg(w_v(i_gp), B_ip, Ni, detJ, dStress_dEpsBar_neg_V, func_upNeg) ! K_UEneg, C2_neg
! C**
! C** 2nd Row
! C**
C K_E,U
C use d(EpsPl_pos)/d(strain)
C derivative of Positive local equivalent plastic strain w.r.t. strain
          CALL get_func_pu(B_ip, w_v, Ni, detJ, dEpsPl_dEps_V, G_TENS, dg_dEps_V, delta_eps_eq_pos, func_pu)  
C
C K_E,E
C use d(EpsPl_pos)/(EpsBar_pos)
C derivative of Positive local equivalent strain w.r.t Positive 'NON-LOCAL' equivalent strain
          CALL getDd_Dr_exp_damage(d_m, kappa_str, dd_dKappa)          
          CALL get_dg_dd(d_m, dg_dd) 
          CALL g_func(d_m, g_val)
C          
          G_TENS = 1.d0 ! already included in UMAT
          CALL get_func_pp(w_v(i_gp), Ni, detJ, lc, B_ip_NL, 
     1      U_phi, dd_dKappa, dg_dd, g_val,     
     2      dEpsPl_dEpsPl_bar, G_TENS, func_pp) ! K_EE
C 
C K_E,Eneg  is zero by default
C use d(EpsPl_pos)/(EpsBar_neg)
C derivative of Positive local equivalent strain w.r.t NEGATIVE 'NON-LOCAL' equivalent strain
          CALL get_func_ppNeg(w_v(i_gp), Ni, detJ, dEpsPl_dEpsPl_neg_bar, func_ppNeg)
! C**
! C** 3rd Row      
! C**
C K_En, U
C use d(EpsPl_neg)/d(strain)
C derivative of NEGATIVE local equivalent plastic strain w.r.t. strain
          CALL get_func_pNegu(B_ip, w_v, Ni, detJ, dEpsPl_neg_dEps_V, G_COMP, dgNeg_dEps_V, delta_eps_eq_neg, func_pNegu)	 
C     
C K_En, E
C use d(EpsPl_neg)/(EpsBar_pos)
C derivative of NEGATIVE local equivalent strain w.r.t Positive 'NON-LOCAL' equivalent strain
         CALL get_func_pNegp( w_v(i_gp), Ni, detJ, dEpsPl_neg_dEpsPl_bar, func_pNegp)
C
C K_En, En
C use d(EpsPl_neg)/(EpsBar_neg)
C derivative of NEGATIVE local equivalent strain w.r.t negative 'NON-LOCAL' equivalent strain
          CALL getDd_Dr_exp_damage_comp(d_m_neg, kappa_str_neg, dd_dKappa_neg)          
          CALL get_dg_dd(d_m_neg, dg_dd_neg)
          CALL g_func(d_m_neg, g_val_neg) ! localizing gradient, interaction function
		  G_COMP = 1.d0 ! already included in UMAT
          CALL get_func_pNegpNeg(w_v(i_gp), Ni, detJ, lc, B_ip_NL,
     1  U_phi_neg, dd_dKappa_neg, dg_dd_neg, g_val_neg, dEpsPl_neg_dEpsPl_neg_bar, 
     2  G_COMP, func_pNegpNeg)   

        enddo! Loop over integration points (i_gp)
C once all integration points in element is completed, element number for
C dummy elements is increased. Hence the integration point data in dummy 
C element is kept track
        noel_k = noel_k + 1 ! related to dummy element numbering          
C ** 
C Construct element internal force column
C **      
C	
		fint_phi=-(f1+f2-fd)
		fint_phi_neg=-(f1_neg+f2_neg-fd_neg)
		fint_tot = 0.d0
        j = 1
        k = 1
        do i=1,36,5
          fint_tot(i) = fint(j)
          fint_tot(i+1) = fint(j+1)
          fint_tot(i+2) = fint(j+2)
          fint_tot(i+3) = fint_phi(k)
          fint_tot(i+4) = fint_phi_neg(k)
          j = j + 3
          k = k + 1
        enddo 
C ** 
C Construct element stiffness matrix
C **           
        CALL construct_K_el(func_uu, func_up, func_pu, func_pp,
     1    func_upNeg, func_ppNeg, func_pNegu, func_pNegp,   
     2    func_pNegpNeg, K_el, DTIME, PNEWDT, KINC)  ! K_ELEMENT 40-by-40 Element stiffness matrix

      end subroutine get_tang_fint