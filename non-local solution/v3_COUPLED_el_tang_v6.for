C
C PARTS OF THE ELEMENT STIFFNESS MATRIX  
C  
C Variables provided by UMAT function are as follows,
C pl_tang:
C dStress_dEpsBar_V:
C dEpsPl_dEps_V: 
C dEpsPl_dEpsPl:
C*** 1st Row
C**************************************
C*** KUU
C**************************************
      subroutine get_func_uu(w_v, B_ip, pl_tang, detJ, func_uu)
        implicit none
! INPUT
        real*8, intent(in) :: w_v,  B_ip(6,24), pl_tang(6,6), detJ
! OUTPUT
        real*8, intent(inout) ::  func_uu(24,24) ! K_UU contribution from i.p.
! LOCAL
        integer i, j, k, l
!C** INPUTS:
!C** pl_tang : dStress/dStrain         	
! CODE:
		func_uu = func_uu + w_v * MATMUL(TRANSPOSE(B_ip), MATMUL(pl_tang, B_ip)) * detJ
      end subroutine get_func_uu  
!c
!c
      subroutine get_func_up(w_v, B_ip, Ni, 
     1    detJ, dStress_dEpsBar_V, func_up)
        implicit none
! INPUT
        real*8, intent(in) :: w_v,  B_ip(6,24), 
     1     Ni(8), detJ, dStress_dEpsBar_V(6)
! OUTPUT
        real*8, intent(inout) ::  func_up(24,8)
! LOCAL
        integer i, j, k
!C** INPUTS:
!C** dStress_dEpsBar_V = dStress/dEpsBar, C2+, D2+        
        do i=1,24
          do j=1,6
            do k=1,8
              func_up(i,k) = func_up(i,k) + 
     1         w_v*B_ip(j,i)*dStress_dEpsBar_V(j)*Ni(k)*detJ
            enddo
          enddo
        enddo    
      end subroutine get_func_up    
C
      subroutine get_func_upNeg(w_v, B_ip, Ni, 
     1    detJ, dStress_dEpsBar_neg_V, func_upNeg)
        implicit none
! INPUT
        real*8, intent(in) :: w_v,  B_ip(6,24), 
     1     Ni(8), detJ, dStress_dEpsBar_neg_V(6)
! OUTPUT
        real*8, intent(inout) ::  func_upNeg(24,8)
! LOCAL
        integer i, j, k
C** INPUTS:
C** dStress_dEpsBar_neg_V = dStress/dEpsNegBar, C2-, D2-         
        do i=1,24
          do j=1,6
            do k=1,8
              func_upNeg(i,k) = func_upNeg(i,k) + 
     1         w_v*B_ip(j,i)*dStress_dEpsBar_neg_V(j)*Ni(k)*detJ
            enddo
          enddo
        enddo    
      end subroutine get_func_upNeg  	  
    
C*** 2nd Row
C**************************************
C***KE+E+
C**************************************     
C
      subroutine get_func_pu(B_ip, w_v, Ni, detJ, 
     3    dEpsPl_dEps_V, g_plus, dg_dEps_V, dEpsPos, func_pu)
        implicit none
! INPUT
        real*8, intent(in) :: B_ip(6,24), w_v, Ni(8), detJ, 
     1  dEpsPl_dEps_V(6), g_plus, dg_dEps_V(6), dEpsPos
! OUTPUT
        real*8, intent(inout) ::  func_pu(8,24)
! LOCAL
	    real*8 d_gPlus_d_eps(6)
        integer i, j, k, l, m
		d_gPlus_d_eps = 0.d0
C NOTE: Simply -d(fd)/d(U)        
C** INPUTS:
C** dEpsPl_dEps_V = d(EpsBar)/d(Eps), C3+, D3+    	   
!		   do i=1,8
!		    do j=1,6
!			 do k=1,24
!		      func_pu(i,k) = func_pu(i,k) -
!        1	  w_v * Ni(i) * (dg_dEps_V(j) *B_ip(j,k) * dEpsPos) * detJ
!	         enddo
!			enddo
!         enddo 	
		   do i=1,8
		    do j=1,6
			 do k=1,24
		      func_pu(i,k) = func_pu(i,k) - 
     1	  w_v * Ni(i) * (dEpsPl_dEps_V(j)*B_ip(j,k)) * detJ
	         enddo
			enddo
           enddo 			   
      end subroutine get_func_pu  
!*
!*	  
      subroutine get_func_pp(w_v, Ni, detJ, lc, B_ip_NL,
     1  U_phi, dd_dKappa, dg_dd, g_val, dEpsPl_dEpsPl_bar, G_TENS, func_pp)
        implicit none
! INPUT
        real*8, intent(in) :: w_v, Ni(8), detJ, lc, B_ip_NL(3,8), 
     1    U_phi(8), dd_dKappa, dg_dd, g_val, dEpsPl_dEpsPl_bar, G_TENS
! OUTPUT
        real*8, intent(inout) ::  func_pp(8,8) ! sifirlama. Uzerine ekliyor
! LOCAL
        integer i, j, k, l, m
        real*8 dF1_dEpsEqPl(8,8), dF2_dEpsEqPl(8,8), dF3_dEpsEqPl(8,8) 
C** INPUTS
C provided as input in LIG_TANG_AND_FINT_V3.for
C dd_dKappa = d(Dt)/d(kappa_t)
C dg_dd = d(g)/d(Dt)
C
C dEpsPl_dEpsPl_bar = d(EpsPl)/d(EpsBar), D4+, C4+
C D4+, C4+ : Derivative of local equivalent plastic strain w.r.t non-local one
		dF1_dEpsEqPl = 0.d0	! 8by8
		dF2_dEpsEqPl = 0.d0	
		dF3_dEpsEqPl = 0.d0	
! C       
        do i=1,8
          do k=1,8
            dF1_dEpsEqPl(i,k) = w_v*Ni(i)*Ni(k)*detJ
          enddo
        enddo
! C             

		dF2_dEpsEqPl = w_v * g_val * (lc**2.d0) * MATMUL(TRANSPOSE(B_ip_NL), B_ip_NL) * detJ
! C below is wv * dg_dEpsPl_plus * lc**2 * B.T B EpsPl_plus detj
        do i=1,8
          do j=1,3
            do k=1,8
              do l= 1,8
                 dF2_dEpsEqPl(i,l) = dF2_dEpsEqPl(i,l) +
     1          w_v * (lc**2.d0) * B_ip_NL(j,i)*B_ip_NL(j,k)*U_phi(k)*
     2         dg_dd * dd_dKappa * Ni(l) *detJ 
              enddo
            enddo
          enddo
        enddo    		
! C   
        do i=1,8
          do j=1,8
            dF3_dEpsEqPl(i,j) = w_v*Ni(i)*dEpsPl_dEpsPl_bar*Ni(j)*detJ
          enddo
        enddo
        func_pp = func_pp + dF1_dEpsEqPl + dF2_dEpsEqPl - dF3_dEpsEqPl !FIX: + idi - yaptim
      end subroutine get_func_pp 
!
!	  
      subroutine get_func_ppNeg(w_v, Ni, detJ, 
     3    dEpsPl_dEpsPl_neg_bar, func_ppNeg)
        implicit none
! INPUT
        real*8, intent(in) :: w_v, Ni(8), detJ, dEpsPl_dEpsPl_neg_bar
! OUTPUT
        real*8, intent(inout) ::  func_ppNeg(8,8) ! sifirlama. Uzerine ekliyor
		integer i, j
C** INPUTS
C** dEpsPl_dEpsPl_neg_bar = d(EpsPl)/d(EpsNegBar), C5+, D5+
C** dEpsPl_dEpsPl_neg_bar : Derivative of local Positive equivalent
C**  plastic strain w.r.t non-local NEGATIVE one    

        do i=1,8
          do j=1,8
            func_ppNeg(i,j) = func_ppNeg(i,j) -  w_v*Ni(i)*
     1        dEpsPl_dEpsPl_neg_bar*Ni(j)*detJ
          enddo
        enddo        
      end subroutine get_func_ppNeg
	  
C*** 3rd Row
C**************************************
C*** KE-E-
C**************************************
C
!
!	  
      subroutine get_func_pNegu(B_ip, w_v, Ni, detJ, 
     3    dEpsPl_neg_dEps_V, g_neg, dgNeg_dEps_V, dEpsNeg, func_pNegu)	 
        implicit none
! INPUT
        real*8, intent(in) :: B_ip(6,24), w_v, Ni(8), detJ, 
     1  dEpsPl_neg_dEps_V(6), g_neg, dgNeg_dEps_V(6), dEpsNeg
! OUTPUT
!C NOTE: Simply -d(fd)/d(U)        
!C** INPUTS:
!C** dEpsPl_neg_dEps_V = d(EpsNegBar)/d(Eps), C3-, D3-   
!C** dgNeg_dEps_V = derivative of g^- (multiplier of incremental eq. pl. strain) w.r.t strain
        real*8, intent(inout) ::  func_pNegu(8,24) ! sifirlama. Uzerine ekliyor
! LOCAL		
        integer i,j,k
		do i=1,8
		 do j=1,6
		  do k=1,24
		   func_pNegu(i,k) = func_pNegu(i,k) - 
     1	  w_v * Ni(i) * (dEpsPl_neg_dEps_V(j)*B_ip(j,k)) * detJ
	      enddo
		 enddo
        enddo 		   		 
      end subroutine get_func_pNegu



      subroutine get_func_pNegpNeg(w_v, Ni, detJ, lc, B_ip_NL,
     1  U_phi_neg, dd_dKappa, dg_dd, g_val, dEpsNegPl_dEpsPlNeg_bar, 
     2  G_COMP, func_pp)
        implicit none
! INPUT
        real*8, intent(in) :: w_v, Ni(8), detJ, lc, B_ip_NL(3,8), 
     1    U_phi_neg(8), dd_dKappa, dg_dd, g_val, dEpsNegPl_dEpsPlNeg_bar, G_COMP
! OUTPUT
        real*8, intent(inout) ::  func_pp(8,8) ! sifirlama. Uzerine ekliyor
! LOCAL
        integer i, j, k, l, m
        real*8 dF1_dEpsEqPl(8,8), dF2_dEpsEqPl(8,8), dF3_dEpsEqPl(8,8) 
C** INPUTS
C provided as input in LIG_TANG_AND_FINT_V3.for
C dd_dKappa = d(Dc)/d(kappa_c)
C dg_dd = d(g)/d(Dc)
C
C dEpsNegPl_dEpsPlNeg_bar = d(EpsNegPl)/d(EpsNegBar), D4-, C4-
C D4+, C4+ : Derivative of NEGATIVE local equivalent plastic strain w.r.t non-local NEGATIVE one        
C
		dF1_dEpsEqPl = 0.d0	!8by8
		dF2_dEpsEqPl = 0.d0	
		dF3_dEpsEqPl = 0.d0
C       
        do i=1,8
         do k=1,8
          dF1_dEpsEqPl(i,k) = w_v*Ni(i)*Ni(k)*detJ
         enddo
        enddo
C             
		dF2_dEpsEqPl = w_v * g_val * (lc**2.d0) * MATMUL(TRANSPOSE(B_ip_NL), B_ip_NL) * detJ		
C        
! C below is wv * dg_dEpsPl_plus * lc**2 * B.T B EpsPl_neg detj
        do i=1,8
          do j=1,3
            do k=1,8
              do l= 1,8
                 dF2_dEpsEqPl(i,l) = dF2_dEpsEqPl(i,l) +
     1          w_v * (lc**2.d0) * B_ip_NL(j,i)*B_ip_NL(j,k)*U_phi_neg(k)*
     2         dg_dd * dd_dKappa * Ni(l) *detJ 
              enddo
            enddo
          enddo
        enddo  

        do i=1,8
          do j=1,8
            dF3_dEpsEqPl(i,j) = w_v*Ni(i)*dEpsNegPl_dEpsPlNeg_bar*Ni(j)*detJ
          enddo
        enddo
		
        func_pp = func_pp + dF1_dEpsEqPl + dF2_dEpsEqPl - dF3_dEpsEqPl
      end subroutine get_func_pNegpNeg   
!
!
      subroutine get_func_pNegp(w_v, Ni, detJ, 
     3    dEpsPl_neg_dEpsPl_bar, func_pNegp)
        implicit none
! INPUT
        real*8, intent(in) :: w_v, Ni(8), detJ, dEpsPl_neg_dEpsPl_bar
! OUTPUT
        real*8, intent(inout) ::  func_pNegp(8,8) ! sifirlama. Uzerine ekliyor
        integer i,j
C** INPUTS
C** dEpsPl_neg_dEpsPl_bar = d(EpsPlNeg)/d(EpsBar), C5+, D5+
C** dEpsPl_neg_dEpsPl_bar : Derivative of NEGATIVE local 
C** equivalent plastic strain w.r.t non-local POSITIVE one    
        do i=1,8
          do j=1,8
            func_pNegp(i,j) = func_pNegp(i,j) - w_v*Ni(i)*
     1        dEpsPl_neg_dEpsPl_bar*Ni(j)*detJ
          enddo
        enddo             
      end subroutine get_func_pNegp   
C
      subroutine construct_K_el(func_uu, func_up, func_pu, func_pp,
     1    func_upNeg, func_ppNeg, func_pNegu, func_pNegp,   
     2    func_pNegpNeg, K_el, DTIME, PNEWDT, KINC)        
        implicit none
! INPUT
        real*8, intent(in) :: func_uu(24,24), func_up(24,8), 
     1    func_pu(8,24), func_pp(8,8)	 
        real*8, intent(in) :: func_upNeg(24,8), func_ppNeg(8,8), 
     1   func_pNegu(8,24),  func_pNegp(8,8), func_pNegpNeg(8,8)
! in &  out
        real*8, intent(inout) :: DTIME, PNEWDT
        integer, intent(in) :: KINC  		
! OUTPUT
        real*8, intent(out) ::  K_el(40,40)
! LOCAL
        real*8 K_dum(40,40)
        integer i, j, k, l, ind_1(40)
C
		K_dum = 0.d0
		K_el = 0.d0
C
C ** 1st Row of KEL
        do i=1,24
          do j=1,24 ! 1:24, 1:24
            K_dum(i,j) = func_uu(i,j)
          enddo
        enddo  
C
        do i=1,24
          do j=1,8 ! 1:24, 25:32
            K_dum(i,j+24)= func_up(i,j)
          enddo
        enddo
C
        do i=1,24
          do j=1,8 ! 1:24, 32:40
            K_dum(i,j+32)= func_upNeg(i,j)
          enddo
        enddo   		
C
C ** 2nd row of KEL
        do i=1,8
          do j=1,24 ! 25:32, 1:24
            K_dum(i+24,j)= func_pu(i,j)
          enddo
        enddo
C
        do i=1,8
          do j=1,8 ! 25:32, 25:32
            K_dum(i+24,j+24)= func_pp(i,j)
          enddo
        enddo     
C        
        do i=1,8
          do j=1,8 ! 25:32, 32:40
            K_dum(i+24,j+32)= func_ppNeg(i,j)
          enddo
        enddo               
C
C ** 3rd row of KEL
C        
        do i=1,8
          do j=1,24 ! 32,40, 1:24
            K_dum(i+32,j)= func_pNegu(i,j)
          enddo
        enddo
C
        do i=1,8
          do j=1,8 ! 32:40, 25:32
            K_dum(i+32,j+24)= func_pNegp(i,j)
          enddo
        enddo     
C        
        do i=1,8
          do j=1,8 ! 32:40, 32:40
            K_dum(i+32,j+32)= func_pNegpNeg(i,j)
          enddo
        enddo  
C
C
C First 3 related to x-y-z count them first
C one they are done. COntinue with non-local quabtitis
C 1 2 3 : x-y-z
C 25: non-local-positive
C 33: non-local-negative
        ind_1 = reshape( (/1,2,3,25,33, 4,5,6,26,34, 7,8,9,27,35,
     1  10,11,12,28,36, 13,14,15,29,37, 16,17,18,30,38, 
     2  19,20,21,31,39, 22,23,24,32,40 /), (/ 40/) )
C          
C
C Be aware not Symmetric
        do i=1,40
          do j=1,40
            K_el(i,j)=K_dum(ind_1(i), ind_1(j))
          enddo
        enddo	
      end subroutine construct_K_el   