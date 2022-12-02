      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

      real*8 E1,E2, G12, G13, nu12, nu13, R_para_t, R_para_c, E_para_f,
     1 nu_perp_para_f,  m_sigma_f, p_perp_perp_t, p_perp_perp_c,
     2 p_perp_para_t, p_perp_para_c, R_perp_t, R_perp_para, R_perp_c,
     3 GG_para_t, GG_para_c, GG_perp_t, GG_perp_c, GG_perp_para, Lc, E3,
     4 G23, nu23, nu21, nu32, nu31, nu_perp_para, E_para  
      
      real*8 del_ft, del_fc, del_m1t, del_m1c, del_m2t, del_m2c, rE_FFt,
     1 rE_FFc, rE_IFFt, rE_IFFc, theta_max, theta_fp
      
      real*8 MYSTRAN(6), strain(3,3), pi
      
      real*8 C_f(6,6), C_m1(6,6), C_m2(6,6), CeSecond(6,6), 
     1 effStress(6), sigma_f(6), sigma_m1(6), sigma_m2(6)
      
      real*8 fE_FFt, fE_FFc, fE_FFt_next, fE_FFc_next
      
      real*8 R_perp_At, R_perp_para_A, R_perp_perp_A

      integer crackFlag, i,j,k,l
            
      real*8 fE_IFF_max, theta, theta_rad, c, s, sigma_n, tau_nt,
     1 tau_nl, p_over_R_t, p_over_R_c, c2, s2, const1, const2,
     2 const3, const4, const5, const6, const7, const8, const9,
     3 const10, const11, const13, const14, const15, const16, const17
     4 fE_IFF_cur, nu_w, fE_IFFc_next, fE_IFFc, fE_IFFt_next, fE_IFFt
      
      integer flag    
      
      real*8 del_ft_next, del_fc_next, A_ft, A_fc      
      integer flag_fc, flag_ft  
      
      real*8 del_m1t_next, del_m1c_next, del_m2t_next, del_m2c_next,
     1 A_m1t, A_m2t, A_m1c, A_m2c
      integer flag_mc, flag_mt    
      
      real*8 C_secant(6,6)      
      
      real*8 Nf1_f(6), Nf2_f(6), Nf3_f(6), C1_f(6), C2_f(6), C3_f(6)
      
      real*8 A_tf, B1_tf, B2_tf, B3_tf
      real*8 del_ft_part(6), C_tg_tf(6,6)  
      
      real*8 A_cf, B1_cf, B2_cf, B3_cf
      real*8 del_fc_part(6), C_tg_cf(6,6)
      
      real*8 NN1(6), NN2(6), NN3(6), CC1(6), CC2(6), CC3(6),
     1 theta_fp_rad   
      
      real*8 A_tm1, A_tm2, BB1, BB2, BB3
      real*8 del_m1t_part(6), C_tg_tm1(6,6), del_m2t_part(6),
     1 C_tg_tm2(6,6)    
      
      real*8 A_cm1, A_cm2, del_m1c_part(6), del_m2c_part(6),
     1 C_tg_cm1(6,6), C_tg_cm2(6,6)   
      
      real*8 tauLimit, feLimit, zero
      
      real*8 GPerpPerp, GPerpPara, G21
      
      real*8 del_ft_v, del_fc_v, del_m1t_v, del_m1c_v, del_m2t_v, 
     1      del_m2c_v, del_ft_v_next, del_fc_v_next, del_m1t_v_next,
     2      del_m1c_v_next, del_m2t_v_next, del_m2c_v_next, eta
      
      
      real*8 reductionFactor, fE
      
      real*8 theta1, theta2, theta3, theta4, TOLERANCE, DIFF, gamma,
     1   a_gs, b_gs, c_gs, fE_gs4, fE_gs3      
      tauLimit=1e-12
      feLimit=1e-12
      zero=1e-12
      
      ! if(NOEL.EQ.556)then
      !   E1 = 0.25d0*PROPS(1)
      !   E2 = 0.25d0*PROPS(2)
      !   G12 = 0.25d0*PROPS(3)
      !   G13 = 0.25d0*PROPS(4)
      !   nu12 = 0.25d0*PROPS(5)
      !   nu13 = 0.25d0*PROPS(6)    
      !   R_para_t = 0.25d0*PROPS(7)
      !   R_para_c = 0.25d0*PROPS(8)
      !   E_para_f = 0.25d0*PROPS(9)
      !   nu_perp_para_f = 0.25d0*PROPS(10)      
      !   m_sigma_f = 0.25d0*PROPS(11)
      !   p_perp_perp_t = 0.25d0*PROPS(12)
      !   p_perp_perp_c = 0.25d0*PROPS(13)
      !   p_perp_para_t = 0.25d0*PROPS(14)
      !   p_perp_para_c = 0.25d0*PROPS(15)   
      !   R_perp_t = 0.25d0*PROPS(16)
      !   R_perp_para = 0.25d0*PROPS(17)
      !   R_perp_c = 0.25d0*PROPS(18)
      !   GG_para_t = 0.25d0*PROPS(19)*1.d0
      !   GG_para_c = 0.25d0*PROPS(20)*1.d0
      !   GG_perp_t = 0.25d0*PROPS(21)*1.d0
      !   GG_perp_c = 0.25d0*PROPS(22)*1.d0
      !   GG_perp_para = 0.25d0*PROPS(23)*1.d0 
      !   Lc = PROPS(24)
      ! else
        E1 = PROPS(1)
        E2 = PROPS(2)
        G12 = PROPS(3)
        G13 = PROPS(4)
        nu12 = PROPS(5)
        nu13 = PROPS(6)    
        R_para_t = PROPS(7)
        R_para_c = PROPS(8)
        E_para_f = PROPS(9)
        nu_perp_para_f = PROPS(10)      
        m_sigma_f = PROPS(11)
        p_perp_perp_t = PROPS(12)
        p_perp_perp_c = PROPS(13)
        p_perp_para_t = PROPS(14)
        p_perp_para_c = PROPS(15)   
        R_perp_t = PROPS(16)
        R_perp_para = PROPS(17)
        R_perp_c = PROPS(18)
        GG_para_t = PROPS(19)*1.d0
        GG_para_c = PROPS(20)*1.d0
        GG_perp_t = PROPS(21)*1.d0
        GG_perp_c = PROPS(22)*1.d0
        GG_perp_para = PROPS(23)*1.d0 
        Lc = PROPS(24)        
      ! endif
C Orthtropic material property, ABAQUS Analysis user manual p.74      
      E3 = E2
      G23 = G13
      nu23 = nu13
      
c     nu_ij/Ei = nu_ji/Ej Poisson's ratios
      nu21 = E2*nu12/E1
      nu32 = E3*nu23/E2
      nu31 = E3*nu13/E1    
      
      nu_perp_para = nu31  
      E_para = E1   

      G21=G12
c###############################################################   
c######################## LOAD SDV #############################
c###############################################################  
      del_ft = STATEV(1)
      del_fc = STATEV(2)
      del_m1t = STATEV(3)
      del_m1c = STATEV(4)
      del_m2t = STATEV(5)
      del_m2c = STATEV(6)
 
      rE_FFt = STATEV(7)
      rE_FFc = STATEV(8)
      rE_IFFt = STATEV(9)
      rE_IFFc = STATEV(10)
      
!     r values must be equal to or larger than 1.0
      rE_FFt = max(1.d0, rE_FFt)
      rE_FFc = max(1.d0, rE_FFc)
      rE_IFFt = max(1.d0, rE_IFFt)
      rE_IFFc = max(1.d0, rE_IFFc)
      
      fE_FFt = STATEV(11)
      fE_FFc = STATEV(12)
      fE_IFFt = STATEV(13)
      fE_IFFc = STATEV(14)
      
      crackFlag = STATEV(15)
c     theta_fp = STATEV(16)
      theta_max = STATEV(17) 
      
c     Viscous damage parameters
c     Regularization is not used.      
c     eta is another problem      
      del_ft_v = STATEV(18)
      del_fc_v = STATEV(19)
      del_m1t_v = STATEV(20)
      del_m1c_v = STATEV(21)
      del_m2t_v = STATEV(22)
      del_m2c_v = STATEV(23)      
      
      ! do i=1,6
      !  MYSTRAN(i) = STRAN(i) + DSTRAN(i)
      ! enddo      
      do i=1,6
       MYSTRAN(i) = STRAN(i) + DSTRAN(i)
      enddo          
      
      strain(1,1) = MYSTRAN(1)
      strain(2,2) = MYSTRAN(2)
      strain(3,3) = MYSTRAN(3)
      strain(1,2) = 0.5d0 * MYSTRAN(4)
      strain(1,3) = 0.5d0 * MYSTRAN(5)
      strain(2,3) = 0.5d0 * MYSTRAN(6)
      strain(2,1) = strain(1,2)
      strain(3,1) = strain(1,3)
      strain(3,2) = strain(2,3)        
c      print*,strain
      
      pi = 4.d0*datan(1.d0)    
      
      call get_C_elastic(CeSecond, E1, E2, E3, nu12, nu13, 
     1 nu23, G12, G13, G23)
      call get_C_slices(C_f, C_m1, C_m2, CeSecond )      
      
      do i=1,6
       effStress(i) = 0.d0
       sigma_f(i) = 0.d0
       sigma_m1(i) = 0.d0
       sigma_m2(i) = 0.d0
      enddo      
      
      do i=1,6 ! calculate effective stresses sigma_eff = C_el:STRAN
        do j=1,6
          sigma_f(i) = sigma_f(i) + C_f(i,j)*MYSTRAN(j)
          sigma_m1(i) = sigma_m1(i) + C_m1(i,j)*MYSTRAN(j)
          sigma_m2(i) = sigma_m2(i) + C_m2(i,j)*MYSTRAN(j)
          effStress(i) = effStress(i) + CeSecond(i,j)*MYSTRAN(j)
        enddo        
      enddo  
      
c###############################################################  
c######################### FF check ############################
c############################################################### 
c### 1. EXPOSURE FACTORS
!     Fiber Fracture (FF)
      if(effStress(1) >= 0.d0) then
        fE_FFt_next = (1.d0/R_para_t)*( effStress(1)-
     1 (nu_perp_para - (E_para/E_para_f)*nu_perp_para_f*m_sigma_f)*
     2 (effStress(2)+effStress(3)))        
        fE_FFc_next = fE_FFc
      else
        fE_FFt_next = fE_FFt
        fE_FFc_next = (-1.d0/R_para_c)*( effStress(1)-
     1 (nu_perp_para - (E_para/E_para_f)*nu_perp_para_f*m_sigma_f)*
     2 (effStress(2)+effStress(3)))
      endif   
     
c###############################################################  
c######################### IFF check ###########################
c###############################################################
      R_perp_At = R_perp_t
      R_perp_para_A = R_perp_para
      R_perp_perp_A = R_perp_c/(2.d0*(1.d0+p_perp_perp_c)) 
      nu_w = 1.d0

      if(abs(rE_IFFt-1.d0)<zero .AND. abs(rE_IFFc-1.d0)<zero)then  
c        real*8 theta1, theta2, theta3, theta4, TOLERANCE, DIFF, gamma,
c     1   a_gs, b_gs, c_gs, fE_gs4, fE_gs3
        theta1 = -90.d0
        theta2 = 90.d0
        TOLERANCE = 1e-5
        DIFF = 1.d0
        gamma = (1.d0 + sqrt(5.d0)) / 2.d0
        
        a_gs = (theta2-theta1)/(1.d0+gamma)
        b_gs = theta2-theta1-a_gs
        c_gs = a_gs/gamma

        theta3 = theta1 + a_gs
        theta4 = theta3 + c_gs
        
        do while(abs(DIFF)>TOLERANCE)
        
          CALL getFE(theta4, effStress, p_perp_perp_t, R_perp_perp_A,
     1  R_perp_At, p_perp_para_t, R_perp_para_A, p_perp_perp_c,
     1  p_perp_para_c, fE_gs4)        
          
          CALL getFE(theta3, effStress, p_perp_perp_t, R_perp_perp_A,
     1  R_perp_At, p_perp_para_t, R_perp_para_A, p_perp_perp_c,
     1  p_perp_para_c, fE_gs3)              
          
          if (fE_gs4 > fE_gs3) then
            theta1 = theta3
          else
            theta2 = theta4
          endif
          
          a_gs = (theta2-theta1)/(1.d0+gamma)
          b_gs = theta2-theta1-a_gs
          c_gs = a_gs/gamma
          theta3 = theta1 + a_gs
          theta4 = theta3 + c_gs  
          
          DIFF = fE_gs3 - fE_gs4
        enddo

        theta_max = theta2
        STATEV(17)=theta_max

        CALL getFE(theta_max, effStress, p_perp_perp_t, R_perp_perp_A,
     1   R_perp_At, p_perp_para_t, R_perp_para_A, p_perp_perp_c,
     1   p_perp_para_c, fE_IFF_max)           
        
        if(fE_IFF_max - 1.d0 > zero)then !check DAMAGE initiation 
          theta_fp = theta_max
          crackFlag = 1 
          STATEV(15) = crackFlag
          STATEV(16) = theta_fp  
        endif
       
        theta_rad=theta_max/pi*180.d0
        c = cos(theta_rad)
        s = sin(theta_rad)
        sigma_n = (c**2.d0)*effStress(2) + (s**2.d0)*effStress(3) 
     1 + 2.d0*c*s*effStress(6)  
         tau_nt = -c*s*effStress(2) + c*s*effStress(3) 
     1            + (c**2.d0-s**2.d0)*effStress(6)
         tau_nl = s*effStress(5) + c*effStress(4)        
        
        if(sigma_n > zero) then
          fE_IFFt_next = fE_IFF_max
          fE_IFFc_next = fE_IFFc
          flag=1
        else
          fE_IFFt_next = fE_IFFt
          fE_IFFc_next = fE_IFF_max
          flag=-1
        endif      
       else ! if damage was previously formed.
         theta_fp = STATEV(16) !failure angle is predefined      
         theta_rad = theta_fp*pi/180.d0
         c = cos(theta_rad)
         s = sin(theta_rad)
         sigma_n = c**2.d0*effStress(2) +s**2.d0 * effStress(3) 
     1             + 2.d0*c*s*effStress(6)
         tau_nt = -c*s*effStress(2) + c*s*effStress(3) 
     1            + (c**2.d0-s**2.d0)*effStress(6)
         tau_nl = s*effStress(5) + c*effStress(4)
         
         if( (tau_nt**2.d0) + (tau_nl**2.d0) < tauLimit) then
           p_over_R_t = 0.d0
           p_over_R_c = 0.d0
           c2=0.d0
           s2=0.d0
         else
           c2 = (tau_nt**2.d0) / ( (tau_nt**2.d0)+(tau_nl**2.d0) )
           s2 = (tau_nl**2.d0) / ( (tau_nt**2.d0)+(tau_nl**2.d0) )
           p_over_R_t = (p_perp_perp_t/R_perp_perp_A)*c2 +
     1  (p_perp_para_t/R_perp_para_A)*s2
           p_over_R_c = (p_perp_perp_c/R_perp_perp_A)*c2 +
     1  (p_perp_para_c/R_perp_para_A)*s2
         endif ! shear stresses check
    
         if(sigma_n > zero) then ! calculate fE due to the sign of sigma_n
           const1=(( (1.d0/R_perp_At)-p_over_R_t)*sigma_n)**2.d0
           const2=(tau_nt/R_perp_perp_A)**2.d0
           const3=(tau_nl/R_perp_para_A)**2.d0
           const4=p_over_R_t*sigma_n
           fE_IFF_cur=sqrt(const1+const2+const3)+ const4
           fE_IFF_cur = fE_IFF_cur/nu_w
           fE_IFFc_next = fE_IFFc
           flag = 1
         else
           const1=(p_over_R_c*sigma_n)**2.d0
           const2=(tau_nt/R_perp_perp_A)**2.d0
           const3=(tau_nl/R_perp_para_A)**2.d0
           const4=(p_over_R_c*sigma_n)
           fE_IFF_cur = sqrt(const1+const2+const3) +const4
           fE_IFF_cur = fE_IFF_cur/nu_w
           fE_IFFt_next = fE_IFFt
           flag = -1
         endif ! sign of sigma_n  
                      
         if(flag == 1) then
           fE_IFFt_next = fE_IFF_cur
         else
           fE_IFFc_next = fE_IFF_cur
         endif       
       endif! check if damage is initiated or not           
c###############################################################   
c#################### FF DAMAGE EVOLUATION #####################
c###############################################################           
       if(effStress(1) > zero)then
         flag_fc=0
         if(fE_FFt_next - rE_FFt < zero) then
           flag_ft=0
           del_ft_next = del_ft
           del_fc_next = del_fc
         else
           flag_ft = 1         
           A_ft = (2.0d0*Lc*(R_para_t**2.0d0))/
     1            (2.0d0*E1*GG_para_t - Lc*R_para_t**2.d0)  
           del_ft_next = 1.0d0 - (1.0d0/fE_FFt_next)*
     1                   exp(A_ft*(1.0d0-fE_FFt_next))
           del_fc_next = del_fc
         endif
      else ! compressive damage in fiber
        flag_ft=0
        if(fE_FFc_next - rE_FFc < zero) then
          flag_fc=0
          del_ft_next = del_ft
          del_fc_next = del_fc
        else
          flag_fc = 1        
          del_ft_next = del_ft
          A_fc = (2.0d0*Lc*(R_para_c**2.0d0))/
     1           (2.0d0*E1*GG_para_c - Lc*R_para_c**2.d0) 
          del_fc_next = 1.0d0 - (1.0d0/fE_FFc_next)*
     1                  exp(A_fc*(1.0d0-fE_FFc_next))
        endif
      endif !check effective stress' sign
! C PREVENT FIBER FAILURE      
      ! del_fc_next=0.d0
      ! del_ft_next=0.d0      
c###############################################################   
c################### IFF DAMAGE EVOLUATION #####################
c############################################################### 
      if(flag==1) then ! tensile damage in matrix
      !print*,'I am in matrix tension damage'
        flag_mc = 0
        if(fE_IFFt_next-rE_IFFt<zero) then
          flag_mt = 0        
          del_m1t_next = del_m1t
          del_m1c_next = del_m1c
          del_m2t_next = del_m2t
          del_m2c_next = del_m2c
        else
          flag_mt = 1
          A_m1t = 2.d0*Lc*(R_perp_t)**2.d0 /
     1     (2.d0*E2*GG_perp_t-Lc*(R_perp_t)**2.d0)
          del_m1t_next = 1.d0-(1.d0/fE_IFFt_next) *
     1     exp(A_m1t*(1.d0-fE_IFFt_next))
          del_m1c_next = del_m1c
          A_m2t = 2.d0*Lc*(R_perp_para)**2.d0 /
     1     (2.d0*G21*GG_perp_para-Lc*(R_perp_para)**2.d0)
          del_m2t_next = 1.d0-(1.d0/fE_IFFt_next) *
     1     exp(A_m2t*(1.d0-fE_IFFt_next))
          del_m2c_next = del_m2c
        endif! matrix damage in tension
      else  ! compression damage test 
       flag_mt = 0
       if(fE_IFFc_next-rE_IFFc<zero) then
          flag_mc = 0      
          del_m1t_next = del_m1t
          del_m1c_next = del_m1c
          del_m2t_next = del_m2t
          del_m2c_next = del_m2c   
        else
          flag_mc = 1  
          A_m1c = (2.d0*Lc*(R_perp_c)**2.d0) /
     1     ((2.d0*E2*GG_perp_c) - (Lc*(R_perp_c**2.d0)))
          del_m1c_next = 1.d0 - (1.d0 / fE_IFFc_next) *
     1     exp(A_m1c*(1.d0-fE_IFFc_next))
          del_m1t_next = del_m1t
          A_m2c = (2.d0*Lc*(R_perp_para**2.d0)) /
     1     ( 2.d0*G21*GG_perp_para - (Lc*(R_perp_para**2.d0)) )
          del_m2c_next = 1.d0 - (1.d0/fE_IFFc_next) *
     1     exp(A_m2c*(1.d0-fE_IFFc_next))
          del_m2t_next = del_m2t          
        endif!matrix damage in compression    
      endif!type of matrix damage     

      ! if(del_m1t_next.GT.0.d0.OR.del_m2t_next.GT.0.d0)then
      !   DTIME = 1e-5
      ! endif
c###############################################################   
c################### UPDATE STRESS #############################
c###############################################################   
c     testing the influence of fiber failure   
c viscous paramters      
      eta = 0.d0
      !eta=5e-3   !then viscous reg. is not used
      !eta=0.0002d0 !abaqus value.
      !eta=0.0003d0 !recommended value
     
      del_ft_v_next=DTIME/(eta+DTIME)*del_ft_next +
     1 eta/(eta+DTIME)*del_ft_v     
      
      del_fc_v_next=DTIME/(eta+DTIME)*del_fc_next +
     1 eta/(eta+DTIME)*del_fc_v  

      del_m1t_v_next=DTIME/(eta+DTIME)*del_m1t_next +
     1 eta/(eta+DTIME)*del_m1t_v  

      del_m1c_v_next=DTIME/(eta+DTIME)*del_m1c_next +
     1 eta/(eta+DTIME)*del_m1c_v  
      
      del_m2t_v_next=DTIME/(eta+DTIME)*del_m2t_next +
     1 eta/(eta+DTIME)*del_m2t_v   
      
      del_m2c_v_next=DTIME/(eta+DTIME)*del_m2c_next +
     1 eta/(eta+DTIME)*del_m2c_v  

      do i=1,6
        STRESS(i)=(1.d0-del_ft_v_next)*(1.d0-del_fc_v_next)*sigma_f(i)+
     1 (1.d0-del_m1t_v_next)*(1.d0-del_m1c_v_next)*sigma_m1(i)+
     2 (1.d0-del_m2t_v_next)*(1.d0-del_m2c_v_next)*sigma_m2(i)
      enddo
      
c###############################################################   
c################ CONSTRUCT C_tg TERMS #########################
c###############################################################
      do i=1,6
        do j=1,6
          C_secant(i,j)=
     1    (1.d0-del_ft_v_next)*(1.d0-del_fc_v_next)*C_f(i,j)+
     2    (1.d0-del_m1t_v_next)*(1.d0-del_m1c_v_next)*C_m1(i,j)+
     3    (1.d0-del_m2t_v_next)*(1.d0-del_m2c_v_next)*C_m2(i,j)
        enddo
      enddo 
      
c###############################################################   
c##################### UPDATE TANGENT ##########################
c###############################################################  
      do i=1,6
        do j=1,6
          DDSDDE(i,j) = C_secant(i,j)
        enddo
      enddo
c####################################################### 
c############### FIBER FAILURE #########################
c####################################################### 
c Shared terms for fiber damage related tangent terms     
      do i=1,6
        Nf1_f(i) = 0.d0
        Nf2_f(i) = 0.d0
        Nf3_f(i) = 0.d0
        C1_f(i) = 0.d0
        C2_f(i) = 0.d0
        C3_f(i) = 0.d0
      enddo
      
      Nf1_f(1)=1.d0
      Nf2_f(2)=1.d0
      Nf3_f(3)=1.d0
        
      do k=1,6
        do i=1,6
          C1_f(k) = C1_f(k) + Nf1_f(i)*CeSecond(i,k)
          C2_f(k) = C2_f(k) + Nf2_f(i)*CeSecond(i,k)
          C3_f(k) = C3_f(k) + Nf3_f(i)*CeSecond(i,k)
        enddo
      enddo!C terms      
      
c Tensile failure     
      if(flag_ft==1) then
        const1=A_ft*(1.d0-fE_FFt_next_next)
        const2=(A_ft*fE_FFt_next)+1.d0
        A_tf = exp(const1)*(const2)/(fE_FFt_next**2.d0)
        
        const3=(-nu_perp_para +
     1   (E_para/E_para_f)*nu_perp_para_f*m_sigma_f)
        B1_tf = 1.d0/R_para_t
        B2_tf =(1.d0/R_para_t)*const3
        B3_tf =(1.d0/R_para_t)*const3
        
        do i=1,6
          del_ft_part(i) = A_tf*( B1_tf*C1_f(i) +
     1     B2_tf*C2_f(i)+B3_tf*C3_f(i) )
        enddo
        
        do i=1,6
          do j=1,6
            C_tg_tf(i,j) = (1.d0-del_fc_v_next)*del_ft_part(i)*sigma_f(j)
          enddo
        enddo          
      endif! if flag_ft is 1
      
c Compressive failure   
      if(flag_fc==1) then
        const1=A_fc*(1.d0-fE_FFc_next)
        const2=(A_fc*fE_FFc_next)+1.d0
        A_cf = exp(const1)*(const2)/(fE_FFc_next**2.d0)
        
        const3= (-nu_perp_para +
     1   (E_para/E_para_f)*nu_perp_para_f*m_sigma_f)
        B1_cf = -1.d0/R_para_c
        B2_cf =-(1.d0/R_para_c)*const3
        B3_cf =-(1.d0/R_para_c)*const3
        
        do i=1,6
          del_fc_part(i) = A_cf*( B1_cf*C1_f(i) +
     1     B2_cf*C2_f(i)+B3_cf*C3_f(i) )
        enddo
        
        do i=1,6
          do j=1,6
            C_tg_cf(i,j) = (1.d0-del_ft_v_next)*del_fc_part(i)*sigma_f(j)
          enddo
        enddo          
      endif! if flag_fc is 1  
      
c####################################################### 
c############ INTER FIBER FAILURE ######################
c####################################################### 
 
c SHARED TERMS      
      if(flag_mt==1 .OR. flag_mc==1)then
        do i=1,6
          NN1(i)=0.d0
          NN2(i)=0.d0
          NN3(i)=0.d0
          CC1(i)=0.d0
          CC2(i)=0.d0
          CC3(i)=0.d0
        enddo  
      
        theta_fp = STATEV(16)
        theta_fp_rad=theta_fp*pi/180.d0
        c = cos(theta_fp_rad)
        s = sin(theta_fp_rad)
        
        NN1(2)=c**(2.d0)
        NN1(3)=s**(2.d0)
        NN1(6)=2.d0*(c*s)
        
        NN2(2)=-(c*s)
        NN2(3)=(c*s)
        NN2(6)=((c**2.d0)-(s**2.d0))
        
        NN3(4)=c
        NN3(5)=s       
  
        do k=1,6
          do i=1,6
            CC1(k) = CC1(k) + NN1(i)*CeSecond(i,k)
            CC2(k) = CC2(k) + NN2(i)*CeSecond(i,k)
            CC3(k) = CC3(k) + NN3(i)*CeSecond(i,k)
          enddo
        enddo     
      endif ! SHARED TERMS IN MATRIX FAILURE
      
c#################################################     
c########### FOR M1&M2 TENSILE DAMAGE ############
c#################################################     
      if(flag_mt==1) then
        const1=A_m1t*(1.d0-fE_IFFt_next)
        const2=((A_m1t*fE_IFFt_next)+1.d0)
        A_tm1 =exp(const1)*const2/(fE_IFFt_next**2.d0)           
        
        if( (tau_nt**2.d0+tau_nl**2.d0) < zero) then
          const1= sqrt( (sigma_n/R_perp_At)**2.d0 + 
     1     (tau_nt/R_perp_perp_A)**2.d0 + (tau_nl/R_perp_para_A)**2.d0 )
          BB1=sigma_n/ (R_perp_At**2.d0 * const1)
          BB2=tau_nt/ (R_perp_perp**2.d0 * const1)
          BB3=tau_nl/ (R_perp_para**2.d0 * const1)
          print*, "matrix tensile shear 0.d0"
        else        
          const1=(tau_nt**2.d0+tau_nl**2.d0)
          const2=(p_perp_para_t*(tau_nl**2.d0)) / 
     1     (R_perp_para_A*const1)
          const3=(p_perp_perp_t*(tau_nt**2.d0)) / 
     1     (R_perp_perp_A*const1)
          const4=( (1.d0/R_perp_At) - const2 - const3 )**2.d0
          const5=( (tau_nl/R_perp_para_A)**2.d0 +
     1     (tau_nt/R_perp_perp_A)**2.d0 )
          const6=sqrt(const5+(sigma_n**2.d0)*const4)
          BB1=const2+const3+(sigma_n*const4) / (const6)        
          
          const1=(tau_nt**2.d0+tau_nl**2.d0)
          const2=(2.d0*p_perp_para_t*(tau_nl**2.d0)*tau_nt) /
     1     (R_perp_para_A*(const1**2.d0))
          const3=(2.d0*p_perp_perp_t*(tau_nt**3.d0)) /
     1     (R_perp_perp_A*(const1**2.d0))
          const4=(2.d0*p_perp_perp_t*tau_nt) / 
     1     (R_perp_perp_A*const1)
          const5=2.d0*tau_nt/(R_perp_perp_A**2.d0)
          const11=2.d0*(sigma_n**2.d0)
          const6=1.d0/(R_perp_At)
          const7=(p_perp_para_t*(tau_nl**2.d0)) /
     1     (R_perp_para_A*const1)
          const8=(p_perp_perp_t*(tau_nt**2.d0)) /
     1     (R_perp_perp_A*const1)
          const9=(tau_nl/R_perp_para_A)**2.d0
          const10=(tau_nt/R_perp_perp_A)**2.d0
          const12=sigma_n*(-const2-const3+const4)
          const13=const5+const11*(const2+const3-const4)*
     1     (const6-const7-const8)
          const14=2.d0*sqrt(const9+const10+sigma_n**2.d0 * 
     1     (const6-const7-const8)**2.d0)
          BB2=const12+(const13/const14)      
          
          const1=(tau_nt**2.d0+tau_nl**2.d0)
          const2=(2.d0*p_perp_para_t*(tau_nl**3.d0)) /
     1     (R_perp_para_A*(const1**2.d0))
          const3=(2.d0*p_perp_perp_t*tau_nl*(tau_nt**2.d0)) /
     1     (R_perp_perp_A*(const1**2.d0))
          const4=(2.d0*p_perp_para_t*tau_nl) / 
     1     (R_perp_para_A*const1)
          const5=2.d0*tau_nl/(R_perp_para_A**2.d0)
          const6=2.d0*(sigma_n**2.d0)
          const7=1.d0/(R_perp_At)
          const8=(p_perp_para_t*(tau_nl**2.d0)) /
     1     (R_perp_para_A*const1)
          const9=(p_perp_perp_t*(tau_nt**2.d0)) /
     1     (R_perp_perp_A*const1)
          const10=(tau_nl/R_perp_para_A)**2.d0
          const11=(tau_nt/R_perp_perp_A)**2.d0
          const12=sigma_n*(-const2-const3+const4)
          const13=const5+const6*(const2+const3-const4)*
     1     (const7-const8-const9)
          const14=2.d0*sqrt(const10+const11+sigma_n**2.d0 * 
     1     (const7-const8-const9)**2.d0)
          BB3=const12+(const13/const14)               
        endif
        do i=1,6
          del_m1t_part(i) = A_tm1*( BB1*CC1(i) +
     1     BB2*CC2(i)+BB3*CC3(i) )
        enddo
        
        do i=1,6
          do j=1,6 ! next olabilir mi ?
            C_tg_tm1(i,j) = (1.d0-del_m1c_v_next)*del_m1t_part(i)*sigma_m1(j)
          enddo
        enddo    
c     m2 part
        const1=A_m2t*(1.d0-fE_IFFt_next)
        const2=((A_m2t*fE_IFFt_next)+1.d0)
        A_tm2=exp(const1)*const2/(fE_IFFt_next**2.d0)  
        
        do i=1,6
          del_m2t_part(i) = A_tm2*( BB1*CC1(i) +
     1     BB2*CC2(i)+BB3*CC3(i) )
        enddo
        
        do i=1,6
          do j=1,6
            C_tg_tm2(i,j) = (1.d0-del_m2c_v_next)*del_m2t_part(i)*sigma_m2(j)
          enddo
        enddo          
      endif! flag_mt is 1
      
c#################################################     
c######### FOR M1&M2 COMPRESSIVE DAMAGE ##########
c#################################################        
      if(flag_mc==1) then 
        const1=exp((1.d0-fE_IFFc_next)*A_m1c)
        const2=(1.d0+(fE_IFFc_next*A_m1c))
        A_cm1 =(const1*const2)/(fE_IFFc_next**2.d0)
        
        if( (tau_nt**2.d0+tau_nl**2.d0) < zero) then
          print*, "matrix compressive shear 0.d0"
        endif  
        const1=(tau_nt**2.d0)+(tau_nl**2.d0)
        const2=p_perp_para_c*(tau_nl**2.d0)  /
     1  (R_perp_para_A * const1)
        const3=p_perp_perp_c*(tau_nt**2.d0)  /
     1  (R_perp_perp_A * const1)
        const4=sigma_n*(const2+const3)**2.d0 ! **2.d0 sonradan eklendi
        const5=(tau_nl/R_perp_para_A)**2.d0 + 
     1  (tau_nt/R_perp_perp_A)**2.d0 
        const6=(sigma_n**2.d0)*(const2+const3)**2.d0
        const7=sqrt(const5+const6)
        BB1=(const2+const3)+(const4/const7)
        
        const1=(tau_nt**2.d0)+(tau_nl**2.d0)
        const2=(2.d0*p_perp_para_c*(tau_nl**2.d0)*tau_nt) /
     1   (R_perp_para_A*(const1**2.d0)) 
        const3=(2.d0*p_perp_perp_c*(tau_nt**3.d0)) /
     1   (R_perp_perp_A*(const1**2.d0))  
        const4=(2.d0*p_perp_perp_c*tau_nt) /
     1   (R_perp_perp_A*const1)   
        const5=(p_perp_para_c*(tau_nl**2.d0)) /
     1   (R_perp_para_A*const1)  
        const6=(p_perp_perp_c*(tau_nt**2.d0)) /
     1   (R_perp_perp_A*const1)  
        const7=(2.d0*tau_nt)/(R_perp_perp_A**2.d0)
        const8=(tau_nl/R_perp_para_A)**2.d0 +
     1   (tau_nt/R_perp_perp_A)**2.d0  
        const9=const7+2.d0*(sigma_n**2.d0)*(-const2-const3+const4) *
     1  (const5+const6)  
        const10=2.d0*sqrt(const8+(sigma_n**2.d0)*(const5+const6)**2.d0)! **2.d0  sonradan eklendi
        const11=sigma_n*(-const2-const3+const4)
        BB2=const11+(const9/const10)
        
        const1=(tau_nt**2.d0)+(tau_nl**2.d0)
        const2=(2.d0*p_perp_para_c*(tau_nl**3.d0)) /
     1   (R_perp_para_A*(const1**2.d0))     
        const3=(2.d0*p_perp_perp_c*tau_nl*(tau_nt**2.d0)) /
     1   (R_perp_perp_A*(const1**2.d0))   
        const4=(2.d0*p_perp_para_c*tau_nl) /
     1   (R_perp_para_A*const1)   
        const5=(p_perp_para_c*(tau_nl**2.d0)) /
     1   (R_perp_para_A*const1)
        const6=(p_perp_perp_c*(tau_nt**2.d0)) /
     1   (R_perp_perp_A*const1)    
        const7=(2.d0*tau_nl)/(R_perp_para_A**2.d0)
        const8=((tau_nl/R_perp_para_A)**2.d0) + 
     1   ((tau_nt/R_perp_perp_A)**2.d0)   
        const9=const7+2.d0*(sigma_n**2.d0)*(-const2-const3+const4) *
     1   (const5+const6)
        const10=2.d0*sqrt(const8+(sigma_n**2.d0)*(const5+const6)**2.d0)! **2.d0 sonradan eklendi
        const11=sigma_n*(-const2-const3+const4)
        BB3=const11+(const9/const10)

        do i=1,6
          del_m1c_part(i) = A_cm1*( BB1*CC1(i) +
     1     BB2*CC2(i)+BB3*CC3(i) )
        enddo
        
        do i=1,6
          do j=1,6
            C_tg_cm1(i,j) = (1.d0-del_m1t_v_next)*del_m1c_part(i)*sigma_m1(j)
          enddo
        enddo   
         
c m2 part        
        const1=exp((1.d0-fE_IFFc_next)*A_m2c)
        const2=(1.d0+(fE_IFFc_next*A_m2c))
        A_cm2 =(const1*const2)/(fE_IFFc_next**2.d0)
        
        do i=1,6
          del_m2c_part(i) = A_cm2*( BB1*CC1(i) +
     1     BB2*CC2(i)+BB3*CC3(i) )
        enddo
        
        do i=1,6
          do j=1,6
            C_tg_cm2(i,j) = (1.d0-del_m2t_v_next)*del_m2c_part(i)*sigma_m2(j)
          enddo
        enddo  
      endif! flag_mc is 1   
      
c###############################################################   
c##################### UPDATE TANGENT ##########################
c###############################################################       
      do i=1,6
        do j=1,6
          DDSDDE(i,j) = C_secant(i,j)
        enddo
      enddo

      if(flag_ft==1) then
       do i=1,6
         do j=1,6
           DDSDDE(i,j) = DDSDDE(i,j)-C_tg_tf(i,j)
         enddo
       enddo
      endif
      
      if(flag_fc==1) then
       do i=1,6
         do j=1,6
           DDSDDE(i,j) = DDSDDE(i,j)-C_tg_cf(i,j)
         enddo
       enddo
      endif     
      
      if(flag_mt==1) then
       do i=1,6
         do j=1,6
           DDSDDE(i,j) = DDSDDE(i,j)-C_tg_tm1(i,j)-C_tg_tm2(i,j)
         enddo
       enddo
      endif     
      
      if(flag_mc==1) then      
        do i=1,6
          do j=1,6
            DDSDDE(i,j) = DDSDDE(i,j)-C_tg_cm1(i,j)-C_tg_cm2(i,j)
          enddo
        enddo
      endif  

c###############################################################   
c##################### UPDATE TRESHOLDS ########################
c###############################################################  
c### Fiber tresholds ###
      if(effStress(1) >= 0.d0)then
        rE_FFt_next = max(1.0d0, fE_FFt_next, rE_FFt)
        rE_FFc_next = rE_FFc
      else
        rE_FFt_next = rE_FFt
        rE_FFc_next = max(1.0d0, fE_FFc_next, rE_FFc)
      endif  
c### Matrix tresholds ###   
      if(flag==1) then  
         rE_IFFt_next = max(1.d0, fE_IFFt_next, rE_IFFt )
         rE_IFFc_next = rE_IFFc
      elseif(flag==-1) then
         rE_IFFc_next = max(1.d0, fE_IFFc_next, rE_IFFc )
         rE_IFFt_next = rE_IFFt
      endif    
c###### SDV######    
      STATEV(1) = del_ft_next
      STATEV(2) = del_fc_next
      STATEV(3) = del_m1t_next 
      STATEV(4) = del_m1c_next
      STATEV(5) = del_m2t_next
      STATEV(6) = del_m2c_next
      
      STATEV(7) = rE_FFt_next
      STATEV(8) = rE_FFc_next
      STATEV(9) = rE_IFFt_next
      STATEV(10) = rE_IFFc_next
      
      STATEV(11) = fE_FFt_next
      STATEV(12) = fE_FFc_next
      STATEV(13) = fE_IFFt_next
      STATEV(14) = fE_IFFc_next
      
      STATEV(15) = crackFlag
      STATEV(16) = abs(theta_fp)    
      STATEV(17) = abs(theta_max)  
      
c     newly added viscously regularized damage variables      
      STATEV(18) = del_ft_v_next
      STATEV(19) = del_fc_v_next
      STATEV(20) = del_m1t_v_next
      STATEV(21) = del_m1c_v_next
      STATEV(22) = del_m2t_v_next
      STATEV(23) = del_m2c_v_next       
c###### SDV ######
      RETURN
      END
      
      subroutine get_C_elastic(CeSecond, E1, E2, E3, nu12, nu13, 
     1 nu23, G12, G13, G23)
        real*8, intent(out) :: CeSecond(6,6)
        real*8, intent(in) :: E1, E2, E3, nu12, nu13, 
     1   nu23, G12, G13, G23              
        real*8 D1111, D2222, D3333, D1122, D1133, 
     1   D2233, D1212, D1313, D2323   
        real*8 nu21, nu31, nu32, delta, const1
        integer i,j
!     nu_ij/Ei = nu_ji/Ej Poisson's ratios
        nu21 = E2*nu12/E1
        nu32 = E3*nu23/E2
        nu31 = E3*nu13/E1          
        
        const1=(1.d0-nu12*nu21-nu32*nu23-nu31*nu13-2.d0*nu21*nu32*nu13)
        delta = 1.d0/const1        
        D1111=E1*(1.d0-nu23*nu32)*delta
        D2222=E2*(1.d0-nu13*nu31)*delta
        D3333=E3*(1.d0-nu12*nu21)*delta
        D1122=E1*(nu21+nu31*nu23)*delta
        D1133=E1*(nu31+nu21*nu32)*delta
        D2233=E2*(nu32+nu12*nu31)*delta
        D1212 = G12
        D1313 = G13
        D2323 = G23
        
        do i=1,6
          do j=1,6
            CeSecond(i,j) = 0.d0
          enddo
        enddo   
        
        CeSecond(1,1)=D1111
        CeSecond(2,2)=D2222
        CeSecond(3,3)=D3333
        CeSecond(1,2)=D1122
        CeSecond(1,3)=D1133
        CeSecond(2,3)=D2233
        CeSecond(4,4)=D1212
        CeSecond(5,5)=D1313
        CeSecond(6,6)=D2323
        CeSecond(2,1)=CeSecond(1,2)
        CeSecond(3,1)=CeSecond(1,3)
        CeSecond(3,2)=CeSecond(2,3)      

      end subroutine get_C_elastic      

      subroutine get_C_slices(C_f, C_m1, C_m2, CeSecond )
        real*8, intent(in) :: CeSecond(6,6) ! input
        real*8, intent(out) :: C_f(6,6), C_m1(6,6), C_m2(6,6) ! output
        integer i, j
        do i=1,6
          do j=1,6
            C_f(i,j)=0.d0
            C_m1(i,j)=0.d0
            C_m2(i,j)=0.d0
          enddo
        enddo
        
        C_f(1,1) = CeSecond(1,1)
      
        C_m1(1,2) = CeSecond(1,2)
        C_m1(1,3) = CeSecond(1,3)
        C_m1(2,1) = CeSecond(2,1)
        C_m1(2,2) = CeSecond(2,2)
        C_m1(2,3) = CeSecond(2,3)
        C_m1(3,1) = CeSecond(3,1)
        C_m1(3,2) = CeSecond(3,2)
        C_m1(3,3) = CeSecond(3,3)
        C_m1(6,6) = CeSecond(6,6)
      
        C_m2(4,4) = CeSecond(4,4)
        C_m2(5,5) = CeSecond(5,5)
      end subroutine get_C_slices 
      
      subroutine getFE(theta, effStress, p_perp_perp_t, R_perp_perp_A,
     1  R_perp_At, p_perp_para_t, R_perp_para_A, p_perp_perp_c,
     1  p_perp_para_c, fE)
      
        real*8, intent(in) :: theta, effStress(6) ! input
        real*8, intent(out) :: fE ! output
        real*8, intent(in) :: p_perp_perp_t, R_perp_perp_A, R_perp_At,
     1  p_perp_para_t, R_perp_para_A, p_perp_perp_c, p_perp_para_c   
        
        real*8 theta_rad, pi, c, s, sigma_n, tau_nt, tau_nl, c2, s2,
     1   p_over_R_t, p_over_R_c
        real*8 TAU_LIMIT, ZERO, const1, const2, const3, const4
 
        
        TAU_LIMIT = 1e-12
        ZERO = 1e-12
        pi = 4.d0*datan(1.d0)
        theta_rad = theta*pi/180.d0
        c = cos(theta_rad)
        s = sin(theta_rad)
        
        sigma_n = (c**2.d0)*effStress(2) + (s**2.d0)*effStress(3) 
     1 + 2.d0*c*s*effStress(6)
        tau_nt = -(c*s*effStress(2)) + c*s*effStress(3) 
     1 + ((c**2.d0)-(s**2.d0))*effStress(6)
        tau_nl = s*effStress(5) + c*effStress(4)        
        
      if( (tau_nt**2.d0)+(tau_nl**2.d0) < TAU_LIMIT) then
        p_over_R_t = 0.d0
        p_over_R_c = 0.d0   
        c2=0.d0
        s2=0.d0
      else
        c2 = (tau_nt**2.d0)/( (tau_nt**2.d0) + (tau_nl**2.d0) )
        s2 = (tau_nl**2.d0)/( (tau_nt**2.d0) + (tau_nl**2.d0) )
        p_over_R_t = (p_perp_perp_t/R_perp_perp_A)*c2 +
     1  (p_perp_para_t/R_perp_para_A)*s2
        p_over_R_c = (p_perp_perp_c/R_perp_perp_A)*c2 +
     1  (p_perp_para_c/R_perp_para_A)*s2
      endif   
      
      if(sigma_n > ZERO) then
        const1=(( (1.d0/R_perp_At)-p_over_R_t)*sigma_n)**2.d0
        const2=(tau_nt/R_perp_perp_A)**2.d0
        const3=(tau_nl/R_perp_para_A)**2.d0
        const4=p_over_R_t*sigma_n
        fE= sqrt(const1+const2+const3)+const4
      else
        const1=(p_over_R_c*sigma_n)**2.d0
        const2=(tau_nt/R_perp_perp_A)**2.d0
        const3=(tau_nl/R_perp_para_A)**2.d0
        const4=p_over_R_c*sigma_n
        fE = sqrt(const1+const2+const3) + const4
      endif ! sign of sigma_n       
        
      end subroutine getFE       
            