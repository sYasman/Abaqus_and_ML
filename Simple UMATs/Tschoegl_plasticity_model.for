      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      ! implicit none
      INCLUDE 'ABA_PARAM.INC'
C Based on the study of Micromechanical analysis of polymer composites reinforced by unidirectional fibres: Part I-Constitutive modelling
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      integer i, j, k, l, m, n
      real*8 eps_pl_n(6), eps_pl_n2D(3,3) !PLASTIC STRAIN FROM PREVIOUS STEP
      real*8 eps_eq_n !EQUIVALENT PLASTIC STRAIN FROM PREVIOUS STEP
      real*8 I_1(3,3), II(3,3,3,3), PP(3,3,3,3) !IDENTITY AND DEVIATORIC PROJECTION TENSORS
      real*8 E, nu, w, sigma_c_0, sigma_t_0, sigma_c_inf, sigma_t_inf, nu_p
      real*8 lame, mu, kappa!LAMEs CONSTANTS
      real*8 total_strain(6), total_strain2D(3,3)
      real*8 eps_el_tr(3,3) ! ELASTIC TRIAL STRAIN
      real*8 eps_el_dev_tr(3,3) !DEVIATORIC PART OF ELASTIC TRIAL STRAIN
      real*8 eps_el_vol_tr !VOLUMETRIC PART OF ELASTIC TRIAL STRAIN
      real*8 stress_dev_tr(3,3), pres_tr !TRIAL DEVIATORIC STRESS AND PRESSURE
      real*8 nr_stres_dev_tr ! NORM OF TRIAL DEVIATORIC STRESS
      real*8 stress_tr(3,3)! TRIAL TOTAL STRESS
      real*8 Q, sigma_c, sigma_t! YIELD STRESSES 
      real*8 I1, J2 !invariants
      real*8 D_el(3,3,3,3)
      integer ind_1(6), ind_2(6)! vectors conints index number to convert 4D DDSDDDE to 2D
      real*8 alpha, kk ! material properties
      real*8 diff, TOL ! local N-R variables
      integer iter_num, iter_max
      real*8 phi_cur, phi_new, phi_pr
      real*8 gamma_cur, gamma_next
      real*8 const1, const2, const3 ! constants for phi_prime
      real*8 sc_pr, st_pr!derivatve of yield stress wrt d_gamma
      real*8 stress_dev_new(3,3), nr_stres_dev_new, pres_new,
     1   stress_new(3,3)
      real*8 delta_eps_p(3,3), delta_eps_eq 
      real*8 Q1, Q2(3,3), Q3(3,3), eps_eq_pr, const4(3,3),const5(3,3), 
     1  const6(3,3), const7(4,4)
      real*8 eps_pl_new2D(3,3) 
      real*8 eps_eq_new, eps_eq_new2
CCCC terms for DDSDDE in Plastic part      
      real*8 t1(3,3), t2(3,3,3,3), t3, t4(3,3), t5(3,3), t6, t7,
     1 t8(3,3), t9(3,3)
      real*8 D_ep_4D(3,3,3,3), D_ep(6,6) 
      real*8 a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c ! fit tparameters
      real*8 phi_left, phi_right
      real*8 phi_tr, eps_eq_newer, phi, phi_pl, phi_dm
        
C#### Parameters for the fit of sigma_t and sigma_c
      a_t=86.29d0
      b_t=3.775d0
      c_t=-61.29d0
      d_t=-241.d0

      a_c=114.2d0
      b_c=3.51d0
      c_c=-48.63d0
      d_c=-462.2d0      

C##### READ SDV FROM PREVIOS STEP #####     
      do i=1,6
        eps_pl_n(i)=0.d0
      enddo
      do i=1, 6
        eps_pl_n(i) = STATEV(i)
      enddo
      eps_eq_n = STATEV(7)  
C##### CONVERT SDV TO 2D FORMAT
      do i=1,3
        do j=1,3
          eps_pl_n2D(i,j)=0.d0
        enddo
      enddo
      eps_pl_n2D(1,1)=eps_pl_n(1)       
      eps_pl_n2D(2,2)=eps_pl_n(2)        
      eps_pl_n2D(3,3)=eps_pl_n(3)        
      eps_pl_n2D(1,2)=0.5d0*eps_pl_n(4)        
      eps_pl_n2D(1,3)=0.5d0*eps_pl_n(5)        
      eps_pl_n2D(2,3)=0.5d0*eps_pl_n(6)  
      eps_pl_n2D(2,1)=eps_pl_n2D(1,2)      
      eps_pl_n2D(3,1)=eps_pl_n2D(1,3)      
      eps_pl_n2D(3,2)=eps_pl_n2D(2,3)  

C##### COMPUTE INDETITIY TENSORS IN FULL SIZE (NOT VOIGT)         
      CALL getIandPterms(I_1, II, PP)   

C##### READ PROPS AND COMPUTE LAMEs CONSTANTS #####         
      E = PROPS(1)
      nu = PROPS(2)
      w = PROPS(3)
      sigma_c_0 = PROPS(4)
      sigma_t_0 = PROPS(5)
      sigma_c_inf = PROPS(6)
      sigma_t_inf = PROPS(7)
      nu_p = PROPS(8)   

      alpha = (9.d0/2.d0)*(1.d0-2.d0*nu_p)/(1.d0+nu_p)
      kk = 1.d0/(1.d0+2.d0*(nu_p**2.d0))      

      lame = E*nu / ((1.d0+nu)*(1.d0-2.d0*nu))
      mu = E / (2.d0*(1.d0+nu))
      kappa = E/(3.d0*(1.d0-2.d0*nu))  

C#### ELASTIC PREDICTOR ####
C#### UPDATE STRAN      
      do i=1,6
        total_strain(i)=0.d0
      enddo
      do i=1,6
        total_strain(i)=STRAN(i)+DSTRAN(i)
      enddo

C#### EXPAND TOTAL STRAIN TO 2D FORM      
      do i=1,3
        do j=1,3
          total_strain2D(i,j)=0.d0
        enddo
      enddo
      total_strain2D(1,1)=total_strain(1)       
      total_strain2D(2,2)=total_strain(2)        
      total_strain2D(3,3)=total_strain(3)        
      total_strain2D(1,2)=0.5d0*total_strain(4)        
      total_strain2D(1,3)=0.5d0*total_strain(5)        
      total_strain2D(2,3)=0.5d0*total_strain(6)  
      total_strain2D(2,1)=total_strain2D(1,2)      
      total_strain2D(3,1)=total_strain2D(1,3)      
      total_strain2D(3,2)=total_strain2D(2,3)

      do i=1,3 !TRIAL ELASTIC STRAIN
        do j=1,3
          eps_el_tr(i,j)=0.d0
        enddo
      enddo
      do i=1,3
        do j=1,3
          eps_el_tr(i,j) = total_strain2D(i,j) - eps_pl_n2D(i,j)
        enddo
      enddo      

!CCCC VOLUMETRIC ELASTIC TRIAL STRAIN I:total_strain2D
      eps_el_vol_tr = eps_el_tr(1,1) +
     1  eps_el_tr(2,2) + eps_el_tr(3,3)       

      do i=1,3 !DEVIATORIC PART OF TRIAL ELASTIC STRAIN, PP:total_strain2D
        do j=1,3
          eps_el_dev_tr(i,j) = 0.d0
        enddo
      enddo  
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
          eps_el_dev_tr(i,j) = eps_el_dev_tr(i,j) +
     1      PP(i,j,k,l) * eps_el_tr(k,l)
            enddo
          enddo
        enddo
      enddo
      
      do i=1,3
        do j=1,3
          stress_dev_tr(i,j)=0.d0
        enddo
      enddo
      do i=1,3
        do j=1,3 !DEVIATORIC TRIAL STRESS
          stress_dev_tr(i,j) = 2.d0*mu*eps_el_dev_tr(i,j)
        enddo
      enddo

    !   nr_stres_dev_tr = 0.d0 ! norm of trial deviatoric stress
    !   do i=1,3
    !     do j=1,3
    !       nr_stres_dev_tr=nr_stres_dev_tr+
    !  1     stress_dev_tr(i,j)*stress_dev_tr(i,j)
    !     enddo
    !   enddo
    !   nr_stres_dev_tr = SQRT(nr_stres_dev_tr)
    !   print*, 'nr_stres_dev_tr:', nr_stres_dev_tr
      
      CALL normOf(stress_dev_tr, nr_stres_dev_tr)  

      pres_tr = kappa*eps_el_vol_tr ! TRIAL PRESSURE

      do i=1,3
        do j=1,3
          stress_tr(i,j)=0.d0
        enddo
      enddo
      do i=1,3
        do j=1,3 ! TOTAL TRIAL STRESS
          stress_tr(i,j) = stress_dev_tr(i,j) + I_1(i,j)*pres_tr
        enddo
      enddo

      ! UPDATE sigma 
      CALL getQandSigma(eps_eq_n, a_t, b_t, c_t, d_t,  
     1   a_c, b_c, c_c, d_c, sigma_c, sigma_t )  

      phi_tr = 3.d0*(nr_stres_dev_tr**2.d0)+6.d0*pres_tr*
     1 (sigma_c-sigma_t)-2.d0*(sigma_c*sigma_t)

      if(phi_tr.LT.0.d0)then
C#### ELASTIC PREDICTION IS CORRECT    
        STRESS(1)=stress_tr(1,1)
        STRESS(2)=stress_tr(2,2)
        STRESS(3)=stress_tr(3,3)
        STRESS(4)=stress_tr(1,2)
        STRESS(5)=stress_tr(1,3)
        STRESS(6)=stress_tr(2,3)     
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                D_el(i,j,k,l)=0.d0
              enddo
            enddo
          enddo
        enddo
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
              D_el(i,j,k,l)=kappa*I_1(i,j)*I_1(k,l)+
     1  2.d0*mu*PP(i,j,k,l)         
              enddo
            enddo
          enddo
        enddo
        ind_1 = reshape( (/1, 2, 3, 1, 1, 2/), (/ 6/) )
        ind_2 = reshape( (/1, 2, 3, 2, 3, 3/), (/ 6/) )

        do i=1,6 !CONVERT DDSDDE TO VOIGT FORM (3,3,3,3)-->(6,6)
          do j=1,6
        DDSDDE(i,j)=D_el(ind_1(i),ind_2(i),ind_1(j),ind_2(j))
          enddo
        enddo   
        
CCCCCC NO CHANGE IN SDV    
        STATEV(1)=eps_pl_n2D(1,1)
        STATEV(2)=eps_pl_n2D(2,2)
        STATEV(3)=eps_pl_n2D(3,3)
        STATEV(4)=2.d0*eps_pl_n2D(1,2)
        STATEV(5)=2.d0*eps_pl_n2D(1,3)
        STATEV(6)=2.d0*eps_pl_n2D(2,3)
        STATEV(7)=eps_eq_n
      else    
! CCCCCC LOCAL N-R TO GET GAMMA_CUR     
        diff = 1.d0
        TOL = 1e-9
        iter_num = 1
        iter_max = 5000
        gamma_cur = 0.d0
        do i=1,3
          do j=1,3
            delta_eps_p(i,j) = 0.d0
          enddo
        enddo
        delta_eps_eq = 0.d0
        eps_eq_new = eps_eq_n
        phi = phi_tr

        do while(diff.GT.TOL.AND.iter_num<iter_max) ! local N-R for d_gamma

CCCCCC COMPUTE D(PHI)/D(GAMMA)     
          CALL getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr)

CCCCCC UPDATE GAMMA_CUR  
        gamma_cur = gamma_cur - phi/phi_pr

CCCCC Regarding gamma value compute new stresses and delta_eps_p
      do i=1,3
        do j=1,3    
          stress_dev_new(i,j) = stress_dev_tr(i,j)/
     1 (1.d0+6.d0*mu*gamma_cur)
        enddo
      enddo

      pres_new = pres_tr/
     1 (1.d0+2.d0*gamma_cur*alpha*kappa)

      do i=1,3
        do j=1,3
      delta_eps_p(i,j)=gamma_cur*
     1 (3.d0*stress_dev_new(i,j)+(2.d0/3.d0)*alpha*pres_new*I_1(i,j)) 
        enddo
      enddo 

      delta_eps_eq=0.d0
      do i=1,3
        do j=1,3
      delta_eps_eq=delta_eps_eq+delta_eps_p(i,j)*delta_eps_p(i,j)
        enddo
      enddo
      delta_eps_eq=SQRT(kk*delta_eps_eq)      
      eps_eq_new = eps_eq_n + delta_eps_eq
  
CCCCC Update sigma values wrt eps_eq_new   
      CALL getQandSigma(eps_eq_new, a_t, b_t, c_t, d_t,  
     1   a_c, b_c, c_c, d_c, sigma_c, sigma_t )  

    !   nr_stres_dev_new = 0.d0
    !   do i=1,3
    !     do j=1,3
    !       nr_stres_dev_new=nr_stres_dev_new+
    !  1      stress_dev_new(i,j)*stress_dev_new(i,j)
    !     enddo
    !   enddo
    !   nr_stres_dev_new = SQRT(nr_stres_dev_new) 

      CALL normOf(stress_dev_new, nr_stres_dev_new)  
      phi = 3.d0*(nr_stres_dev_new**2.d0)+6.d0*pres_new*
     1 (sigma_c-sigma_t)-2.d0*(sigma_c*sigma_t)

      diff = ABS(phi)
      iter_num = iter_num + 1 
      ! print*,'This is diff:', diff

      enddo !end of local N-R
CCCC END OF LOCAL N-R

CCCCCC COMPUTE D(PHI)/D(GAMMA)    
          CALL getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr)  

        if(ABS(phi_pr).LT.1e-12)then
          print*, 'phi_pr is almost zero...'
          CALL XIT()
        endif

        if(gamma_cur.LT.0.d0)then
          print*, 'gamma_cur is negative', gamma_cur
          CALL XIT()
        endif
      
        if(iter_num.EQ.iter_max)then
          print*, '*** EROOR ***'
          print*,'Too much iteration..'
          print*, '*** EROOR ***'
          CALL XIT()
        endif

CCCCC Regarding gamma value compute new stresses and delta_eps_p
        do i=1,3
          do j=1,3
        eps_pl_new2D(i,j)=eps_pl_n2D(i,j)+delta_eps_p(i,j)
          enddo
        enddo

        do i=1,3
          do j=1,3
        stress_new(i,j)=stress_dev_new(i,j)+I_1(i,j)*pres_new
          enddo
        enddo
        CALL getInvariants(stress_new, stress_dev_new, I1, J2) 

        STRESS(1)=stress_new(1,1)
        STRESS(2)=stress_new(2,2)
        STRESS(3)=stress_new(3,3)
        STRESS(4)=stress_new(1,2)
        STRESS(5)=stress_new(1,3)
        STRESS(6)=stress_new(2,3)   

CCCCC BEGINNING OF CONSISTENT TANGENT MODULUS
        do i=1,3
          do j=1,3 !d(p)/d(strain)
            t1(i,j)=kappa*I_1(i,j) /
     1       (1.d0+2.d0*kappa*alpha*gamma_cur)
          enddo
        enddo

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3 !d(S)/d(strain)
              t2(i,j,k,l)=2.d0*mu*PP(i,j,k,l)/
     1    (1.d0+6.d0*mu*gamma_cur)                
              enddo
            enddo
          enddo
        enddo

        !d(p)/d(d_gamma)
        t3 = pres_tr*(-2.d0*alpha*kappa)/
     1   ((1.d0+2.d0*gamma_cur*alpha*kappa)**2.d0)

        do i=1,3
          do j=1,3 !d(S)/d(d_gamma)
          t4(i,j)=stress_dev_tr(i,j)*(-6.d0*mu)/
     1     ((1.d0+6.d0*mu*gamma_cur)**2.d0)
          enddo
        enddo

        do i=1,3
          do j=1,3
            t5(i,j)=0.d0
          enddo
        enddo
        do i=1,3 ! denoted as 5.1
          do j=1,3
            do k=1,3
              do l=1,3!d(||S_tr||)/d(d_gamma)  ! norm of TRIAL STRESS
                t5(k,l)=t5(k,l) +
     1    (2.d0*mu*PP(k,l,i,j)*stress_dev_tr(i,j) / nr_stres_dev_tr)     
              enddo
            enddo
          enddo
        enddo

        t6 = 6.d0*nr_stres_dev_tr /
     1   ((1.d0+6.d0*mu*gamma_cur)**2.d0)

        t7 = 6.d0*(sigma_c-sigma_t) /
     1   (1.d0+2.d0*gamma_cur*alpha*kappa)
        ! print*, 'this is t7', t7

        do i=1,3
          do j=1,3
          t8(i,j) = t6*t5(i,j)+t7*(kappa*I_1(i,j))
          enddo
        enddo

        do i=1,3
          do j=1,3 !d(gamma_cur)/d(eps_new)
          t9(i,j)=(-1.d0/phi_pr)*t8(i,j)
          enddo
        enddo

CCCCC d(gamma_cur)/d(eps_new)
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
        D_ep_4D(i,j,k,l)=(t1(k,l)+t3*t9(k,l))*I_1(i,j) + 
     1    t2(i,j,k,l)+t4(i,j)*t9(k,l)  
              enddo
            enddo
          enddo
        enddo 
        ind_1 = reshape( (/1, 2, 3, 1, 1, 2/), (/ 6/) )
        ind_2 = reshape( (/1, 2, 3, 2, 3, 3/), (/ 6/) )

        do i=1,6 !CONVERT DDSDDE TO VOIGT FORM (3,3,3,3)-->(6,6)
          do j=1,6
    !     D_ep(i,j)=0.5d0*(D_ep_4D(ind_1(i),ind_2(i),ind_1(j),ind_2(j))+
    !  1     D_ep_4D(ind_1(j),ind_2(j),ind_1(i),ind_2(i))). Since D is not symmetric below is used.
        D_ep(i,j)=D_ep_4D(ind_1(i),ind_2(i),ind_1(j),ind_2(j))
          enddo
        enddo  
               
        do i=1,6 !CONVERT DDSDDE TO VOIGT FORM (3,3,3,3)-->(6,6)
          do j=1,6
            DDSDDE(i,j)=D_ep(i,j)
          enddo
        enddo   
CCCCC END OF CONSISTENT TANGENT MODULUS

        STATEV(1)=eps_pl_new2D(1,1)
        STATEV(2)=eps_pl_new2D(2,2)
        STATEV(3)=eps_pl_new2D(3,3)
        STATEV(4)=2.d0*eps_pl_new2D(1,2)
        STATEV(5)=2.d0*eps_pl_new2D(1,3)
        STATEV(6)=2.d0*eps_pl_new2D(2,3)
        STATEV(7)=eps_eq_new
      endif
      RETURN
      END
C###################################################
C################### SUBROUTINES ################### 
C###################################################   
      subroutine getIandPterms(I_1, II, PP)     
!CCCC SUBROTINE RETURNS IDENTITY TENSORS (2D AND 4D) AND DEVIATORIC PROJECTION TENSOR(4D)
!CCCC I_1 -> SECOND ORDER IDENTITY TENSOR        
!CCCC II -> FOURTH ORDER IDENTITY TENSOR        
!CCCC PP -> FOURTH ORDER DEVIATORIC PROJECTION TENSOR       
        implicit none
        real*8, intent(out) :: I_1(3,3), II(3,3,3,3), PP(3,3,3,3)
        integer i, j, k, l, m, n
        do i=1,3
          do j=1,3
            I_1(i,j)=0.d0
            do k=1,3
              do l=1,3
                II(i,j,k,l)=0.d0
                PP(i,j,k,l)=0.d0
              enddo
            enddo
          enddo
        enddo
        I_1(1,1)=1.d0
        I_1(2,2)=1.d0
        I_1(3,3)=1.d0

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                II(i,j,k,l)=0.5d0*(I_1(i,k)*I_1(j,l)+I_1(i,l)*I_1(j,k))
                ! PP(i,j,k,l)=II(i,j,k,l)-(1.d0/3.d0)*I_1(i,j)*I_1(k,l)
              enddo
            enddo
          enddo
        enddo

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                ! II(i,j,k,l)=0.5d0*(I_1(i,k)*I_1(j,l)+I_1(i,l)*I_1(j,k))
                PP(i,j,k,l)=II(i,j,k,l)-(1.d0/3.d0)*I_1(i,j)*I_1(k,l)
              enddo
            enddo
          enddo
        enddo        
      endsubroutine getIandPterms

      subroutine getQandSigma(eps_eq_n, a_t, b_t, c_t, d_t, 
     1   a_c, b_c, c_c, d_c, sigma_c, sigma_t )  
!CCCC SUBROUTINE RETURNS YIELD STRESSES IN TENSION AND COMPRESSION
!CCCC INPUTS : sigma_c_0,sigma_c_inf, sigma_t_0,sigma_t_inf--> STERNGTH VALUES      
!CCCC INPUTS : w-->FITTING PARAMETER
!CCCC INPUTS : eps_eq_n -->EQUIVALENT PLASTIC STRAIN FROM PREVIOUS STEP    
!CCCC OUPUTS : sigma_c, sigma_t --> NEW YIELD STRESSES     
!CCCC OUPUTS : Q --> MULTIPLER OF STRESS DIFFERENCE     
      implicit none
      real*8, intent(in) :: eps_eq_n, a_t, b_t, c_t, d_t, a_c, 
     1 b_c, c_c, d_c
      real*8, intent(out) ::sigma_c, sigma_t 
    
      sigma_t = a_t*EXP(b_t*eps_eq_n)+c_t*EXP(d_t*eps_eq_n)       
      sigma_c = a_c*EXP(b_c*eps_eq_n)+c_c*EXP(d_c*eps_eq_n)         

      endsubroutine getQandSigma

      subroutine getInvariants(A, A_dev, I1, J2)  
!CCCC SUBROUTINE RETURNS INVARIANTS OF 3by3 STRESS MATRIX
!CCCC I1: 1ST INVARIANT OF TOTAL STRESS MATRIX
!CCCC J2: 2ND INVARIANT OF DEVIATORIC STRESS MATRIX                 
!CCCC based on 
!CCCC https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Invariants_of_the_stress_deviator_tensor             
      implicit none
      real*8, intent(in) :: A(3,3), A_dev(3,3)
      real*8, intent(out) :: I1, J2  
      integer i, j
      I1 = 0.d0
      J2 = 0.d0
      do i=1,3
        I1 = I1 + A(i,i)
        do j=1,3  
          J2 = J2+0.5d0*(A_dev(i,j)*A_dev(i,j))
        enddo
      enddo
      endsubroutine getInvariants     

      subroutine normOf(A, nrA)  
!CCCC SUBROUTINE RETURNS NORM OF 3by3 MATRIX    
!CCCC nrT = sqrt(T_ij T_ij)   
      implicit none
      real*8, intent(in) :: A(3,3)
      real*8, intent(out) :: nrA  
      integer i, j

      nrA = 0.d0
      do i=1,3
        do j=1,3
          nrA = nrA + A(i,j)*A(i,j)
        enddo
      enddo
      nrA = SQRT(nrA)
      end subroutine normOf

      subroutine getPhi(d_gamma, stress_dev_tr, pres_tr, mu, kappa, alpha,
     1 kk, eps_eq_n, phi)  
!CCCC SUBROUTINE RETURNS PHI VALUE  
      implicit none
      real*8, intent(in) :: d_gamma, stress_dev_tr(3,3), pres_tr, mu, kappa,
     1   alpha, kk, eps_eq_n
      real*8, intent(out) :: phi
      integer i, j, k, l, m, n

      real*8 stress_dev_new(3,3), pres_new, stress_new(3,3) ! stress terms
      real*8 nr_stres_dev_new ! norm of deviatoric stress
      real*8 I_1(3,3), II(3,3,3,3), PP(3,3,3,3) ! identity and projection tensors
      real*8 I1, J2 !Stress invariants
      real*8 delta_eps_p(3,3) !increment in plasstic strain
      real*8 delta_eps_eq !increment in equivalent plastic strain
      real*8 Q, sigma_c, sigma_t !yield strengths parameters
      real*8 eps_eq_new

CCCC CODING STARTS HERE CCCC  
      do i=1,3
        do j=1,3    
          stress_dev_new(i,j) = stress_dev_tr(i,j)/
     1 (1.d0+6.d0*mu*d_gamma)
        enddo
      enddo
      CALL normOf(stress_dev_new, nr_stres_dev_new)

      pres_new = pres_tr/
     1 (1.d0+2.d0*d_gamma*alpha*kappa)

      CALL getIandPterms(I_1, II, PP)  
      do i=1,3
        do j=1,3
          stress_new(i,j)=stress_dev_new(i,j)+pres_new*I_1(i,j)
        enddo
      enddo

      CALL getInvariants(stress_new, stress_dev_new, I1, J2) 

      CALL getDeltaEpsP(d_gamma, stress_dev_tr, pres_tr,kk, 
     1 mu, alpha, kappa, I_1, eps_eq_n, delta_eps_p,  
     2 delta_eps_eq, eps_eq_new) 

      CALL getQandSigma(eps_eq_new, sigma_c, sigma_t ) 

      phi = 3.d0*(nr_stres_dev_new**2.d0)+6.d0*pres_new*
     1 (sigma_c-sigma_t)-2.d0*(sigma_c*sigma_t)
      end subroutine getPhi
  
      subroutine getDeltaEpsP(d_gamma, stress_dev_tr, pres_tr,kk, 
     1 mu, alpha, kappa, I_1, eps_eq_n, delta_eps_p,  
     2 delta_eps_eq, eps_eq_new)  
!CCCC SUBROUTINE RETURNS PHI VALUE  
      implicit none
      real*8, intent(in) :: d_gamma, stress_dev_tr(3,3), pres_tr, 
     1   kk, mu, alpha, kappa, I_1(3,3), eps_eq_n
      real*8, intent(out) :: delta_eps_p(3,3), delta_eps_eq, eps_eq_new
      integer i, j
      real*8 stress_dev_new(3,3), pres_new, I1, J2, stress_new(3,3)
CCCCC Regarding gamma value compute new stresses and delta_eps_p
      do i=1,3
        do j=1,3    
          stress_dev_new(i,j) = stress_dev_tr(i,j)/
     1 (1.d0+6.d0*mu*d_gamma)
        enddo
      enddo

      pres_new = pres_tr/
     1 (1.d0+2.d0*d_gamma*alpha*kappa)

      do i=1,3
        do j=1,3
          stress_new(i,j)=stress_dev_new(i,j)+I_1(i,j)*pres_new
        enddo
      enddo

      CALL getInvariants(stress_new, stress_dev_new, I1, J2) 

      do i=1,3
        do j=1,3
      delta_eps_p(i,j)=d_gamma*
     1 (3.d0*stress_dev_new(i,j)+(2.d0/3.d0)*alpha*pres_new*I_1(i,j)) 
        enddo
      enddo 

      delta_eps_eq=0.d0
      do i=1,3
        do j=1,3
      delta_eps_eq=delta_eps_eq+delta_eps_p(i,j)*delta_eps_p(i,j)
        enddo
      enddo
      delta_eps_eq=SQRT(kk*delta_eps_eq)      
      eps_eq_new = eps_eq_n + delta_eps_eq
      end subroutine getDeltaEpsP

      subroutine getPhiPr(nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1 alpha, kappa, kk, delta_eps_p, delta_eps_eq, stress_dev_tr,
     2 I_1, a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c, eps_eq_new,
     3 sigma_c, sigma_t, iter_num, phi_pr)  
    
!CCCC SUBROUTINE RETURNS PHI_PR VALUE  
      implicit none
      real*8, intent(in) :: nr_stres_dev_tr, mu, gamma_cur, pres_tr,
     1  alpha, kappa, kk, delta_eps_p(3,3), delta_eps_eq, 
     2  stress_dev_tr(3,3), I_1(3,3), a_t, b_t, c_t, d_t, a_c, 
     3  b_c, c_c, d_c, eps_eq_new, sigma_c, sigma_t
      integer, intent(in) :: iter_num
      real*8 eps_eq_pr       
      real*8, intent(out) :: phi_pr
      real*8 const1, const2, const3, const4(3,3), const5(3,3),
     1  const6(3,3), const7(3,3), Q1, Q2(3,3), Q3(3,3),
     2  st_pr, sc_pr
      integer i, j    
      
CCCCCC COMPUTE D(PHI)/D(GAMMA)     
        const1 = -36.d0*(nr_stres_dev_tr**2.d0)*mu/
     1    ((1.d0+6.d0*mu*gamma_cur)**3.d0) 
        const2 = -12.d0*pres_tr*(alpha*kappa)/
     1   ((1.d0+2.d0*gamma_cur*alpha*kappa)**2.d0)
        const3 = 6.d0*pres_tr/
     1    (1.d0+2.d0*gamma_cur*alpha*kappa)  
        Q1 = 1.d0
        do i=1,3
          do j=1,3
            if(iter_num.EQ.1)then
              Q2(i,j) = 0.d0
            else
              Q2(i,j) = kk * delta_eps_p(i,j)/delta_eps_eq 
            endif
          const4(i,j)=3.d0*stress_dev_tr(i,j)/(1.d0+6.d0*mu*gamma_cur)
          const5(i,j)=(2.d0/3.d0)*alpha*
     1      (pres_tr/(1.d0+2.d0*gamma_cur*alpha*kappa))*I_1(i,j)
          const6(i,j)=(-18.d0)*stress_dev_tr(i,j)*mu/
     1      ((1.d0+6.d0*mu*gamma_cur)**2.d0)
          const7(i,j)=(-4.d0/3.d0)*(alpha**2.d0)*pres_tr*kappa /
     1      ((1.d0+2.d0*gamma_cur*alpha*kappa)**2.d0)*I_1(i,j)    
          enddo
        enddo
        do i=1,3
          do j=1,3
          Q3(i,j)=(const4(i,j)+const5(i,j)) + 
     1      gamma_cur*(const6(i,j)+const7(i,j))
          enddo
        enddo
        eps_eq_pr = 0.d0
        do i=1,3
          do j=1,3 ! d(eps_eq_cur)/d(gamma_cur)
            eps_eq_pr = eps_eq_pr + Q1*Q2(i,j)*Q3(i,j)
          enddo
        enddo

        st_pr = eps_eq_pr*
     1   (a_t*b_t*EXP(b_t*eps_eq_new)+c_t*d_t*EXP(d_t*eps_eq_new))
        sc_pr = eps_eq_pr*
     1   (a_c*b_c*EXP(b_c*eps_eq_new)+c_c*d_c*EXP(d_c*eps_eq_new))   

        phi_pr = const1 + const2*(sigma_c-sigma_t) +
     1    const3*(sc_pr-st_pr) - 
     2    2.d0*(sc_pr*sigma_t + sigma_c*st_pr)        
      end subroutine getPhiPr      