      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'KARC.blc'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      

      integer i, j, k, ii, jj, kk
      real*8 U_k(60)
      real*8 D_el(6,6)
      real*8 x(20), y(20), z(20)
    
      real*8 pl_props(10), strain_ip(6), pl_stress(6)
      real*8 fint(60), K_el(60,60)
      real*8 hist_n(56), hist_new(56)
      real*8 tangent_puck(6,6)

!CCC MATERIAL PROPERTIES      
      real*8 e1p(3)

      pi = 4.d0*datan(1.d0) 

      e1p = [PROPS(8), PROPS(9), PROPS(10)] ! NO ORIENTATION

      if(JELEM.EQ.1)then
        noel_k=1
      endif

      do i=1,NPROPS ! load properties
        pl_props(i)=PROPS(i)
      enddo

      do i=1, NSVARS ! load history_n variables
        hist_n(i) = SVARS(i)
      enddo
      
      do i=1,56
        hist_new(i)=0.d0
      enddo

      do i=1,60 ! Initiation of RHS and AMATRX
        RHS(i,1) = 0.d0
        do j=1,60
          AMATRX(i,j)=0.d0
        enddo!j
      enddo!i

      do i=1,20 ! load coordinates
        x(i) = COORDS(1,i)
        y(i) = COORDS(2,i)
        z(i) = COORDS(3,i)
      enddo

      if(LFLAGS(2).EQ.1)then
        print*,'This is not LARGE strain analysis.'
        print*,'Analysis is termianted'
        CALL XIT()
      endif
      if(LFLAGS(3).EQ.1)then

        CALL get_fint(x, y, z, U, hist_n, pl_props,
     1 hist_new, fint, e1p)
        do i=1,NSVARS
          SVARS(i)=hist_new(i)
        enddo     
        CALL get_tang(x, y, z, D_el, hist_n, hist_new, pl_props, U, 
     1   K_el, e1p)

        do i=1,60
          RHS(i,1)=fint(i)
          do j=1,60
            AMATRX(i,j)=K_el(i,j)
          enddo
        enddo
        ! CALL XIT()
      endif
      
      if(LFLAGS(3).EQ.2)then           
C     PUCK BLOCK
        do i=1,n_pl_props
          pl_props(i)=PROPS(i)
        enddo       
        do i=1,NSVARS
          hist_n(i) = SVARS(i)
        enddo
        CALL get_tang(x, y, z, D_el, hist_n, hist_new, pl_props, U, 
     1   AMATRX, e1p)
C     END OF PUCK BLOCK    
      endif
      
      if(LFLAGS(3).EQ.3)then
c     define current Damping Matrix     
        do i=1,60
          do j=1,60
            AMATRX(i,j)=0.d0
          enddo!i
        enddo!j  
      endif 
      
      if(LFLAGS(3).EQ.4)then  
c     define current Mass Matrix, firstly this step is checked by solver     
        do i=1,60
          do j=1,60
            AMATRX(i,j)=0.d0
          enddo!j
        enddo!i 
      endif  
      
      if(LFLAGS(3).EQ.5)then   
c     define current Residual Vector   
        do i=1,60
          RHS(i,1)=0.d0
        enddo  
        do i=1,n_pl_props
          pl_props(i)=PROPS(i)
        enddo      
        do i=1,NSVARS
          hist_n(i) = SVARS(i)
        enddo
        CALL get_fint(x, y, z, U, hist_n, pl_props,
     1 hist_new, RHS, e1p)
      endif    
      
      if(LFLAGS(3).EQ.6)then
c     define current Mass Matrix and Residual Vector
        do i=1,60
          RHS(i,1)=0.d0
          do j=1,60
            AMATRX(i,j)=0.d0
          enddo!i
        enddo!j  

        do i=1,n_pl_props
          pl_props(i)=PROPS(i)
        enddo       
        do i=1,NSVARS
          hist_n(i) = SVARS(i)
        enddo
        CALL get_fint(x, y, z, U, hist_n, pl_props,
     1 hist_new, RHS, e1p)        
      endif
      
      RETURN
      END 
        
      subroutine getJacobian(x, y, z, xi, eta, tau, Jac, B_el, detJ)
      implicit none
      real*8, intent(out) :: Jac(3,3), B_el(6,60), detJ
      real*8, intent(in) :: x(20), y(20), z(20), xi, eta, tau
      integer i, j, k
      
      real*8 J_inv(3,3), N_X(3,20), Pmat(3,20), XYZ(20,3)
      real*8 dummy
      do i=1,3
        do j=1,3
          Jac(i,j) = 0.d0
          J_inv(i,j) = 0.d0
        enddo!j
      enddo!i
      
      do i=1,6
        do j=1,60
          B_el(i,j)=0.d0
        enddo!j
      enddo!i
      
      do i=1,3
        do j=1,20
          N_X(i,j) = 0.d0
          Pmat(i,j) = 0.d0
        enddo!j
      enddo!i
      
      do i=1,20
        do j=1,3
          XYZ(i,j) = 0.d0
        enddo!j
      enddo!j  
      
      do i=1,20
        XYZ(i,1) = x(i)
        XYZ(i,2) = y(i)
        XYZ(i,3) = z(i)
      enddo
      do i=1,3
        do j=1,20
          Pmat(i,j)=0.d0
        enddo
      enddo  
c     Pi : Ni,xi Ni,eta, Ni,tau, obtained from MATLAB
        Pmat(1,1) = (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0) +
     2        ((eta - 1.d0)*(tau - 1.d0)*(eta + tau + xi + 2.d0))/8.d0
        Pmat(1,2) = (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0) - 
     2  ((eta - 1.d0)*(tau - 1.d0)*(eta + tau - xi + 2.d0))/8.d0
        Pmat(1,3) = -((eta+1.d0)*(tau - 1.d0)*(eta - tau + xi - 2.d0))
     2  /8.d0 - (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0)
        Pmat(1,4) = - ((eta + 1.d0)*(tau - 1.d0)*(tau - eta + xi + 2.d0))
     2  /8.d0 - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0)
        Pmat(1,5) = - ((eta - 1.d0)*(tau + 1.d0)*(eta - tau + xi + 2.d0))
     2  /8.d0 - (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0)
        Pmat(1,6) = ((eta - 1.d0)*(tau + 1.d0)*(eta - tau - xi + 2.d0))
     2  /8.d0 - (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0)
        Pmat(1,7) = (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0) +
     2  ((eta + 1.d0)*(tau + 1.d0)*(eta + tau + xi - 2.d0))/8.d0
        Pmat(1,8) = (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0) -
     2  ((eta + 1.d0)*(tau + 1.d0)*(eta + tau - xi - 2.d0))/8.d0
        Pmat(1,9) = - (xi/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(tau - 1.d0) -
     2  ((eta - 1.d0)*(tau - 1.d0)*(xi + 1.d0))/4.d0
        Pmat(1,10) = (eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau - 1.d0)
        Pmat(1,11) = (xi/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau - 1.d0) +
     2  ((eta + 1.d0)*(tau - 1.d0)*(xi + 1.d0))/4.d0
        Pmat(1,12) = -(eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau - 1.d0)
        Pmat(1,13) = (xi/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(tau + 1.d0) +
     2  ((eta - 1.d0)*(tau + 1.d0)*(xi + 1.d0))/4.d0
        Pmat(1,14) = -(eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau + 1.d0)
        Pmat(1,15) = - (xi/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau + 1.d0)
     2  - ((eta + 1.d0)*(tau + 1.d0)*(xi + 1.d0))/4.d0
        Pmat(1,16) = (eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau + 1.d0)
        Pmat(1,17) = -(tau/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(tau + 1.d0)
        Pmat(1,18) = (tau/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(tau + 1.d0)
        Pmat(1,19) = -(tau/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau + 1.d0)
        Pmat(1,20) = (tau/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(tau + 1.d0)

        Pmat(2,1) = (xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)*(eta + tau + xi +
     2  2.d0) + (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0)
        Pmat(2,2) = - (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0) -
     2  (xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)*(eta + tau - xi + 2.d0)
        Pmat(2,3) = - (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0) -
     2  (xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)*(eta - tau + xi - 2.d0)
        Pmat(2,4) = (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0) -
     2  (xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)*(tau - eta + xi + 2.d0)
        Pmat(2,5) = - (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0) -
     2  (xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)*(eta - tau + xi + 2.d0)
        Pmat(2,6) = (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0) +
     2  (xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)*(eta - tau - xi + 2.d0)
        Pmat(2,7) = (xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)*(eta + tau +
     2  xi - 2.d0) + (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0)
        Pmat(2,8) = - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0)
     2  - (xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)*(eta + tau - xi - 2.d0)
        Pmat(2,9) = -(xi/4.d0 - 1.d0/4.d0)*(tau - 1.d0)*(xi + 1.d0)
        Pmat(2,10) = (eta/4.d0 - 1.d0/4.d0)*(tau - 1.d0)*(xi + 1.d0)
     2  + ((eta + 1.d0)*(tau - 1.d0)*(xi + 1.d0))/4.d0
        Pmat(2,11) = (xi/4.d0 - 1.d0/4.d0)*(tau - 1.d0)*(xi + 1.d0)
        Pmat(2,12) = - (eta/4.d0 - 1.d0/4.d0)*(tau - 1.d0)*(xi - 1.d0)
     2  - ((eta + 1.d0)*(tau - 1.d0)*(xi - 1.d0))/4.d0
        Pmat(2,13) = (xi/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi + 1.d0)
        Pmat(2,14) = - (eta/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi + 1.d0)
     2  - ((eta + 1.d0)*(tau + 1.d0)*(xi + 1.d0))/4.d0
        Pmat(2,15) = -(xi/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi + 1.d0)
        Pmat(2,16) = (eta/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi - 1.d0) 
     2  + ((eta + 1.d0)*(tau + 1.d0)*(xi - 1.d0))/4.d0
        Pmat(2,17) = -(tau/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi - 1.d0)
        Pmat(2,18) = (tau/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi + 1.d0)
        Pmat(2,19) = -(tau/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi + 1.d0)
        Pmat(2,20) = (tau/4.d0 - 1.d0/4.d0)*(tau + 1.d0)*(xi - 1.d0)

        Pmat(3,1) = (xi/8.d0-1.d0/8.d0)*(eta-1.d0)*(eta+tau+xi+2.d0) 
     2  + (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0)
        Pmat(3,2) = - (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau - 1.d0) 
     2  - (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(eta + tau - xi + 2.d0)
        Pmat(3,3) = (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0) 
     2  - (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(eta - tau + xi - 2.d0)
        Pmat(3,4) = - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau - 1.d0) 
     2  - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau - eta + xi + 2.d0)
        Pmat(3,5) = (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0) 
     2  - (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)*(eta - tau + xi + 2.d0)
        Pmat(3,6) = (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(eta - tau 
     2  - xi + 2.d0) - (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)*(tau + 1.d0)
        Pmat(3,7) = (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(eta + tau + xi 
     2  - 2.d0) + (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0)
        Pmat(3,8) = - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(tau + 1.d0) 
     2  - (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)*(eta + tau - xi - 2.d0)
        Pmat(3,9) = -(xi/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(xi + 1.d0)
        Pmat(3,10) = (eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi + 1.d0)
        Pmat(3,11) = (xi/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi + 1.d0)
        Pmat(3,12) = -(eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi - 1.d0)
        Pmat(3,13) = (xi/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(xi + 1.d0)
        Pmat(3,14) = -(eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi + 1.d0)
        Pmat(3,15) = -(xi/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi + 1.d0)
        Pmat(3,16) = (eta/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi - 1.d0)
        Pmat(3,17) = - (tau/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(xi - 1.d0) 
     2  - ((eta - 1.d0)*(tau + 1.d0)*(xi - 1.d0))/4.d0
        Pmat(3,18) = (tau/4.d0 - 1.d0/4.d0)*(eta - 1.d0)*(xi + 1.d0) 
     2  + ((eta - 1.d0)*(tau + 1.d0)*(xi + 1.d0))/4.d0
        Pmat(3,19) = - (tau/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi + 1.d0) 
     2  - ((eta + 1.d0)*(tau + 1.d0)*(xi + 1.d0))/4.d0
        Pmat(3,20) = (tau/4.d0 - 1.d0/4.d0)*(eta + 1.d0)*(xi - 1.d0) 
     2  + ((eta + 1.d0)*(tau + 1.d0)*(xi - 1.d0))/4.d0 
        
c       P_mat(i,j), XYZ(j,k)
        do i=1,3
          do j=1,3
            Jac(i,j)=0.d0
          enddo
        enddo  
        do i=1,3
          do k=1,3
            do j=1,20
              Jac(i,k) = Jac(i,k) + Pmat(i,j)*XYZ(j,k)
            enddo
          enddo
        enddo  

        detJ = Jac(1,1)*( Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2) ) +
     1    Jac(1,2)*( Jac(2,3)*Jac(3,1) - Jac(2,1)*Jac(3,3) ) +
     2    Jac(1,3)*( Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1) )
        
        do i=1,3
          do j=1,3
            J_inv(i,j)=0.d0
          enddo
        enddo  
        J_inv(1,1)= ( Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2) )/detJ
        J_inv(1,2)= ( Jac(1,3)*Jac(3,2) - Jac(1,2)*Jac(3,3) )/detJ
        J_inv(1,3)= ( Jac(1,2)*Jac(2,3) - Jac(1,3)*Jac(2,2) )/detJ
        J_inv(2,1)= ( Jac(2,3)*Jac(3,1) - Jac(2,1)*Jac(3,3) )/detJ
        J_inv(2,2)= ( Jac(1,1)*Jac(3,3) - Jac(1,3)*Jac(3,1) )/detJ
        J_inv(2,3)= ( Jac(1,3)*Jac(2,1) - Jac(1,1)*Jac(2,3) )/detJ
        J_inv(3,1)= ( Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1) )/detJ
        J_inv(3,2)= ( Jac(1,2)*Jac(3,1) - Jac(1,1)*Jac(3,2) )/detJ
        J_inv(3,3)= ( Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1) )/detJ 
        
c       J_inv(i,j)*Pmat(j,k)
        do i=1,3
          do j=1,20
            N_X(i,j)=0.d0
          enddo
        enddo  
        do i=1,3
          do k=1,20
            do j=1,3
              N_X(i,k) = N_X(i,k) + J_inv(i,j)*Pmat(j,k)
            enddo
          enddo
        enddo  
        
        k=1
        do i=1, 58, 3
          B_el(1,i)=N_X(1,k)
          
          B_el(2,i+1)=N_X(2,k)
          
          B_el(3,i+2)=N_X(3,k)

          B_el(4,i)=N_X(2,k)
          B_el(4,i+1)=N_X(1,k)
          
          B_el(5,i)=N_X(3,k)
          B_el(5,i+2)=N_X(1,k)
          
          B_el(6,i+1)=N_X(3,k)
          B_el(6,i+2)=N_X(2,k)
          k=k+1
        enddo  
      end subroutine getJacobian   
      
      subroutine getXi(xi, eta, tau, w)
        implicit none
        real*8, intent(out)  :: xi(8), eta(8), tau(8), w(8)
        
        real*8 val
        val = dsqrt(1.d0/3.d0)
        xi(1)=-val
        xi(2)=val
        xi(3)=-val
        xi(4)=val
        xi(5)=-val
        xi(6)=val
        xi(7)=-val
        xi(8)=val  

        eta(1)=-val
        eta(2)=-val
        eta(3)=val
        eta(4)=val
        eta(5)=-val
        eta(6)=-val
        eta(7)=val
        eta(8)=val
        
        tau(1)=-val
        tau(2)=-val
        tau(3)=-val
        tau(4)=-val
        tau(5)=val
        tau(6)=val
        tau(7)=val
        tau(8)=val  
        
        w(1)=1.d0
        w(2)=1.d0
        w(3)=1.d0
        w(4)=1.d0
        w(5)=1.d0
        w(6)=1.d0
        w(7)=1.d0
        w(8)=1.d0        
      end subroutine getXi      

      subroutine get_tang(x, y, z, D_el, hist_n, hist_new, pl_props, U_v,
     1    K_el, e1p)
        implicit none
        real*8, intent(out) :: K_el(60, 60)
        real*8, intent(in) :: x(20), y(20), z(20), D_el(6,6),
     1  hist_n(56), pl_props(7), U_v(60), hist_new(56), e1p(3)        
        
        integer i_gp, i, j, k, l, start_, end_
        real*8 func(60, 60), pl_tang(6,6), pl_stress(3,3)
        real*8 xi_v(8), eta_v(8), tau_v(8), w_v(8)!v:vector
        real*8 Jac(3,3), B_ip(6,60), detJ, dummy
        real*8 ip_strain(6)

C Below is used for the transformation DDSDDE from local to global
        real*8 e2p(3), e3p(3), Qt(3,3)
        real*8 strain(3,3), strainL(3,3), ip_strain_L(6)
        real*8 l1, l2, l3, m1, m2, m3, n1, n2, n3, Q_bar(6,6), p_tangent_G(6,6),
     1   Q2(6,6) 

        CALL getQ(e1p, e2p, e3p, Qt)

        do i=1,60
          do j=1,60
            K_el(i,j) = 0.d0
            func(i,j) = 0.d0
          enddo
        enddo
        
        CALL getXi(xi_v, eta_v, tau_v, w_v)
        
        do i_gp=1,8
          CALL getJacobian(x, y, z, xi_v(i_gp), 
     1         eta_v(i_gp), tau_v(i_gp), Jac, B_ip, detJ)
          
c         el_strain=B_el(i,j)*U_v(j,1)
          do i=1,6
            ip_strain(i)=0.d0
          enddo  
          do i=1,6
            do j=1,60
              ip_strain(i) = ip_strain(i)+B_ip(i,j)*U_v(j) ! Global strain in Voigt Notation
            enddo
          enddo   

          start_ = 7*(i_gp-1)+1
          end_ = start_+6   
          
        CALL get_stress( pl_props, hist_n(start_:end_), 
     1     ip_strain, hist_new(start_:end_), pl_stress, pl_tang)   
   
c         BT(i,j)D(j,k),B(k,l)    
          do i=1,60
            do l=1,60
              do j=1,6
                do k=1,6
                  func(i,l) = func(i,l) + w_v(i_gp) *
     1             B_ip(j,i)*pl_tang(j,k)*B_ip(k,l)*detJ
                enddo!k
              enddo!j
            enddo!l
          enddo!i 
        enddo!i_gp 
        
        do i=1,60
          do j=1,60
            K_el(i,j) = func(i,j)
          enddo
        enddo       
      end subroutine get_tang       
      
      subroutine get_fint(x, y, z, U_v, hist_n, pl_props,
     1 hist_new, fint, e1p)
!CCC TO DO nhist yap ona gore history var. ata     
        implicit none
        INCLUDE 'KARC.blc'
        
        real*8, intent(out) :: fint(60), hist_new(56)
        real*8, intent(in) :: x(20), y(20), z(20),
     1          U_v(60), hist_n(56), pl_props(7), e1p(3) !v:vector
        
        integer i_gp, i, j, k, l, start_, end_
        
        real*8 func(60), xi_v(8), eta_v(8), tau_v(8), w_v(8), 
     1      Jac(3,3), B_ip(6,60), detJ, ip_strain(6), pl_stress(6)
        real*8 e2p(3), e3p(3), Qt(3,3)
        real*8 strain(3,3), strainL(3,3), ip_strain_L(6)
        real*8 stress_L(3,3), stressG(3,3), ip_stress_G(6)
        real*8 pl_tang(6,6)
        CALL getQ(e1p, e2p, e3p, Qt) 
   
        do i=1,60
          func(i) = 0.d0
          fint(i) = 0.d0
        enddo
        
        CALL getXi(xi_v, eta_v, tau_v, w_v)
        do i_gp=1,8
          CALL getJacobian(x, y, z, xi_v(i_gp), eta_v(i_gp), 
     1    tau_v(i_gp), Jac, B_ip, detJ)
          
c         el_strain=B_el(i,j)*U_v(j,1)
          do i=1,6
            ip_strain(i)=0.d0
          enddo  
          do i=1,6
            do j=1,60
              ip_strain(i) = ip_strain(i)+B_ip(i,j)*U_v(j) ! ip_strain:Global strain in Voigt notation
            enddo
          enddo
     
          start_ = 7*(i_gp-1)+1
          end_ = start_+6

          CALL get_stress( pl_props, hist_n(start_:end_), 
     1     ip_strain, hist_new(start_:end_), pl_stress, pl_tang)

          do i=1,6
            STRESS_k(noel_k, i_gp, i) = pl_stress(i)
            STRAIN_k(noel_k, i_gp, i) = ip_strain(i)
          enddo  

          EPS_PL_K(noel_k, i_gp,1)=hist_new(start_)
          EPS_PL_K(noel_k, i_gp,2)=hist_new(start_+1)
          EPS_PL_K(noel_k, i_gp,3)=hist_new(start_+2)
          EPS_PL_K(noel_k, i_gp,4)=hist_new(start_+3)
          EPS_PL_K(noel_k, i_gp,5)=hist_new(start_+4)
          EPS_PL_K(noel_k, i_gp,6)=hist_new(start_+5)

          EPS_EQ_K(noel_k, i_gp)=hist_new(start_+5)

c         func=B'(i,j)el_stress(j,1)
          do i=1,60
            do j=1,6
              func(i)=func(i) + w_v(i_gp)*
     1                  B_ip(j,i)*pl_stress(j)*detJ
            enddo
          enddo       
        enddo! i_gp  
        noel_k = noel_k + 1
        do i=1,60
          fint(i) = -func(i)
        enddo  
      end subroutine get_fint   
      

      subroutine getQ(e1p, e2p, e3p, Qt)
    ! #################################################################    
    !   Subrotine, returns Transformatio nmatrix wrt given unit normal  
    !   Assumed GLOBAL unit normals
    !      |e3
    !      |
    !      |
    !      |__________e2
    !     /
    !    /
    !   /
    !  /e1             
    !  Regarding those LOCAL unit normals and Transformation matrix Q
    !  are computed.
    !  NEW LOCAL orientations unit normal wrt GLOBAL axis
    !  e2p, e3p LOCAL axis for 2' and 3' axis
    !  Q transformation matrix                
    ! #################################################################  
      implicit none
      real*8, intent(in) :: e1p(3)
      real*8, intent(out) :: e2p(3), e3p(3), Qt(3,3)
      real*8 e1(3), e2(3), e3(3), c1(3), c2
      integer i, j 
      e1 = [ 1.d0, 0.d0, 0.d0 ] !Unit vectors
      e2 = [ 0.d0, 1.d0, 0.d0 ]
      e3 = [ 0.d0, 0.d0, 1.d0 ]
      if(e1(1)-e1p(1) < 1e-6)then
        e2p = [ 0.d0, 1.d0, 0.d0 ]
        e3p = [ 0.d0, 0.d0, 1.d0 ]               
      else
      ! e'_3 = c1/c2,
      ! c1 = e1 x e'_1
      ! c2 = ||e1 x e'_1||     
        c1(1) = e1(2)*e1p(3) - e1(3)*e1p(2)    
        c1(2) = e1(3)*e1p(1) - e1(1)*e1p(3)    
        c1(3) = e1(1)*e1p(2) - e1(2)*e1p(1) 
        c2 = sqrt(c1(1)**2  + c1(2)**2 + c1(3)**2)  
        do i=1,3
          e3p(i)=c1(i)/c2
        enddo  

      ! e'_2 = e'_3 x e'_1
        e2p(1) = e3p(2)*e1p(3) - e3p(3)*e1p(2)    
        e2p(2) = e3p(3)*e1p(1) - e3p(1)*e1p(3)    
        e2p(3) = e3p(1)*e1p(2) - e3p(2)*e1p(1) 
      endif
      do i=1,3
        do j=1,3
          Qt(i,j)=0.d0
        enddo  
      enddo  
      Qt(1,1) = e1p(1)
      Qt(1,2) = e1p(2)
      Qt(1,3) = e1p(3)

      Qt(2,1) = e2p(1)
      Qt(2,2) = e2p(2)
      Qt(2,3) = e2p(3)      

      Qt(3,1) = e3p(1)
      Qt(3,2) = e3p(2)
      Qt(3,3) = e3p(3)      
      end subroutine getQ      



      subroutine get_stress( pl_props, hist_n, ip_strain_L,
     1         hist_new, pl_stress, pl_tang) !(start_:end_) 
      implicit none
!CCC Input Output      
      real*8, intent(in) :: pl_props(10), hist_n(7), ip_strain_L(6)
      real*8, intent(out) :: hist_new(7), pl_stress(6), pl_tang(6,6)

!CCC Local variables     
      real*8 a_t, b_t, c_t, d_t, a_c, b_c, c_c, d_c ! fit tparameters
      real*8 eps_pl_n2D(3,3), eps_pl_n(6), eps_eq_n
      real*8 I_1(3,3), II(3,3,3,3), PP(3,3,3,3) !IDENTITY AND DEVIATORIC PROJECTION TENSORS
      real*8 E, nu, sigma_c_0, sigma_t_0, sigma_c_inf, sigma_t_inf, nu_p ! material parameters
      real*8 lame, mu, kappa!LAMEs CONSTANTS
      real*8 alpha, kk ! material const.s
      real*8 e1p(3) ! unit normal of fiber axis 
      real*8 total_strain(6), total_strain2D(3,3), eps_el_tr(3,3)  
      real*8 eps_el_vol_tr, eps_el_dev_tr(3,3)
      real*8 stress_dev_tr(3,3), nr_stres_dev_tr, pres_tr, stress_tr(3,3) 
      real*8 sigma_c, sigma_t 
      real*8 phi_tr
      real*8 D_el(3,3,3,3)
      integer ind_1(6), ind_2(6)
      real*8 diff, TOL, gamma_cur, delta_eps_p(3,3), delta_eps_eq,
     1  eps_eq_new, phi, phi_pr    
      real*8 stress_dev_new(3,3), pres_new, nr_stres_dev_new
      integer iter_num, iter_max   
      real*8 eps_pl_new2D(3,3), stress_new(3,3)    
      real*8 t1(3,3), t2(3,3,3,3), t3, t4(3,3), t5(3,3), t6, t7,
     1 t8(3,3), t9(3,3)   
      real*8 D_ep_4D(3,3,3,3), D_ep(6,6)     
      integer i,j,k,l,m,n

!CCC START CODING 
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
        eps_pl_n(i) = hist_n(i)
      enddo
      eps_eq_n = hist_n(7)  
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
      
C##### READ PROPS AND COMPUTE LAMEs CONSTANTS #####         
      E = pl_props(1)
      nu = pl_props(2)
      sigma_c_0 = pl_props(3)
      sigma_t_0 = pl_props(4)
      sigma_c_inf = pl_props(5)
      sigma_t_inf = pl_props(6)
      nu_p = pl_props(7)  
      e1p(1) = pl_props(8)  
      e1p(2) = pl_props(9)  
      e1p(3) = pl_props(10)  

      alpha = (9.d0/2.d0)*(1.d0-2.d0*nu_p)/(1.d0+nu_p)
      kk = 1.d0/(1.d0+2.d0*(nu_p**2.d0))      

      lame = E*nu / ((1.d0+nu)*(1.d0-2.d0*nu))
      mu = E / (2.d0*(1.d0+nu))
      kappa = E/(3.d0*(1.d0-2.d0*nu))        
     
C##### COMPUTE INDETITIY TENSORS IN FULL SIZE (NOT VOIGT)         
      CALL getIandPterms(I_1, II, PP)   

C#### ELASTIC PREDICTOR ####
C#### UPDATE STRAN      
      do i=1,6
        total_strain(i)=0.d0
      enddo
      do i=1,6
        total_strain(i)=ip_strain_L(i)
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
      ! print*, 'phi_tr:', phi_tr
      if(phi_tr.LT.0.d0)then
C#### ELASTIC PREDICTION IS CORRECT    
        pl_stress(1)=stress_tr(1,1)
        pl_stress(2)=stress_tr(2,2)
        pl_stress(3)=stress_tr(3,3)
        pl_stress(4)=stress_tr(1,2)
        pl_stress(5)=stress_tr(1,3)
        pl_stress(6)=stress_tr(2,3)  
        
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
        pl_tang(i,j)=D_el(ind_1(i),ind_2(i),ind_1(j),ind_2(j))
          enddo
        enddo   
        
CCCCCC NO CHANGE IN SDV    
        hist_new(1)=eps_pl_n2D(1,1)
        hist_new(2)=eps_pl_n2D(2,2)
        hist_new(3)=eps_pl_n2D(3,3)
        hist_new(4)=2.d0*eps_pl_n2D(1,2)
        hist_new(5)=2.d0*eps_pl_n2D(1,3)
        hist_new(6)=2.d0*eps_pl_n2D(2,3)
        hist_new(7)=eps_eq_n
        ! print*, 'pl_stress:', pl_stress
      else  
        ! print*, 'PLASTIC PART IS NOT IMPLEMENTED YET' 
        ! CALL XIT()
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

      CALL normOf(stress_dev_new, nr_stres_dev_new)  
      phi = 3.d0*(nr_stres_dev_new**2.d0)+6.d0*pres_new*
     1 (sigma_c-sigma_t)-2.d0*(sigma_c*sigma_t)

      diff = ABS(phi)
      iter_num = iter_num + 1 

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
        ! CALL getInvariants(stress_new, stress_dev_new, I1, J2) 

        pl_stress(1)=stress_new(1,1)
        pl_stress(2)=stress_new(2,2)
        pl_stress(3)=stress_new(3,3)
        pl_stress(4)=stress_new(1,2)
        pl_stress(5)=stress_new(1,3)
        pl_stress(6)=stress_new(2,3)  
        
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
            pl_tang(i,j)=D_ep(i,j)
          enddo
        enddo   
CCCCC END OF CONSISTENT TANGENT MODULUS

        hist_new(1)=eps_pl_new2D(1,1)
        hist_new(2)=eps_pl_new2D(2,2)
        hist_new(3)=eps_pl_new2D(3,3)
        hist_new(4)=2.d0*eps_pl_new2D(1,2)
        hist_new(5)=2.d0*eps_pl_new2D(1,3)
        hist_new(6)=2.d0*eps_pl_new2D(2,3)
        hist_new(7)=eps_eq_new        
      endif
      end subroutine get_stress
      
    !   subroutine get_tangent(pl_props, hist_n, 
    !  1   hist_new, ip_strain_L, p_tangent, e1p)! material stiffness in LOCAL orient     
    !   implicit none
    !   real*8, intent(in) :: pl_props, hist_n, 
    !  1   hist_new, ip_strain_L, p_tangent, e1p
    !   ! real*8, intent(out) :: e2p(3), e3p(3), Qt(3,3)
      
    !   end subroutine get_tangent


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