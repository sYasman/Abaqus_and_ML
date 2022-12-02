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
      
      real*8 EE,nu
      integer i, j, k, ii, jj, kk
      real*8 U_k(NDOFEL,1)
      real*8 D_el(6,6)
      real*8 x(20,1), y(20,1), z(20,1), STRAIN(6,1), STRESS(6,1)
      
      EE = PROPS(1)
      nu = PROPS(2)
      
      if(JELEM.EQ.1)then
        noel_k=1
      endif
      
      print*,'This is NDOFEL', NDOFEL
      print*,'This is NNODE', NNODE
      print*,'This is MCRD', MCRD
      do i=1,NDOFEL
        U_k(i,1)=0.d0
      enddo  
      
      CALL getDel(EE, nu, D_el)
      
      do i=1,NNODE
        x(i,1) = COORDS(1,i)
        y(i,1) = COORDS(2,i)
        z(i,1) = COORDS(3,i)
      enddo
      
      do i=1,NDOFEL
        U_k(i,1) = U(i)
      enddo  
      if(LFLAGS(2).EQ.1)then
        print*,'This is not LARGE strain analysis.'
        print*,'Analysis is termianted'
        CALL XIT
      endif
      
      if(LFLAGS(3).EQ.1)then
c       Define RHS and AMATRX    
        do i=1,NDOFEL
          RHS(i,1) = 0.d0
          do j=1,NDOFEL
            AMATRX(i,j)=0.d0
          enddo!j
        enddo!i
        
        CALL getK(x,y,z,D_el, AMATRX)
        CALL getFint(x, y, z, U_k, D_el, RHS, STRAIN, STRESS)
        
        do i=1,6
          SVARS(i)=STRAIN(i,1)
          SVARS(i+6)=STRESS(i,1)
        enddo 
        
      endif
      
      if(LFLAGS(3).EQ.2)then   
c     define current Stiffness Matrix   
        do i=1,NDOFEL
          do j=1,NDOFEL
            AMATRX(i,j)=0.d0
          enddo!i
        enddo!j  
        CALL getK(x,y,z,D_el, AMATRX)
      endif
      
      if(LFLAGS(3).EQ.3)then
c     define current Damping Matrix     
        do i=1,NDOFEL
          do j=1,NDOFEL
            AMATRX(i,j)=0.d0
          enddo!i
        enddo!j  
      endif 
      
      if(LFLAGS(3).EQ.4)then  
c     define current Mass Matrix     
        do i=1,NDOFEL
          do j=1,NDOFEL
            AMATRX(i,j)=0.d0
          enddo!j
        enddo!i 
      endif  
      
      if(LFLAGS(3).EQ.5)then   
c     define current Residual Vector   
        do i=1,NDOFEL
          RHS(i,1)=0.d0
        enddo  
        CALL getFint(x, y, z, U_k, D_el, RHS, STRAIN, STRESS)
      endif    
      
      if(LFLAGS(3).EQ.6)then
c     define current Mass Matrix and Residual Vector
      
        do i=1,NDOFEL
          RHS(i,1)=0.d0
          do j=1,NDOFEL
            AMATRX(i,j)=0.d0
          enddo!i
        enddo!j  
        CALL getFint(x, y, z, U_k, D_el, RHS, STRAIN, STRESS)         
      endif
      
      RETURN
      END
      
      subroutine getDel(EE, nu, D_el)
        implicit none      
c     input and output      
        real*8, intent(out) :: D_el(6,6)
        real*8, intent(in) :: EE, nu   
      
c     local variables
        real*8 c1
        integer i, j

c     INITILIZATION     c        
        do i=1,6
          do j=1,6
            D_el(i,j)=0.d0
          enddo
        enddo        
c     INITILIZATION     c   
        c1 = EE / ( (1.d0+nu)*(1.d0-2.d0*nu) )
        
        D_el(1,1)=c1*(1.d0-nu)
        D_el(1,2)=c1*nu
        D_el(1,3)=c1*nu
        
        D_el(2,1)=c1*nu
        D_el(2,2)=c1*(1.d0-nu)
        D_el(2,3)=c1*nu
        
        D_el(3,1)=c1*nu
        D_el(3,2)=c1*nu
        D_el(3,3)=c1*(1.d0-nu)      
        
        D_el(4,4)=c1*(1.d0-2.d0*nu)/2.d0
        D_el(5,5)=c1*(1.d0-2.d0*nu)/2.d0
        D_el(6,6)=c1*(1.d0-2.d0*nu)/2.d0
      end subroutine getDel      
        
      subroutine getJacobian(x, y, z, xi, eta, tau, Jac, B_el, detJ)
      implicit none
      real*8, intent(out) :: Jac(3,3), B_el(6,60), detJ
      real*8, intent(in) :: x(20,1), y(20,1), z(20,1), xi, eta, tau
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
        XYZ(i,1) = x(i,1)
        XYZ(i,2) = y(i,1)
        XYZ(i,3) = z(i,1)
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
      
      subroutine getGaussData(n_gp, w, gp)
        real*8, intent(out) :: w(n_gp,1), gp(n_gp,1)
        integer, intent(in) :: n_gp
        integer i
        
        do i=1,n_gp
          w(i,1)=0.d0
          gp(i,1)=0.d0
        enddo 
        if(n_gp.EQ.1)then
          w(1,1)=2.d0
          gp(1,1)=0.d0
        elseif(n_gp.EQ.2)then
          w(1,1) = 1.d0
          w(2,1) = 1.d0
          gp(1,1) = -(1.d0/dsqrt(3.d0))
          gp(2,1) = (1.d0/dsqrt(3.d0))
        elseif(n_gp.EQ.3)then
          w(1,1)=(5.d0/9.d0)
          w(2,1)=(8.d0/9.d0)
          w(3,1)=(5.d0/9.d0)
          gp(1,1)=-dsqrt(3.d0/5.d0)
          gp(2,1)=dsqrt(0.d0)
          gp(3,1)=dsqrt(3.d0/5.d0)  
        endif
      end subroutine getGaussData
      
      subroutine getXi(xi, eta, tau, w)
        implicit none
        real*8, intent(out)  :: xi(8,1), eta(8,1), tau(8,1), w(8,1)
        
        real*8 val
        val = dsqrt(1.d0/3.d0)
        xi(1,1)=-val
        xi(2,1)=val
        xi(3,1)=-val
        xi(4,1)=val
        xi(5,1)=-val
        xi(6,1)=val
        xi(7,1)=-val
        xi(8,1)=val  

        eta(1,1)=-val
        eta(2,1)=-val
        eta(3,1)=val
        eta(4,1)=val
        eta(5,1)=-val
        eta(6,1)=-val
        eta(7,1)=val
        eta(8,1)=val
        
        tau(1,1)=-val
        tau(2,1)=-val
        tau(3,1)=-val
        tau(4,1)=-val
        tau(5,1)=val
        tau(6,1)=val
        tau(7,1)=val
        tau(8,1)=val  
        
        w(1,1)=1.d0
        w(2,1)=1.d0
        w(3,1)=1.d0
        w(4,1)=1.d0
        w(5,1)=1.d0
        w(6,1)=1.d0
        w(7,1)=1.d0
        w(8,1)=1.d0        
      end subroutine getXi      
      
      subroutine getK(x, y, z, D_el, K_el)
        implicit none
        real*8, intent(out) :: K_el(60, 60)
        real*8, intent(in) :: x(20,1), y(20,1), z(20,1), D_el(6,6)
        
        integer i_gp, j_gp, k_gp, i, j, k, l
        
        real*8 func1(60, 60), func(60, 60)
        real*8 xi_v(8,1), eta_v(8,1), tau_v(8,1), w_v(8,1)!v:vector
        real*8 Jac(3,3), B_el(6,60), detJ, dummy
        real*8 w(3,1), gp(3,1)!so far hard coded
        integer n_gp
        n_gp=3
        do i=1,60
          do j=1,60
            K_el(i,j) = 0.d0
            func(i,j) = 0.d0
          enddo
        enddo
        
        CALL getGaussData(n_gp, w, gp)
        
        do i_gp=1,n_gp
         do j_gp=1,n_gp
          do k_gp=1,n_gp
           CALL getJacobian(x, y, z, gp(i_gp,1), 
     1         gp(j_gp,1), gp(k_gp,1), Jac, B_el, detJ)
           do i=1,60
            do l=1,60
             do j=1,6
              do k=1,6
               func(i,l) = func(i,l) + w(i_gp,1)*w(j_gp,1)*w(k_gp,1)*
     1                      B_el(j,i)*D_el(j,k)*B_el(k,l)*detJ
              enddo!k
             enddo!j
            enddo!l
           enddo!i            
          enddo!k_gp
         enddo!j_gp
        enddo!i_gp 
    
        do i=1,60
          do j=1,60
            K_el(i,j)= func(i,j)
          enddo
        enddo          
      end subroutine getK        
      
      subroutine getFint(x, y, z, U_v, D_el,fint,el_strain,el_stress)
        implicit none
        INCLUDE 'KARC.blc'
        
        real*8, intent(out) :: fint(60,1), el_strain(6,1), 
     1          el_stress(6,1)
        real*8, intent(in) :: x(20,1), y(20,1), z(20,1), D_el(6,6), 
     1          U_v(60,1) !v:vector
        
        integer i_gp, j_gp, k_gp, i, j, k
        integer n_gp

        real*8 func(60,1), xi_v(8,1), eta_v(8,1), tau_v(8,1), w_v(8,1)
        real*8 Jac(3,3), B_el(6,60), detJ, dummy
        real*8 w(3,1), gp(3,1)
        n_gp=3
        CALL getGaussData(n_gp, w, gp)        
        do i=1,60
          func(i,1) = 0.d0
          fint(i,1) = 0.d0
        enddo
        
        k=1
        do i_gp=1,n_gp
         do j_gp=1,n_gp
          do k_gp=1,n_gp
           CALL getJacobian(x, y, z, gp(i_gp,1), gp(j_gp,1), 
     1    gp(k_gp,1), Jac, B_el, detJ)
           do i=1,6
            el_strain(i,1)=0.d0
           enddo  
           do i=1,6
            do j=1,60
             el_strain(i,1) = el_strain(i,1)+B_el(i,j)*U_v(j,1)
            enddo
           enddo  

c        el_stress=D_el(i,j)el_strain(j,1)          
           do i=1,6
            el_stress(i,1)=0.d0
           enddo  
           do i=1,6
            do j=1,6
             el_stress(i,1) = el_stress(i,1) +
     1                        D_el(i,j)*el_strain(j,1) 
            enddo
            STRESS_k(noel_k,k,i) = el_stress(i,1) 
            STRAIN_k(noel_k,k,i) = el_strain(i,1) 
           enddo
          k=k+1
          print*,'This is k', k
c         func=B'(i,j)el_stress(j,1)
           do i=1,60
            do j=1,6
             func(i,1)=func(i,1) + w(i_gp,1)*w(j_gp,1)*w(k_gp,1)*
     1                 B_el(j,i)*el_stress(j,1)*detJ
            enddo
           enddo

          enddo
         enddo
        enddo 
        
        do i=1,60
          fint(i,1) = -func(i,1)
        enddo  
        noel_k = noel_k+1
      end subroutine getFint   