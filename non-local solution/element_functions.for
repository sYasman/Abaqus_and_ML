      subroutine getJacobian(x, y, z, xi, eta, tau, Jac, B_el, detJ)
!CCCC JACOBIAN AND B MATRIX FOR STRESS RELATED PARTS FOR PRIMARY VARIABLE
        implicit none
C OUTPUT        
        real*8, intent(out) :: Jac(3,3), B_el(6,24), detJ
C INPUT        
        real*8, intent(in) :: x(8), y(8), z(8), xi, eta, tau
C LOCAL VARIABLES        
        integer i, j, k
        real*8 J_inv(3,3), N_X(3,8), Pmat(3,8), XYZ(8,3)
C
        do i=1,3
          do j=1,3
            Jac(i,j) = 0.d0
            J_inv(i,j) = 0.d0
          enddo!j
        enddo!i
C
        do i=1,6
          do j=1,24
            B_el(i,j)=0.d0
          enddo!j
        enddo!i
C
        do i=1,3
          do j=1,8
            N_X(i,j) = 0.d0
            Pmat(i,j) = 0.d0
          enddo!j
        enddo!i
C
        do i=1,8
          do j=1,3
            XYZ(i,j) = 0.d0
          enddo!j
        enddo!j
        do i=1,8
          XYZ(i,1) = x(i)
          XYZ(i,2) = y(i)
          XYZ(i,3) = z(i)
        enddo
c     Pi : Ni,xi Ni,eta, Ni,tau, obtained from MATLAB
        Pmat(1,1) = -((eta - 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,2) = ((eta - 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,3) = -((eta + 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,4) = ((eta + 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,5) = ((eta - 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,6) = -((eta - 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,7) = ((eta + 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,8) = -((eta + 1.d0)*(tau + 1.d0))/8.d0
C
        Pmat(2,1) = -(xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,2) = (xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,3) = -(xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,4) = (xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,5) = (xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,6) = -(xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,7) = (xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,8) = -(xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)
C
        Pmat(3,1) = -(xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,2) = (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,3) = -(xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,4) = (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,5) = (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,6) = -(xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,7) = (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,8) = -(xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)
c       P_mat(i,j), XYZ(j,k)
        do i=1,3
          do k=1,3
            do j=1,8
              Jac(i,k) = Jac(i,k) + Pmat(i,j)*XYZ(j,k)
            enddo
          enddo
        enddo
C
        detJ = Jac(1,1)*( Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2) ) +
     1    Jac(1,2)*( Jac(2,3)*Jac(3,1) - Jac(2,1)*Jac(3,3) ) +
     2    Jac(1,3)*( Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1) )
C
        ! if(abs(detJ).LT.1e-12)then
        !   print*, '*** Error ***'
        !   print*, 'Jacobian is almost zero...'
        !   print*, '*** Error ***'
        !   CALL XIT()
        ! endif
C
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
          do k=1,8
            do j=1,3
              N_X(i,k) = N_X(i,k) + J_inv(i,j)*Pmat(j,k)
            enddo
          enddo
        enddo
C
        k=1
        do i=1, 22, 3
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
!CCCC GAUSS POINTS FOR STRESS RELATED ELEMENTS FOR PRIMARY VARIABLE
        implicit none
C OUTPUT        
        real*8, intent(out)  :: xi(8), eta(8), tau(8), w(8)
C LOCAL
        real*8 val
C        
        val = dsqrt(1.d0/3.d0)
        xi(1)=-val
        xi(2)=val
        xi(3)=-val
        xi(4)=val
        xi(5)=-val
        xi(6)=val
        xi(7)=-val
        xi(8)=val
C
        eta(1)=-val
        eta(2)=-val
        eta(3)=val
        eta(4)=val
        eta(5)=-val
        eta(6)=-val
        eta(7)=val
        eta(8)=val
C
        tau(1)=-val
        tau(2)=-val
        tau(3)=-val
        tau(4)=-val
        tau(5)=val
        tau(6)=val
        tau(7)=val
        tau(8)=val
C
        w(1)=1.d0
        w(2)=1.d0
        w(3)=1.d0
        w(4)=1.d0
        w(5)=1.d0
        w(6)=1.d0
        w(7)=1.d0
        w(8)=1.d0
      end subroutine getXi

      subroutine getJacobian_NL(x, y, z, xi, eta, tau, Jac, B_el, 
     1  detJ, Ni)
!CCCCFOR NON-LOCAL SOLUTION RELATED ELEMENTS FOR NON-LOCAL PART
        implicit none
C OUTPUT        
        real*8, intent(out) :: Jac(3,3), B_el(3,8), detJ, Ni(8)
C INPUT        
        real*8, intent(in) :: x(8), y(8), z(8), xi, eta, tau
C LOCAL        
        integer i, j, k
        real*8 J_inv(3,3), N_X(3,8), Pmat(3,8), XYZ(8,3)
C
        do i=1,3
          do j=1,3
            Jac(i,j) = 0.d0
            J_inv(i,j) = 0.d0
          enddo!j
        enddo!i
C
        do i=1,3
          do j=1,8
            B_el(i,j)=0.d0
          enddo!j
        enddo!i
C
        do i=1,3
          do j=1,8
            N_X(i,j) = 0.d0
            Pmat(i,j) = 0.d0
          enddo!j
        enddo!i
C
        do i=1,8
          do j=1,3
            XYZ(i,j) = 0.d0
          enddo!j
        enddo!j
C        
        do i=1,8
          XYZ(i,1) = x(i)
          XYZ(i,2) = y(i)
          XYZ(i,3) = z(i)
        enddo
C        
        do i=1,8
          Ni(i)=0.d0
        enddo
C        
        Ni(1)=(1.d0/8.d0)*(1.d0-xi)*(1.d0-eta)*(1.d0-tau)
        Ni(2)=(1.d0/8.d0)*(1.d0+xi)*(1.d0-eta)*(1.d0-tau)
        Ni(3)=(1.d0/8.d0)*(1.d0+xi)*(1.d0+eta)*(1.d0-tau)
        Ni(4)=(1.d0/8.d0)*(1.d0-xi)*(1.d0+eta)*(1.d0-tau)
        Ni(5)=(1.d0/8.d0)*(1.d0-xi)*(1.d0-eta)*(1.d0+tau)
        Ni(6)=(1.d0/8.d0)*(1.d0+xi)*(1.d0-eta)*(1.d0+tau)
        Ni(7)=(1.d0/8.d0)*(1.d0+xi)*(1.d0+eta)*(1.d0+tau)
        Ni(8)=(1.d0/8.d0)*(1.d0-xi)*(1.d0+eta)*(1.d0+tau)
c     Pi : Ni,xi Ni,eta, Ni,tau, obtained from MATLAB
        Pmat(1,1) = -((eta - 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,2) = ((eta - 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,3) = -((eta + 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,4) = ((eta + 1.d0)*(tau - 1.d0))/8.d0
        Pmat(1,5) = ((eta - 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,6) = -((eta - 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,7) = ((eta + 1.d0)*(tau + 1.d0))/8.d0
        Pmat(1,8) = -((eta + 1.d0)*(tau + 1.d0))/8.d0
C
        Pmat(2,1) = -(xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,2) = (xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,3) = -(xi/8.d0 + 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,4) = (xi/8.d0 - 1.d0/8.d0)*(tau - 1.d0)
        Pmat(2,5) = (xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,6) = -(xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,7) = (xi/8.d0 + 1.d0/8.d0)*(tau + 1.d0)
        Pmat(2,8) = -(xi/8.d0 - 1.d0/8.d0)*(tau + 1.d0)
C
        Pmat(3,1) = -(xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,2) = (xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,3) = -(xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,4) = (xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,5) = (xi/8.d0 - 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,6) = -(xi/8.d0 + 1.d0/8.d0)*(eta - 1.d0)
        Pmat(3,7) = (xi/8.d0 + 1.d0/8.d0)*(eta + 1.d0)
        Pmat(3,8) = -(xi/8.d0 - 1.d0/8.d0)*(eta + 1.d0)
C Jacobian
        do i=1,3
          do k=1,3
            do j=1,8
              Jac(i,k) = Jac(i,k) + Pmat(i,j)*XYZ(j,k)
            enddo
          enddo
        enddo
C
        detJ = Jac(1,1)*( Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2) ) +
     1    Jac(1,2)*( Jac(2,3)*Jac(3,1) - Jac(2,1)*Jac(3,3) ) +
     2    Jac(1,3)*( Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1) )
C     
        ! if(ABS(detJ).LT.1e-12)then
        !   print*, '*** Error ***'
        !   print*, 'Jacobian is almost zero...'
        !   print*, '*** Error ***'
        !   CALL XIT()
        ! endif
C        
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
          do k=1,8
            do j=1,3
              N_X(i,k) = N_X(i,k) + J_inv(i,j)*Pmat(j,k)
            enddo
          enddo
        enddo
C
        B_el(1,1) = N_X(1,1)
        B_el(1,2) = N_X(1,2)
        B_el(1,3) = N_X(1,3)
        B_el(1,4) = N_X(1,4)
        B_el(1,5) = N_X(1,5)
        B_el(1,6) = N_X(1,6)
        B_el(1,7) = N_X(1,7)
        B_el(1,8) = N_X(1,8)
C
        B_el(2,1) = N_X(2,1)
        B_el(2,2) = N_X(2,2)
        B_el(2,3) = N_X(2,3)
        B_el(2,4) = N_X(2,4)
        B_el(2,5) = N_X(2,5)
        B_el(2,6) = N_X(2,6)
        B_el(2,7) = N_X(2,7)
        B_el(2,8) = N_X(2,8)
C
        B_el(3,1) = N_X(3,1)
        B_el(3,2) = N_X(3,2)
        B_el(3,3) = N_X(3,3)
        B_el(3,4) = N_X(3,4)
        B_el(3,5) = N_X(3,5)
        B_el(3,6) = N_X(3,6)
        B_el(3,7) = N_X(3,7)
        B_el(3,8) = N_X(3,8)
      end subroutine getJacobian_NL

      subroutine getXi_NL(xi, eta, tau, w)
!CCCC GAUSS POINTS FOR NON-LOCAL SOLUTION RELATED ELEMENTS
        implicit none
        real*8, intent(out)  :: xi(8), eta(8), tau(8), w(8)
C
        real*8 val
C        
        val = dsqrt(1.d0/3.d0)
        xi(1)=-val
        xi(2)=val
        xi(3)=-val
        xi(4)=val
        xi(5)=-val
        xi(6)=val
        xi(7)=-val
        xi(8)=val
C
        eta(1)=-val
        eta(2)=-val
        eta(3)=val
        eta(4)=val
        eta(5)=-val
        eta(6)=-val
        eta(7)=val
        eta(8)=val
C
        tau(1)=-val
        tau(2)=-val
        tau(3)=-val
        tau(4)=-val
        tau(5)=val
        tau(6)=val
        tau(7)=val
        tau(8)=val
C
        w(1)=1.d0
        w(2)=1.d0
        w(3)=1.d0
        w(4)=1.d0
        w(5)=1.d0
        w(6)=1.d0
        w(7)=1.d0
        w(8)=1.d0
      end subroutine getXi_NL