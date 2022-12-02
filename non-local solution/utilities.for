      subroutine getIandPterms(I_1, II, PP)
        implicit none
C OUTPUT        
        real*8, intent(out) :: I_1(3,3), II(3,3,3,3), PP(3,3,3,3)
C LOCAL        
        integer i, j, k, l, m, n
C        
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
C
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
              II(i,j,k,l)=0.5d0*(I_1(i,k)*I_1(j,l)+I_1(i,l)*I_1(j,k))
              PP(i,j,k,l)=II(i,j,k,l)-(1.d0/3.d0)*I_1(i,j)*I_1(k,l)
              enddo
            enddo
          enddo
        enddo
      endsubroutine getIandPterms

      subroutine normOf(A, nrA)
!CCCC SUBROUTINE RETURNS NORM OF 3by3 MATRIX
        implicit none
C INPUT        
        real*8, intent(in) :: A(3,3)
C OUTPUT        
        real*8, intent(out) :: nrA
C LOCAL        
        integer i, j
C
        nrA = 0.d0
        do i=1,3
          do j=1,3
            nrA = nrA + A(i,j)*A(i,j)
          enddo
        enddo
        nrA = SQRT(nrA)
      end subroutine normOf

      subroutine getInvariants(A, A_dev, I1, J2)
!CCCC SUBROUTINE RETURNS INVARIANTS OF 3by3 STRESS MATRIX
        implicit none
C INPUT      
        real*8, intent(in) :: A(3,3), A_dev(3,3)
C OUTPUT      
        real*8, intent(out) :: I1, J2
C LOCAL      
        integer i, j
C      
        I1 = 0.d0
        J2 = 0.d0
        do i=1,3
          I1 = I1 + A(i,i)
          do j=1,3
            J2 = J2+0.5d0*(A_dev(i,j)*A_dev(i,j))
          enddo
        enddo
      endsubroutine getInvariants

      subroutine getIandPtermsV(I_1v, I_dy_Iv, IIv, PPv)
C Definig sencond and fourth order idenitiy matrices in voigt notation
        implicit none
C OUTPUT        
        real*8, intent(out) :: I_1v(6), IIv(6,6),I_dy_Iv(6,6),PPv(6,6)
C LOCAL        
        integer i, j
C        
        do i=1,6
          I_1v(i)=0.0
          do j=1,6
            IIv(i,j)=0.d0
            I_dy_Iv(i,j)=0.d0
          enddo
        enddo
C
        do i=1,3
          I_1v(i)=1.d0
          IIv(i,i)=1.d0
          IIv(i+3,i+3)=0.5d0
          do j=1,3
            I_dy_Iv(i,j)=1.d0
          enddo
        enddo
C
        do i=1,6
          do j=1,6
            PPv(i,j)=IIv(i,j)-(1.d0/3.d0)*I_1v(i)*I_1v(j)
          enddo
        enddo
C
! PP_MATRIX_BURADA.PDF --> p15
        do i=4,6
          PPv(i,i) = 2.d0*PPv(i,i)
        enddo ! diagonal terimler 0.50 idi (4,4 5,5 6,6) 1.00 yapildi
      endsubroutine getIandPtermsV      

      subroutine getDet(A, det)
        implicit none
C
        real*8, intent(in)::A(3,3)      
        real*8, intent(out)::det     
        real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33
        a11=A(1,1)
        a12=A(1,2)
        a13=A(1,3)
        a21=A(2,1)
        a22=A(2,2)
        a23=A(2,3)
        a31=A(3,1)
        a32=A(3,2)
        a33=A(3,3)
C        
        det = a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + 
     @ a13*a21*a32 - a13*a22*a31
      end subroutine getDet    

      subroutine getInv(A, invA)
        implicit none
C
        real*8, intent(in)::A(3,3)      
        real*8, intent(out)::invA(3,3)     
        real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33, detA
        integer i, j
        a11=A(1,1)
        a12=A(1,2)
        a13=A(1,3)
        a21=A(2,1)
        a22=A(2,2)
        a23=A(2,3)
        a31=A(3,1)
        a32=A(3,2)
        a33=A(3,3)
C        
        CALL getDet(A, detA)
        invA(1,1)=a22*a33 - a23*a32
        invA(1,2)=-a12*a33 + a13*a32
        invA(1,3)=a12*a23 - a13*a22
        invA(2,1)=-a21*a33 + a23*a31
        invA(2,2)=a11*a33 - a13*a31
        invA(2,3)=-a11*a23 + a13*a21
        invA(3,1)=a21*a32 - a22*a31
        invA(3,2)=-a11*a32 + a12*a31
        invA(3,3)=a11*a22 - a12*a21

        do i=1,3
          do j=1,3
            invA(i,j) = invA(i,j) / detA
          enddo
        enddo
      end subroutine getInv           