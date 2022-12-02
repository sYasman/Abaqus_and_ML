      subroutine get_R_n(R, n)
C ************
C n: Speed of the g function to reach zero. Bigger faster
C R: The last value of g
C R=1 --> CONVENTIONAL IMPLICIT GRADIENT DAMAGE MODEL
C ************       
C
        ! n = 7.d0 ! Çok iyi sonuc verdi
        ! R = 0.005d0 ! ORG   
C        
        implicit none
! INPUT
! OUTPUT
        real*8, intent(out) ::  R, n
! LOCAL
C Valeus from paper Sakar et al.   
        ! n = 5.d0 ! Torsion test icin bu secilecek.
        n = 1.d0 ! ORG
        R = 0.005d0 ! ORG
		! R = 1.d0 ! CIGD olur..
      end subroutine get_R_n

      subroutine g_func(d_m, g_val)
        implicit none
! INPUT
        real*8, intent(in) :: d_m
! OUTPUT
        real*8, intent(out) ::  g_val
! LOCAL
        real*8 R, n, t1, t2

        CALL get_R_n(R, n)
        t1 = (1.d0-R)*EXP(-n*d_m)+R-EXP(-n)
        t2 = 1.d0 - EXP(-n)

        if(ABS(t2).LT.1e-6)then
          print*, 'Almost zero error in g_func'
          print*, 't2 is zero'          
          CALL XIT()
        endif

        g_val = t1/t2    
      end subroutine g_func       
      
      subroutine get_dg_dd(d_m, dg_dd)
C**
C** DERIVATIVE OF INTERACTION FUNCTION W.R.T. DAMAGE
C**        
        implicit none
! INPUT
        real*8, intent(in) :: d_m
! OUTPUT
        real*8, intent(out) ::  dg_dd
! LOCAL
        real*8 R, n, t1, t2

        CALL get_R_n(R, n)
        t1 = n*(R-1.d0)*EXP(-n*d_m)
        t2 = 1.d0-EXP(-n)

        if(ABS(t2).LT.1e-6)then
          print*, 'Almost zero error in get_dg_dd'
          CALL XIT()
        endif      
        dg_dd = t1/t2
      end subroutine get_dg_dd    