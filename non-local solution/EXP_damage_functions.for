C**************************************************************
C*********** DAMAGE FUNCTIONS FOR EQ. PLASTIC STRAIN **********  
C**************************************************************     
C*** The maximum value of damae is limited to bo 0.990d0
C*** dd_dr after d_m reaches 0.990d0 becomes 0.d0
      subroutine getRcEta(kappa_0, d_max, aa, bb)   
C aa := alpha             
C bb := beta    
C d_max: maximum damage value. To prevent instabilities   
C kappa_0: damage initiation, non-local equivalent plastic strain    
C regarding the sign of the largest principal_equivalent_strain, 
C it is decided if the loading is tension or compression. Then inspecting load_flag
C kappa_0 is chosen.
C
        implicit none
! OUTPUT
        real*8, intent(out) :: kappa_0, d_max, aa, bb
C        
C damage initiation
		kappa_0 = 2.00e-2
C
C damage parameters
        aa = 0.95d0
	    bb = 7.d0 
C   
C limiting the max. damage     
        d_max=0.99d0                                 
      endsubroutine getRcEta 
C
C
C      
      subroutine getRcEta_comp(kappa_0, d_max, aa, bb)   
C aa := alpha             
C bb := beta    
C d_max: maximum damage value. To prevent instabilities   
C kappa_0: damage initiation, non-local equivalent plastic strain    
C regarding the sign of the largest principal_equivalent_strain, 
C it is decided if the loading is tension or compression. Then inspecting load_flag
C kappa_0 is chosen.
C
        implicit none
! OUTPUT
        real*8, intent(out) :: kappa_0, d_max, aa, bb
C        
C damage initiation
        kappa_0 = 15.00e-2		
	
C damage parameters
        aa = 0.95d0
        bb = 7.d0				
C
        d_max=0.99d0                                 
      endsubroutine getRcEta_comp             
      
      subroutine exp_damage(r, d_m)   
        implicit none
C INPUT
        real*8, intent(in) :: r ! (r = kappa_str)
C OUTPUT
        real*8, intent(out) :: d_m
C LOCAL
        real*8 kappa_0, rc, eta, t1, d_max, k_i, k_c, beta
        real*8 aa, bb
C
        CALL getRcEta(k_i, d_max, aa, bb) 
C
        if(r.LT.k_i)then
          d_m = 0
        elseif(r.GE.k_i)then
          d_m = 1.0d0 - EXP(-bb*(r-k_i))
        endif
C
        if(d_m.GT.d_max)then
          d_m=d_max
        endif
C
        if(d_m.LT.0.d0)then
          print*, 'Negative damage error'
          CALL XIT()
        endif
C
        if(d_m.GT.1.d0)then
          print*, 'Damage is larger than 1.0'
          CALL XIT()          
        endif
      end subroutine exp_damage    
      
      subroutine exp_damage_comp(r, d_m)   
        implicit none
C INPUT
        real*8, intent(in) :: r ! (r = kappa_str)
C OUTPUT
        real*8, intent(out) :: d_m
C LOCAL
        real*8 kappa_0, rc, eta, t1, d_max, k_i, k_c, beta
        real*8 aa, bb
C
        CALL getRcEta_comp(k_i, d_max, aa, bb) 
C
        if(r.LT.k_i)then
          d_m = 0
        elseif(r.GE.k_i)then
          d_m = 1.0d0 - EXP(-bb*(r-k_i))
        endif
C
        if(d_m.GT.d_max)then
          d_m=d_max
        endif	
C
        if(d_m.LT.0.d0)then
          print*, 'Negative damage error'
          CALL XIT()
        endif
C
        if(d_m.GT.1.d0)then
          print*, 'Damage is larger than 1.0'
          CALL XIT()          
        endif
      end subroutine exp_damage_comp         

      
      subroutine getDd_Dr_exp_damage(d_m, r, dd_dr)
        implicit none
! INPUT
        real*8, intent(in) :: r, d_m
C r: kappa_str, the most severe equivalent plastic strain, non-local        
! OUTPUT
        real*8, intent(out) ::  dd_dr
! LOCAL
        real*8 kappa_0, t1, t2, t3, d_max, k_i, k_c
        real*8 aa, bb        
C
        CALL getRcEta(k_i, d_max, aa, bb) 
C
        if(d_m.GT.k_i)then
          dd_dr = bb*exp(-bb*(r-k_i))   
        else
          dd_dr = 0.0d0
        endif
C
        if(d_m.GT.d_max)then
          dd_dr = 0.d0
        endif        
C                
      end subroutine getDd_Dr_exp_damage   

      subroutine getDd_Dr_exp_damage_comp(d_m, r, dd_dr)
        implicit none
! INPUT
        real*8, intent(in) :: r, d_m
C r: kappa_str, the most severe equivalent plastic strain, non-local        
! OUTPUT
        real*8, intent(out) ::  dd_dr
! LOCAL
        real*8 kappa_0, t1, t2, t3, d_max, k_i, k_c
        real*8 aa, bb
C
        CALL getRcEta_comp(k_i, d_max, aa, bb) 
C
        if(d_m.GT.k_i)then
          dd_dr = bb*exp(-bb*(r-k_i))   
        else
          dd_dr = 0.0d0
        endif
C
        if(d_m.GT.d_max)then
          dd_dr = 0.d0
        endif        
C                
      end subroutine getDd_Dr_exp_damage_comp         