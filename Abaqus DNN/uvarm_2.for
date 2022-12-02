      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'KARC.blc'      
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      
        integer i, j, numUEL
        real*8 stressFull(3,3),dev_stress(3,3),eye(3,3),sp, 
     1   stress_vm

        real*8 strainFull(3,3), strain_tr, strain_dev(3,3), str_eq,
     1 num_nor_el, INCR       
C The dimensions of the variables FLGRAY, ARRAY and JARRAY
C must be set equal to or greater than 15.      
      
c     numUEL : number of UEL elements 
c user elemanların sonuncusu
		INCR = 7232 ! C3D8R_son - C3D8R_ilk 
        numUEL = 27336  ! *Element, type=C3D8T  ^^ numara dummy den önce son eleman
		
		numUEL = 87514 ! TEST_DEPTH_EFFECT
		numUEL = 5750 ! TEST_MODEL
		numUEL = 149948!NO_PBC_FINE
		! numUEL = 2875!PBC_DEPTH_ONE_EL
		numUEL = 74974 ! SYMM_MODEL
		numUEL = 78186 ! SYMMP_COMP
		numUEL = 11238
		numUEL = 9896 ! SMALL_RVE
		numUEL = 1000
		! numUEL = 16160
		! numUEL = 8000 ! SIMPLE CUBE FINE
C** ---------------- ONLY REQUIRED UVAR PROPS ----------        
        do i=1,3
          do j=1,3
            stressFull(i,j)=0.d0
            dev_stress(i,j)=0.d0
            eye(i,j)=0.d0
          enddo
        enddo   

        eye(1,1)=1.d0
        eye(2,2)=1.d0
        eye(3,3)=1.d0

        stressFull(1,1)=STRESS_k(NOEL-numUEL,NPT,1)
        stressFull(2,2)=STRESS_k(NOEL-numUEL,NPT,2)
        stressFull(3,3)=STRESS_k(NOEL-numUEL,NPT,3)
        stressFull(1,2)=STRESS_k(NOEL-numUEL,NPT,4)
        stressFull(1,3)=STRESS_k(NOEL-numUEL,NPT,5)
        stressFull(2,3)=STRESS_k(NOEL-numUEL,NPT,6)
        stressFull(2,1)=stressFull(1,2)
        stressFull(3,1)=stressFull(1,3)
        stressFull(3,2)=stressFull(2,3)      

        sp = (1.d0/3.d0)*(stressFull(1,1)+stressFull(2,2)+
     1  stressFull(3,3))

        do i=1,3
          do j=1,3
            dev_stress(i,j)=stressFull(i,j)-sp*eye(i,j)
          enddo
        enddo
        stress_vm = 0.d0

        do i=1,3
          do j=1,3
            stress_vm = stress_vm + dev_stress(i,j)*dev_stress(i,j)
          enddo
        enddo

        ! stress_vm = SQRT(2.d0/3.d0 * stress_vm)
        stress_vm = SQRT(3.d0/2.d0 * stress_vm)

C equivalent strain
        strainFull(1,1)=STRAIN_k(NOEL-numUEL,NPT,1)
        strainFull(2,2)=STRAIN_k(NOEL-numUEL,NPT,2)
        strainFull(3,3)=STRAIN_k(NOEL-numUEL,NPT,3)
        strainFull(1,2)=STRAIN_k(NOEL-numUEL,NPT,4)
        strainFull(1,3)=STRAIN_k(NOEL-numUEL,NPT,5)
        strainFull(2,3)=STRAIN_k(NOEL-numUEL,NPT,6)
        strainFull(2,1)=0.50d0*strainFull(1,2)
        strainFull(3,1)=0.50d0*strainFull(1,3)
        strainFull(3,2)=0.50d0*strainFull(2,3)   

        strain_tr = strainFull(1,1) + strainFull(2,2) + strainFull(3,3)
        
        do i=1,3
          do j=1,3
            strain_dev(i,j) = strainFull(i,j) -
     @         (1.0d0/3.0d0)*strain_tr*eye(i,j)
          enddo
        enddo

        str_eq = 0.d0

        do i=1,3
          do j=1,3
            str_eq = str_eq + 
     @        (2.d0/3.d0)*strain_dev(i,j)*strain_dev(i,j)
          enddo
        enddo

        str_eq = SQRT(str_eq)
C**
C** Store data in UVARs
C**
        UVAR(1)=stress_vm ! von misses stress   
        UVAR(2)=str_eq ! Eq. strain
        UVAR(3)=EPS_EQ_K(NOEL-numUEL,NPT) ! equivalent plastic strain
        UVAR(4)=EPS_EQ_N_K(NOEL-numUEL,NPT) ! equivalent plastic strain
        UVAR(5)=DM_K(NOEL-numUEL,NPT) ! positive damage        
        UVAR(6)=DM_N_K(NOEL-numUEL,NPT) ! negative damage        
        UVAR(7)=STRESS_k(NOEL-numUEL,NPT,1)   ! S11
        UVAR(8)=STRESS_k(NOEL-numUEL,NPT,2)   ! S22           
        UVAR(9)=STRESS_k(NOEL-numUEL,NPT,3)   ! S33           
        UVAR(10)=STRESS_k(NOEL-numUEL,NPT,4)  ! S12           
        UVAR(11)=STRESS_k(NOEL-numUEL,NPT,5)  ! S13           
        UVAR(12)=STRESS_k(NOEL-numUEL,NPT,6)  ! S23    
		! if(KINC.EQ.2.and.NPT.EQ.3)then
		! print*, DS_DE(NOEL-numUEL,NPT, 1, 1)
		! CALL XIT()
		! endif
        UVAR(15) = DS_DE(NOEL-numUEL,NPT, 1, 1)
        UVAR(16) = DS_DE(NOEL-numUEL,NPT, 1, 2)
        UVAR(17) = DS_DE(NOEL-numUEL,NPT, 1, 3)
        UVAR(18) = DS_DE(NOEL-numUEL,NPT, 1, 4)
        UVAR(19) = DS_DE(NOEL-numUEL,NPT, 1, 5)
        UVAR(20) = DS_DE(NOEL-numUEL,NPT, 1, 6)
		
        UVAR(21) = DS_DE(NOEL-numUEL,NPT, 2, 1)
        UVAR(22) = DS_DE(NOEL-numUEL,NPT, 2, 2)
        UVAR(23) = DS_DE(NOEL-numUEL,NPT, 2, 3)
        UVAR(24) = DS_DE(NOEL-numUEL,NPT, 2, 4)
        UVAR(25) = DS_DE(NOEL-numUEL,NPT, 2, 5)
        UVAR(26) = DS_DE(NOEL-numUEL,NPT, 2, 6)		
		
        UVAR(27) = DS_DE(NOEL-numUEL,NPT, 3, 1)
        UVAR(28) = DS_DE(NOEL-numUEL,NPT, 3, 2)
        UVAR(29) = DS_DE(NOEL-numUEL,NPT, 3, 3)
        UVAR(30) = DS_DE(NOEL-numUEL,NPT, 3, 4)
        UVAR(31) = DS_DE(NOEL-numUEL,NPT, 3, 5)
        UVAR(32) = DS_DE(NOEL-numUEL,NPT, 3, 6)			
		
        UVAR(33) = DS_DE(NOEL-numUEL,NPT, 4, 1)
        UVAR(34) = DS_DE(NOEL-numUEL,NPT, 4, 2)
        UVAR(35) = DS_DE(NOEL-numUEL,NPT, 4, 3)
        UVAR(36) = DS_DE(NOEL-numUEL,NPT, 4, 4)
        UVAR(37) = DS_DE(NOEL-numUEL,NPT, 4, 5)
        UVAR(38) = DS_DE(NOEL-numUEL,NPT, 4, 6)		

        UVAR(39) = DS_DE(NOEL-numUEL,NPT, 5, 1)
        UVAR(40) = DS_DE(NOEL-numUEL,NPT, 5, 2)
        UVAR(41) = DS_DE(NOEL-numUEL,NPT, 5, 3)
        UVAR(42) = DS_DE(NOEL-numUEL,NPT, 5, 4)
        UVAR(43) = DS_DE(NOEL-numUEL,NPT, 5, 5)
        UVAR(44) = DS_DE(NOEL-numUEL,NPT, 5, 6)		

        UVAR(45) = DS_DE(NOEL-numUEL,NPT, 6, 1)
        UVAR(46) = DS_DE(NOEL-numUEL,NPT, 6, 2)
        UVAR(47) = DS_DE(NOEL-numUEL,NPT, 6, 3)
        UVAR(48) = DS_DE(NOEL-numUEL,NPT, 6, 4)
        UVAR(49) = DS_DE(NOEL-numUEL,NPT, 6, 5)
        UVAR(50) = DS_DE(NOEL-numUEL,NPT, 6, 6)		
		
        UVAR(51) = dEpsPl_dEpsPl_bar_pos(NOEL-numUEL,NPT)		
        UVAR(52) = dEpsPl_dEpsPl_bar_neg(NOEL-numUEL,NPT)		
		
        UVAR(53)=STRAIN_k(NOEL-numUEL,NPT,1)   ! E11
        UVAR(54)=STRAIN_k(NOEL-numUEL,NPT,2)   ! E22           
        UVAR(55)=STRAIN_k(NOEL-numUEL,NPT,3)   ! E33           
        UVAR(56)=STRAIN_k(NOEL-numUEL,NPT,4)  ! E12           
        UVAR(57)=STRAIN_k(NOEL-numUEL,NPT,5)  ! E13           
        UVAR(58)=STRAIN_k(NOEL-numUEL,NPT,6)  ! E23   	

        UVAR(59)=EPS_PL_K(NOEL-numUEL,NPT,1)   ! EPS_PL_K_11
        UVAR(60)=EPS_PL_K(NOEL-numUEL,NPT,2)   ! EPS_PL_K_22           
        UVAR(61)=EPS_PL_K(NOEL-numUEL,NPT,3)   ! EPS_PL_K_33           
        UVAR(62)=EPS_PL_K(NOEL-numUEL,NPT,4)  ! EPS_PL_K_12           
        UVAR(63)=EPS_PL_K(NOEL-numUEL,NPT,5)  ! EPS_PL_K_13           
        UVAR(64)=EPS_PL_K(NOEL-numUEL,NPT,6)  ! EPS_PL_K_23   
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
      RETURN
      END