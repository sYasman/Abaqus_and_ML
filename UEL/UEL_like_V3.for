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
C
      integer i, j, k, u_first
      real*8 x(8), y(8), z(8), pl_props(10),
     1  fint(40), K_el(40,40), hist_n(112), hist_new(112), lc
C     
      lc = PROPS(3)
C MODIFIED HERE.
	  u_first = 1
C **  
C ** Must be the first user element's number
C ** *Element, type=U1 numara 
      ! if(JELEM.EQ.1)then
      if(JELEM.EQ.u_first)then  
        noel_k=1
      endif
C ** NOTICE:
C IF UEL is in the second assembly, then change 1 to first number of USET      
C
      do i=1,NPROPS ! load properties
        pl_props(i)=PROPS(i)
      enddo
C
      do i=1, NSVARS ! load history_n variables
        hist_n(i)=SVARS(i)
        hist_new(i)=0.d0        
      enddo
C
      do i=1,NDOFEL ! Initiation of RHS and AMATRX
        RHS(i,1)=0.d0
        fint(i)=0.d0
        do j=1,NDOFEL
          AMATRX(i,j)=0.d0
          K_el(i,j)=0.d0
        enddo!j
      enddo!i
C
      do i=1,NNODE ! load coordinates
        x(i) = COORDS(1,i)
        y(i) = COORDS(2,i)
        z(i) = COORDS(3,i)
      enddo  
C      
      if(LFLAGS(2).EQ.1)then
        print*,'This is not LARGE strain analysis.'
        print*,'Analysis is termianted'
        CALL XIT()
      endif
C      
      if(LFLAGS(3).EQ.1)then
        if(KINC.LT.1.AND.JELEM.EQ.u_first)then
	      print*, 'Analysis has been started with the following parameters'
!		  print*, 'pl_props:', pl_props
! E, nu, lc, aa, bb, sigma_t_inf(N_U), nu_p, e1p(N_U), e2p(N_U), e3p(N_U)
		  print*, 'E: ', pl_props(1)	
		  print*, 'nu: ', pl_props(2) 	
		  print*, 'lc: ', pl_props(3) 	
		  print*, 'nu_p: ', pl_props(7) 			  
		endif	  
	  
        CALL get_tang_fint(x, y, z, hist_n, hist_new, pl_props, 
     1  U, K_el, fint, lc, KINC, DTIME, PNEWDT)     
C
        do i=1,NDOFEL
          RHS(i,1)=fint(i)
          do j=1,NDOFEL
            AMATRX(i,j)=K_el(i,j)
          enddo
        enddo
C     
        do i=1,NSVARS
          SVARS(i)=hist_new(i)
        enddo        
      endif
C
      if(LFLAGS(3).EQ.2)then
        print*, 'flag 2'
        CALL XIT()        
      endif
C
      if(LFLAGS(3).EQ.3)then
        print*, 'flag 3'
        CALL XIT()        
c     define current Damping Matrix
      endif
C
      if(LFLAGS(3).EQ.4)then ! ilk buna giriyor
        ! print*, 'flag 4'
        ! CALL XIT()        
c     define current Mass Matrix, firstly this step is checked by solver
      endif
C
      if(LFLAGS(3).EQ.5)then
        print*, 'flag 5'
        CALL XIT()        
C     define current Residual Vector
      endif
C
      if(LFLAGS(3).EQ.6)then
        print*, 'flag 6'
        CALL XIT()
C     define current Mass Matrix and Residual Vector
      endif
C
      RETURN
      END    

   
