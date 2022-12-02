C for assigning ip variables to dummy elements
      include 'uvarm_2.for' 

C Jacobian, gauss points, B_operator
      include 'element_functions.for' 

C get_fint and get_tang functions
	  include 'v3_COUPLED_LIG_tang_and_fint_V9.for' ! vectorized	  
	  
C Constructs element stiffness matrix the one 40-by-40
	  include 'v3_COUPLED_el_tang_v6.for' ! extra lines are removed  

C Utility functions for localizing gradient damage model      
      include 'localizing_gradient_functions.for'

C Damage functions
      ! include 'new_EXP_damage_functions.for' ! alpha is 1.0d0 version of TEST_damage_functions
      include 'EXP_damage_functions.for' ! alpha is 1.0d0 version of TEST_damage_functions	
  

C Material responsee subroutines
C Similar to UMAT subroutinnes        
	  include 'v9_COUPLED_MOD_v8_17_UMAT_like_v34.for' ! NR wo G_tens and G_comp
	  ! include 'v9_COUPLED_MOD_v8_19_UMAT_like_v35.for' ! bi-section 	  

C Helper functions      
C Identitiy, Projection, their voight notations
      include 'utilities.for' 
      
C main subroutine
C entry point of Abaqus solver  
      include 'UEL_like_V3.for' ! three-field model
      
      



      