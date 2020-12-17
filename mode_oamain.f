      MODULE MODE_OAMAIN

!Purpose: main subroutine for calling SOA equilibrium routines
!Similar programs recieved from B. Pun and R. Griffin, 2005
!Rewritten by Alf Grini, alf.grini@cnrm.meteo.fr, 2005

      USE modd_glo
      USE modd_glodef

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE MPMPO(          
     &   TEMPK,                    !I [K] Temperature
     &   RH,                       !I [-] Relative humidity
     &   CB,                       !O [ug/m3] total (aerosol+gas) conc
     &   CPT,                      !I [ug/m3] primary aerosol concentration
     &   GAS,                      !O [ug/m3] gas phase concentrations
     &   AERO_ORG,                 !O [ug/m3] organic phase concentrations
     &   AERO_AQ,                  !I [ug/m3] aquous phase concentration of acid and ions
     &   PARTORG,                  !O [ug/m3] aerosol phase concentration (aq+org)
     &   LWC,                      !I [ug/m3] liquid water content available for partitioning
     &   acHP,                     !I [mol/kg_{water}] H+ concentrations
     &   DELTALWC,                 !O [ug/m3] change in LWC
     &   ORGANION,                 !O [mole/m3] anion charge
     &   Gamma_aq_Raoult, Gamma_aq_Henry !ddw -- output the gammas
     &  )

       USE mode_soatinit           !Module which calculates parameters only dependent on temperature
       USE mode_firstguess         !Module which calculates first guesses of components
       USE mode_soaeql             !Module which does the equilibrium between aquous phase, organic phase and gas phase
       USE mode_zsrpun             !Module to get water associated with aquous organics
      USE modd_aunifacparam       !Module with unifac coefficients for aquous phase
      USE modd_bunifacparam       !Module with unifac coefficients for organic phase

       IMPLICIT NONE

!INPUTS/OUTPUTS
      REAL, DIMENSION(:), INTENT(IN)       :: TEMPK      !I [K] Temperature
      REAL, DIMENSION(:), INTENT(IN)       :: RH         !I [0-1] Relative humidity
      REAL, DIMENSION(:,:), INTENT(IN)     :: CB         !I [ug/m3] total (g+p) aerosol organic species
      REAL, DIMENSION(:,:), INTENT(IN)     :: CPT        !I [ug/m3]
      REAL, DIMENSION(:,:), INTENT(OUT)    :: GAS        !I [ug/m3] gas of organic species
      REAL, DIMENSION(:,:), INTENT(OUT)    :: AERO_AQ    !I [ug/m3] aqueous species
      REAL, DIMENSION(:,:), INTENT(OUT)    :: AERO_ORG   !I [ug/m3] organic phase species
      REAL, DIMENSION(:,:), INTENT(OUT)    :: PARTORG    !I [ug/m3] aerosol phase organic species
      REAL, DIMENSION(:), INTENT(IN)       :: LWC        !I [ug/m3] liquid water content of aerosols
      REAL, DIMENSION(:), INTENT(IN)       :: acHP       !I [mol_{H+}/kg_{water}] concentration of H+ ions
      REAL, DIMENSION(:), INTENT(OUT)      :: DELTALWC   !I [ug/m3] LWC assiciated with organics
      REAL, DIMENSION(:), INTENT(OUT)      :: ORGANION   !I [mole/m3] organic anion concentratios

!Local variables
      REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2)) :: VP  ![torr] temp. dependent vapor pressures of organic phase
      REAL, DIMENSION(SIZE(AERO_AQ,1),SIZE(AERO_AQ,2))   :: K   ![m3/ug and mole/kg] temp. dependent Henry's law and diss. const
      REAL, DIMENSION(:,:,:), ALLOCATABLE                :: SI_ORG  ![-] temp. dependent unifac coefficient
      REAL, DIMENSION(:,:,:), ALLOCATABLE                :: SI_AQ   ![-] temp. dependent unifac coefficient

!Small variables
      INTEGER                              :: I         ![idx] counter for main type A components
      INTEGER                              :: J         ![idx] counter for sub type A components
      INTEGER                              :: COMP_IDX  ![idx] index for right type A comp (sub or main)
      INTEGER                              :: COMP_IDX2 ![idx] help counter for idx of main components
      REAL                                 :: XX        ![idx] account for difference in molecular weight
! ddw -- output gamma
      REAL, DIMENSION(:,:), ALLOCATABLE    :: GAMMA_AQ_HENRY   ![-] activity coefficients for aquous phase (Henry's law std. state)
      REAL, DIMENSION(:,:), ALLOCATABLE    :: GAMMA_AQ_RAOULT  ![-] activity coefficients for aquous phase (Raoult's law std. state)

      ALLOCATE(SI_ORG(SIZE(AERO_ORG,1),NFUNC_ORG,NFUNC_ORG))
      ALLOCATE(SI_AQ(SIZE(AERO_ORG,1),NFUNC_AQ,NFUNC_AQ))

! ddw- initilize gamma_aq
      GAMMA_AQ_HENRY(:,:) = 0.0
      GAMMA_AQ_RAOULT(:,:) = 0.0

!Initialize the negative charge associated with type AQ organics (moles/m3)
      ORGANION(:)=0.d0

!Initialize the liquid water content associated with AQ organics (ug/m3)
      deltaLWC(:)=0.d0

!Initialize AERO_AQ
      AERO_AQ(:,:)=0.d0

!Initialize AERO_ORG
      AERO_ORG(:,:)=0.d0

!Initialize PARTORG
      PARTORG(:,:)=0.d0

!Initialize GAS
      GAS(:,:)=0.d0

!Get termperature corrected vapor pressures 
       CALL VP_GET(           
     &     TEMPK,             !I [K] Temperature
     &     VP                 !I [torr] saturation vapor pressures
     &    )

!Get temperature corrected henry's law constants 
        CALL AQCONST_GET(            
     &    TEMPK,                  !I [K] temperature
     &    K                       !O [m3/ug and mole/kg] Henry's law and dissociation constants
     &    )

!Get temperature dependent coefficients needed in UNIFAC parameterization
       CALL SI_GET(                 
     &     TEMPK,                     !I [K] temperature
     &     A_ORG,                     !I [units?] term in UNIFAC parameterization
     &     SI_ORG,                    !O [units?] term in UNIFAC parameterization
     &     NFUNC_ORG                 !I [nbr] number of functional group in mixture
     &    )
    
!Do the same thing, but now get the coefficients for aquous phase
       CALL SI_GET(                 
     &    TEMPK,                     !I [K] temperature
     &    A_AQ,                     !I [units?] term in UNIFAC parameterization
     &    SI_AQ,                    !O [units?] term in UNIFAC parameterization
     &    NFUNC_AQ                 !I [nbr] number of functional group in mixture
     &    )

!Get first guesses for aerosols
       CALL GUESS_AERO(   
     &    CB,              !I [ug/m3] total (gas + aerosol + ions) of component
     &    AERO_ORG,        !O [ug/m3] liquid aerosol concentrations
     &    AERO_AQ,         !O [ug/m3] solid aerosol concentration
     &    GAS,             !O [ug/m3] gaseous concentration
     &    LWC,             !I [ug/m3] liquid water content already available for partitioning
     &    acHP,            !I [umol/kg_{water}] proton concentration
     &    VP,              !I [torr] saturation vapor pressure of organic precursors
     &    K                !I [m3/ug and mole/kg] Henry's law and dissociation constants
     &    )

!Use the guesses to iterate for the right solution
      CALL SOAEQL(             
     &    acHP,                 !I [mol/kg_{water}] proton concentration
     &    LWC,                 !I [ug/m3] liquid water content
     &    CPT,                 !I [ug/m3] primary organic concentration
     &    TEMPK,               !I [K] temperature
     &    GAS,                 !I/O [ug/m3] gas phase concentrations of SOA
     &    AERO_AQ,             !I/O [ug/m3] aquous phase concentrations of SOA and ions
     &    AERO_ORG,            !I/O [ug/m3] organic phase concentrations of SOA
     &    CB,                  !I [ug/m3] total SOA concentration (GAS + AQ + ORG)
     &    K,                  !I [m3/ug and mole/kg] Henry's law and dissocation constants at right T
     &    SI_ORG,             !I [?] T-dependent coefficient in UNIFAC parameterization
     &    SI_AQ,              !I [?] T-dependent coefficient in UNIFAC parameterization
     &    VP,                 !I [torr] saturation vapor pressure
     &    GAMMA_AQ_RAOULT, GAMMA_AQ_HENRY   ! ddw --- output the coefficient 
     &)

!Get the wataer associated with the extra organics using ZSR
      CALL ZSRPUN(          
     &    AERO_AQ,              !I [ug/m3] guess for aerosol concentrations
     &    RH,                   !I [0-1] relative humidity
     &    DELTALWC              !O [ug/m3] liquid water content
     &    )
   

      COMP_IDX=1
      DO I=1,NBSP

       !Start with summing the organic part
       PARTORG(:,I)=AERO_ORG(:,I)

       !Conserve COMP_IDX of main component
       COMP_IDX2 = COMP_IDX
       
       !Initiate partorg from this component
       PARTORG(:,I)=PARTORG(:,I) + AERO_AQ(:,COMP_IDX)
       COMP_IDX=COMP_IDX+1
       
       !Add the contribution of the ions
       DO J=2,NK(I)
          XX=MW_SOA(COMP_IDX2)    !MW of main component
     &          /MW_SOA(COMP_IDX)   !MW of this component
          
          !Sum PARTORG
          PARTORG(:,I) = PARTORG(:,I) + AERO_AQ(:,COMP_IDX)*XX

          !Get anion concentration in mole/m3
          organion(:)=organion(:)  
     &         +AERO_AQ(:,COMP_IDX)/MW_SOA(COMP_IDX)  !umole/m3
     &          *dble(J-1)                            !number of negative charges
     &          *1.d-6                                  !==> mole/m3

          !Prepare for next component
          COMP_IDX=COMP_IDX+1
       ENDDO

      ENDDO

      DEALLOCATE(SI_ORG)
      DEALLOCATE(SI_AQ)

      END SUBROUTINE MPMPO
    
      END MODULE MODE_OAMAIN
