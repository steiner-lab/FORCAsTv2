      MODULE MODE_UNIFAC
  
!*******************************************************************
!Theory: Basics of UNIFAC group combination method to get activity 
!coefficients for mixtures. This is nicely outlined in 
!Marcolli and Peter, ACP, vol 5, pp 1501-1527, 2005
!http://overview.sref.org/1680-7324/acp/2005-5-1545
!In the following I follow notations of their eqn 1-9
  
!Approximately following original code recieved from Pierre Tulet
!who got it from Betty Pun who (it seems) got it from Pradeep Saxena
  
!Adopted to MESONH by Alf Grini, CNRM, 2005
!*******************************************************************
!MAIN CHANGES:
!1) Doing calculations on a vector instead of on a point
!2) A lot of calcualtions taken out and done in unifac_ini routines

      USE modd_bunifacparam
      USE modd_glodef

      IMPLICIT NONE

      PUBLIC
      PRIVATE :: LNGAMMA_RESIDUAL_GET, LNGAMMA_COMBINATORY_GET
  
      CONTAINS
  
      SUBROUTINE ACT_COEFF_GET(      
     &  NU,                         !I [nbr] number of functional groups in molecule I
     &  X,                          !I [frc] mole fraction of molecule I
     &  QG,                         !I [m2(?)] group surface area parameter
     &  GAMA,                       !O [-] activity coefficient
     &  THTAGP,                     !I [-] surface area ratio of groups (j) in pure component (i)
     &  Q,                          !I [m2] surface area of pure component
     &  R,                          !I [m3] total volume of pure component
     &  L,                          !I [?] unifac parameter pure component
     &  SI,                         !I [?] temperature dependent term
     &  NMOL,                       !I [nbr] total number of molecules
     &  NFUNC                       !I [nbr] total number of functional groups (e.g. CH2, NO3 ..)
     &  )

! Purpose: Get the activity coefficients for any mixture using the UNIFAC method
! Book (Eqn. 3.8 & 3.9): Fredenslund et al. 1977. Vapor-liqud Equilibriua Using UNIFAC.  
      IMPLICIT NONE
    
!INPUT
      INTEGER, INTENT(IN)                     :: NFUNC     ![nbr] number of functional groups
      INTEGER, INTENT(IN)                     :: NMOL      ![nbr] number of molecules in question
      INTEGER, DIMENSION(:,:)                 :: NU        ![nbr] number of func. groups (j) in molec (i)
      REAL, DIMENSION(:,:)                    :: X         ![-] molar fraction of components
      REAL, DIMENSION(:), INTENT(IN)          :: R         ![m3] volume of molecules
      REAL, DIMENSION(:), INTENT(IN)          :: Q         ![m2] surface of molecule
      REAL, DIMENSION(:), INTENT(IN)          :: QG        ![m2] surface of groups
      REAL, DIMENSION(:), INTENT(IN)          :: L         ![??] UNIFAC parameter
      REAL, DIMENSION(:,:), INTENT(IN)        :: THTAGP    ![-] sfc area frc of group j in comp i
      REAL, DIMENSION(:,:,:), INTENT(IN)      :: SI        ![?] temperature dependent term

!OUTPUT
      REAL, DIMENSION(:,:), INTENT(OUT)       :: GAMA      ![-] activity coefficient of comp i in mix 

!LOCAL
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2))    :: LNGAMMA_C ![-] log of combinatory part of gamma
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2))    :: LNGAMMA_R ![-] log of residual part of gamma
      INTEGER                                 :: I         ![idx] counter for molecules

!Get the combinatorial part of gamma
      CALL LNGAMMA_COMBINATORY_GET (       
     &    Q,                                 !I [m2] surface of one molecule
     &    R,                                 !I [m3] volume of one molecule
     &    L,                                 !L [?] unifac parameter
     &    X,                                 ![-] molar fraction of components
     &    LNGAMMA_C,                         ![-] log of combinatory part of act. coeff
     &    NMOL                               ![nbr] number of molecules
     &    )

!Get the residual part of gamma
      CALL LNGAMMA_RESIDUAL_GET(          
     &    QG,                              !I [m2] surface of groups
     &    NU,                              !I [nbr] number of functional groups in components
     &    X,                               !I [-] molar fraction of species
     &    THTAGP,                          !I [-] area fraction of groups in pure components
     &    SI,                              !I [?] temperature dependent term
     &    LNGAMMA_R,                       !O [-] log of residual part of gamma
     &    NMOL,                            !I [nbr] number of molecules 
     &    NFUNC                            !I [nbr] number of functional groups
     &    )
!Get the total activity coefificent
      DO I=1,NMOL
       GAMA(:,I)=EXP(LNGAMMA_C(:,I) + LNGAMMA_R(:,I))
      ENDDO

      END SUBROUTINE ACT_COEFF_GET

!***************************************************************************

      SUBROUTINE LNGAMMA_COMBINATORY_GET(   
     &   Q,                                 !I [m2] surface of one molecule
     &   R,                                 !I [m3] volume of one molecule
     &   L,                                 !L [?] unifac parameter
     &   X,                                 ![-] molar fraction of components
     &   LNGAMMA_C,                         ![-] log of combinatory part of act. coeff
     &   NMOL                               ![nbr] number of molecules
     &  )
    
!Purpose: Find the combinatory part of the activity coefficient
      IMPLICIT NONE
    
!INPUT/OUTPUT
      REAL, DIMENSION(:), INTENT(IN)       :: R         ![m3] volume of molecules
      REAL, DIMENSION(:), INTENT(IN)       :: Q         ![m2] surface of molecule
      REAL, DIMENSION(:), INTENT(IN)       :: L         ![??] UNIFAC parameter
      REAL, DIMENSION(:,:), INTENT(IN)     :: X         ![-] molar fraction of components
      REAL, DIMENSION(:,:), INTENT(OUT)    :: LNGAMMA_C ![-] log of combinatory part of gamma
      INTEGER, INTENT(IN)                  :: NMOL      ![nbr] number of molecules

!Local variables
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2)) :: RX         ![m3] volume of one component
      REAL, DIMENSION(SIZE(X,1))           :: SUMRX      ![m3] sum of volume 
      REAL, DIMENSION(SIZE(X,1))           :: SUMRXINV   ![m-3] inverse sum of volume
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2)) :: QX         ![m2] surface of one component
      REAL, DIMENSION(SIZE(X,1))           :: SUMQX      ![m2] sum of area 
      REAL, DIMENSION(SIZE(X,1))           :: SUMQXINV   ![m-2] inverse sum of area
      REAL, DIMENSION(SIZE(X,1))           :: SUMXL      ![??] sum of "L"
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2)) :: PHI        ![??] parameter in eqn 4
      REAL, DIMENSION(SIZE(X,1),SIZE(X,2)) :: THETA      ![??] parameter in eqn 4
      INTEGER                              :: I          ![idx] counter for molecules
      REAL, PARAMETER                      :: ZHALF=5.d0 ![-] constant in UNIFAC parameterization


      SUMRX(:)=0.d0
      SUMQX(:)=0.d0
      SUMXL(:)=0.d0
      DO I=1,NMOL
       RX(:,I)=X(:,I)*R(I)          !Needed for numerator in eqn 4
       QX(:,I)=X(:,I)*Q(I)          !Needed for numerator in eqn 4
       SUMRX(:)=SUMRX(:) + RX(:,I)  !Needed in denominator in eqn 4
       SUMQX(:)=SUMQX(:) + QX(:,I)  !Needed in denominator in eqn 4
       SUMXL(:)=SUMXL(:) + X(:,I)*L(I)  !Needed at the end of eqn 3
      ENDDO
      SUMRXINV(:) = 1.d0/SUMRX(:)
      SUMQXINV(:) = 1.d0/SUMQX(:)

!Solve eqn 4 (for volume ==> phi and surface ==> theta)
      DO I=1,NMOL
       PHI(:,I)   = RX(:,I)*SUMRXINV(:)
       THETA(:,I) = QX(:,I)*SUMQXINV(:)
      ENDDO
    
!Get the combinatorial part of the activity coefficient
!This is obtained through solving eqn 3
      DO I=1,NMOL
       WHERE(X(:,I).gt.0.)
          LNGAMMA_C(:,I) =                             
     &          log(PHI(:,I)/X(:,I))                   
     &          + ZHALF*Q(I)*log(THETA(:,I)/PHI(:,I))  
     &          + L(I)                                 
     &          - PHI(:,I)/X(:,I)*SUMXL(:)
       ELSEWHERE
          LNGAMMA_C(:,I)=0.d0
       ENDWHERE
      ENDDO
    
      END SUBROUTINE LNGAMMA_COMBINATORY_GET

!*********************************************************

      SUBROUTINE LNGAMMA_RESIDUAL_GET(   
     &  QG,                              !I [m2] surface of groups
     &  NU,                             !I [nbr] number of functional groups in components
     &  X,                              !I [-] molar fraction of species
     &  THTAGP,                         !I [-] area fraction of groups in pure components
     &  SI,                             !I [-] temperature dependent term
     &  LNGAMMA_R,                      !O [-] log of residual part of gamma
     &  NMOL,                           !I [nbr] number of molecules 
     &  NFUNC                          !I [nbr] number of functional groups
     &   )

!***********************************************
!Purpose:
!Get the residual part of the activity coefficients
!basically, solve eqn 7-9 in MaP05
! ************************************************
    
      IMPLICIT NONE
!INPUTS
      REAL, DIMENSION(:), INTENT(IN)        :: QG      ![m2] surface of groups
      INTEGER, DIMENSION(:,:), INTENT(IN)   :: NU      ![nbr] number of groups (j)  in molecule (i)
      REAL, DIMENSION(:,:), INTENT(IN)      :: X       ![frc] molar fraction of components
      REAL, DIMENSION(:,:), INTENT(IN)      :: THTAGP  ![frc] area fraction of groups in pure comp
      REAL, DIMENSION(:,:,:)                :: SI      ![?] temperature dependent term
      INTEGER, INTENT(IN)                   :: NMOL    ![nbr] number of molecules
      INTEGER, INTENT(IN)                   :: NFUNC   ![nbr] number of functional groups


!OUTPUTS
      REAL, DIMENSION(:,:)                  :: LNGAMMA_R  ![-] log of residual part of act. coeff
    
!LOCAL
      REAL, DIMENSION(SIZE(X,1))            :: SUMXGM     ![nbr] number of groups in mixtures
      REAL, DIMENSION(SIZE(X,1))            :: SUMXGMINV  ![1/nbr] (number of groups in mix)**(-1)
      REAL, DIMENSION(SIZE(X,1))            :: SUMJX      ![nbr] total number of one group
      REAL, DIMENSION(SIZE(X,1),SIZE(QG))   :: XGM        ![-] fraction of groups (of total)
      REAL, DIMENSION(SIZE(X,1))            :: SUMTGM     ![m2] sum of functional group area
      REAL, DIMENSION(SIZE(X,1))            :: SUMTGMINV  ![m-2] 1/"sum of functional group area"
      REAL, DIMENSION(SIZE(X,1),SIZE(QG))   :: THTAGM     ![frc] surface fraction of groups
      REAL, DIMENSION(SIZE(X,1),SIZE(QG))   :: TTSIM      ![-] mix: part of term 2 of eq 8 (sum{theta x psi})
      REAL, DIMENSION(SIZE(X,1))            :: TERM3      ![-] term 3 of equation 8
      REAL, DIMENSION(SIZE(X,1), SIZE(QG))  :: TTSIP      ![-] pure: part of term 2 of eq 8 (sum{theta x psi})
      REAL, DIMENSION(SIZE(X,1), SIZE(X,2),SIZE(QG))  :: GAMPLN     ![-] log of large gamma, pure comp
      REAL, DIMENSION(SIZE(X,1), SIZE(QG))  :: GAMMLN     ![-] log of large gamma, mix
      INTEGER                               :: I          ![idx] counter for molecules
      INTEGER                               :: J, J1,J2   ![idx] counter for functional groups

!Get sum of "molar fractions" of functional groups 
!(i.e. NOT molar fractions of components)
      SUMXGM(:)=0.d0  !sum of molar fraction of functional groups in pure state
      DO J=1,NFUNC
       DO I=1,NMOL
          SUMXGM(:)=SUMXGM(:) + DBLE(NU(I,J))*X(:,I)
       ENDDO
      ENDDO
      SUMXGMINV(:)=1.d0/SUMXGM(:) !ZTOT7 in orig code

      DO J=1,NFUNC
       SUMJX(:)=0.d0  !Sum of one particular functional group
       DO I=1,NMOL
          SUMJX(:)=DBLE(NU(I,J))*X(:,I) + SUMJX(:)
       ENDDO
       XGM(:,J)=SUMJX(:)*SUMXGMINV(:)  !Mole fraction of one group in mixture
      ENDDO

!Get total group area in mixture
      SUMTGM(:)=0.d0
      DO J=1,NFUNC
       SUMTGM(:)=SUMTGM(:)      !Old sum
     &       + XGM(:,J)*QG(J)   !Sum of functional group area of component j
      ENDDO
      SUMTGMINV(:)=1.d0/SUMTGM(:)
    
!Get theta (eqn 9) for the mixture
      DO J=1,NFUNC
       THTAGM(:,J)=QG(J)*XGM(:,J)*SUMTGMINV(:)
      ENDDO

!Get values of SI (temperature dependent), don't need to be in the iteration!
!eqn 9 in MaP05
!Sent in as input to save computer time
!TEMPKINV(:) = 1.d0/TEMPK(:)
!DO J1=1,NFUNC
!   DO J2=1,NFUNC
!      SI(:,J1,J2) = EXP(-A(J1,J2)*TEMPKINV(:))
!   ENDDO
!ENDDO
    
!Get capital gammas, (eqn 8) for pure and mix
      DO J1=1,NFUNC
      TTSIM(:,J1)=0.d0
       DO J2=1,NFUNC
          !TTSIM is part of second term in eqn 8
          TTSIM(:,J1)=TTSIM(:,J1)+THTAGM(:,J2)*SI(:,J2,J1)
       ENDDO
      ENDDO
    
!DONE ALL INTERMEDIATE EQNATIONS ==> ATTACK EQN 8
      DO J1=1,NFUNC
       TERM3(:)=0.d0   !Term 3 in eqn 8 (ZOT10 in Saxena's code)
       !Now: first get the third term, and then solve all eqn 8 for all func. groups
       DO J2=1,NFUNC
          TERM3(:)=TERM3(:)+THTAGM(:,J2)*SI(:,J1,J2)/TTSIM(:,J2)
       ENDDO
       !Now solve the whole eqn 8 for mixture
       GAMMLN(:,J1)=QG(J1)*(1. - LOG(TTSIM(:,J1)) - TERM3(:))
      ENDDO

!DO THE SAME FOR PURE COMPONENTS
      DO I=1,NMOL
       DO J1=1,NFUNC
          TTSIP(:,J1)=0.d0
          DO J2=1,NFUNC
             !Sum of seconnd term in eqn 8
             TTSIP(:,J1)=TTSIP(:,J1)          !old value
     &             + THTAGP(I,J2)*SI(:,J2,J1)
          ENDDO
       ENDDO

       DO J1=1,NFUNC
          TERM3(:)=0.d0  !ZOT11 in original code
          DO J2=1,NFUNC
             !Get term 3 in equation 8 for pure components
             TERM3(:)=TERM3(:) + THTAGP(I,J2)*SI(:,J1,J2)/TTSIP(:,J2)
          ENDDO
          
          !Get capital gamma (eqn 8) for pure components
          GAMPLN(:,I,J1)=QG(J1)*(1.d0 - LOG(TTSIP(:,J1)) -TERM3(:))
          
       ENDDO
       
      ENDDO  !Loop on moles

!Get the log of the residual part of the activity coefficient 
!Equation 7 in MaP05
      DO I=1,NMOL
       LNGAMMA_R(:,I)=0.d0
       DO J=1,NFUNC
          LNGAMMA_R(:,I)=LNGAMMA_R(:,I)     
     &          + DBLE(NU(I,J))*(GAMMLN(:,J)-GAMPLN(:,I,J))
       ENDDO
      ENDDO


      END SUBROUTINE LNGAMMA_RESIDUAL_GET

!******************************************************************

      SUBROUTINE UNIFAC_INI(            
     &  QG,                              !I [m2] surface of functional groups
     &  RG,                              !I [m3] volume of functional groups
     &  NU,                              !I [nbr] number of functional groups in molec
     &  THTAGP,                          !O [frc] surface fraction of group (j) in molecule (i)
     &  Q,                               !O [m2] surface of molecule
     &  R,                               !O [m3] volume of molecule
     &  L,                               !O [?] UNIFAC parameter for molecule
     &  NMOL,                            !I [nbr] number of molecules used
     &  NFUNC                            !I [nbr] number of functional groups used
     &  )
    
!Purpose: Set all non-time varying variables needed for UNIFAC formulation
!Following the equations given in 
!Marcolli and Peter, ACP, vol 5, pp 1501-1527, 2005
!http://overview.sref.org/1680-7324/acp/2005-5-1545
!In the following I follow notations of their eqn 1-9
    
      IMPLICIT NONE
    
!IN
      INTEGER                                    :: NMOL  ![nbr] number of molecules
      INTEGER                                    :: NFUNC ![nbr] number of functional groups
      REAL, DIMENSION(NFUNC), INTENT(IN)         :: QG    ![m2] group surface parameter
      REAL, DIMENSION(NFUNC), INTENT(IN)         :: RG    ![m3] group volume parameter
      INTEGER, DIMENSION(NMOL,NFUNC), INTENT(IN) :: NU    ![nbr] number of groups (j) in comp (i)
      REAL, DIMENSION(NMOL,NFUNC), INTENT(OUT)   :: THTAGP![-] sfc fraction of group (j) in comp (i)
    
!OUT
      REAL, DIMENSION(NMOL), INTENT(OUT)         :: Q     ![m2] total surface area of molec
      REAL, DIMENSION(NMOL), INTENT(OUT)         :: R     ![m3] total volume of molec
      REAL, DIMENSION(NMOL), INTENT(OUT)         :: L     ![?] IUPAC parameter
    
!LOCAL VARIABLES
      REAL                                       :: SUMXGP    ![-] sum of groups in one comp
      REAL                                       :: SUMXGPINV ![-] 1/"sum of groups in one comp"
      REAL, DIMENSION(NMOL,NFUNC)                :: XGP       ![-] fraction of groups in one comp
      REAL                                       :: SUMTGP    ![m2] sum sfc area of one group in one comp
      REAL                                       :: SUMTGPINV ![m-2] 1/"sum sfc area of one group in one comp"
      INTEGER                                    :: I         ![idx] counter for components
      INTEGER                                    :: J         ![idx] counter for functional group
    
      REAL, PARAMETER          :: ZHALF=5.d0               ! parameter in UNIFAC 
  
!First: get total surface (Q) and volume (R) for all molecules
      Q(:)=0.d0
      R(:)=0.d0
    
!This can be moved to an initialization routine, cause it will always
!give same answer for Q and R
!Solve eqn 6:
      DO J=1,NFUNC
       DO I=1,NMOL
          Q(I) = Q(I) + dble(NU(I,J))*QG(J) !==> total surface of molec
          R(I) = R(I) + dble(NU(I,J))*RG(J) !==> total volume of molec
       ENDDO
      ENDDO
    
!Get the parameter L (eqn 5)
      DO I=1,NMOL
       L(I)=ZHALF*(R(I) - Q(I)) - (R(I) -1.d0)
      ENDDO
    
!==> all the below can be done in initialization
!Get molar fraction of group n in pure component
      DO I=1,NMOL
       SUMXGP = 0.d0
       DO J=1,NFUNC
          SUMXGP = SUMXGP + DBLE(NU(I,J)) !sum of molar groups in ONE component
       ENDDO
       SUMXGPINV = 1.d0/SUMXGP
       DO J=1,NFUNC
          XGP(I,J)=DBLE(NU(I,J))*SUMXGPINV  !Molar fraction of group j in pure component i
       ENDDO
       ENDDO
    
!Get the area fraction
      DO I=1,NMOL
       SUMTGP = 0.d0
       DO J=1,NFUNC
          SUMTGP = SUMTGP + XGP(I,J)*QG(J)  !Sum of surface area of all groups in one comp
       ENDDO
       SUMTGPINV = 1.d0/SUMTGP   !Inverse sum of surface area of all groups in one comp
       DO J=1,NFUNC
          !Get theta, (eqn 9) for pure components
          THTAGP(I,J) = QG(J)*XGP(I,J)*SUMTGPINV
       ENDDO
      ENDDO
    
    
      END SUBROUTINE UNIFAC_INI

!***********************************************************************************
      SUBROUTINE ORG_UNIFAC_INI()

!Purpose: Set the values of declared in the bunifacparam module

      USE MODD_BUNIFACPARAM 

      IMPLICIT NONE

!Notes: 8 primary compounds: 
!             nC29
!             C4 diacid
!             naphthalene - 2,6-diacid
!             benzo(ghi) perylene
!             biomarker compound
!             4-carboxybenzoic acid
!             C18 acid
!             1,5-ditertbutyl-4,8-di(dimethylpropyl)-decalin

!order for functional groups: CH3, CH2, CH, C, C=C, C=CH
!aroC, pahC, tol, ethylbenz, OH, phen, ketone, alde, COOH, NO2,
!CH3C=O, CHOC(oxide), CHnONO2(nitrate,n=0), CHnONO2(n=1),CHnONO2(n=2), 
!CHnOOH(hydroperoxy,n=0),CHnOOH(n=1),CH2=CH -----by ddw-mpmpo

!Group volume parameters
      RG_ORG(1:NFUNC_ORG) = (/0.9011, 0.6744, 0.4469, 0.2195,   
     &    0.6605, 0.8886,  
     &    0.5313, 0.3652, 1.2663, 1.0396, 1.0000, 0.8952,         
     &    1.6724, 0.9980, 1.3013, 1.4199,           
     &    1.6724, 0.9103, 1.6697, 1.8971, 2.1246, 1.1320, !ddw 
     &    1.3594, 1.3454/) !ddw-mpmpo 
    
!Group surface parameters
      QG_ORG(1:NFUNC_ORG)= (/0.8480, 0.5400, 0.2280, 0.0, 0.4850,0.6760,
     &    0.4000, 0.1200, 0.9680, 0.6600, 1.2000, 0.6800,        
     &    1.4880, 0.9480, 1.2440, 1.1040,
     &    1.4480, 0.4680, 1.3282, 1.5562, 1.8682, 0.897, !ddw 
     &    1.1250, 1.1760/) !ddw-mpmpo 

!Set number of functional groups in each molecule. 
!Each column corresponds to a functional group (ddw-mpmpo)
!ddw-mpmpo: the first 8 rows in NU_ORG correspond to POA species
!ddw-mpmpo: the 9:end rows correspond to the SOA surrogates
      NU_ORG(:,1)=(/ 2, 0, 0, 0, 8, 2, 1, 12, 0, 0, 0, 2,  
     &              3, 0, 1, 2, 2, 2, 2, !CH3
     &              0, 1, 1, 1, 0, 1/)   !ddw-mpmopo
      NU_ORG(:,2)=(/27, 2, 0, 0, 11, 0, 16, 6, 0, 0, 0, 12,
     &              3, 0, 0, 0, 0, 2, 4,
     &              0, 2, 1, 0, 1, 1/)   !CH2
      NU_ORG(:,3)=(/0, 0, 0, 0, 6, 2, 0, 6, 0, 0, 0, 1, 2, 
     &              0, 0, 0, 2, 3, 2,
     &              0, 0, 0, 0, 0, 0 /)   !CH
      NU_ORG(:,4)=(/0, 0, 0, 0, 5, 0, 0, 4, 0, 0, 0, 0, 1, 
     &              0, 0, 0, 0, 0, 2,
     &              0, 0, 0, 0, 0, 0 /)   !C
      NU_ORG(:,5)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     &              0, 0, 1,0, 0, 0,   !C==C
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,6)=(/0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
     &              0, 2, 1, 1, 0, 0,   !C==CH
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,7)=(/ 0, 0, 8, 12, 0, 0, 0, 0, 2, 3, 7, 0, 
     &              0, 0, 0, 0, 0, 0, 0, !aromaticC
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,8)=(/0, 0, 2, 10, 0, 0, 0, 0, 0, 0, 2, 0, 0,
     &              0, 0, 0, 0, 0,0,  !pahC(AC)
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,9)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0,
     &              0, 0, 0, 0, 0,0, !tol(ACCH3)
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,10)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     &              0, 0, 0, 0, 0, 0,0, !EB(ACCH2)
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,11)=(/0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 
     &              0, 0, 1, 1, 1,0,
     &              0, 2, 1, 0, 1, 1/)  !OH
      NU_ORG(:,12)=(/0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
     &              0, 0, 0, 0, 0 ,0, 
     &              0, 0, 0, 0, 0, 0 /) !phenol 
      NU_ORG(:,13)=(/ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
     &              0, 0, 0, 1, 1,0, !ket
     &              0, 0, 0, 0, 0, 0 /)  
      NU_ORG(:,14)=(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
     &              0, 1, 2, 0, 1,0 , !ald
     &              1, 0, 0, 0, 0, 1/)
      NU_ORG(:,15)=(/0, 2, 2, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 
     &              2, 2, 0, 1, 0,0 ,
     &              0, 0, 0, 0, 0, 0 /)  !COOH
      NU_ORG(:,16)=(/0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1,
     &              0, 0, 0, 0, 0,2,
     &              0, 0, 0, 0, 0, 0/) !CNO2
      NU_ORG(:,17)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              1, 0, 0, 0, 1, 0/) !CH3C=O
      NU_ORG(:,18)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 1, 0, 0, 0, 0/) !CHOC
      NU_ORG(:,19)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 1, 0, 0, 1/) !CONO2
      NU_ORG(:,20)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 1, 0/) !CHONO2
      NU_ORG(:,21)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 1, 0, 0/) !CH2ONO2
      NU_ORG(:,22)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 1, 0, 0/) !C-OOH
      NU_ORG(:,23)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 1/) !CHOOH
      NU_ORG(:,24)=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &              0, 0, 0, 0, 0, 0,
     &              0, 0, 1, 1, 0, 0/) !CH2=CH

!order for functional groups: CH3, CH2, CH, C, C=C, C=CH
!aroC, pahC, tol, ethylbenz, OH, phen, ketone, alde, COOH, NO2,
!CH3C=O, CHOC(oxide), CHnONO2(nitrate,n=0), CHnONO2(n=1),CHnONO2(n=2), 
!CHnOOH(hydroperoxy,n=0),CHnOOH(n=1),CH2=CH -----by ddw-mpmpo

!CH3
      A_ORG(:,1)=(/ 0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &    -11.1200, -69.700, -69.700, 156.400, 10000.00,             
     &    26.7600, 505.700, 315.300, 5541.0,
     &   26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,-35.36/) !ddw-mpmpo
!CH2
      A_ORG(:,2)=(/ 0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &    -11.1200, -69.700, -69.700, 156.400, 10000.00,             
     &    26.7600, 505.700, 315.300, 5541.0,
     &   26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,-35.36/) !ddw-mpmpo
!CH
      A_ORG(:,3)=(/ 0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &    -11.1200, -69.700, -69.700, 156.400, 10000.00,             
     &    26.7600, 505.700, 315.300, 5541.0,
     &   26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,-35.36/) !ddw-mpmpo
!C
      A_ORG(:,4)=(/ 0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &    -11.1200, -69.700, -69.700, 156.400, 10000.00,             
     &    26.7600, 505.700, 315.300, 5541.0,
     &   26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,-35.36/) !ddw-mpmpo

!C==C
      A_ORG(:,5)=(/-200.00, -200.00, -200.00, -200.00, 0.0, 0.0, 
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,   
     &    732.200, -82.9200, 0.0, 349.200, 0.0,
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)
!C==CH
      A_ORG(:,6)=(/-200.00, -200.00, -200.00, -200.00, 0.0, 0.0, 
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,   
     &    732.200, -82.9200, 0.0, 349.200, 0.0,
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)
!aromaticC
      A_ORG(:,7)=(/61.1300, 61.1300,                            
     &    61.1300, 61.1300, 340.700, 340.700, 0.0, 0.0,        
     &    -146.800, -146.800, 89.600, 270.200, 140.100, 0.0,   
     &       62.3200, 1824.00,
     &    140.10,  344.42,  0.0, 0.0, 0.0,0.0, 0.0, 38.810/)
!pah
      A_ORG(:,8)=(/61.1300, 61.1300,                            
     &    61.1300, 61.1300, 340.700, 340.700, 0.0, 0.0,        
     &    -146.800, -146.800, 89.600, 270.200, 140.100, 0.0,   
     &       62.3200, 1824.00,
     &    140.10,  344.42,  0.0, 0.0, 0.0,0.0, 0.0, 38.810/)
!tol
      A_ORG(:,9)= (/ 76.500, 76.500, 76.500, 76.500, 4102.00, 4102.00, 
     &    167.00, 167.00, 0.0, 0.0, 25.8200, 10000.00, 365.800, 
     &    0.0,268.200, -127.800,
     &    365.80, 510.32,0.0, 0.0, 0.0, 0.0, 0.0, 74.15/)
!EB
      A_ORG(:,10)= (/ 76.500, 76.500, 76.500, 76.500, 4102.00, 4102.00, 
     &    167.00, 167.00, 0.0, 0.0, 25.8200, 10000.00, 365.800, 
     &    0.0,268.200, -127.800,
     &    365.80, 510.32,0.0, 0.0, 0.0, 0.0, 0.0, 74.15/)
!OH
      A_ORG(:,11)=(/986.500, 986.500, 986.500, 986.500,693.900,693.900,
     &    636.100, 636.100, 803.200, 803.200, 0.0, -274.500,     
     &    164.500, -404.800, -151.00, 561.600,
     &    164.500, 244.67, 818.97, 818.97, 818.97,342.92, 342.92,524.1/)
!phenol
      A_ORG(:,12)=(/912.200, 912.200,                                
     &    912.200, 912.200, 926.300, 926.300, 1174.00, 1174.00,  
     &    674.300, 674.300, -442.100, 0.0, -246.800, 0.0, 0.0, 0.0,
     &    -133.1,0.0,0.0,0.0,0.0,0.0,0.0,526.10 /)
!ket
      A_ORG(:,13)=(/476.400, 476.400, 476.400, 476.400, 524.500,   
     &    524.500, 25.7700, 25.7700, -52.100, -52.100, 84.00,  
     &    -158.800, 0.0, 128.00, -297.800, 0.0,
     &    0.0, 569.18,188.72, 188.72,188.72,380.94,380.94,182.60 /) 
!ald
      A_ORG(:,14)=(/677.00, 677.00,                                  
     &    677.00, 677.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 441.800, 
     &    0.0, -37.3600, 0.0, 0.0, 0.0,
     &    -37.36,-1.25,-179.38,-179.38,-179.38,408.88,408.88,448.8 /)
!COOH
      A_ORG(:,15)=(/663.500, 663.500, 663.500,                       
     &    663.500, 730.400, 730.400, 537.400, 537.400, 603.800, 
     &    603.800, 119.00, 0.0, 669.400, 0.0, 0.0, 0.0 ,
     &    669.4, 600.78,1173.3,1173.3,1173.3,1479.0,1479.0,318.90 /)
!CNO2
      A_ORG(:,16)=(/543.00,                                           
     &    543.00, 543.00, 543.00, 0.0, 0.0, 194.900, 194.900,    
     &    4448.00, 4448.00, 157.100, -413.48, 0.0, 0.0, 0.0, 0.0,
     &     548.5, 0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
!CH3C=O (carbonyl)
      A_ORG(:,17)=(/ 476.40, 476.40,476.40,476.40,182.60,182.60,
     &      25.77, 25.77, 
     &      -52.1, -52.1, 84.0,-356.1,0.0, 128.0, -297.80,-101.5, 
     &      0.0,569.18, 188.72,188.72,188.72, 380.94, 380.94,182.60/) !ddw-mpmpo
!CHOC (oxide)
      A_ORG(:,18)=(/408.3,408.3,408.3,408.3, 219.9, 219.9,171.49, 
     &      171.49, -184.68, -184.68, 6.39, 0.0, 0.0, 79.71,12.55,0.0,
     &      -288.93, 0.0, 0.0,0.0,0.0,0.0,0.0, 219.9/) !ddw-mpmpo
!CONO2
      A_ORG(:,19)=(/500.95,500.95,500.95,500.95,10326.0,10326.0,
     &     0.0, 0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0,-100.17,0.0,
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!CHONO2
      A_ORG(:,20)=(/500.95,500.95,500.95,500.95,10326.0,10326.0, 0.0,
     &      0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0, -100.17,0.0,
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!CH2ONO2
      A_ORG(:,21)=(/500.95,500.95,500.95,500.95,10326.0,10326.0, 0.0,
     &      0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0, -100.17,0.0,
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!C-OOH
      A_ORG(:,22)=(/977.56,977.56,977.56,977.56,475.91,475.91,0.0,
     &      0.0,0.0,0.0, -330.28, 0.0, 0.0, -387.63, -501.23, 0.0, 
     &      -350.58, 0.0, 545.66, 545.66,545.66,0.0, 0.0, 475.91/) !ddw-mpmpo
!CHOOH
      A_ORG(:,23)=(/977.56,977.56,977.56,977.56,475.91,475.91,0.0,
     &      0.0,0.0,0.0, -330.28, 0.0, 0.0, -387.63, -501.23, 0.0, 
     &      -350.58, 0.0, 545.66, 545.66,545.66,0.0, 0.0, 475.91/) !ddw-mpmpo
!CH2=CH
      A_ORG(:,24)=(/-200.00, -200.00, -200.00, -200.00, 0.0, 0.0, 
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,   
     &    732.200, -82.9200, 0.0, 349.200, 0.0,
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)

      END SUBROUTINE ORG_UNIFAC_INI


!***************************************************************************************
      SUBROUTINE AQ_UNIFAC_INI()

!Purpose: Set the values of declared in the bunifacparam module

      USE modd_aunifacparam

      IMPLICIT NONE

!Functional groups (17+8 in total):
!CH3, CH2, CH, C=C, C=CH, ACH(aromatic), AC(aromatic), tol, EB, OH
!phenol, ketone, aldehyde, COOH, aroNO2, water, C,
!CH3C=O, CHOC(oxide), CHnONO2(nitrate,n=0), CHnONO2(n=1),CHnONO2(n=2), 
!CHnOOH(hydroperoxy,n=0),CHnOOH(n=1),CH2=CH -----by ddw-mpmpo

!Group volume parameters
      RG_AQ(1:NFUNC_AQ) = (/0.9011, 0.6744, 0.4469, 0.6605, 0.8886, 
     &    0.5313, 0.3652, 1.2663, 1.0396, 1.0000, 0.8952,  
     &    1.6724, 0.9980, 1.3013, 1.4199, 0.9200, 0.2195,
     &    1.6724, 0.9103, 1.6697, 1.8971, 2.1246, 1.1320, !ddw 
     &    1.3594, 1.3454/) !ddw-mpmpo 

!Group surface parameters
      QG_AQ(1:NFUNC_AQ)= (/0.8480, 0.5400, 0.2280, 0.4850, 0.6760, 
     &    0.4000, 0.1200, 0.9680, 0.6600, 1.2000, 0.6800, 
     &    1.4880, 0.9480, 1.2440, 1.1040, 1.4000, 0.0, 
     &    1.4480, 0.4680, 1.3282, 1.5562, 1.8682, 0.897, !ddw 
     &    1.1250, 1.1760/) !ddw-mpmpo 

!Set number of functional groups in each molecule. 
!Each row corresponds to a molecule ==> 
!Functional groups:
                !CH3, CH2, CH, C=C, C=CH, aromC, pah, tol, EB, OH(10), phenol, ketone, aldehyde(13), COOH, aroNO2, water, C
                !CH3C=O(18), CHOC(oxide,19), CHnONO2(20,nitrate,n=0), CHnONO2(n=1),CHnONO2(n=2), 
                !CHnOOH(hydroperoxy,n=0),CHnOOH(n=1),CH2=CH -----by ddw-mpmpo
      NU_AQ( 1,:)=(/0   ,0   ,0   ,0   ,0    ,2    ,0   ,2   ,0  
     &        ,0   ,1      ,0      ,0        ,1    ,1     ,0   ,0, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 2,:)=(/0   ,0   ,0   ,0   ,0    ,3    ,0   ,2   ,0  
     &        ,0   ,0      ,0      ,1        ,1    ,0     ,0   ,0, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 3,:)=(/0   ,0   ,0   ,0   ,0    ,7    ,2   ,1   ,0  
     &         ,0   ,0      ,0      ,0        ,0    ,1     ,0   ,0, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 4,:)=(/2  ,12   ,1   ,0   ,0    ,0    ,0   ,0   ,0  
     &        ,1   ,0      ,0      ,0        ,0    ,1     ,0   ,0, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 5,:)=(/3   ,3   ,2   ,0   ,0    ,0    ,0   ,0   ,0  , 
     &        1   ,0      ,0      ,0        ,0    ,1     ,0   ,1, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 6,:)=(/0   ,0   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,0  
     &         ,0      ,0      ,0        ,2    ,0     ,0   ,0, 
     &          0   ,   0      ,0      ,0        ,0    ,0     ,0   ,0 /)
      NU_AQ( 7,:)=(/1   ,0   ,0   ,0   ,2    ,0    ,0   ,0   ,0  ,0  
     &         ,0      ,0      ,1        ,2    ,0     ,0   ,0, 
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 8,:)=(/2   ,0   ,0   ,1   ,1    ,0    ,0   ,0   ,0  ,1 
     &         ,0      ,0      ,2        ,0    ,0     ,0   ,0,
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ( 9,:)=(/2   ,0   ,2   ,0   ,1    ,0    ,0   ,0   ,0  ,1  
     &         ,0      ,1      ,0        ,1    ,0     ,0   ,0,
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ(10,:)=(/2   ,2   ,3   ,0   ,0    ,0    ,0   ,0   ,0  ,1  
     &         ,0      ,1      ,1        ,0    ,0     ,0   ,0,
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ(11,:)=(/2   ,4   ,2   ,0   ,0    ,0    ,0   ,0   ,0  ,0  
     &         ,0      ,0      ,0        ,0    ,2     ,0   ,2,
     &         0   ,0      ,0      ,0        ,0    ,0     ,0   ,0 /) 
      NU_AQ(12,:)=(/0   ,0   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,0  
     &         ,    0   ,0   ,1   ,0   ,0    ,0    ,0, 
     &              1   ,0   ,0   ,0   ,0    ,0    ,0,  0 /)!S12:MGLY 
      NU_AQ(13,:)=(/1   ,2   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,2  
     &         ,    0   ,0   ,0   ,0   ,0    ,0    ,0, 
     &              0   ,1   ,0   ,0   ,0    ,0    ,0,  0 /)!S13:IEPOX 
      NU_AQ(14,:)=(/1   ,1   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,1  
     &         ,    0   ,0   ,0   ,0   ,0    ,0    ,0, 
     &              0   ,0   ,1   ,0   ,0    ,0    ,0,  1 /)!S14:IHN2 
      NU_AQ(15,:)=(/1   ,0   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,0  
     &         ,    0   ,0   ,0   ,0   ,0    ,0    ,0, 
     &              0   ,0   ,0   ,0   ,1    ,1    ,0,  1 /)!S15:INPB 
      NU_AQ(16,:)=(/0   ,1   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,1  
     &         ,    0   ,0   ,0   ,0   ,0    ,0    ,0, 
     &              1   ,0   ,0   ,1   ,0    ,0    ,0,  0 /)!S16:MVK 
      NU_AQ(17,:)=(/1   ,1   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,1  
     &         ,    0   ,0   ,1   ,0   ,0    ,0    ,0, 
     &              0   ,0   ,1   ,0   ,0    ,0    ,1,  0/)!S17:dinitrates 
      NU_AQ(18,:)=(/0   ,0   ,0   ,0   ,0    ,0    ,0   ,0   ,0  ,0  
     &         ,    0   ,0   ,0   ,0   ,0    ,1    ,0, 
     &              0   ,0   ,0   ,0   ,0    ,0    ,0,  0/)!Old S12 (should be water) ddw 

!Copy and paste values from Griffin
!interaction parameters group-group

!Functional groups (17+8 in total):
!CH3, CH2, CH, C=C, C=CH, aromatic carbon, pah, tol, EB, OH
!phenol, ketone, aldehyde, COOH, aroNO2, water, C,
!CH3C=O, CHOC(oxide), CHnONO2(nitrate,n=0), CHnONO2(n=1),CHnONO2(n=2), 
!CHnOOH(hydroperoxy,n=0),CHnOOH(n=1),CH2=CH -----by ddw-mpmpo

!Copy and paste values from Griffin
!interaction parameters group-group

!CH3
      A_AQ(:,1)=(/ 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &      -11.1200, -69.700, -69.700, 156.400, 10000.00,     
     &      26.7600, 505.700, 315.300, 5541.00, 300.0, 0.0,
     &      26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,
     &     -35.36/) !ddw-mpmpo
!CH2
      A_AQ(:,2)=(/ 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &      -11.1200, -69.700, -69.700, 156.400, 10000.00,     
     &      26.7600, 505.700, 315.300, 5541.00, 300.0, 0.0,
     &      26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,
     &     -35.36/) !ddw-mpmpo
!CH
      A_AQ(:,3)=(/ 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &      -11.1200, -69.700, -69.700, 156.400, 10000.00,     
     &      26.7600, 505.700, 315.300, 5541.00, 300.0, 0.0,
     &      26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,
     &     -35.36/) !ddw-mpmpo
!C=C
      A_AQ(:,4) =(/-200.00, -200.00, -200.00, 0.0, 0.0,      
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,       
     &    732.200, -82.9200, 0.0, 349.200, 0.0, 692.7, -200.00, 
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)
!C=CH
      A_AQ(:,5) =(/-200.00, -200.00, -200.00, 0.0, 0.0,      
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,       
     &    732.200, -82.9200, 0.0, 349.200, 0.0, 692.7, -200.00, 
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)

!aro(ACH)
      A_AQ(:,6) = (/61.1300, 61.1300, 61.1300, 340.700, 340.700, 0.0, 
     &    0.0, -146.800, -146.800, 89.600, 270.200, 140.100, 0.0,     
     &    62.3200, 1824.00, 362.3, 61.1300,
     &    140.10,  344.42,  0.0, 0.0, 0.0,0.0, 0.0, 38.810/)
!pah(AC:aromatic)
      A_AQ(:,6) = (/61.1300, 61.1300, 61.1300, 340.700, 340.700, 0.0, 
     &    0.0, -146.800, -146.800, 89.600, 270.200, 140.100, 0.0,     
     &    62.3200, 1824.00, 362.3, 61.1300,
     &    140.10,  344.42,  0.0, 0.0, 0.0,0.0, 0.0, 38.810/)

!tol(ACCH3)
      A_AQ(:,8) = (/                                                 
     &    76.500, 76.500, 76.500, 4102.00, 4102.00, 167.00,      
     &    167.00, 0.0, 0.0, 25.8200, 10000.00, 365.800, 0.0,     
     &    268.200, -127.800, 377.6, 76.500,
     &    365.80, 510.32,0.0, 0.0, 0.0, 0.0, 0.0, 74.15/)
!EB(ACCH2)
      A_AQ(:,9) = (/                                                 
     &    76.500, 76.500, 76.500, 4102.00, 4102.00, 167.00,      
     &    167.00, 0.0, 0.0, 25.8200, 10000.00, 365.800, 0.0,     
     &    268.200, -127.800, 377.6, 76.500,
     &    365.80, 510.32,0.0, 0.0, 0.0, 0.0, 0.0, 74.15/)
!OH
      A_AQ(:,10) = (/                                                
     &    986.500, 986.500, 986.500, 693.900, 693.900,           
     &    636.100, 636.100, 803.200, 803.200, 0.0, -274.500,     
     &    164.500, -404.800, -151.00, 561.600, -229.1, 986.500,
     &    164.500, 244.67, 818.97, 818.97, 818.97,342.92, 342.92,524.1/)

!phenol(ACOH:aromatic carbon-alcohol)
      A_AQ(:,11) = (/                                                 
     &    912.200, 912.200, 912.200, 926.300, 926.300, 1174.00, 1174.00,      
     &    674.300, 674.300, -442.100, 0.0, -246.800, 0.0, 0.0,       
     &    0.0, 324.5, 912.200,
     &    -133.1,0.0,0.0,0.0,0.0,0.0,0.0,526.10 /)
!ketone
!ddw-mpmpo: treat ketone as carbonyl
      A_AQ(:,12) = (/                                             
     &    476.400, 476.400, 476.400, 524.500,                     
     &    524.500, 25.7700, 25.7700, -52.100, -52.100, 84.00,     
     &    -158.800, 0.0, 128.00, -297.800, 0.0, -195.4, 476.400,
     &    0.0, 569.18,188.72, 188.72,188.72,380.94,380.94,182.60 /) 
!aldehyde
      A_AQ(:,13) = (/                                             
     &    677.00,                                                 
     &    677.00, 677.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 441.800,  
     &    0.0, -37.3600, 0.0, 0.0, 0.0, -257.3, 677.00,
     &    -37.36,-1.25,-179.38,-179.38,-179.38,408.88,408.88,448.8 /)
!COOH
      A_AQ(:,14) = (/                                               
     &    663.500, 663.500,                                      
     &    663.500, 730.400, 730.400, 537.400, 537.400, 603.800,  
     &    603.800, 119.00, 0.0, 669.400, 0.0, 0.0, 0.0, -14.090,663.500,
     &    669.4, 600.78,1173.3,1173.3,1173.3,1479.0,1479.0,318.90 /)
!aro NO2
      A_AQ(:,15) = (/                                                
     &       543.00, 543.00, 543.00, 0.0, 0.0, 194.900, 194.900,    
     &       4448.00, 4448.00, 157.100, -413.48, 0.0, 0.0, 0.0, 0.0,    
     &       399.50, 543.00,
     &     548.5, 0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
!water
      A_AQ(:,16) = (/                                               
     &       1318.00, 1318.00, 1318.00, 634.2, 634.2, 903.8,        
     &       903.8, 5695.0, 5695.0, 353.5, -601.8, 472.5, 232.7,    
     &       -66.17, 360.7, 0.0, 1318.00,
     &   472.5,833.21,681.78,681.78,681.78,795.55,795.55,270.60/)
!C
      A_AQ(:,17)=(/ 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200,  
     &      -11.1200, -69.700, -69.700, 156.400, 10000.00,     
     &      26.7600, 505.700, 315.300, 5541.00, 300.0, 0.0,
     &      26.76, 21.49, -75.718,-75.718,-75.718,-23.233,-23.233,
     &     -35.36/) !ddw-mpmpo
!CH3C=O (carbonyl)
      A_AQ(:,18)=(/ 476.40, 476.40,476.40,182.60,182.60,25.77, 25.77, 
     &      -52.1, -52.1, 84.0,-356.1,0.0, 128.0, -297.80,-101.5,
     &      -195.4,476.40, 
     &      0.0,569.18, 188.72,188.72,188.72, 380.94, 380.94,182.60/) !ddw-mpmpo
!CHOC (oxide)
      A_AQ(:,19)=(/408.3,408.3,408.3, 219.9, 219.9,171.49,171.49, 
     &      -184.68, -184.68, 6.39, 0.0, 0.0, 79.71,12.55,0.0,
     &      -144.77,408.3,
     &      -288.93, 0.0, 0.0,0.0,0.0,0.0,0.0, 219.9/) !ddw-mpmpo
!CONO2
      A_AQ(:,20)=(/500.95,500.95,500.95,10326.0,10326.0, 0.0,
     &      0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0, -100.17, 0.0, 
     &      142.65, 500.95, 
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!CHONO2
      A_AQ(:,21)=(/500.95,500.95,500.95,10326.0,10326.0, 0.0,
     &      0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0, -100.17, 0.0, 
     &      142.65, 500.95, 
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!CH2ONO2
      A_AQ(:,22)=(/500.95,500.95,500.95,10326.0,10326.0, 0.0,
     &      0.0,0.0,0.0,37.631, 0.0, 0.0, 402.0, -100.17, 0.0, 
     &      142.65, 500.95, 
     &      -197.93,0.0,0.0,0.0,0.0,-86.279,-86.279,10326./) !ddw-mpmpo
!C-OOH
      A_AQ(:,23)=(/977.56,977.56,977.56,475.91,475.91, 0.0,
     &      0.0,0.0,0.0, -330.28, 0.0, 0.0, -387.63, -501.23,
     &      0.0, -341.18, 977.56, 
     &      -350.58, 0.0, 545.66, 545.66,545.66,0.0, 0.0, 475.91/) !ddw-mpmpo
!CHOOH
      A_AQ(:,24)=(/977.56,977.56,977.56,475.91,475.91, 0.0,
     &      0.0,0.0,0.0, -330.28, 0.0, 0.0, -387.63, -501.23,
     &      0.0, -341.18, 977.56, 
     &      -350.58, 0.0, 545.66, 545.66,545.66,0.0, 0.0, 475.91/) !ddw-mpmpo
!CH2=CH
      A_AQ(:,25) =(/-200.00, -200.00, -200.00, 0.0, 0.0,      
     &    -97.7800, -97.7800, -269.700, -269.700, 8694.00,       
     &    732.200, -82.9200, 0.0, 349.200, 0.0, 692.7, -200.00, 
     &    42.92, -2.80, -294.43, -294.43, -294.43,-57.949,-57.949,0.0 /)

      END SUBROUTINE AQ_UNIFAC_INI

!*************************************************************************

      SUBROUTINE ZSR_INI()
      use modd_binsolu
      implicit none
      molalbin(:,1) = (/ 555.600, 0.65000, 0.30762, 0.19305, 
     &   .13537, 0.10038, 0.07666, 0.05927, 0.04567, 0.03431, 0.00000 /)
      molalbin(:,2) = (/ 555.600, 1.98908, 0.95252, 0.60644, 
     &  0.43292, 0.32835, 0.25820, 0.20764, 0.16923, 0.13830, 0.00000 /)
      molalbin(:,3) = (/ 555.600,63.80219,31.65389,20.93034, 
     & 15.57215,12.35612,10.21267, 8.68083, 7.53173, 6.63831, 0.00000 /)
      molalbin(:,4) = (/ 555.600, 2.81667, 1.37478, 0.89336, 
     &  0.65205, 0.50667, 0.40951, 0.33965, 0.28692, 0.24558, 0.00000 /)
      molalbin(:,5) = (/ 555.600, 1.37882, 0.63792, 0.39248, 
     &  0.27079, 0.19847, 0.15067, 0.11672, 0.09124, 0.07117, 0.00000 /)
      molalbin(:,6) = (/ 555.600, 0.49396, 0.22574, 0.13569, 
     &  0.09018, 0.06232, 0.04319, 0.02890, 0.01740, 0.00755, 0.00000 /)
      molalbin(:,7) = (/ 555.600, 0.45742, 0.21651, 0.13549, 
     &  0.09439, 0.06918, 0.05184, 0.03891, 0.02848, 0.01921, 0.00000 /)
      molalbin(:,8) = (/ 555.600, 0.86272, 0.40057, 0.24605, 
     &  0.16855, 0.12160, 0.08988, 0.06671, 0.04860, 0.03329, 0.00000 /)
      molalbin(:,9) = (/ 555.600, 0.55832, 0.26300, 0.16412, 
     &  0.11418, 0.08375, 0.06308, 0.04780, 0.03574, 0.02543, 0.00000 /)
      molalbin(:,10)= (/ 555.600, 0.92947, 0.42461, 0.25726, 
     &  0.17392, 0.12435, 0.09149, 0.06809, 0.05035, 0.03601, 0.00000 /)
      molalbin(:,11) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
      ! ddw-mpmpo
      molalbin(:,12) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
      molalbin(:,13) = (/ 555.600,63.80219,31.65389,20.93034, 
     & 15.57215,12.35612,10.21267, 8.68083, 7.53173, 6.63831, 0.00000 /)
      molalbin(:,14) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
      molalbin(:,15) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
      molalbin(:,16) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
      molalbin(:,17) =(/ 555.600, 0.54200, 0.25400, 0.15800, 
     &  0.11000, 0.08000, 0.06000, 0.04600, 0.03400, 0.02400, 0.00000 /)
    
      END SUBROUTINE ZSR_INI

!*************************************************************************************  
  
      END MODULE MODE_UNIFAC


