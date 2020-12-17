!program MAIN

!  MODULE soa_main_mpmpo  
!Purpose: main program do drive the MPMPO routines
!Refernces: Griffin et al, J. atm. Chem, vol 44, pp 171-190, 2003

      SUBROUTINE MAIN_MPMPO(jjcpt, jjcpart,jjtemp,jjrh,jjlwc,jjhplus,
     &                  jjpar,jjorganion,jjgamma_aq_r,jjgamma_aq_h,
     &                   LOGDEV, LAYER)

      USE mode_oamain
      USE modd_aunifacparam
      USE modd_bunifacparam
      USE mode_unifac
      USE modd_glo, only: LBOX, LPRINT, NAAERO
      USE parameters, only:nmpmpo    !ddw-mpmpo

!MAIN PROGRAM FOR BETTY PUN ORGANIC AEROSOL EQUILIBRIUM SOLVER
!WRITTEN SO THAT IT WORKS ON ONE VECTOR OF ATMOSPHERIC POINTS

!CONTAINS
!  SUBROUTINE MAIN_MPMPO(jjcpt,jjcpart,jjtemp,jjrh,jjlwc,jjhplus, 
!     &            jjpar,jjorganion,jjh2oorg)
      IMPLICIT NONE
  
! these variables are added to exchange the input and output to MAIN_MPMPO
! ddw-mpmpo modified  jjpar(nmpmpo),jjcpart(nmpmpo)
       REAL jjcpt(8),jjcpart(nmpmpo),jjtemp,jjrh,jjlwc,jjplu
       REAL jjpar(nmpmpo),jjorganion(nmpmpo),jjgamma_aq_r(1, NAAERO),
     & jjgamma_aq_h(1, NAAERO),jjhplus
       integer LOGDEV, LAYER
! end


      INTEGER, PARAMETER  :: NPOINTS = 1    !Number of physical points in a vector

      REAL, ALLOCATABLE, DIMENSION(:)   :: TEMPK          ![K] temperature 
      REAL, ALLOCATABLE, DIMENSION(:)   :: RH             ![-] relative humidity 

      REAL, ALLOCATABLE, DIMENSION(:,:)  :: CPT           ![ug/m3] primary organic aerosol
      REAL, ALLOCATABLE, DIMENSION(:,:)  :: AERO_ORG      ![ug/m3] organic phase SOA
      REAL, ALLOCATABLE, DIMENSION(:,:)  :: AERO_AQ       ![ug/m3] aquous phase SOA
      REAL, ALLOCATABLE, DIMENSION(:,:)  :: WORG(:,:)     ![ug/m3] organic aersol concentrations
      REAL, ALLOCATABLE, DIMENSION(:,:)  :: GASORG(:,:)   ![ug/m3] organic gas concentrations
      REAL, ALLOCATABLE, DIMENSION(:,:)  :: PARTORG(:,:)  ![ug/m3] organic particle concentrations
      REAL, ALLOCATABLE, DIMENSION(:)    :: RWATER(:)     ![ug/m3] LWC available for partitioning
      REAL, ALLOCATABLE, DIMENSION(:)    :: VARTMP(:)     ![mol/kg_{water}] proton concentration
      REAL, ALLOCATABLE, DIMENSION(:)    :: ORGANION(:)   ![mol/m3] anion concentration
      REAL, ALLOCATABLE, DIMENSION(:)    :: LWCORG(:)     ![ug/m3] LWC associated with org. aerosols
      INTEGER                            :: II            ![idx] a random counter
      INTEGER                            :: JJ            ![idx] another random counter

!ddw -- output the gammas
      REAL, DIMENSION(:,:), ALLOCATABLE    :: GAMMA_AQ_HENRY   ![-] activity coefficients for aquous phase (Henry's law std. state)
      REAL, DIMENSION(:,:), ALLOCATABLE    :: GAMMA_AQ_RAOULT  ![-] activity coefficients for aquous phase (Raoult's law std. state)

!ALLOCATE THE MEMORY NEEDED:
       ALLOCATE(TEMPK(NPOINTS))
       ALLOCATE(RH(NPOINTS))
       ALLOCATE(WORG(NPOINTS,NBSP))  
       ALLOCATE(CPT(NPOINTS,NBSPOA))
       ALLOCATE(AERO_ORG(NPOINTS,NBSP))
       ALLOCATE(AERO_AQ(NPOINTS,NAAERO))
       ALLOCATE(GASORG(NPOINTS,NBSP))
       ALLOCATE(PARTORG(NPOINTS,NBSP))
       ALLOCATE(RWATER(NPOINTS))
       ALLOCATE(VARTMP(NPOINTS))
       ALLOCATE(ORGANION(NPOINTS))
       ALLOCATE(LWCORG(NPOINTS))
! ddw -- output gamma
      ALLOCATE(GAMMA_AQ_RAOULT(SIZE(AERO_ORG,1),NAAERO))
      ALLOCATE(GAMMA_AQ_HENRY(SIZE(AERO_ORG,1),NAAERO))

!Check if box model, allow for prints
        IF(NPOINTS.eq.1)LBOX=.TRUE.

!Set LPRINT
      LPRINT=.TRUE.
      LPRINT=.FALSE.
      IF(LPRINT.AND..NOT.LBOX)THEN
       stop"cannot activate LPRINT except in box version"
       ENDIF
      LBOX=.TRUE.

!Initialize the temperature
!TEMPK(:)=20.0 + 283.15d0
       TEMPK(:) = jjtemp  


!Initialze relative humidity
       !RH(:)=0.80
        RH(:) = jjrh

!Initialize the starting values: total organic concentrations
       DO II =1,NPOINTS
      !WORG (II,:) = (/ 3.7924456E-03, 2.5479761E-3, 4.1438680E-3, 1.6882788E-1 &
      !     , 6.8781158E-2 ,4.8853220E-2, 8.2842094E-4, 2.7897930E-3  &
      !     , 5.9604807E-4, 4.6155168E-7 /)
      !
       DO JJ=1,NBSP
        WORG(II,JJ)=0.01d0
      !WORG(II,:)= (/ 3.44003041E-07, 6.20606691E-12, 9.99999997E-07, 9.99999997E-07 &
      !     , 9.99999997E-07, 2.37479409E-15, 4.41267002E-06, 1.16432034E-06          &
      !     ,9.99999997E-07 ,9.99999997E-07 /)
       ENDDO
          ENDDO

! added by james
        WORG(1,1) = jjcpart(1)
        WORG(1,2) = jjcpart(2)
        WORG(1,3) = jjcpart(3)
        WORG(1,4) = jjcpart(4)
        WORG(1,5) = jjcpart(5)
        WORG(1,6) = jjcpart(6)
        WORG(1,7) = jjcpart(7)
        WORG(1,8) = jjcpart(8)
        WORG(1,9) = jjcpart(9)
        WORG(1,10) =jjcpart(10)
        WORG(1,11) =jjcpart(11)
! added by ddw-mpmpo
        WORG(1,12) = jjcpart(12)
        WORG(1,13) = jjcpart(13)
        WORG(1,14) = jjcpart(14)
        WORG(1,15) = jjcpart(15)
        WORG(1,16) = jjcpart(16)
        WORG(1,17) = jjcpart(17)

        DO JJ =1, NBSP
        IF ( WORG(1,JJ) .lt. 1.0e-12 ) then
           WORG(1,JJ) = 1.0e-12
        ENDIF
          ENDDO



!Initialize the starting values for POA
       DO II=1,NPOINTS
        DO JJ=1,NBSPOA
         CPT(II,JJ)=0.002d0
        ENDDO
       !cpt(II,:)= (/ 9.99999997E-07, 9.99999997E-07, 9.99999997E-07, 9.99999997E-07 &
       !     ,9.99999997E-07, 9.99999997E-07, 9.99999997E-07, 9.99999997E-07 /)
        ENDDO

! added by james
        CPT(1,1) = jjcpt(1)
        CPT(1,2) = jjcpt(2)
        CPT(1,3) = jjcpt(3)
        CPT(1,4) = jjcpt(4)
        CPT(1,5) = jjcpt(5)
        CPT(1,6) = jjcpt(6)
        CPT(1,7) = jjcpt(7)
        CPT(1,8) = jjcpt(8)

        DO JJ=1, NBSPOA
          IF ( CPT(1,JJ) .lt. 1.0e-12) then
              CPT(1,JJ) = 1.0e-12 
         ENDIF
          ENDDO

!Initialize the starting values for the gas phase concentrations
       DO II = 1,NPOINTS
       GASORG(II,:)=(/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 
     &      0.d0, 0.d0, 0.d0, 0.d0, 0.d0 , 0.d0,
     &      0.d0, 0.d0, 0.d0, 0.d0, 0.d0 , 0.d0 /) !ddw-mpmpo

!GASORG(II,:) = (/3.7924456E-03, 2.5479750E-3, 4.1438677E-3, 1.6878345E-1  &
!     , 6.8775396E-2, 4.8833221E-2, 8.2842015E-4, 2.7897911E-3           &
!     , 5.9604729E-4 ,4.6155168E-7/)
       ENDDO

!Initialize the aerosol phase concentrations
!ddw-mpmpo
       DO II = 1, NPOINTS
        PARTORG(II,:) = (/ 0., 0., 0., 0.  
     &     ,0., 0., 0., 0.            
     &     ,0., 0. ,0.               
     &     ,0., 0., 0, 0., 0., 0/)
        ENDDO
  
!LWC available for partitioning (ug/m3)
! RWATER=  18.8d0

       if (jjlwc .gt. 0.0001) then
        RWATER(:) = jjlwc
        else
        RWATER(:) = 0.0001
       endif
          

!Initialize the proton concentration for all points
         VARTMP(:) = 1.d-10
         VARTMP(:) = jjhplus

!        if ( LAYER .eq. 1) then
!        WRITE(LOGDEV,*) 'MPMPO_CNRM INPUT'
!        WRITE(LOGDEV,*) TEMPK(:), RWATER(:), VARTMP(:)
!        WRITE(LOGDEV,*) WORG(1,1), WORG(1,2), WORG(1,3),WORG(1,4)
!        WRITE(LOGDEV,*) WORG(1,5), WORG(1,6), WORG(1,7),WORG(1,8)
!        WRITE(LOGDEV,*) WORG(1,9), WORG(1,10), WORG(1,11)
!        WRITE(LOGDEV,*) CPT(1,1), CPT(1,2), CPT(1,3), CPT(1,4)
!        WRITE(LOGDEV,*) CPT(1,5), CPT(1,6), CPT(1,7), CPT(1,8)
!        endif


!Initialize the organic anion concentration (why??? this is output!!)
           ORGANION(:) = 1.000000
  
!Initialize the LWC associated with organics
           LWCORG(:)  =    0.000000  

!Set unifac coefficients for group a
!including RG, QG, NU, A (ddw-mpmpo)
      CALL AQ_UNIFAC_INI()

!Set unifac coefficients for group b
!including RG, QG, NU, A (ddw-mpmpo)
      CALL ORG_UNIFAC_INI()

!Calculate non time varying unifac stuff for aquous phase
!i.e. THTAGP (ddw-mpmpo)
      CALL UNIFAC_INI(  
     &  QG_AQ           !I [m2] surface of functional groups, input
     &  ,RG_AQ          !I [m3] volume of functional groups,input
     &  ,NU_AQ          !I [nbr] number of functional groups in molec, input
     &  ,THTAGP_AQ      !O [frc] surface fraction of group (j) in molecule (i), input
     &  ,Q_AQ           !O [m2] surface of molecule, output
     &  ,R_AQ           !O [m3] volume of molecule, output
     &  ,L_AQ           !O [?] UNIFAC parameter for molecule, output
     &  ,NMOL_AQ        !I [nbr] number of molecules used, input
     &  ,NFUNC_AQ       !I [nbr] number of functional groups used, input
     &  )
!Calculate non time varying unifac stuff for group organic phase
!i.e. THTAGP (ddw-mpmpo)
       CALL UNIFAC_INI(  
     &  QG_ORG           !I [m2] surface of functional groups, input
     &  ,RG_ORG          !I [m3] volume of functional groups, input
     &  ,NU_ORG          !I [nbr] number of functional groups in molec, input
     &  ,THTAGP_ORG      !O [frc] surface fraction of group (j) in molecule (i), input
     &  ,Q_ORG           !O [m2] surface of molecule, output
     &  ,R_ORG           !O [m3] volume of molecule, output
     &  ,L_ORG           !O [?] UNIFAC parameter for molecule, output
     &  ,NMOL_ORG        !I [nbr] number of molecules used, input
     &  ,NFUNC_ORG       !I [nbr] number of functional groups used, input
     &  )

!Set molality of solvent in binary mix with water at several RH
!including molalbin
        CALL ZSR_INI()

!DO THE ACTUAL CALCULATIONS
        CALL MPMPO(      
     &  TEMPK           !I [K] Temperature, input
     &  ,RH             !I [-] relative humidity, input
     &  ,WORG           !I [ug/m3] organic total conconentrations, input
     &  ,CPT            !I [ug/m3] primary organic aerosol, input
     &  ,GASORG         !O [ug/m3] organic gas concentrations, output
     &  ,AERO_ORG       !O [ug/m3] organic phase SOA, output
     &  ,AERO_AQ        !O [ug/m3] aquous phase SOA, output
     &  ,PARTORG        !O [ug/m3] total organic particle concentrations, output
     &  ,RWATER         !I [ug/m3] LWC available for partitioning, input
     &  ,VARTMP         !I [mol/kg_{water}] proton concentrations, input
     &  ,LWCORG         !O [ug/m3] liquid water concent associated with organics, output
     &  ,ORGANION       !O [mol/m3] organic anion concentrations, output
     &  ,GAMMA_AQ_Raoult, Gamma_aq_Henry   !ddw -- output the gammas
     &  ) 

! output
      DO II=1,NBSP
       jjpar(II) = PARTORG(1,II)
       jjorganion(II) = AERO_ORG(1, II)   !ddw -- output the organic part
      ENDDO
!ddw commented out this line       jjorganion = ORGANION(1)
!ddw commented out this line       jjh2oorg = LWCORG(1)
      jjgamma_aq_r(1,:) = GAMMA_AQ_Raoult(1, :)  !ddw 
      jjgamma_aq_h(1,:) = GAMMA_AQ_Henry(1, :)  !ddw 

!      print*, 'raoult in soa =', GAMMA_AQ_Raoult(1, 24)
!ddw      jjh2oorg = GAMMA_AQ_Henry(1, :)   !ddw
! end of output

! Comment out the print
!  write(6,*)"                     GAS CONCENTRATION        %DEVIATION FROM MB"
!  DO II=1,NBSP
!     write(6,*)"AFTER: GASORG",II,GASORG(:,II), (WORG(:,II)-PARTORG(:,II)-GASORG(:,II))/WORG(:,II)*100.d0,"%"
!  ENDDO
!  write(6,*)"                     PARTICLE CONCENTRATION     %OF TOTAL       "
!  DO II=1,NBSP
!     write(6,*)"AFTER: PARTORG",II,PARTORG(:,II), PARTORG(:,II)/WORG(:,II)*100,"%"
!  ENDDO
!  write(6,*)"AFTER LWCORG", LWCORG
!  write(6,*)"AFTER RWATER", RWATER
!  write(6,*)"AFTER ORGANION", ORGANION

      DEALLOCATE(TEMPK)
      DEALLOCATE(RH)
      DEALLOCATE(WORG)
      DEALLOCATE(CPT)
      DEALLOCATE(AERO_ORG)
      DEALLOCATE(AERO_AQ)
      DEALLOCATE(GASORG)
      DEALLOCATE(PARTORG)
      DEALLOCATE(RWATER)
      DEALLOCATE(VARTMP)
      DEALLOCATE(ORGANION)
      DEALLOCATE(LWCORG)

          END SUBROUTINE MAIN_MPMPO
!  END MODULE soa_main_mpmpo 
!  END PROGRAM MAIN


  
