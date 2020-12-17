       MODULE MODD_GLO

       use modd_glodef

       IMPLICIT NONE

! Change SOA surrogate number from 10 to 11 by James

! flag to use binary solution and zsr if = 1 */
! zsrflag = 0 calls unifac and solves implicit
! equation for a.c. water = RH using newt1 */ 
       INTEGER, PARAMETER :: zsrflag = 1 
     
! flag to use NEWT in Type B module 
! and type A module with absorption when = 1, 
! if = 0; don't use NEWT 
       INTEGER, PARAMETER :: Newtflag = 0

! saturationflag = 1 means to 
! use saturation to determine particulate-phase
! concentration when inorganic particle is dry. 
! If = 0, use absorption */
       INTEGER, PARAMETER  :: saturationflag = 0

! Criterion for setting intital particle 
! conc, compared to VP in torr or VP 
! rated by PAOM */
       REAL, PARAMETER :: VPCrit = 1.d-7

!Critical value for LWC above which much aerosol can dissolve
       REAL, PARAMETER :: LWCCRIT = 2.0d0 ![ug/m3] 

!Critical value for Henry's law constant above which much aerosol will dissolve
       REAL, PARAMETER :: HLCRIT = 1.d-2

       REAL, PARAMETER :: TINY = 1.d-8

!Set the universal gas constant in non SI units
       REAL, PARAMETER       :: R_UNIV = 8.206e-5  !gas constant in units of (m3*atm)/(mol*K) 
       INTEGER, PARAMETER     :: NBSPOA = 8 !Number of POA species
! changed NBSP from 10 to 11, James, and NAAERO from 17 to 18
! ddw-mpmpo: changed the NBSP from 11 to 17
       INTEGER, PARAMETER     :: NBSP = 11+6  !Number of main SOA species
! ddw-mpmpo: for now, set isoprene ions to zero
! ddw-mpmpo: so NAAERO is 18+6(new main surrogates)+2(new ions)
       INTEGER, PARAMETER     :: NAAERO= 18+6 !Number of SOA species and their ions

!Molar weight of all parameters including ions
!ddw-mpmpo: added MW for the 6 new surrogates
       REAL, PARAMETER, DIMENSION(NAAERO) ::  MW_SOA = (/  
     &  211.0, 210.0,         !Comp #1 (NK=2)
     &  178.0, 177.0,        !Comp #2 (NK=2)
     &  217.0,               !Comp #3 (NK=1)
     &  303.0,               !Comp #4 (NK=1)
     &  217.0,               !Comp #5 (NK=1)
     &  90.0, 89.0, 88.0,    !Comp #6 (NK=3)
     &  184.0, 183.0, 182.0,  !Comp #7 (NK=3)
     &  154.0,                !Comp #8 (NK=1)
     &  186.0, 185.0,         !Comp #9 (NK=2)
     &  186.0,                !Comp #10 (NK=1)
     &  260.0,                !Comp #11 (NK=1)
     &  72.0 ,                !Comp #12 (NK=1)
     &  118.0,                !Comp #13 (NK=1)
     &  147.0,                !Comp #14 (NK=1)
     &  163.0,                !Comp #15 (NK=1)
     &  149.0,                !Comp #16 (NK=1)
     &  208.0                 !Comp #17 (NK=1)
     &  /)

!Molecular weight of primary organic paerosols
       REAL, PARAMETER, DIMENSION(NBSPOA) :: 
     & MW_POA = (/ 408., 118., 216., 276., 412., 166., 284., 390. /)

!Molecular weight of water
       REAL, PARAMETER                    :: MW_WATER = 18.0  ![g/mol] molar weight of water

!Parameter needed to get the saturation vapor pressures of organics
!The below is cut and pasted from aerodriv.f recieved from Griffin
!ddw-mpmpo: need to modify HBN, TAUVP, TBOIL
!ddw-mpmpo: set to 0.0 for now
       REAL, PARAMETER, DIMENSION(NBSP)  :: HBN = (/6.7e-03,
     &   5.61e-03,0.0,3.3e-03,4.61e-03, 
     &   1.57e-02,7.68e-03,6.49e-03,7.6e-03,5.37e-03, 0.0, 
     &   1.57e-02, 0.0, 0.0, 0.0, 0.0, 0.0/) ! ddw-mpmpo: added by ddw

       REAL, PARAMETER, DIMENSION(NBSP) :: TAUVP=(/3.5,2.5,6.0,13.5,1.0,
     &    0.0,2.5,2.0,3.0,5.0, 4.0,
     &    0.0, 0.0, 0.0, 0.0, 0.0, 0.0/) ! ddw-mpmpo: added by ddw

       REAL, PARAMETER, DIMENSION(NBSP)  :: TBOIL = (/685.3,634.0,645.5,
     &  672.5,566.3,560.0,698.0,
     &  575.0,679.0,615.0, 573.9,
     &  560.0, 345., 500., 500., 500., 500.0/) ! ddw-mpmpo: assume Tboil is the boiling point in Kelvin 

!partition parameters H and K
!H is in units of m3/ug estimated based on Suzuki et al., 1992
!K is in units of [mol/kg water] (same as {H+}) with; 
!ddw-mpmpo: K is in units of [m3 air/L water] (Griffin 2003) 
!concentrations of molecules in ions in the same molar units
       REAL, PARAMETER, DIMENSION(NAAERO) :: K_298 = (/ 
     &   1.82e-2, 1.7e-3,     !Comp 1
     &   1.202e-4, 7.33e-5,   !Comp 2
     &   1.38e-6,             !Comp 3
     &   2.455e-7,            !Comp 4
     &   2.884e-5,            !Comp 5 
     &   2.512e-2, 5.4e-2, 5.2e-5,   !Comp 6,
     &   22.01, 3.7e-5, 3.9e-6,      !Comp 7
     &   1.47e-4,                    !Comp 8
     &   0.0489, 6.52e-4,            !Comp 9 
     &   9.55e-4,                    !Comp 10
     &   7.586e-6,                   !Comp 11
     &   2.512E-2, !ddw-mpmpo:Comp 12(3.24E4 M atm-1 for MGLY; 4.15E5 for GLYX) 
     &   1.952,  !ddw-mpmpo:Comp 13 (8.0E7 M atm-1)
     &   4.9E-2, !ddw-mpmpo:Comp 14(1.7E4 M atm-1;1.97E4*101.325 M/atm)
     &   4.9E-2, !ddw-mpmpo:Comp 15(1.7E4 M atm-1;1.97E4*101.325 M/atm)
     &   4.9E-2, !ddw-mpmpo:Comp 16(1.7E4 M atm-1;1.97E4*101.325 M/atm)
     &   2.44    !ddw-mpmpo:Comp 17(1.0E8 M atm-1)
     &   /)  ! ddw-mpmpo: the values for the new surrogates(from GEOS-Chem)

!Total number of (main+sub-components (ions)) per main components
!NK=3 means main component can dissociate twice (acid=H2A)
!NK=2 means main component can dissociate once (acid=HA)
!NK=1 means main component can not dissociate to ions
       INTEGER, PARAMETER, DIMENSION(NBSP)  :: NK=(/ 2, 2, 1, 1,
     &         1, 3, 3, 1, 2, 1, 1,
     &         1, 1, 1, 1, 1, 1 /)  !ddw-mpmpo: added for the new surrogates

!Are we running box model?
       LOGICAL                               :: LBOX=.FALSE.

!Do we want a lot of prints           
       LOGICAL                               :: LPRINT=.FALSE.
       END MODULE MODD_GLO
