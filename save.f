      subroutine print_header(nlev,nspec,z,levcpy,levhtr,ichem)
      USE parameters, ONLY: nrect_j
      USE cacm3_Parameters, ONLY: nreact,nvar
c 
c Prints header row(s) to "CACHE" output files for 
c meteorology, canopy structure, and soil data.

      implicit none !shc
      INTEGER             LEVCPY,      LEVHTR,      NLEV,        NSPEC !shc
      dimension z(nlev)
      double precision z
      integer i,ichem

cxc      if (ichem.eq.1) then	!ka - RACM chemistry scheme selected
c
c Set up output files for RACM species concentrations
c file names are fort.1xx where xx is species #
c
      write(318,'(x,"Time",10x,"Level",2x,411(i10,3x))') (i, i=1,nspec)
      write(319,'(x,"Time",10x,"Level",2x,51(i10,3x))')
     & (i, i=1,nrect_j)               ! photolysis reaction rates
      write(320,'(x,"Time",10x,"Level",2x,756(i10,3x))') (i, i=1,nreact)
      write(321,'(x,"Time",10x,"Level",2x,411(i10,3x))') (i, i=1,nspec)
      write(322,'(x,"Time",10x,"Level",2x,407(i10,3x))') (i, i=1,nvar)
c
c And for RACM source and sink rates
c file names are fort.7xx where xx is species # 
c
cx      do i = 1, nspec  
cx         write(i+700, 500) z
cx      enddo
c 
cxc      endif
c
c Write headers to output files for CACHE data
c First cache_met.out
c  
      write(41,971) 'Time','Level','PAR','Tair','U','Kh',
     &   'HeatFlux','Patm','rhoair','RH','Theta'
      write(41,971) 's','m','W m-2','K','m s-1','m2 s-1',
     &   'K m s-1','Pa','kg m-3','%','K'
c
c then cache_misc.out
c 
      write(42,972) 'Time','Radiation','(a)at_canopy_top',
     &    '(b)at_soil_surface'
      write(42,973) 'Time','(a)Solar','(a)Thermal','(a)Net',
     &    '(b)Solar','(b)Thermal','(b)Net','Zenith'
      write(42,973) 's','Wm-2','Wm-2','Wm-2','Wm-2','Wm-2','Wm-2','-'
c
c then cache_cpy.out
c 
      write(51,975) 'Time','Level','Tleaf1','Tleaf2','Tleaf3','Tleaf4',
     &'Tleaf5','Tleaf6','Tleaf7','Tleaf8','Tleaf9','Tleaf10','Rstom1',
     &'Rstom2','Rstom3','Rstom4','Rstom5','Rstom6','Rstom7','Rstom8',
     &'Rstom9','Rstom10','frarea1','frarea2','frarea3','frarea4',
     &'frarea5','frarea6','frarea7','frarea8','frarea9','frarea10',
     &'PAR Extinction','Sens. heat flux','Lat. heat flux'
      write(51,975) 's','m','K','K','K','K','K','K','K','K','K','K',
     &'umol-1 m2 s','umol-1 m2 s','umol-1 m2 s','umol-1 m2 s',
     &'umol-1 m2 s','umol-1 m2 s','umol-1 m2 s','umol-1 m2 s',
     &'umol-1 m2 s','umol-1 m2 s','-','-','-','-','-','-','-','-','-',
     &'-','-','K m s-1','K m s-1'

  500 format(1x, 60(:f12.2))
  971 format(x,a10,12(x,a10))
  972 format(x,2(x,a10),2(x,a20))
  973 format(x,a10,7(x,a12))
  975 format(2(x,a10),30(x,a10),3(x,a15))

      end

c ###########################################################
c Write out data directly to output files set up above in 
c print_header subroutine at time intervals set in main.f
      subroutine print_data(nlev,nspec,tstamp,ichem,averad,T,RelH, 
     &    wind,P,theta,rhoair,akh,vcp,photout,levcpy,z,coszen, 
     &    source,ustara,levhtr,rnlam,rnet,radabv,frarea,templf,rst
     &    , depn,chemr8_updated,rrat_updated)
c ###########################################################
      USE parameters, ONLY: nrect_j
      USE cacm3_Parameters, ONLY: nreact, nvar
      USE cacm3_Precision, ONLY: dp
      implicit none !shc
      !shc added the following declarations to allow "implicit none"
      INTEGER             NSPEC,       LEV
      INTEGER             I,           K
      DOUBLE PRECISION    DSTNET,            DSTRAD
      Double PRECISION    TSOIL
      DOUBLE PRECISION    POTVIS
      REAL(KIND=dp), DIMENSION(nlev,NREACT):: rrat_updated 
      REAL(KIND=dp), DIMENSION(nlev,NVAR):: chemr8_updated 
      !shc end adding new declarations

!ka      dimension T(nlev),RelH(nlev),wind(nlev),P(nlev),theta(nlev),	!cxio
!ka     &    rhoair(nlev),akh(nlev),vcp(nlev,nspec),photout(23,nlev),	!cxio
!ka     &    z(nlev),source(nlev,nspec),QHflx(nlev),akm(nlev),		!cxio
!ka     &    radabv(3),rnlam(3,99),rnet(99),rst(10,99)			!cxio
!ka      double precision tstamp,averad,T,RelH,wind,P,theta,rhoair,	!cxio
!ka     &    akh,akm,vcp,photout,z,source,ustara,QHflx,			!cxio
!ka     &    radabv,rnlam,rnet,coszen,rst				!cxio

      dimension T(nlev),RelH(nlev),wind(nlev),P(nlev),theta(nlev),	!cxio
     &    rhoair(nlev),akh(nlev),vcp(nlev,nspec),photout(nrect_j,nlev),	!cxio
     &    z(nlev),source(nlev,nspec),QHflx(nlev),radabv(3),		!cxio
     &    rnlam(3,99),rnet(99),frarea(10,99),templf(10,99),rst(10,99)	!cxio
      double precision tstamp,averad,T,RelH,wind,P,theta,rhoair,	!cxio
     &    akh,akm,vcp,photout,z,source,ustara,QHflx,			!cxio
     &    radabv,rnlam,rnet,coszen,frarea,templf,rst			!cxio
      integer nlev,levcpy,levhtr,j,ichem

      dimension nair(nlev)
      double precision nair,mslP,Av,RD,RUC
      double precision AKHtow,QHtow,dthdz,dakhdz

      common /rad7/ facte(99),extf(99),potvis
      common/cpysrc/ evaps(99),heats(99)
      double precision depn(nlev,nspec)
!cxio Next 4 lines to be removed for altered inputn
!cxio      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)	!cxio
!cxio     &,templf(10,99),tsoil						!cxio

!cxio      double precision extf,facte,frarea,templf			!cxio
!cxio      double precision heats,evaps					!cxio

      double precision extf,facte
      double precision heats,evaps
      integer jtot1
      jtot1=levcpy-levhtr+1
      do lev = 1, nlev
        write(318,991) tstamp,z(lev),(vcp(lev, i), i=1,nspec)
        write(319,993) tstamp,z(lev),(photout(i, lev), i=1,nrect_j)
        write(320,992) tstamp,z(lev),(rrat_updated(lev,i), i=1,nreact)
        write(321,991) tstamp,z(lev),(depn(lev,i), i=1,nspec)
        write(322,994) tstamp,z(lev),(chemr8_updated(lev,i), i=1,nvar)
      enddo
991   format(f10.1,x,f10.4,2x,411(E13.4E3)) 
992   format(f10.1,x,f10.4,2x,756(E13.4E3)) 
993   format(f10.1,x,f10.4,2x,51(E13.4E3)) 
994   format(f10.1,x,f10.4,2x,407(E13.4E3)) 
c
c Write out RACM source/sink rates (mol mol-1 s-1)
c 
cx      do i = 1, nspec 
cx         write (i+700, 510) tstamp, source(:, i)
cx      enddo
cxc      endif
c 
c Calculate and output kinematic heat flux at all levels
c 
      do k=2,nlev-1                                                    
        dthdz = (theta(k+1)-theta(k-1))/(z(k+1)-z(k-1))            
        QHflx(k) = -akh(k)*dthdz ! kinematic heat flux (K m s-1)          
      enddo                                                            
      QHflx(1)=0.                                                 
      QHflx(nlev)=0. 
c
c Or calculate kinematic heat flux at tower-top measurement 
c height only - levels are site-dependent.
c 
      dakhdz = (akh(20)-akh(19))/(z(20)-z(19))
      AKHtow = akh(19)+(dakhdz*(29.-z(19)))
      dthdz = (theta(20)-theta(19))/(Z(20)-Z(19))
      QHtow = -AKHtow*dthdz
c
c     Calculate new number density of air
      mslP = 101325.0
      RUC = 82.05746
      Av = 6.0221E+23
      RD = RUC/Av
      do k = 1,nlev
         nair(k) = P(k)/mslP/(RD*T(k))
      enddo

       do j=1,nlev     !cxio
         write(41,981) tstamp,z(j),radabv(1)*extf(j),T(j),wind(j), !cxio
     &     akh(j),QHflx(j),P(j),rhoair(j),RelH(j),theta(j)         !cxio
       enddo           !cxio

      write(42,982)    !cxio
     &  tstamp,-(rnlam(1,jtot1)+rnlam(2,jtot1)),!cxio
     &   -rnlam(3,jtot1),-rnet(jtot1),          !cxio
     &   -(rnlam(1,1)+rnlam(2,1)),-rnlam(3,1),-rnet(1),coszen !cxio

c 
c Canopy variables only valid in crown space.
c 
      do j=2,jtot1     !cxio
       write(51,985) tstamp,z(j+levhtr-1),(templf(i,j),i=1,10),   !cxio
     &   (rst(i,j),i=1,10),(frarea(i,j),i=1,10),                  !cxio
     &   extf(j+levhtr-1),heats(j+levhtr-1),evaps(j+levhtr-1)     !cxio
      enddo            !cxio

!comment out the following statements for new i/o
!cxio      do j=1,nlev			!cxio
!cxio        write(41,981) tstamp,z(j),radabv(1)*extf(j),T(j),wind(j),	!cxio
!cxio     &     akh(j),akm(j),QHflx(j),P(j),rhoair(j),RelH(j),theta(j)	!cxio
!cxio      enddo			!cxio

!cxio      write(42,982) 		!cxio
!cxio     &  tstamp,-(rnlam(1,levcpy)+rnlam(2,levcpy)),			!cxio
!cxio     &   -rnlam(3,levcpy),-rnet(levcpy),radabv(2),			!cxio
!cxio     &   -(rnlam(1,1)+rnlam(2,1)),-rnlam(3,1),-rnet(1),coszen	!cxio

!      write(42,*) tstamp,coszen,averad,radabv(2),-rnlam(2,levcpy),
!     & -rnet(levcpy)

!cxio      do j=1,levcpy		!cxio
!cxio       write(51,985) tstamp,z(j),(templf(i,j),i=1,10),!(rst(i,j),i=1,10),	!cxio
!cxio     &   (frarea(i,j),i=1,10),extf(j),heats(j),evaps(j)			!cxio
!cxio      enddo			!cxio

  510 format(1x,f12.0,60(:e12.4))		!cxio
  981 format(x,f10.1,x,f10.4,x,9e12.4)		!cxio
  982 format(x,f10.1,x,7e12.4)			!cxio
c  982 format(x,f10.1,x,8e12.4)		!cxio
  985 format(x,f10.1,x,f10.4,x,33e12.4)	!cxio
      end
