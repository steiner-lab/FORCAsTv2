!======================================================================================================================!
!                                                                                                                      !
!     Program:      FORCAST                                                                                            !
!                   Forest Canopy Atmosphere Transfer                                                                  !
!                   RCIM (the Reduced CalTech Isoprene mechanism) version                                              !
!                                                                                                                      !
!     Version:      2.0                                                                                                !
!                                                                                                                      !
!     Last Update:  Feb 2019                                                                                           !
!                                                                                                                      !
!     Contact:                                                                                                         !
!                                                                                                                      !
!                                                                                                                      !
!                                                                                                                      !
!                                                                                                                      !
!                                                                                                                      !
!======================================================================================================================!
!**********************************************************************************************************************!
!************************************************ LEGAL NOTICE ********************************************************!
!**********************************************************************************************************************!
!                                                                                                                      !
!      Originally was Canopy-chemistry model CACHE (Renate Forkel, IMK-IFU)                                            !
!      Major parts of the model (radiation canopy structure) are from CUPID by John Norman                             !
!      RCIM chemical mechanism by Wennberg 2018 et al. (Chemical Reviews)                                              !
!                                                                                                                      !
!======================================================================================================================!
!**********************************************************************************************************************!
!******************************************** STILL TO BE IMPROVED ****************************************************!
!**********************************************************************************************************************!
!                                                                                                                      !
!      Water extraction by roots                                                                                       !
!      Photolysis                                                                                                      !
!                                                                                                                      !
!**********************************************************************************************************************!
!======================================================================================================================!
      PROGRAM FORCAST
       
      use parameters
      use cacm3_parameters
      use module_interf
      USE cacm3_Precision, ONLY:dp 

      implicit none 

!      INTEGER             ndspec !ka - added deposition of meoh and etoh; added ndspec
      INTEGER             iph
      INTEGER             MAXEMI
      INTEGER             MAXRECT
      INTEGER             MAXPHR
      INTEGER             MAXPHPT
      real hours ! simulation length in hours      
      real pp, tt, rho_phy,  moist
      DOUBLE PRECISION    TEMPS, ETA
      DOUBLE PRECISION    A,           AIRMAS,      AJNN
      DOUBLE PRECISION    AKH,   AKHK,        AKHTOW
      DOUBLE PRECISION    ALPHAB,      AMX
      DOUBLE PRECISION    AXLOG,       B,           CKHCPY,      CLUMP
      DOUBLE PRECISION    CSIGW,       CSIGW0,      CSIGW1
      DOUBLE PRECISION    CSIGWCPY,    CSIGWTOW,    CUSTAR,      CUSTAR0
      DOUBLE PRECISION    CUSTAR1,     CUSTARCPY,   CUSTARTOW,   DAKH
      DOUBLE PRECISION    DCKHDZ,      DSIGWDZ,     DTL
      DOUBLE PRECISION    DTL2,        DTLAST,      DTLF,        DUM1
      DOUBLE PRECISION    DUM2,        DUSTARDZ,    EMISS,       EVAL
      DOUBLE PRECISION    P1,          P1A,         P1B
      DOUBLE PRECISION    PHOT1, PHOTOUT
      DOUBLE PRECISION    POTBM1,      POTBM2,      POTDIF,      POTNIR
      DOUBLE PRECISION    PSAT,        PSIDIF,      PSIMN
      DOUBLE PRECISION    PSISV1,      PSITP1,      PTIME,       QVAP1
      DOUBLE PRECISION    QVAP2,       RATDSAV,     RATNSAV
      DOUBLE PRECISION    RFACTORCPY,  RFACTORTOW,  RLAT,        T0
      DOUBLE PRECISION    T1,          T30,         TA,          TAU
      DOUBLE PRECISION    TEMLF1,            TIMIN,       TL
      DOUBLE PRECISION    TLFSUM,      TRSAV1 !TSTART
      DOUBLE PRECISION    USTARA,      UU
      DOUBLE PRECISION    USTARB  
      DOUBLE PRECISION    VCP,          VCP1,        VG
      DOUBLE PRECISION    VG0,         VP,    VPA,         WATABS
      DOUBLE PRECISION    WIND,  XLA,         XROOT
      INTEGER             I,           IDAY,        IGO,         IHR
      INTEGER             IPASS,       ITERH,       ITIME
      INTEGER             IZO,         IZU,         J,           JSPEC
      INTEGER             JTIME,       K,           KSTRT,       LEV
      INTEGER             NOITER,  NTIMES
      DOUBLE PRECISION    COSDEC,      COSLAT,      DECL,        DECMAX
      DOUBLE PRECISION    DLONG,       EQTM,        SINDEC,      SINLAT
      DOUBLE PRECISION    TANLAT
      DOUBLE PRECISION    P,     RELH,  RHOAIR
      DOUBLE PRECISION    T,     THETA
      DOUBLE PRECISION    ALBDRY,      ALBWET,      CP3,         ETAMAX
      DOUBLE PRECISION    ETAMIN,      FELDKA,      POR,         RHO3
      DOUBLE PRECISION    TREFK,       TREFP,       XKR
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    EAVE,    ETOTWT,      EVAPG
      DOUBLE PRECISION    HAVE,    HEATG,   HTOT,        RPLNT
      DOUBLE PRECISION    WATERG
      DOUBLE PRECISION    CONTOT,  CPESTR,      CPHSTR,      ECPYS
      DOUBLE PRECISION    ETOTW,   EVSMIC,      EVTOT,       HCPYS
      DOUBLE PRECISION    HSOIL,   SCOND,             WCPYS
      INTEGER             IHRWET
      DOUBLE PRECISION    EVAPS,   HEATS
      DOUBLE PRECISION    SUNAZM,      VIEWAZ,  VIEWZN
      DOUBLE PRECISION    XINTV
      INTEGER             NOAZMV,      NOZENV
      DOUBLE PRECISION    ADVEC, ADVEMS
      DOUBLE PRECISION    BVOC, BVOCEMS, BEMIS
      DOUBLE PRECISION    SOVOCE  
      DOUBLE PRECISION    FAPI,        FAPI2,       FISO,        FISO2
      DOUBLE PRECISION    FLIM,        FLIM2,       FMACR,       FMACR2
      DOUBLE PRECISION    FACET,       FACET2,      FACALD,      FACALD2
      DOUBLE PRECISION    SUMSAPI,     SUMSISO,     SUMSLIM
      DOUBLE PRECISION    AXA,  BXB,  CXC
      DOUBLE PRECISION    D1ETA,             D1T
      DOUBLE PRECISION    D2ETA,             D2T
      DOUBLE PRECISION    DATM, DXD
      DOUBLE PRECISION    H,    P21,  P21T
      DOUBLE PRECISION    PSI,  PSIETA
      DOUBLE PRECISION    PSIT, PSITET
      DOUBLE PRECISION    RHO0, RHO1, RHO22
      DOUBLE PRECISION    XK,   XL11
      DOUBLE PRECISION    XL210,             XL22
      DOUBLE PRECISION    XLAM
      REAL(KIND=dp)    :: Z, ZF, DZ, DZF
      INTEGER             IITIME,      LLEV
      INTEGER             INRAD,       IWPM2
      DOUBLE PRECISION    AFL1,        AFLH,        BFL1,        BFL2
      DOUBLE PRECISION    BFLH,        FNET,        GDOWN,       GFNET
      DOUBLE PRECISION    PATM,        QGRND,       QSATG,       TAUSUM
      DOUBLE PRECISION    TSUM,        TSURF,       VSUM
      INTEGER             IWRITE,             KMAX
      DOUBLE PRECISION    CLAI,    CT,      DF
      DOUBLE PRECISION    DISTLS,  FR,      TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    PARABV
      DOUBLE PRECISION    DEPN,       DEPR
      DOUBLE PRECISION    DEPRAT,     SOURCE
      DOUBLE PRECISION    ALEAF,    EMIS,        EMISOL
      DOUBLE PRECISION    EXPDIF,  RLAYR, RLEAF
      DOUBLE PRECISION    RSOIL,    TLAYR, TLEAF
      DOUBLE PRECISION    DSTNET,            DSTRAD
      DOUBLE PRECISION    FRAREA,            TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    COSZEN,      FBEAM,    RADTOP,   RATIO
      DOUBLE PRECISION    RATIOD,      RATION,      ZENANG
      DOUBLE PRECISION    BMFLX, D,     RNDIV
      DOUBLE PRECISION    TSFC,        U
      DOUBLE PRECISION    RLFDIF,   RLFDIR,   TLFDIF
      DOUBLE PRECISION    TLFDIR
      DOUBLE PRECISION    FACTE,   POTVIS
      DOUBLE PRECISION    HPSI,        RHLEAF,  RSLEAF
      DOUBLE PRECISION    RSNOVP, RST
      DOUBLE PRECISION    ANSTOM,      PSI1,        PSI2,        RADN
      DOUBLE PRECISION    RCUT20,      RSEXP,       RSM,         RSMIN
      DOUBLE PRECISION    TRSMAX,      TRSMIN,      TRSOPT
      DOUBLE PRECISION    AROOT,       CPYTR,       FROOT,   PSISUM
      DOUBLE PRECISION    PSITOP,      PSIXY,       RESROT,  ROOTSM
      DOUBLE PRECISION    ROOTUP,  RROOT
      DOUBLE PRECISION    DELS,        SL,          ZS
      DOUBLE PRECISION    ALPHADAY,    ALPHAN,      VGDAY,       VGN
      DOUBLE PRECISION    ZRDAY,       ZRN
      REAL  TIMLOC, deltim, tstart
      REAL(KIND=dp), DIMENSION(NREACT):: RCONST 
      REAL(KIND=dp), DIMENSION(nlev,NREACT):: rrat_updated 
      REAL(KIND=dp), DIMENSION(nlev,NVAR):: chemr8_updated 
!**********************************************************************
!        MPMPO (the Model to Predict the Multiphase Partitioning of
!        Organics)
!**********************************************************************
      DOUBLE PRECISION boxcacm(nspec)     
      LOGICAL is_call_mpmpo
      ! 53 photolysis reactions in WRF/KPP coupler
      REAL ::  
     &  ph_o31d,     ph_o33p,     ph_no2,      ph_no3o2,     ph_no3o,   
     &  ph_hno2,     ph_hno3,     ph_hno4,     ph_h2o2,      ph_ch2or,  
     &  ph_ch2om,    ph_ch3cho,   ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, 
     &  ph_ch3cocho, ph_hcochest, ph_ch3o2h,   ph_ch3coo2h,  ph_ch3ono2,
     &  ph_hcochob,  ph_macr,     ph_hobr,     ph_br2,       ph_bro,    
     &  ph_n2o5,     ph_o2,       ph_pan,      ph_acet,      ph_mglo,   
     &  ph_hno4_2,   ph_c6h5cho,  ph_unald,    ph_n2o,       ph_pooh,   
     &  ph_mpan,     ph_mvk,      ph_etooh,    ph_prooh,     ph_onitr,  
     &  ph_acetol,   ph_glyald,   ph_hyac,     ph_mek,       ph_open,   
     &  ph_gly,      ph_acetp,    ph_xooh,     ph_isooh,     ph_alkooh, 
     &  ph_mekooh,   ph_tolooh,   ph_terpooh, 
     &  ph_ooh1, ph_ooh2, ph_rpr2301, ph_ald2, ph_ooh6101,    ! new photolysis
     &  ph_ooh6102, ph_hac, ph_ap6101, ph_ap6102, ph_nald, 
     &  ph_hpald, ph_pacald, ph_rp6301, ph_ooh6301, ph_rp7102, 
     &  ph_pina, ph_nrpa, ph_edlm, ph_lmkt, ph_rpr8401, 
     &  ph_rpr8501, ph_ooh8601, ph_ooh8602,  
     &  ph_1co2ooh, ph_3ooh4co,  ph_3ooh4oh, ph_1oh2ooh, ph_mvk3oh4co, 
     &  ph_mvk3co4oh, ph_1co4ooh, ph_1ooh4co, ph_mcrenol, ph_mvkenol,
     &  ph_mvk3ooh4co, ph_1n4co, ph_1co4n, ph_3co4n
      
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      INTEGER             ISORT
      integer ncacm !shc 2015-04
      parameter (ncacm=351) !shc ncacm === nspec in cacm_Parameters.f90; updated for cacm2.0 ka
      double precision   vset !shc settling velocity for aerosol
      double precision   alwc !shc aerosol liquid water content

!      parameter (nlev=40,nsoil=15)    ! number of levels  
!      parameter (nspec=397, nrect=884)!ka - added deposition of meoh and etoh; added ndspec
      parameter(maxemi=nspec) ! maximal number of emissions
      parameter(maxrect=500)
      parameter(maxphr=100, maxphpt=100)

      external f, jac, rate
      logical lprin

      dimension vcp(nlev,nspec), phot1(nrect_j,nlev),
     &vp(nlev),wind(nlev),akh(nlev),photout(nrect_j,nlev),ajnn(nlev)
      dimension temlf1(10,99)
      double precision photoj(maxphr)
      real vc(nspec)
      common /info/ iitime,llev
      common /treespec/ isort
      common /inputr/ inrad,iwpm2
      double precision emivalue(maxemi)
      common /emissions/ emivalue
      common/quellen/source(nlev,nspec),depn(nlev,nspec),
     & deprat(nlev,nspec),depr(nspec)  !ddw deleted ndspec
       !ka - added deposition of meoh and etoh; altered nspec to ndspec
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
      common/height/z(nlev),zf(nlev),dz(nlev),dzf(nlev)
!ka - added for sonics data; 4 lines inserted
      common/oheight/isonics,zobstow,zmodtow,levobstow
     &,zobscpy,zmodcpy,levobscpy,zmxd,levmxd     !ka - added for sonics data
      integer levobstow,levobscpy,levmxd,isonics !ka - added for sonics data
      double precision zobstow,zmodtow,zobscpy,zmodcpy,zmxd !ka - added for sonics data
      logical zflgt,zflgc,zflgm                  !ka - added for sonics data
!ka - added for sonics data; end of inserted lines
      common/atmo1/ t(nlev),theta(nlev),p(nlev),rhoair(nlev)
     & ,relh(nlev)
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/cpy2/hsoil,hcpys,evtot,etotw(99),contot(99),scond(10,99)
     &,ecpys,cphstr,cpestr,wcpys,evsmic,ihrwet(99)
      common/misc1/iwrite(9,99),kmax
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot

      common/astron/eqtm,decl,sindec,cosdec,decmax,sinlat,coslat,
     & tanlat,dlong
      common/deg/sunazm,viewzn(99),viewaz(99),nozenv,noazmv
     & ,xintv(20,90)
      common /rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3)
     &,expdif(99),rlayr(3,30),tlayr(3,30)
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common /rad3/radtop(3),fbeam(3),coszen,zenang
     & ,ratiod,ration,ratio
      common/rad4/d(3,99),u(3,99),bmflx(3,99),rnet(99),rndiv(99),
     &tsfc,rnlam(3,99)
      common/rad6/rlfdif(3),tlfdif(3),rlfdir(3),tlfdir(3)
      common /rad7/ facte(99),extf(99),potvis
      common/resis1/rhleaf(99),rsleaf(10,99),rsnovp(10,99),
     & rst(10,99),hpsi
      common/resis2/radn,anstom,rcut20,rsmin,trsopt,trsmin,trsmax,rsm
     &,rsexp,psi1,psi2
      common/root1/froot(99),resrot(99),rootsm,psixy,psitop,rroot
     &,rootup(99),cpytr,psisum,aroot
      common/cpy1/etotwt,htot,rplnt,evapg(99),heatg(99),eave(99)
     &,have(99),waterg(99)
      common/cpysrc/ evaps(99),heats(99)
      common/soil1/zs(nsoil),sl,dels
      common/feld/ eta(nsoil),psi(nsoil),psit(nsoil), temps(nsoil),
     & xk(nsoil),psieta(nsoil),
     & psitet(nsoil),xlam(nsoil),axa(nsoil),bxb(nsoil),
     & cxc(nsoil),dxd(nsoil),d1t(nsoil),
     & d2t(nsoil),d1eta(nsoil),d2eta(nsoil),
     & xl11(nsoil),xl22(nsoil),rho0(nsoil),
     & rho1(nsoil),rho22(nsoil),xl210(nsoil),
     & p21(nsoil),p21t(nsoil),h(nsoil),datm(nsoil)
      common/bodk/ etamin,etamax,xkr,por,trefp,trefk,cp3,rho3,
     &feldka,albdry,albwet
      common/intf/ qgrnd,qsatg,tausum,tsum,vsum,tsurf,patm,
     &fnet,gdown,gfnet,afl1,aflh,bfl1,bfl2,bflh
      data akh/nlev*0./
      common/stratpar/ vgday,vgn,zrday,zrn,alphaday,alphan,stable
      logical stable
! ka - added beta (bexp array)
      common/efact/ EFno,EFsyn(14),EFpl(14),sovoce(9),bvoc(nlev,14),
     & bemis(nlev,13),bvocems(14),advec(nlev,10),advems(10),bexp(13)
      double precision EFno,EFsyn,EFpl,bexp
! ka - added beta;  end of added section

      common /par/parabv
      dimension akm(nlev)
      double precision akm
      double precision tstamp,extf

      common/etc/fiso,fapi,flim,fmacr,sumsiso,sumsapi,sumslim,fapi2
     &,fiso2,flim2,fmacr2

      logical l30
      integer ivar,ichem,impmpo !ka - added MPMPO switch to inputn
      double precision UV,avekh
      real  dpt
      dimension inconc(nlev,nspec),vcnc(nlev,nmpmpo),flxkh(nlev)  
! ddw - change vcnc to aerosol array at all heights
! ka - added to calculate fluxes at all heights
      double precision inconc,vcnc,flxkh  !ka - added to calculate fluxes at all heights
      dimension gascno(ncacm),gascni(ncacm)
      dimension parcno(nmpmpo),parcni(nmpmpo)
      double precision gascno,gascni,parcno,parcni
      double precision rnlam,rnet,radabv

      integer ii

      do k=1,3
       fbeam(k)=0.
      enddo

      do ivar=1,ncacm
       gascno(ivar)=0.0
       gascni(ivar)=0.0
      enddo

      do ivar=1,nmpmpo
       parcno(ivar)=0.0
       parcni(ivar)=0.0
      enddo

      do k=1,nlev
       do ivar=1,nmpmpo
        vcnc(k,ivar)=0.0
       enddo
       flxkh=0.0  !ka - added to calculate fluxes at all heights
      enddo

!**********************************************************************************************************************!
!                           
! Simulation length and number of time steps
!                           
!**********************************************************************************************************************!
      timloc=0.                          ! Local time in hours 
      tstart= timloc*3600.               ! Start time in seconds

      hours = 48.0     
      deltim= 1.*60.                       ! Time step in seconds
      ntimes= 60.*60./deltim * hours     ! Number of time steps (48 hours)

      ptime = timloc+dpt                 ! Next print time step in hours
      dpt= 30./60.                       ! print frequency in hours

!**********************************************************************************************************************!
!
! Start the run-time log
!
!**********************************************************************************************************************!
      write(06,*) '48 hour simulation of canopy exchange processes'
      write(06,*) 'using the FORCAsT 1-D model'

!**********************************************************************************************************************!
!
! Initialize profiles of T, RH, and chemical species
! ka - added MPMPO switch to inputn
!
!**********************************************************************************************************************!
      call init(vcp,phot1,emiss,wind,rlat,clump,ustara,dpt,ichem,impmpo)
      do i=1, nlev
         akh(i)=0.1
      enddo
      ratdsav=ratiod
      ratnsav=ration
      iday=1

!**********************************************************************************************************************!
!
! Print out header rows of output files
!
!**********************************************************************************************************************!
      ! Initialise the time stamp (for printing)
      tstamp = 0.0
      t30=0.                      ! 30-min reads

      call print_header(nlev,nspec,z,levcpy,levhtr,ichem)
       
!**********************************************************************************************************************!
!
! Now create an array to pass initial conditions to 
! chemistry box model - will be mapped to appropriate
! species in the chemistry scheme
!
!**********************************************************************************************************************!
      do k=1,nlev
       do i=1,nspec
        inconc(k,i)=vcp(k,i)
       enddo
      enddo
!**********************************************************************************************************************!
!
! Print initial data to output files at first time step
!
!**********************************************************************************************************************!
      write(06,'(/," Model time: ",a5,"s; ",a5)') '0','0.00'
      write(06,'(" Simulation is ",f6.1," % complete")')100.0*0.0/ntimes
     
      averad = 0.0 
      call print_data(nlev,nspec,tstamp,ichem,averad,T,RelH,    !cxio
     &    wind,P,theta,rhoair,akh,vcp,photout,levcpy,z,coszen, !cxio
     &    source,ustara,levhtr,rnlam,rnet,radtop,frarea,templf,rst
     &    ,depn,chemr8_updated,rrat_updated)
!
! Read PAR from fort.500 input file
!
!**********************************************************************************************************************!
      if(inrad.eq.1) then
         dtl2=0.
         open(unit=500,file='./data/par.dat', status='unknown')
         read (500, *)   ! read header
         timin=0.
      endif

!**********************************************************************************************************************!
!
! Start the time loop here. 
! Radiation profiles, emissions, deposition, vertical exchange, turbulent transport, 
! advection and chemistry are calculated for each time step
!
!**********************************************************************************************************************!
      do 100 itime=1,ntimes

       if(timloc.ge.24.) then  ! New day?
       timloc=timloc-24.
       tstart=timloc*3600.
       jday=jday+1
       iday=iday+1
       icumdy=icumdy+1
       endif
       
       ! print 
       lprin=.false.
       if(timloc.ge.ptime-deltim/3600.) lprin=.true.
       if(itime.eq.ntimes) lprin=.true.
       ! print by ak
       l30=.false.
       if(tstamp.eq.t30) l30=.true.  

!**********************************************************************************************************************!
!
! Consider cloud effects
! 
!**********************************************************************************************************************!
      ratiod=ratdsav
      ration=ratnsav

      call declin

      call zenith
!**********************************************************************************************************************!
! Read incoming radiation from external input file for calulating ratiod (cloud effect)
!**********************************************************************************************************************!
      if(inrad.eq.1) then
      if(l30) read (500, *) timin,averad !Read every half-hour

      potvis=0.
      potnir=0.
      if (coszen.gt.0.01) then
      !  correct for curvature of atmos in airmas
      airmas=(sqrt(coszen**2+.0025)-coszen)/.00125
      !  correct for refraction(good to 89.5 deg.)
      airmas=airmas-2.8/(90.-zenang/pid180)**2
      potbm1=600.*exp(-.160*airmas)
      potvis=(potbm1+(600.-potbm1)*.4)*coszen
      potdif=(600.-potbm1)*.4*coszen
      uu=1/coszen
      axlog=dlog10(uu)
      a=10**(-1.195+.4459*axlog-.0345*axlog*axlog)
      watabs=1320.*a
      potbm2=720.*exp(-.05*airmas)-watabs
      if(potbm2.lt.0.)potbm2=0.
      eval=(720.-potbm2-watabs)*.54*coszen
      potnir=eval+potbm2*coszen
        if(iwpm2.eq.1) then
        ratiod=averad/(potvis+potnir) !  averad includes total solar spectrum (SWdown) in W/m**2
        else
        ratiod=(averad/4.6)/potvis    !  averad is PAR in micromoles per (m2 seconds)
        endif
      endif
       else
!  Modify if you want to introduce cloud cover without specifying PAR on extra file
!  Example:
      if(iday.eq.2) then
       if(timloc.ge.8..and.timloc.lt.12.5)  then
       if(timloc.ge.11..and.timloc.lt.12.)  ratiod=0.1
       ratiod=0.6
       ration=0.6
       endif
       if(timloc.gt.19.) then
       ratiod=0.4
       endif
      endif
      endif

!**********************************************************************************************************************!
c Radiation for all layers and 10 angle classes
!**********************************************************************************************************************!
      call declin
      
      call zenith
      
      call dstlit(zenang,sunazm,kmax,rlfdir,tlfdir)
      
      if (coszen.gt.0.01) call radin4

      call difint(clump,sunazm)
      
      do k=1,nlev
       vp(k)=0.01*vcp(k,lh2o)*p(k)
       psat=610.7*exp(17.1536*(t(k)-273.15)/(t(k)-38.33))
       p1=vcp(k,lh2o)*p(k)
       relh(k)=p1/psat
      enddo
      

      ta=t(30)
      vpa=vp(30) 

      call skyir(ta,vpa) 

!**********************************************************************************************************************!
!
! Iteration over stress 
! Iterate leaf temperature 
!
!**********************************************************************************************************************!
      ipass=0
      noiter=0
      iterh=0

      kstrt=3
      if (coszen.gt.0.01) kstrt=1 ! always longwave, shortwave only if sun is up
 1060 if(ipass.eq.1) kstrt=3
      
      ! No idea what this is
      call radiat(kstrt,coszen,radtop,fbeam,zenang,sunazm,clump,lprin)


      ! Calculate stomatal resistance
      call stoma(lprin,t,relh)
      
      ! Leaf surface temperature
      call lfebal(t,vp,wind,lprin)

      igo=0
      tlfsum=0.
c----------------------------------
      do 1200 i=1,itotp1
      do 1200 j=2,jtot
      if(abs(templf(i,j)-temlf1(i,j)).gt.0.2) igo=1
      tlfsum=tlfsum+abs(templf(i,j)-temlf1(i,j))
 1200 continue
c----------------------------------
      dtlf=tlfsum/(itotp1*(jtot-1))
 1218 format(' mean abs tleaf_new-tleaf= ',f6.3,' noiter=',i3)
      if(noiter.gt.50) goto 1219
      if(igo.eq.1) then      
      ipass=1
      noiter=noiter+1
      do 1100 i=1,itotp1
      do 1100 j=2,jtot
      ! tried weighting .5 and .5 but convergence is only half as fast.
      templf(i,j)=.9*templf(i,j)+.1*temlf1(i,j)
 1100 temlf1(i,j)=templf(i,j)
      goto 1060
      endif

      ! Water extraction by roots
 1219 call rootex
      
      call stress
      psidif=psitop-psitp1
      if(iterh.gt.50)write(06,1265)ihr,iterh,psitop/100.,cpytr
     &,psidif/100.
 1265 format(' ihr=',i3,' iterh=',i3,' psitop=',f6.2,' cpytr=',f8.2
     &,' psidif=',f6.2)


      if(iterh.gt.50)goto 1300
     
      if(psitp1.lt.psi1.and.iterh.eq.0) goto 1255
      
      if(psitop.ge.psi1.and.iterh.eq.0) goto 1300
 1255 psimn=(psi1+psi2)/2.
      if(iterh.gt.0) goto 1270
      ! psitp1 is previous psitop.   psisv1 is psitot before calc of cpytr
      ! from previous iteration. it goes with trsav1 for first coord. of
      ! slope estm of cpytr vs psitop relation.
      trsav1=cpytr
      psisv1=psitp1
      if(abs(psitop-psitp1).gt.0.) psitop=(psitop+psitp1)/2.
      if(abs(psitop-psi2).le.0.) psitop=psitop+10.
      if(psitop.lt.psi2) psitop=psimn-sqrt((psi2-psitop)/(-psitop))*
     &(psi1-psi2)
      psitp1=psitop
      iterh=iterh+1
      goto 1060
 1270 if(abs(psitop-psitp1).lt.10.)goto 1300
      amx=(trsav1-cpytr)/(psisv1-psitp1)
      b=trsav1-amx*psisv1
      xroot=0.66667*rroot+1./rootsm
      psitop=(-b*xroot+psisum)/(1.+amx*xroot)
      ! psisum in above stm was originally psisum/rootsm, chen, 02/19/90.
      if(iterh.gt.10) psitop=(psitop+psitp1)/2.
      if(psitop.lt.psi2) psitop=psimn-sqrt((psi2-psitop)/(-psitop))*
     &(psi1-psi2)/2.
      trsav1=cpytr
      psisv1=psitp1
      psitp1=psitop
      iterh=iterh+1
      goto 1060
 1300 continue
!**********************************************************************************************************************!
! end iteration over stress 
!**********************************************************************************************************************!


      do k=1,nlev
       evaps(k)=0.
       heats(k)=0.
      enddo
      
      do k=2,jtot
       evaps(levhtr-2+k)=evapg(k)
       heats(levhtr-2+k)=heatg(k)
      enddo
!**********************************************************************************************************************!
!
! Vertical exchange 
!
! Wind (fixed specified profile) 
! 
c Calculation of vertical exchange coefficients (KH and KM) - based on 
c parameterisations of Stroud and Wolfe. Assume linear interpolations
c between sonic heights within the canopy, and then smooth the profile
c to match the above canopy parameterisation of Forkel at 1km.
!
!**********************************************************************************************************************!
      ! Read wind data from external CABINEX data file
       if (itime.eq.1) vg0=wind(21) ! in case of error on first read
       if (l30.and.(itime.ne.1)) read(5034, *) timin, vg0
       if (vg0 .ne. -999.) vg = vg0
      
      ! Parameters for wind profiles
      if(rnet(jtot).gt.0.) then
c       vg = vgday
       alphab=alphaday
       zrough=zrday
      else  
c       vg = vgn  
       alphab=alphan
       zrough=zrn
      endif
      
      ! Read observed ustar (ustarb) to estimate rbl and Kh; added by ddw
      if(l30) read (5081, *) timin,CUSTAR1
      ustarb = custar1
      ! Wind profile

      call wprofil(vg,alphab,wind,ustarb,tstamp)

      ! Eddy diffusivity
      call atk(vg, rlat,z,wind,akh,ajnn,lprin,tstamp,ustarb)
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
      ! ka - added for sonics data to allow switching between sites       !
      ! start of inserted section                                         !
      ! Re-calculation of within canopy KH profile (Stroud et al. 2005)   ! 
      ! isonics is read from inputn canopy sonics are available (e.g. HF) !
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
      select case (isonics)           
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
      ! First calculation is for sonics at two heights (e.g. UMBS)        !
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
     
       case (2)                
      
      ! lev35 is 1043 m; 18:00 the heat flux at the top of the canopy becomes negative 
      levmxd = 29                     ! daytime; lev=29 is 215 m 
      if(timloc.ge.18.) levmxd = 35    ! nighttime

      ! Read USTAR and SIGW every 30 min for in-canopy (bottom) anemometer
      if (l30) read(5082, *) timin, CUSTAR0
      if (l30) read(5092, *) timin, CSIGW0
      if (CUSTAR0 .ne. -999.) CUSTARcpy = CUSTAR0
      if (CSIGW0 .ne. -999.) CSIGWcpy = CSIGW0
      ! solve for KH at this level by calculating TL
      TL = 0.3*hcpy/CUSTARcpy
      CKHcpy = CSIGWcpy**2.*TL
      ! linearly interpolate KH from trunk height to in-canopy observation height
      DCKHDZ = (CKHcpy-AKH(levhtr))/(zobscpy-Z(levhtr))
      do k=levhtr,levobscpy
         AKH(k) = AKH(levhtr)+DCKHDZ*(Z(k)-Z(levhtr))
      enddo
      ! save modeled KH at tower observation height for adjusting ABL to remove discontinuity
      AKHtow = AKH(levobstow)
      ! read USTAR and SIGW every 30 min for above canopy (tower) anemometer
      ! the following line was commented out by ddw the 5081 file was read above before atk
      ! if (l30) read (5081, *) timin, CUSTAR1 
      if(l30) read (5091, *) timin,CSIGW1
      if (CUSTAR1 .ne. -999.) CUSTARtow = CUSTAR1
      if (CSIGW1 .ne. -999.) CSIGWtow = CSIGW1
      ! linearly interpolate SIGW, as in Stroud et al. 2005, 
      ! do the same for USTAR (though USTAR at both heights are about
      ! equal and Stroud only linear interpolated SIGW) using the gradient
      ! between the measurements at the two heights
      DSIGWDZ = (CSIGWtow-CSIGWcpy)/(zobstow-zobscpy)
      DUSTARDZ = (CUSTARtow-CUSTARcpy)/(zobstow-zobscpy)
      do k=levobscpy+1,levobstow
        CSIGW = CSIGWcpy+DSIGWDZ*(Z(k)-zobscpy)
        CUSTAR = CUSTARcpy+DUSTARDZ*(Z(k)-zobscpy)
        TL = 0.3*hcpy/CUSTAR
        AKH(k) = CSIGW**2.*TL
      enddo
      ! And finally adjust ABL to remove discontinuity up to ~1 km
         DAKH = AKH(levobstow)-AKHtow
      !ka         DCKHDZ = (AKH(35)-AKH(levobstow))/(Z(35)-Z(levobstow))
      !ka         do k=levobstow+1,35
      DCKHDZ = (AKH(levmxd)-AKH(levobstow))/(Z(levmxd)-Z(levobstow)) !ka - 1km is level levmxd
      do k=levobstow+1,levmxd                                        !ka - 1km is level levmxd
         AKHK = AKH(levobstow)+DCKHDZ*(Z(k)-Z(levobstow))
         if (akh(k).lt.akhk) then
            AKH(k) = AKHK
         else
            AKH(k) = AKH(k)
         endif
      enddo

      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
      ! Second case is for sonics only available above the canopy
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
     
      case(1)    
      
      !c read USTAR and SIGW every 30 min for above canopy (tower) anemometer 
c      if(l30) read (5081, *) timin, CUSTAR1   !ddw custar1 was read before atk
      if(l30) read (5091, *) timin, CSIGW1
      if (CUSTAR1 .ne. -999.) CUSTARtow = CUSTAR1
      if (CSIGW1 .ne. -999.) CSIGWtow = CSIGW1
      ! solve for KH at this level by calculating TL    
      TL =  0.3*hcpy/CUSTARtow
      CKHcpy=CSIGWtow**2.*TL
      ! linearly interpolate KH from trunk height to the tower top
      DCKHDZ = (CKHcpy-AKH(levhtr))/(zobstow-Z(levhtr))     
      do k=levhtr,levobstow   
         AKH(k) = AKH(levhtr)+DCKHDZ*(Z(k)-Z(levhtr))
      enddo
      ! And finally adjust ABL to remove discontinuity up to ~1 km
      DAKH = AKH(levobstow)-AKHtow
      !ka         DCKHDZ = (AKH(35)-AKH(levobstow))/(Z(35)-Z(levobstow))
      !ka         do k=levobstow+1,35
      DCKHDZ = (AKH(levmxd)-AKH(levobstow))/(Z(levmxd)-Z(levobstow)) !ka - 1km is level levmxd
      do k=levobstow+1,levmxd                                        !ka - 1km is level levmxd
         AKHK = AKH(levobstow)+DCKHDZ*(Z(k)-Z(levobstow))
         if (akh(k).lt.akhk) then
            AKH(k) = AKHK
         else
            AKH(k) = AKH(k)
         endif
      enddo
      
      end select
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!
      ! End of inserted section
      !-------------------------------------------------------------------!
      !-------------------------------------------------------------------!

!**********************************************************************************************************************!
!
! Sources due to emissions
! "Steinbrecher" emissions - Guenther algorithms adapted to include
! pool and synthesis emissions for monoterpenes and other VOCs
!
! ka - isort no longer required as emission factors read from inputn
!**********************************************************************************************************************!
      call sourcest(vcp,lprin,deltim,tstamp)
!**********************************************************************************************************************!
! Sinks due to deposition
!**********************************************************************************************************************!
      call sinks(vcp,akh,lprin,tstamp,ustarb)
!**********************************************************************************************************************!
!
! Advection
! Modify this section to incorporate advection of energy (heat)
! and mass (concentrations). For most sites this will be wind 
! direction dependent, and will predominantly consist of key
! anthropogenic pollutants (NOx, CO, SO2, O3, VOCs)
! This is the default CACHE advection scheme - will ultimately
! remove and replace with SC's scheme or version thereof
!
!**********************************************************************************************************************!
       call advect(l30,itime,vg,vcp,tstamp)
      
!**********************************************************************************************************************!
!
! The R factor correction by Makar et al. 1999
! Turbulent transport of chemical species & water vapor modified
! after Stroud et al canopy model scheme. Some of scheme relies
! on within canopy sonics data (currently not available at HF).
! These sections of scheme commented out but will be put back in
! when/if data becomes available.
!**********************************************************************************************************************!
      
      !-----------------------------------------------------! 
      !-----------------------------------------------------! 
      ! ka - added for sonics data to allow for switching between sites; 
      ! start of inserted section
      !-----------------------------------------------------! 

      select case (isonics)
     
      !-----------------------------------------------------! 
      !mid-canopy and above canopy sonics available
      !-----------------------------------------------------! 
      case(2)
      TL = 0.3*hcpy/CUSTARcpy
      TAU = 4.*TL
      RFACTORcpy = ((1.-exp(-TAU/TL))*(TAU-TL)**(3./2.))/
     &     (TAU-TL+TL*exp(-TAU/TL))**(3./2.)
      do k = 1,levobscpy 
         AKH(k) = AKH(k)*RFACTORcpy
      enddo
      ! Interpolate u* from canopy top to top of tower, and 
      ! repeat adjustment between canopy and tower tops. 
      do k=levobscpy+1,levobstow
         CUSTAR = CUSTARcpy+DUSTARDZ*(Z(K)-zobscpy)
         TL = 0.3*hcpy/CUSTAR
         TAU = 4.*TL
         RFACTORtow = ((1.-exp(-TAU/TL))*(TAU-TL)**(3./2.))/
     &      (TAU-TL+TL*exp(-TAU/TL))**(3./2.)
         AKH(k) = AKH(k)*RFACTORtow
      enddo
      ! calculate an average KH at observation height for calculating fluxes
      ! ka - fluxes calculated at all heights
      ! ka avekh=0.5*(akh(levobstow-1)+akh(levobstow))	
      
      !-----------------------------------------------------! 
      ! Only above canopy sonics available
      !-----------------------------------------------------! 
      case(1)
      TL = 0.3*hcpy/CUSTARtow 
      TAU = 4.*TL
      RFACTORtow = ((1.-exp(-TAU/TL))*(TAU-TL)**(3./2.))/
     &     (TAU-TL+TL*exp(-TAU/TL))**(3./2.)
      do k = 1,levobstow 
         AKH(k) = AKH(k)*RFACTORtow
      enddo
      ! calculate an average KH at observation height for calculating fluxes
      ! ka avekh=0.5*(akh(levobstow-1)+akh(levobstow))	!ka - added to calculate fluxes
      
      end select
      !-----------------------------------------------------! 
      ! End of inserted section 
      !-----------------------------------------------------! 
      !-----------------------------------------------------! 
     
      ! ka - added to calculate fluxes at all heights
      do k=2,nlev
         flxkh(k)=0.5*(akh(k-1)+akh(k))
      enddo
      
!**********************************************************************************************************************!
!
! Update concentrations and temperature after turbulent transport
!
!**********************************************************************************************************************!
      if (l30) read (5093, *) timin, alwc !shc aerosol water content [=] microgram/m3
      if (l30) read (5094, *) timin, vset !shc setting velocity [=] m/s

      call newc(akh,vcp,vcnc,vset) 

      ! Now move on to calculate the new temperature profile 
      call newt(akh)
c!**********************************************************************************************************************!
c!
c! Here the previously calculated KH is back-tracked to give KM
c! for the vertical exchange of mass (chemical species) instead
c! of heat, in advance of calling newc (CACHE's turbulent 
c! transport subroutine)
c!
c!**********************************************************************************************************************!
c      ! And now re-readjust KH
c      select case (isonics)
c      case (2)
c      ! Readjust in-canopy KH
c      do k = 1,levobscpy
c         AKH(K) = AKH(K)/RFACTORcpy
c      enddo
c      ! And then for mid-canopy to tower
c      do k=levobscpy+1,levobstow
c         CUSTAR = CUSTARcpy+DUSTARDZ*(Z(K)-zobscpy)
c         TL = 0.3*hcpy/CUSTAR
c         TAU = 4.*TL
c         RFACTORtow = ((1.-exp(-TAU/TL))*(TAU-TL)**(3./2.))/
c     &       (TAU-TL+TL*exp(-TAU/TL))**(3./2.)
c         AKH(K) = AKH(K)/RFACTORtow
c      enddo
c      
c      case (1)
c
c      ! Readjust KH to the tower
c      do k = 1,levobstow
c         TL = 0.3*hcpy/CUSTARtow
c         TAU = 4.*TL
c         RFACTORtow = ((1.-exp(-TAU/TL))*(TAU-TL)**(3./2.))/
c     &       (TAU-TL+TL*exp(-TAU/TL))**(3./2.)
c         AKH(K) = AKH(K)/RFACTORtow
c      enddo
c      end select

!**********************************************************************************************************************!
!
! Chemistry
!
!**********************************************************************************************************************!
       ! Level loop for chemistry
       do 10 k=1,nlev    
         lev=k
!         lev = z(k)
         t0=tstart
         t1=tstart+deltim

         ! dtlast=0.
         dtlast=dtl2
         iitime=itime
         llev=lev

         ! species concentrations
         ! emivalue in 1/s , set to zero, emissions are handled within newc
         do j=1,nspec
            vc(j)=vcp(k,j)
            emivalue(j)=0.0 
         enddo  

         ! biogenic emission rates
         do j=1,14
            bvocems(j)=bvoc(k,j)
         enddo
         ! advection rates
         do j=1,10
            advems(j)=advec(k,j)
         enddo
         
!         ! deposition rates
!         ! ka - added deposition of meoh and etoh; altered nspec to ndspec
!         do j=1,nspec 
!            depr(j)=deprat(k,j)
!         enddo
!
         ! concentrations

         tt=t(lev)
         pp=p(lev)
         moist = vc(lh2o)*18/28.8 
         rho_phy = pp/(tt*287.058)

         avekh=flxkh(lev)  !ka - added for calculating fluxes at all levels
         ! humidity in volume mixing ration of water vapor
         ! s=qva(i,j,k)/psb(i,j)

         dum1=vc(lmacr)

!**********************************************************************************************************************!
!                 RACM-MIM from WRF  
!**********************************************************************************************************************!
c Photolysis frequencies (preliminary)
         photoj( 1)= phot1( 1,lev)           !no2
         photoj( 2)= phot1( 2,lev)            !no3 -> no
         photoj( 3)= phot1( 3,lev)             !no3 -> no2 + o3p
         photoj( 4)= phot1( 4,lev)       !o3 -> o3p
         photoj( 5)= phot1( 5,lev)        !o3 -> o1d
         photoj( 6)= phot1( 6,lev)           !hono
         photoj( 7)= phot1( 7,lev)        !h2o2
         photoj( 8)= phot1( 8,lev)        !ooh1 (cacm)
         photoj( 9)= phot1( 9,lev)        !ooh2 (cacm)
         photoj(10)= phot1(10,lev)        !rpr2301 (cacm)
         photoj(11)= phot1(11,lev)        !hcho -> h2 + co
         photoj(12)= phot1(12,lev)        !hcho -> 2 ho2 + co
         photoj(13)= phot1(13,lev)        !ald
         photoj(14)= phot1(14,lev)        !ald2
         photoj(15)= phot1(15,lev)        !ketl
         photoj(16)= phot1(16,lev)        !ooh6101
         photoj(17)= phot1(17,lev)        !ooh6102
         photoj(18)= phot1(18,lev)        !hac
         photoj(19)= phot1(19,lev)        !ap6101
         photoj(20)= phot1(20,lev)        !ap6102
         photoj(21)= phot1(21,lev)        !nald
         photoj(22)= phot1(22,lev)        !hpald
         photoj(23)= phot1(23,lev)        !pacald
         photoj(24)= phot1(24,lev)        !macr
         photoj(25)= phot1(25,lev)        !rp6301
         photoj(26)= phot1(26,lev)        !ooh6301
         photoj(27)= phot1(27,lev)        !gly -> 0.45 hcho + 1.55 co + 0.80 ho2
         photoj(28)= phot1(28,lev)        !rp7102
         photoj(29)= phot1(29,lev)        !mgly (ch3cocho)
         photoj(30)= phot1(30,lev)        !pina
         photoj(31)= phot1(31,lev)        !nrpa
         photoj(32)= phot1(32,lev)        !edlm
         photoj(33)= phot1(33,lev)        !lmkt
         photoj(34)= phot1(34,lev)        !rpr8401
         photoj(35)= phot1(35,lev)        !rpr8501
         photoj(36)= phot1(36,lev)        !ooh8601
         photoj(37)= phot1(37,lev)        !ooh8602
         photoj(38)= phot1(38,lev)        !ISOP-POW
         photoj(39)= phot1(39,lev)        !ISOP-POW
         photoj(40)= phot1(40,lev)        !ISOP-POW
         photoj(41)= phot1(41,lev)        !ISOP-POW
         photoj(42)= phot1(42,lev)        !ISOP-POW
         photoj(43)= phot1(43,lev)        !ISOP-POW
         photoj(44)= phot1(44,lev)        !ISOP-POW
         photoj(45)= phot1(45,lev)        !ISOP-POW
         photoj(46)= phot1(46,lev)        !ISOP-POW
         photoj(47)= phot1(47,lev)        !ISOP-POW
         photoj(48)= phot1(48,lev)        !ISOP-POW
         photoj(49)= phot1(49,lev)        !ISOP-POW
         photoj(50)= phot1(50,lev)        !ISOP-POW
         photoj(51)= phot1(51,lev)        !ISOP-POW

c Simple parameterization of photolysis frequencies.
c Preliminary !!!!!!!
c photolysis profile follows PAR profile
c levhtr=3, and extf = 0 for the bottom three layers

      if(lev.le.levhtr) then
        extf(lev)=facte(1)
      else
        if((lev-levhtr).le.jtot) extf(lev)=facte(lev-levhtr+1) ! extf between layer 4 and 15 is set to zero (why??)
        if(lev.gt.levcpy) extf(lev)=1.
      endif
c      extf=1. ! Test: no photolysis profile

c scale to local time 
      do iph=1,nrect_j
         if(coszen.gt.0.1) then
           photoj(iph)=photoj(iph)*ratiod
           photoj(iph)=photoj(iph)*coszen*extf(lev)
         else
           photoj(iph)=0.0
         endif
         photout(iph,lev)=photoj(iph)
      enddo
            ph_no2          = photoj( 1)  
            ph_no3o2        = photoj( 2)   
            ph_no3o         = photoj( 3) 
            ph_o33p         = photoj( 4)
            ph_o31d         = photoj( 5) 
            ph_hno2         = photoj( 6)
            ph_h2o2         = photoj( 7)  
            ph_ooh1         = photoj( 8) ! no3-->no+o2 
            ph_ooh2         = photoj( 9)
            ph_rpr2301        = photoj(10) ! hcho-->h2
            ph_ch2om        = photoj(11) ! hcho-->2ho2
            ph_ch2or       = photoj(12) ! ald 
            ph_ch3cho       = photoj(13) ! op1
            ph_ald2         = photoj(14) ! op2
            ph_ch3coc2h5     = photoj(15) ! paa
            ph_ooh6101    = photoj(16) ! ket (same as hket)
            ph_ooh6102       = photoj(17)  
            ph_hac      = photoj(18) 
            ph_ap6101     = photoj(19) ! mgly
            ph_ap6102     = photoj(20) ! dcb
            ph_nald      = photoj(21) ! onit
            ph_hpald         = photoj(22) ! macr 
            ph_pacald    = photoj(23) ! hket 
            ph_macr    = photoj(24) ! hket 
            ph_rp6301    = photoj(25) ! hket 
            ph_ooh6301    = photoj(26) ! hket 
            ph_hcochob    = photoj(27) ! hket 
            ph_rp7102    = photoj(28) ! hket 
            ph_ch3cocho    = photoj(29) ! hket 
            ph_pina    = photoj(30) ! hket 
            ph_nrpa    = photoj(31) ! hket 
            ph_edlm    = photoj(32) ! hket 
            ph_lmkt    = photoj(33) ! hket 
            ph_rpr8401    = photoj(34) ! hket 
            ph_rpr8501    = photoj(35) ! hket 
            ph_ooh8601    = photoj(36) ! hket 
            ph_ooh8602    = photoj(37) ! hket 
            ph_1co2ooh = photoj(38) ! ISOP-POW
            ph_3ooh4co = photoj(39)  
            ph_3ooh4oh = photoj(40) 
            ph_1oh2ooh = photoj(41) 
            ph_mvk3oh4co = photoj(42) 
            ph_mvk3co4oh = photoj(43)
            ph_1co4ooh = photoj(44)
            ph_1ooh4co  = photoj(45)
            ph_mcrenol = photoj(46)
            ph_mvkenol= photoj(47)
            ph_mvk3ooh4co = photoj(48)
            ph_1n4co = photoj(49)
            ph_1co4n = photoj(50)
            ph_3co4n = photoj(51)
            
c          if(vc(lno).le.1.0E-11.and.lev.lt.18.and.
c     & timloc.lt.5.or.timloc.gt.16)
c            vc(lno) = 1.2*vc(lno)
            CALL cacm_interface(timloc, 
     &     nspec, vc,  deltim,         
     &     pp, tt, rho_phy, moist,   
     &     ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o,  
     &     ph_hno2, ph_hno3, ph_hno4, ph_h2o2, ph_ch2or,  
     &     ph_ch2om, ph_ch3cho, ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, 
     &     ph_ch3cocho, ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2,
     &     ph_hcochob, ph_macr, ph_hobr, ph_br2, ph_bro,  
     &     ph_n2o5, ph_o2, ph_pan, ph_acet, ph_mglo, 
     &     ph_hno4_2, ph_c6h5cho, ph_unald, ph_n2o, ph_pooh,  
     &     ph_mpan, ph_mvk, ph_etooh, ph_prooh, ph_onitr,
     &     ph_acetol, ph_glyald, ph_hyac, ph_mek, ph_open, 
     &     ph_gly, ph_acetp, ph_xooh, ph_isooh, ph_alkooh, 
     &     ph_mekooh, ph_tolooh, ph_terpooh,
     &     ph_ooh1, ph_ooh2, ph_rpr2301, ph_ald2, ph_ooh6101,    ! new photolysis
     &     ph_ooh6102, ph_hac, ph_ap6101, ph_ap6102, ph_nald, 
     &     ph_hpald, ph_pacald, ph_rp6301, ph_ooh6301, ph_rp7102, 
     &     ph_pina, ph_nrpa, ph_edlm, ph_lmkt, ph_rpr8401, 
     &     ph_rpr8501, ph_ooh8601, ph_ooh8602,
     &    ph_1co2ooh, ph_3ooh4co,  ph_3ooh4oh, ph_1oh2ooh, ph_mvk3oh4co,
     &     ph_mvk3co4oh, ph_1co4ooh, ph_1ooh4co, ph_mcrenol, ph_mvkenol,
     &     ph_mvk3ooh4co, ph_1n4co, ph_1co4n, ph_3co4n, 
     &     RCONST)

         chemr8_updated(k,:) = chemr8 
         rrat_updated(k,:) = rrat 

         dum2=vc(lmacr)-dum1
         if(lev.eq.levcpy+1) fmacr2=dum2/deltim*1e9

         if(lev.eq.1) then
           dtl=dmin1(dble(deltim),dtlast)
         else
           dtl=dmin1(dtl,dtlast)
         endif

         do jspec=1,nspec
            vcp(lev,jspec)=max(dble(vc(jspec)),0.d0) !shc "0.0"->"0.d0
         enddo

           
           
!**********************************************************************
!        MPMPO (the Model to Predict the Multiphase Partitioning of
!        Organics)
!**********************************************************************
         is_call_mpmpo = .true. !shc
         if (impmpo.eq.0) is_call_mpmpo = .false. !ka added impmpo switch via inputn

         do jspec=1,nspec
            boxcacm(jspec) = max(dble(vc(jspec)),0.d0)
            depr(jspec)=deprat(k,jspec)
         enddo

         do ivar=1,nmpmpo
           parcni(ivar)=vcnc(lev,ivar)
         enddo

         call call_mpmpo(tstamp,deltim,lev, z(lev),boxcacm, tt,
     &       parcni,parcno, relh(k), depr, alwc, is_call_mpmpo) 
         ! pass updated boxcacm to CACHE
         do jspec=1,nspec
            vcp(lev,jspec) = boxcacm(jspec)
         enddo
         
         ! pass updated particle conc. to CACHE
         do jspec=1,nmpmpo
            vcnc(lev, jspec)=max(parcno(jspec), 0.d0)
         enddo

  10     continue    
         ! End of level loop
      
!**********************************************************************************************************************!
! Fluxes at canopy top  (convert to ug/m^2/s) 
!**********************************************************************************************************************!

      fiso=-akh(levcpy)*
     &           (vcp(levcpy+1,lisop)*1.e9-vcp(levcpy,lisop)*1.e9)/
     &               (z(levcpy+1)-z(levcpy))*40./1000.*68.
      fapi=-akh(levcpy)*
     &       (vcp(levcpy+1,lapin)*1.e9-vcp(levcpy,lapin)*1.e9)/
     &               (z(levcpy+1)-z(levcpy))*40./1000.*136.
! cacm doens't have limonene
!      flim=-akh(levcpy)*
!     &     (vcp(levcpy+1,llim)*1.e9-vcp(levcpy,llim)*1.e9)/
!     &               (z(levcpy+1)-z(levcpy))*40./1000.*136.
      fmacr=-akh(levcpy)*
     &      (vcp(levcpy+1,lmacr)*1.e9-vcp(levcpy,lmacr)*1.e9)/
     &               (z(levcpy+1)-z(levcpy))*40./1000.*70.
      facet=-akh(levcpy)*
     &      (vcp(levcpy+1,lketl)*1.e9-vcp(levcpy,lketl)*1.e9)/
     &               (z(levcpy+1)-z(levcpy))*40./1000.*70.
      facald=-akh(levcpy)*
     &      (vcp(levcpy+1,lald1)*1.e9-vcp(levcpy,lald1)*1.e9)/
     &               (z(levcpy+1)-z(levcpy))*40./1000.*70.

!**********************************************************************************************************************!
! Fluxes above canopy top, height depends on grid spacing and choice of izo
!**********************************************************************************************************************!
      izo=levcpy+3
      izu=izo-1
      fapi2=-akh(izu)*
     &       (vcp(izo,lapin)*1.e9-vcp(izu,lapin)*1.e9)/
     &               (z(izo)-z(izu))*40./1000.*136.
      fiso2=-akh(izu)*
     &       (vcp(izo,lisop)*1.e9-vcp(izu,lisop)*1.e9)/
     &               (z(izo)-z(izu))*40./1000.*68.
!      flim2=-akh(izu)*
!     &       (vcp(izo,llim)*1.e9-vcp(izu,llim)*1.e9)/
!     &               (z(izo)-z(izu))*40./1000.*136.

      fmacr2=-akh(izu)*
     &       (vcp(izo,lmacr)*1.e9-vcp(izu,lmacr)*1.e9)/
     &               (z(izo)-z(izu))*40./1000.*70.

!**********************************************************************************************************************!
! New soil temperature and moisture
!**********************************************************************************************************************!
      call erdbod

      do 222 i=1,nsoil
      if (eta(i).lt.etamin) goto 223
      goto 224

 223  write(06,3900)i,eta(i)
 3900 format(1x,i5,e20.10,'!!!')

      eta(i)=etamin
 224  if(eta(i).ge.etamax) goto 225
      goto 222
 225  write(06,3900)i,eta(i)

      eta(i)=etamax
 222  continue

!**********************************************************************************************************************!
! Boundary condition for temperature and humidity
!**********************************************************************************************************************!
 
      p1b=p21(2)*h(2)
      p1a=0.01*vcp(1,lh2o)*p(1)
      xla=(xlam(1)+xlam(2))*.5
      fnet=-rnet(1)
      qvap1=vcp(1,lh2o)*0.622  !approximately
      qvap2=vcp(2,lh2o)*0.622
      patm=p(1)
!**********************************************************************************************************************!
! New surface temperature and moisture
!**********************************************************************************************************************!
       call grenzn(fnet,patm,p21(1),h(1),rhoair(1),psi(1),
     &p21t(1),psieta(1),xl11(1),xl11(2),xl22(1),xl22(2),psit(1),
     &d1t(1),d1t(2),d1eta(1),d1eta(2),d2t(1),d2t(2),d2eta(1),d2eta(2),
     &xl210(1),xlam(1),xlam(2),akh(1),z(1),z(2),zs(1),zs(2),
     &theta(1),theta(2),t(1),vcp1,qvap2,
     &temps(1),temps(2),eta(1),eta(2),
     &qgrnd,qsatg,afl1,aflh,bflh,bfl1,bfl2,deltim,tausum,jtime,
     &p1a,p1b)

      tsurf=temps(1)
      tsfc=tsurf-273.15

      if(afl1.ge.0.0) vsum=vsum+afl1*deltim
      if(afl1.lt.0.0) tsum=tsum-afl1*deltim
      vcp(1,lh2o)=vcp1      
 
!**********************************************************************************************************************!
! Update time step; t1 is the next time step of timloc
!**********************************************************************************************************************!

      timloc=t1/3600.

      if(lprin) then
       if(timloc.ge.24.) then
        ptime=ptime-24.
       endif
        ptime=ptime+dpt-1.e-11
      endif

      tstart=tstart+deltim
      dtl2= dtl

!**********************************************************************************************************************!
! Print output data to files at user-set time intervals
!**********************************************************************************************************************!

      if (l30 .and. itime .ne. 1) then                   
         write(06,'(/," Model time: ",f10.2," s; ",f8.2)') 
     &   timin,timloc                                    
         write(06,'(" Simulation is ",f6.1," % complete")')
     &   100.0*(itime-1)/ntimes                          

         CALL print_data(nlev,nspec,tstamp,ichem,averad,T,RelH,
     &    wind,P,theta,rhoair,akh,vcp,photout,levcpy,z,coszen, 
     &    source,ustara,levhtr,rnlam,rnet,radtop,frarea,templf,rst,
     &    depn, chemr8_updated,rrat_updated)
      endif                                        
      if (l30) t30 = t30+1800.

      ! Advance simulation time step
      tstamp = tstamp + deltim

 100  continue
      ! End of time loop
     
      stop 'end of program FORCAST'
      end

!! ****************************************************************
      block data
!! ****************************************************************
! 
!      INTEGER NLEVEL, NSOIL
      DOUBLE PRECISION ALPHA, BETA
      DOUBLE PRECISION AAA, AK,  BBB,  CCC
      DOUBLE PRECISION CP,  DDD, GRAV, P0
      DOUBLE PRECISION PI,  R0,  R1,   TCONV
      DOUBLE PRECISION XL
!      parameter (nlevel=40,nsoil=15)
!      common/imp/alpha,beta
      common/konst/p0,aaa,bbb,ccc,ddd,ak,r0,r1,grav,cp,xl,tconv
!
!c       parameters for the crank-nicolson scheme (alpha+beta=2.)
      data alpha/2.0/,beta/0.0/
!c       miscellaneous parameters
      data p0/1000.00/,aaa/1.809567/,
     & bbb/17.269388/,ccc/4717.3061/,ddd/35.86/,
     &tconv/273.15/,ak/0.4/,grav/9.81/
      data r0/287.05/r1/461.51/,cp/1004.4/,xl/2499000./
      end
c ##############################################################
c ##############################################################
c Initialise the model with profiles of all prognostic variables 
c and many parameters
c ##############################################################
      subroutine init(vcp,phot1,emiss,wind,rlat,clump,ustara,dpt,
     &  ichem,impmpo)  !ka - added MPMPO switch to inputn
      USE cacm3_Parameters, ONLY: NSPEC, NREACT, NVAR  
      USE cacm3_Precision,  ONLY:dp 
      use parameters  
      implicit none 
!      INTEGER             NLEV
!      INTEGER             NSOIL
!      INTEGER             NSPEC
!      INTEGER             NRECT
      INTEGER             MAXEMI
      INTEGER             MAXRECT
      INTEGER             MAXPHR
      INTEGER             MAXPHPT
      DOUBLE PRECISION    CLUMP, DPT, EMISS !shc emiss not used
      DOUBLE PRECISION    PHOT1, RLAT, USTARA
      DOUBLE PRECISION    VCP, WIND
      DOUBLE PRECISION    ADUM,  ALPHAB, BDUM, CIDCAI
      DOUBLE PRECISION    CLAYFC, DFMIN, DGEOM,DMAX
      DOUBLE PRECISION    DTLEAF, E, EFISOIN
      DOUBLE PRECISION    EFOVOCIN, EFTERPIN, EPS, ESAT
      DOUBLE PRECISION    FACTIR, GR, HA, HH
      DOUBLE PRECISION    HPSI, P1, PES, PLSPC
      DOUBLE PRECISION    PSIANF, PSITP1, PSOLID
      DOUBLE PRECISION    REFLAIIN, RGAS, ROWSPC, SANDFC
      DOUBLE PRECISION    SGEOM, SILTFC, SOVOCIN, STDLNG
      DOUBLE PRECISION    STERPIN, SUMD, TIMIN, TSTAMP
      DOUBLE PRECISION    TT, VG, VG0, XLAT
      DOUBLE PRECISION    XLONG, XX, ZM, ZMM
      INTEGER             I,   IBT,   IER,   IGUENTH
      INTEGER             II,  IMUNU, IO,    IST
      INTEGER             J,   JJ,    JMAX,  JMIN
      INTEGER             JZ,  K,     LEV,   M
      INTEGER             N
      DOUBLE PRECISION    COSDEC, COSLAT, DECL, DECMAX
      DOUBLE PRECISION    DLONG, EQTM, SINDEC, SINLAT
      DOUBLE PRECISION    TANLAT
      DOUBLE PRECISION    P, RELH, RHOAIR
      DOUBLE PRECISION    T, THETA, PRES
      DOUBLE PRECISION    ALBDRY, ALBWET, CP3, ETAMAX
      DOUBLE PRECISION    ETAMIN, FELDKA, POR, RHO3
      DOUBLE PRECISION    TREFK,  TREFP,  XKR
      DOUBLE PRECISION    HCPY, HTR,  SIZELF, TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY, LEVHTR
      DOUBLE PRECISION    CONTOT, CPESTR, CPHSTR, ECPYS
      DOUBLE PRECISION    ETOTW, EVSMIC, EVTOT, HCPYS
      DOUBLE PRECISION    HSOIL, SCOND, WCPYS
      INTEGER             IHRWET
      DOUBLE PRECISION    SUNAZM, VIEWAZ,  VIEWZN
      DOUBLE PRECISION    XINTV
      INTEGER             NOAZMV, NOZENV
      DOUBLE PRECISION    ADVEC, ADVEMS
      DOUBLE PRECISION    BVOC, BVOCEMS, BEMIS
      DOUBLE PRECISION    SOVOCE

      DOUBLE PRECISION    AXA, BXB,  CXC
      DOUBLE PRECISION    D1ETA, D1T
      DOUBLE PRECISION    D2ETA, D2T
      DOUBLE PRECISION    DATM, DXD,  ETA
      DOUBLE PRECISION    H, P21,  P21T
      DOUBLE PRECISION    PSI,  PSIETA
      DOUBLE PRECISION    PSIT, PSITET
      DOUBLE PRECISION    RHO0, RHO1, RHO22
      DOUBLE PRECISION     XK,   XL11
      DOUBLE PRECISION    XL210, XL22
      DOUBLE PRECISION    XLAM
      DOUBLE PRECISION    A0, B0, C0, DEL
      DOUBLE PRECISION    ZMIN1, ZMIN2
      INTEGER             INRAD, IWPM2
      DOUBLE PRECISION    FRWET, FRWTMX, PILAST
      DOUBLE PRECISION    PINT, PINT1, TWATER, WTP
      DOUBLE PRECISION    AFL1, AFLH, BFL1, BFL2
      DOUBLE PRECISION    BFLH, FNET,  GDOWN, GFNET
      DOUBLE PRECISION    PATM, QGRND, QSATG, TAUSUM
      DOUBLE PRECISION    TSUM, TSURF, VSUM
      DOUBLE PRECISION    AAA, AK,  BBB, CCC
      DOUBLE PRECISION    CP,  DDD, GRAV, P0
      DOUBLE PRECISION    PII, R0, R1, TCONV
      DOUBLE PRECISION    XL
!      DOUBLE PRECISION    PI, PID180, PID2, SIGMA
      INTEGER             IWRITE, KMAX
      DOUBLE PRECISION    CLAI,    CT, DF
      DOUBLE PRECISION    DISTLS,  FR, TOTLAI
      INTEGER             ITOT, ITOTP1,JTOT
      DOUBLE PRECISION    XINT,        XINTZ
      INTEGER             ISPHER,      NALPHA,      NXINTZ
      DOUBLE PRECISION    BKV,   CAIR,   CONMIN,  CUCOND
      DOUBLE PRECISION    D1,    FACJ,   OX,      RASTOM
      DOUBLE PRECISION    RSFAC, RXCHAM, ZNON
      INTEGER             IC3C4
      DOUBLE PRECISION    ALEAF, EMIS, EMISOL
      DOUBLE PRECISION    EXPDIF, RLAYR, RLEAF
      DOUBLE PRECISION    RSOIL,  TLAYR, TLEAF
      DOUBLE PRECISION    DSTNET,  DSTRAD
      DOUBLE PRECISION    FRAREA,  TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    COSZEN, FBEAM1, RADABV, RATIO
      DOUBLE PRECISION    RATIOD, RATION, ZENANG
      DOUBLE PRECISION    BMFLX, D,     RNDIV
      DOUBLE PRECISION    RNET,    RNLAM, TSFC, U
      DOUBLE PRECISION    ANSTOM, PSI1,  PSI2, RADN
      DOUBLE PRECISION    RCUT20, RSEXP, RSM, RSMIN
      DOUBLE PRECISION    TRSMAX, TRSMIN, TRSOPT
      DOUBLE PRECISION    AROOT,       CPYTR,  FROOT,   PSISUM
      DOUBLE PRECISION    PSITOP,      PSIXY,  RESROT,  ROOTSM
      DOUBLE PRECISION    ROOTUP,  RROOT
      DOUBLE PRECISION    DELS, SL,  ZS
      DOUBLE PRECISION    AKS,   AN, ASOIL, BD
      DOUBLE PRECISION    BSOIL, BX, CSOIL, DSOIL
      DOUBLE PRECISION    ESOIL, PE, WS
      DOUBLE PRECISION    ALPHADAY,    ALPHAN,      VGDAY,       VGN
      DOUBLE PRECISION    ZRDAY,       ZRN
      DOUBLE PRECISION    CK,    CLAM, CPSI
      DOUBLE PRECISION    CPSIT,  SEK, SEL, SEP
      DOUBLE PRECISION    SEPT,  SK, SPP,  SPT
      DOUBLE PRECISION    SXL
      DOUBLE PRECISION    temps
      INTEGER             ISK,         ISL,         ISP,         ISPT
      REAL                TIMLOC, deltim
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      INTEGER             ISORT
      DOUBLE PRECISION    sorg

!      parameter (nlev=40,nsoil=15)    ! number of levels  
!      parameter (nspec=84, nrect=249)
      parameter(maxemi=nspec) ! maximal number of emissions
      parameter(maxrect=500)
      parameter(maxphr=100, maxphpt=100)

      double precision emivalue(maxemi)
      dimension vcp(nlev,nspec), phot1(nrect_j,nlev),
     &wind(nlev)
      dimension gr(99)
      integer ichem,impmpo  !ka - added MPMPO switch to inputn
      double precision grndT0,towT0,lapse !ka - moved input of initial temperature to inputn

      common/konst/p0,aaa,bbb,ccc,ddd,ak,r0,r1,grav,cp,xl,tconv
      common/astron/eqtm,decl,sindec,cosdec,decmax,sinlat,coslat,
     & tanlat,dlong
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
      common/atmo1/ t(nlev),theta(nlev),p(nlev),rhoair(nlev)
     & ,relh(nlev)
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/cpy2/hsoil,hcpys,evtot,etotw(99),contot(99),scond(10,99)
     &,ecpys,cphstr,cpestr,wcpys,evsmic,ihrwet(99)
      common/deg/sunazm,viewzn(99),viewaz(99),nozenv,noazmv
     & ,xintv(20,90)
      common /emissions/ emivalue
      REAL(KIND=dp) :: z, zf, dz, dzf
      common/height/z(nlev),zf(nlev),dz(nlev),dzf(nlev)

      common/oheight/isonics,zobstow,zmodtow,levobstow
     &,zobscpy,zmodcpy,levobscpy,zmxd,levmxd			!ka - added for sonics data
      integer levobstow,levobscpy,levmxd,isonics		!ka - added for sonics data
      dimension zbot(nlev)					!ka - added for sonics data
      double precision zobstow,zmodtow,zobscpy,zmodcpy,zmxd,zbot!ka - added for sonics data
      logical zflgt,zflgc,zflgm					!ka - added for sonics data

      common/grid/ del,zmin1,zmin2,a0(nlev),b0(nlev),c0(nlev)
      common/inter1/wtp(99),frwet(99),frwtmx,pint(99),pilast(99)
     &,pint1(99),twater
      common/misc1/iwrite(9,99),kmax
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/misc6/xint,xintz(99),nalpha,ispher,nxintz
      common/photo2/ox,cair,rastom,cucond,d1,bkv,facj,rxcham,rsfac
     &             ,znon,conmin,ic3c4
      common /rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3)
     &,expdif(99),rlayr(3,30),tlayr(3,30)
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common /rad3/radabv(3),fbeam1(3),coszen,zenang
     &,ratiod,ration,ratio
      common/rad4/d(3,99),u(3,99),bmflx(3,99),rnet(99),rndiv(99),
     &tsfc,rnlam(3,99)
      common/resis2/radn,anstom,rcut20,rsmin,trsopt,trsmin,trsmax,rsm
     &,rsexp,psi1,psi2
      common/root1/froot(99),resrot(99),rootsm,psixy,psitop,rroot
     &,rootup(99),cpytr,psisum,aroot
      common/soil1/zs(nsoil),sl,dels
      common/soil4/pe,bx,bd,aks,an,ws,asoil,bsoil,csoil,dsoil,esoil
      common/feld/ eta(nsoil),psi(nsoil),psit(nsoil),temps(nsoil),
     &xk(nsoil),psieta(nsoil),
     &psitet(nsoil),xlam(nsoil),axa(nsoil),bxb(nsoil),
     &cxc(nsoil),dxd(nsoil),d1t(nsoil),
     &d2t(nsoil),d1eta(nsoil),d2eta(nsoil),
     &xl11(nsoil),xl22(nsoil),rho0(nsoil),
     &rho1(nsoil),rho22(nsoil),xl210(nsoil),
     &p21(nsoil),p21t(nsoil),h(nsoil),datm(nsoil)
      common/stuetzw/sep(18),spp(18),sek(17), sk(17),
     &sept(10),spt(10),sel(10),sxl(10),
     &cpsi(17,3),ck(16,3),cpsit(9,3),clam(9,3),
     &isp,isk,ispt,isl
      common/intf/ qgrnd,qsatg,tausum,tsum,vsum,tsurf,patm,
     &fnet,gdown,gfnet,afl1,aflh,bfl1,bfl2,bflh
      common/bodk/ etamin,etamax,xkr,por,trefp,trefk,cp3,rho3,
     &feldka,albdry,albwet

      common/stratpar/ vgday,vgn,zrday,zrn,alphaday,alphan,stable !shc
      logical stable

      common /treespec/ isort
      common /inputr/ inrad,iwpm2
!ka - added beta (bexp array)
      common/efact/ EFno,EFsyn(14),EFpl(14),sovoce(9),bvoc(nlev,14),
     & bemis(nlev,13),bvocems(14),advec(nlev,10),advems(10),bexp(13)
      double precision EFno,EFsyn,EFpl,bexp,fu
!ka - added beta; end of added section
      dimension sterpin(3),sovocin(9)

      pii=4.*atan(1.)      ! pi
      rgas=8.314
      decmax=sin(23.44*pid180)

c ##############################################################
c Set up the direct write output data files
c ##############################################################
      !------------------------------!
      !     Chemistry output data    !
      !------------------------------!
      open(unit=318,file='out/conc.out',status='unknown')        !ddw - chemical concentration
      open(unit=319,file='out/photo.out',status='unknown')       !ddw - photolysis rates
      open(unit=320,file='out/rrat.out',status='unknown')        !ddw - reaction rates
      open(unit=321,file='out/depn.out',status='unknown')        !ddw - reaction rates
      open(unit=322, file='out/gprd.out',status='unknown')

!      open(unit=315,file='out/cacm_flux.out',status='unknown')   !ka - added canopy-top fluxes
!      open(unit=317,file='out/cacm_zmix.out',status='unknown')   !ka - added vertical transport
!      open(unit=32, file='out/cacm_emis.out',status='unknown')
!      open(unit=33, file='out/cacm_depn.out',status='unknown')
!      open(unit=335,file='cacm_vdep.out',status='unknown')
!      open(unit=34, file='out/cacm_phot.out',status='unknown')

!ka - altered production and loss output files; added reaction rate output file
!      open(unit=36,file='out/cacm_rr8s.out',status='unknown')
!ka - altered production and loss output files; end of addition


      !------------------------------!
      !     Met output data          !
      !------------------------------!
      open(unit=41,file='out/cache_met.out',status='unknown')
      open(unit=42,file='out/cache_misc.out',status='unknown')

      !------------------------------!
      ! Canopy and vegetation data   !
      !------------------------------!
      open(unit=51,file='out/cache_cpy.out',status='unknown')

      !------------------------------!
      !          Soil data           !
      !------------------------------!
c      open(unit=61,file='soil.out',status='unknown')

      !------------------------------!
      !      MPMPO output data       !
      !------------------------------!
      open(unit=85,file='out/mpmpo_aer.out',status='unknown') !shc
      open(unit=86,file='out/mpmpo_gas.out',status='unknown') !shc
      open(unit=87,file='out/mpmpo_diff.out',status='unknown') !shc
      open(unit=88,file='out/mpmpo_gamma_aq_r.out',status='unknown') !ddw
      open(unit=89,file='out/mpmpo_gamma_aq_h.out',status='unknown') !ddw
      open(unit=90,file='out/mpmpo_org.out',status='unknown') !ddw

c ##############################################################
c Read parameters from file inputn - modify to be site-specific
c ##############################################################

      open(unit=15,file='inputn_07222016',status='old')
 
      read(15,*) jday,iyear                      ! ##### Start day and year
      icumdy=jday
      read(15,*) xlat,xlong,stdlng               ! ##### geographical location
      read(15,*) ichem,impmpo                    ! ##### chemistry mechanism switches; ka added MPMPO switch
      read(15,*) hcpy,htr,levhtr,levcpy,del      ! ##### canopy height, trunk space height, top level of trunk spece, top level of canopy, parameter for grid spacing
      read(15,*) isonics,zobstow,zobscpy,zmxd      ! ##### number of sonics measurement heights,above canopy observation height,in-canopy observation height,mixed layer height	!ka - added for sonics data
      read(15,*) tlai,sizelf                     ! ##### total projected LAI, leaf size
      read(15,*) (df(n),n=1,levcpy-levhtr)       ! ##### lai fractions per layer
      read(15,*) stable
      read(15,*) vgday,vgn                       ! ##### geostrophic wind (day, night) 
      read(15,*) zrday,zrn,alphaday,alphan       ! ##### roughness height (day, night), wind profile parameters (day, night)
      read(15,*) itot,jmax,jmin,dfmin,clump,kmax 
      read(15,*) inrad, iwpm2                    ! ##### read incoming solar radation (inrad=1), input is in W/m**2 (iwpm2=1) or micromoles per (m2 seconds) (iwpm2=0) 
      read(15,*) ratiod,ration                   ! modifiers for incoming radiation (only important when incoming radiation is not read from file (i.e. inrad=0)
      read(15,*) emis,emisol
      read(15,*) (rsoil(m),m=1,kmax),(rleaf(m),m=1,kmax), !reflectivities, transmissions
     &(tleaf(m),m=1,kmax)
      nxintz=9

!!!!!!!!!!!!!!!!!!!!!!!!
! echo inputn 
!!!!!!!!!!!!!!!!!!!!!!!!
      if (ichem.eq.0) then
!shc   write(06,'x,a') 'CACM chemistry scheme selected'
       write(06,'(x,a)') 'CACM chemistry scheme selected'
      elseif (ichem.eq.1) then
!shc   write(06,'x,a') 'RACM chemistry scheme selected'
       write(06,'(x,a)') 'RACM chemistry scheme selected'
      else
!shc   write(06,'x,a') 'Choice of chemistry scheme invalid'
!shc   write(06,'x,a') 'ichem must be 0 (CACM) or 1 (RACM-MIM)'
!shc   write(06,'x,a') 'CACM chemistry scheme will be used'
       write(06,'(x,a)') 'Choice of chemistry scheme invalid'
       write(06,'(x,a)') 'ichem must be 0 (CACM) or 1 (RACM-MIM)'
       write(06,'(x,a)') 'CACM chemistry scheme will be used'
       ichem=0
      endif

c ##############################################################
c                       Run-time outputs
c ##############################################################
      write(06,'(/," Model timestep",f5.1,"s; output at",f5.2,"hr")')
     &deltim,dpt
      write(06,'(" Simulation begins: 00:00 on Day",i4," of ",i4)')
     & jday,iyear
      write(06,'(" Location: ",f4.1,"N, ",f4.1,"W")') xlat,xlong
      write(06,'(" Number of model levels: ",i3)') nlev
      write(06,'(/,x,"Canopy levels ",i3," to ",i3)')
     & levhtr,levcpy
      write(06,'(x,"Canopy between ",f4.1," and ",f4.1," m")') 
     & htr,hcpy
      write(06,'(/,x,"Total LAI",f4.1)') tlai
      write(06,'(/,x,"Geostrophic wind (day, night)",x,2f5.1,"m")') 
     & vgday,vgn
      write(06,'(x,"Surface roughness (day, night)",x,2f5.2)') 
     & zrday,zrn
      write(06,'(x,"Wind attenuation parameter (day, night)",x,2f5.2)')
     & alphaday,alphan
      
      ! Canopy structure
      write(06,'(/,x,a45)')
     & 'Canopy types and leaf angle distributions:' 
      write(06, '(x,a16,4x,a20,3a8)') 
     & 'canopy type','f(theta)','theta','meu','neu'
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'planophile','2(1+cos(2*theta))/pi',26.76,2.770,1.172
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'erectophile','2(1-cos(2*theta))/pi',63.24,1.172,2.770
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'plagiophile','2(1-cos(4*theta))/pi',45.00,3.326,3.326
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'extremophile','2(1+cos(4*theta))/pi',45.00,0.433,0.433
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'uniform','2/pi',45.00,1.000,1.000
      write(06, '(x,a16,4x,a20,3f8.3)') 
     & 'spherical','sin(theta)',57.30,1.101,1.930
      write(06, '(x,a)') 
     & 'Source: Naren Goel & Don Strebel (1984), Agron. J., 76, 800-802'
c#####################################################
      call cpstruct(imunu,gr)
c#####################################################
c ##############################################################
c And biogenic emissions
c Read in synthesis and pool emission factors for isoprene, 
c individual monoterpenes, oxyVOCs and other: iso,api,bpi,lim,
c other MTs,sqt,meoh,acetaldehyde,acetone,mvk-mcr,mbo,ovoc
c Then soil NO emission factor
c ##############################################################
      read(15,*) (EFsyn(i),i=1,14)
      read(15,*) (EFpl(i),i=1,14)
      read(15,*) (bexp(i),i=1,13)      ! ka - added beta (bexp array); temperature activity factor for pool emissions
      read(15,*) EFno
      
c Emission factors for Isoprene and Monoterpene ef***(ii)
c where ii=1: Pool, ii=2: Synthesis term
c Numbers are in nmol per m**2 proj. leaf area per sec
c  
      fu=1./40.9*1.e-9
      do i=1,14
       EFsyn(i)=EFsyn(i)*fu
       EFpl(i)=EFpl(i)*fu
      enddo
c 
c  this is where stom resis, photosynthesis, resp. etc. parameters
c  are read in, they are to remain fixed.
c  read in rs vs light
777   read(15,*)rcut20,rsmin,anstom,radn          ! ##### Cuticular resistance, min. stomatal resistance, 
c  read rs vs temp
778   read(15,*)trsopt,trsmax,trsmin,rsexp,rsm    ! #####
c  read rs vs leaf psi
780   read(15,*)psi1,psi2
c  convert potentials from bars to j/kg
      psi1=psi1*100.
      psi2=psi2*100.
c
c  additions to subroutine in plant for leaf photosynthesis by farquhar eqs.
      read(15,*)rastom,d1,bkv
c  input root resis in m4 kg-1 s-1 (corn=3.e6) assume .6 of plant
c    resistance is in root
      read(15,*)rroot
      read(15,*)aroot

c  rxcham blr of one side of leaf in chamber used to meas p.s. (s/m)
c  rastom ratio of stom resis on top and bot of leaf.
c  rsfac  factor to adjust blr by for stom cond or resis ratio(li-6200)
      rsfac=(rastom**2+1.)/(rastom+1.)**2
c  cucond cuticular cond to water (mol/m2/s)
c  d1    vpd for start of stom closure (mbar)
c  d2    vpd for end of stom closure (mbar)
c  bkv   parameter for vpd response of stomata
c 
      hpsi=1.
      psitp1=psi1
      psitop=0.
      
      dlong=(stdlng-xlong)/15.
      sinlat=sin(xlat*pid180)
      coslat=cos(xlat*pid180)
      tanlat=sinlat/coslat
      rlat=xlat*pi/180.

c Solar declination
c#####################################################
      call declin 

c#####################################################
      itotp1=itot+1

c ##############################################################
c  Height levels 
c  To change the height levels modify variables in file inputn:
c  canopy height 'hcpy', trunk height 'htr', level of canopy top
c 'levcpy' and top of trunk space 'levhtr'. Select 'del' to get 
c the desired model domain height (note that not all combinations
c of hcpy, htr, levcpy, levhtr, & del are reasonable or possible)
c ##############################################################
      tt=exp((levhtr-1)*del)
      hh=exp((levcpy-1)*del)
      xx=((levcpy-1)*del)/((levhtr-1)*del)
      zmin1=(hcpy-htr*xx)/(hh+xx-tt*xx-1)
      zmin2=(zmin1+htr-zmin1*tt)/((levhtr-1)*del)

      do k=1,nlev
       z(k)=zmin1*exp((k-1)*del)+zmin2*(k-1)*del-zmin1
c ka - added for sonics data; one line in order to identify the correct level for the observations
       if (k.gt.1) zbot(k)=0.5*(z(k-1)+z(K))   !ka - only roughly true, but only used to identify observation levels
      enddo

! ka - added for sonics data; start of inserted lines
      zbot(1)=0.d0       !ka - added for observation levels
      zflgt=.true.
      zflgc=.true.
c ka - added next lines to find correct levels for the observations
      do k=40,1,-1
       if ((zobstow.gt.zbot(k)).and.(zflgt)) then
        levobstow=k
        zflgt=.false.
       endif
       if ((zobscpy.gt.zbot(k)).and.(zflgc)) then
        levobscpy=k
        zflgc=.false.
       endif
       if ((zmxd.gt.zbot(k)).and.(zflgm)) then   !ka - added for setting level to smooth to
        levmxd=k
        zflgm=.false.
       endif
      enddo
      zmodtow=z(levobstow)      !find mid-point height of above canopy obs level
      zmodcpy=z(levobscpy)      !find mid-point height of in-canopy obs level
!ka - added for sonics data; end of inserted lines

      write(06,'(/,x,"Model level mid-point heights (m):")')
      write(06,'(x,40i8)') (k,k=1,nlev) 
      write(06,'(x,40f8.2)') (z(k), k=1,nlev)
!ka - added for sonics data; start of inserted lines
      write(06,'(/,x,"Above canopy observation height (m):",f8.2)') 
     &   zobstow
      write(06,'(,x,"Model level and mid-point height (m):",i8,f8.2)') 
     &   levobstow,zmodtow
      write(06,'(/,x,"In-canopy observation height (m):",f8.2)')
     &   zobscpy
      write(06,'(,x,"Model level and mid-point height (m):",i8,f8.2)') 
     &   levobscpy,zmodcpy
      write(06,'(/,x,"Model level for KH smoothing:",i8)')
     &   levmxd
!ka - added for sonics data; end of inserted lines

      write(85,'(x,a)') 
     &'Conc. of condensible surrogate (ug m-3) in aerosol phase'
      write(86,'(x,a)') 
     &'Conc. of condensible surrogate (ug m-3) in gas phase'
      write(87,'(x,a)') 
     &'Condensation rate of surrogate (ug m-3 per time step)'

      zf(nlev)=z(nlev)+0.5*(z(nlev)-z(nlev-1))
      do k=1,nlev-1
       zf(k)=0.5*(z(k+1)+z(k))
       if(k.gt.1) dzf(k)=zf(k)-zf(k-1)
       dz(k)=z(k+1)-z(k)
      enddo
      dzf(1)=zf(1)
      dzf(nlev)=dzf(nlev-1)
      dz(nlev)=dz(nlev-1)

      jtot=levcpy-levhtr+1

c ##############################################################
c Specify leaf area index per layer in inputn (df)
c - ultimately would be best if df(z) = fn of z (default for 
c deciduous and for coniferous, following Wolfe et al and refs 
c therein); but currently user-specified. Should be modified
c by site using obs or generic values for ecosystem type.
c ##############################################################

      sumd=0.
      do ii=levhtr+1,levcpy
       jj=ii-levhtr
       sumd=sumd+df(jj)
      enddo

      do ii=1,levcpy-levhtr
       df(ii)=tlai/sumd*df(ii)
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!
! report lai vertical distribution
!!!!!!!!!!!!!!!!!!!!!!!!
      write(06,'(/,x,"Leaf area distribution:")')
      write(06,'(x,10i7)') (ii,ii=1,levcpy-levhtr)
      write(06,'(x,10f7.3)') (df(ii),ii=1,levcpy-levhtr)
!!!!!!!!!!!!!!!!!!!!!!!!
      
c ##############################################################
c Initialize atmospheric variables - use site-specific obs where
c possible => will need modifying for every new site and every
c new simulation time period 
c ##############################################################
      if(nlev.gt.60) 
     & write(06,'("Canopy radiation model permits only 60 layers,",
     & " but nlev is ",i4)') nlev 
c ##############################################################
c Nudge model with site obs (amb, 7/29/11) - firstly tower top
c sonic data (ie above the canopy) and then within canopy obs.
c If not available for a particular site use above canopy data only.
c ##############################################################
c Read headers from input files
      open(unit=5034,file='./data/wspabove.dat', status='unknown')
      open(unit=5044,file='./data/fort.5044', status='unknown') !???
      open(unit=5081,file='./data/ustarabove.dat', status='unknown')
      open(unit=5091,file='./data/sigwabove.dat', status='unknown')
      open(unit=5093,file='./data/alwc.dat', status='unknown') !???
      open(unit=5094,file='./data/vset.dat', status='unknown') !???

      read (5034, *)  ! U @ 29m (Harvard Forest) @ 36.4m  (UMBS CABINEX)
      read (5044, *)  ! UD @ 29m (Harvard Forest) @ 36.4m  (UMBS CABINEX)
      read (5081, *)  ! USTAR @ 29m (Harvard Forest) @ 34 m (UMBS CABINEX)
      read (5091, *)  ! SIGW @ 29m (Harvard Forest) @ 34 m (UMBS CABINEX)
      read (5093, *)  ! alwc !shc
      read (5094, *)  ! vset !shc

      if (isonics.gt.1) then
       open(unit=5082,file='./data/ustarbelow.dat', status='unknown')
       open(unit=5092,file='./data/sigwbelow.dat', status='unknown')
       read (5082, *)  ! USTAR @ 20.6 m ie in-canopy !shc for UMBS CABINEX setup
       read (5092, *)  ! SIGW @ 20.6 m ie in-canopy  !shc for UMBS CABINEX setup
      endif
c Read first line for initial conditions
      read (5034, *) timin,vg0

c check for error (amb, 7/29/11)
      if (vg0 .ne. -999.) vg = vg0
      
c wind
      if (vg0 .eq. -999.) vg = vgn

      zrough=zrn
      ustara = vg/30.


c ##############################################################
c (results in too small values for wind speed in trunk space)
c altered for UMBS cabinex observations =>
c likely to require modification for HF conditions
c #############################################################
      alphab=alphan

c Wind profile
c #############################################################
      call wprofil(vg,alphab,wind,ustara,tstamp)
c #############################################################
c Initialise atmospheric conditions
c 
c Temperature profile based on ground level temperature measurement 
c and adiabatic lapse rate calculated from radiosonde soundings.
c Site and time period specific values set in inputn
c 
      read(15,*) grndT0,towT0,lapse   !ka - added initial temperature at ground level and adiabatic lapse rate to inputn
c 
      do lev=1,nlev
       k=lev
c pressure
       pres=1013.25*100.  ! Surface presssure, pascal
       p(lev)=pres*((1.0-z(lev)/44308.0)**(1.0/0.19023))
c temperature
!shc   t(lev)=towT0-(z(lev)-29.)*lapse	!ka - ABL temperature profile based on initial temperature at 29m and lapse rate (from inputn) !shc Harvard Forest setup
       t(lev)=293.45-z(lev)*.0046 !shc for UMBS CABINEX, based on FORECaSTv1.0

!shc   if (z(lev).lt.29.) t(lev)=grndT0-z(lev)*(towT0-grndT0)/29. !shc Harvard Forest setup
       t(lev)=grndT0-z(lev)*lapse !ka - ABL temperature profile based on initial temperature at ground-level and lapse rate (from inputn)
c 
!ka       t(lev)=towT0-(z(lev)-z(levobstow))*lapse	!ka - ABL temperature profile based on initial temperature at tower and lapse rate (from inputn)
!ka       if (z(lev).lt.z(levobstow)) t(lev)=grndT0-z(lev)*(towT0-grndT0)/(z(levobstow)-2.0) !ka - assuming "ground-level" T at 2.0m
       theta(lev)= t(lev)*((100000.0/p(lev))**.286)
       dtleaf=0.

c air density in kg m-3 (see transp.f) (amb, 8/9/11)
       RHOAIR(lev) = P(lev)/287./T(lev)

c Rel. humidity -> water vapor (volume mixing ratio [ppmv*10**6], 
c corresponding to mole fraction)
c Modify for each site #########################################
       relh(lev)=.4
       if(z(lev).lt.500.) 
     &  relh(lev)=0.80-z(lev)*.4/500.          
       p1=relh(lev)*610.7*exp(17.1536*(t(lev)-273.15)
     &    /(t(lev)-38.33))
       vcp(lev,lh2o)=p1/p(lev) 

      enddo  ! end level loop
c#############################################################
c Initialise concentrations of chemical species at all levels
c and set photolysis frequencies for a zenith angle of zero
c#############################################################
      xx=0.01    ! parameter for height "correction" - default
      
      do i=1,nspec
         vcp(lev,i)=0.
      enddo
      
      do lev=1,nlev
         k=lev

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         ! Inital values are for zenith angle=0
         ! time units for rate constants : 1/min 
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         vcp(lev,lco2)=395*1.0E-6                 !    co2
         vcp(lev,lh2so4)=2.8e5/(p(lev)/101325./((82.05746/6.0221e23)
     & *t(lev)))*exp(-z(lev)*xx)*1.0E-6!SULF
         vcp(lev,lso2)=0.1e-3*exp(-z(lev)*xx)*1.0E-6 !SO2 - set from UMBS
         vcp(lev,lch4)=1.76*1.0E-6                  !CH4
!ka      VAR(110)=10.11e-3*exp(-(lev-29.)*xx*xx) !MEOH
!ka       if (lev.lt.29.) VAR(110)=10.11e-3	!obs @29m
         vcp(lev,laroh)=0.06e-3*exp(-z(lev)*xx)*1.0E-6 !AROH (toluene/xylene)
         vcp(lev,lpan2301)=1.0e-3*exp(-z(lev)*xx)*1.0E-6     !PAN(2301)
!cxc      VAR(114)=0.03e-3*exp(-lev*xx)              !PAN !ka - added for testing
         vcp(lev,letoh)=1.5e-3*.1*exp(-z(lev)*xx)*1.0E-6     !ETOH
         vcp(lev, larol)=0.2e-3*.1*exp(-z(lev)*xx)*1.0E-6  !AROL (toluene/xylene)
!cxc      VAR(144)=0.125e-3*exp(-lev*xx)!N2O5	!ka - for testing
!          vcp(lev,ln2o5)=1.50e-9*exp(-z(lev)*xx)!N2O5	!ka - set to give NOx:NOy=0.4 at HF
          vcp(lev, lhono)=154.2e-12*exp(-z(lev)*xx)!HONO
          vcp(lev,lalkm)=0.3e-9*exp(-z(lev)*xx) !ALKM (alkanes ~C8)
          vcp(lev,lh2o2)=0.5e-9*exp(-z(lev)*xx*xx)!H2O2
          vcp(lev,lhno3)=0.1e-9*exp(-z(lev)*xx) !HNO3
          vcp(lev, lketl)=2.5e-9*exp(-z(lev)*xx*xx)  !KETL (acetone)
          vcp(lev,lethe)=0.02e-9*.1*exp(-z(lev)*xx) !ETHE 
!!cxc      VAR(300)=0.04e-3*exp(-lev*xx)*0.2!BPIN
          vcp(lev,lbpin)=0.0        !BPIN	!ka - changed apr 2015
          vcp(lev, lco)=100.0e-9                 !CO
          vcp(lev,ldlmn)=0.04e-3*exp(-z(lev)*xx)*1.0E-6  !DLMN
!!cxc      VAR(316)=0.04e-3*exp(-lev*xx)*0.8!APIN !ka - changed apr 2015
          vcp(lev, lapin)=0.04e-3*exp(-z(lev)*xx)*1.0E-6 !APIN
!          vcp(lev,lhcho)=0.5e-3*exp(-z(lev)*xx)*1.0E-6    !HCHO 
          vcp(lev,lhcho)=5.0*exp(-z(lev)*xx)*1.0E-9    !HCHO 
          vcp(lev,lisop)=0.2e-3*exp(-z(lev)*xx)*1.0E-6   !ISO
          vcp(lev,lald1)=0.39e-3*exp(-z(lev)*xx)*1.0E-6     !ALD1 (acetaldehyde)
          vcp(lev,lmgly)=20.0e-6*exp(-z(lev)*xx)*1.0E-6     !MGLY 
          vcp(lev, lmvk)=0.50*exp(-z(lev)*xx)*1.0E-9 !MVK
          vcp(lev, lno3)=2.1e-6*exp(-z(lev)*xx)*1.0E-6      !NO3
          vcp(lev, lo3)= 42.0e-3*1.0E-6+z(lev)*0.10e-3*1.0E-6   !O3
          if(vcp(lev,lo3).gt.60.0e-9) vcp(lev,lo3)=60.0e-9
          if(z(lev).gt.1150.) vcp(lev,lo3)=60.0e-9
          vcp(lev,lmacr)=0.40e-3*exp(-z(lev)*xx)*0.5*1.0E-6 !MCR(+MVK)
          vcp(lev, lno2)= 0.834e-3*exp(-z(lev)*0.1*xx)*1.0E-6!NO2
          if(z(lev).gt.20.) vcp(lev, lno2)= 1.048e-9*exp(-z(lev)*0.1*xx)
          if(z(lev).gt.34.) vcp(lev, lno2)= 1.253e-9*exp(-z(lev)*0.1*xx)
          if(z(lev).gt.300.) vcp(lev, lno2)=200.0e-12
          if(z(lev).gt.600.) vcp(lev, lno2)=100.0e-12
          if(z(lev).gt.1500.) vcp(lev, lno2)=50.0e-12
          vcp(lev,lno)=3.0e-6*exp(-z(lev)*0.1*xx)*1.0E-6 !NO
          vcp(lev,lh2o)=0.0198*exp(-z(lev)*xx)  ! dependent on pressure
!
          sorg = 2.*2.44e-2/real(ncondensable)*1.0E-6 !shc ppm*1.0E-6 ---> g/mol
          do i=1, nvar
             if ( mw_of_icacm(i) > 0. ) then
                vcp(lev,i) = sorg/mw_of_icacm(i)*exp(-z(lev)*0.125*xx) 
             end if
          end do

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         ! Inital values are for zenith angle=0
         ! time units for rate constants : 1/s 
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         phot1( 1,lev)=8.29E-3            !no2
         phot1( 2,lev)=1.87E-2             !no3 -> no
         phot1( 3,lev)=1.703E-1             !no3 -> no2 + o3p
         phot1( 4,lev)=4.57E-4        !o3 -> o3p
         phot1( 5,lev)=3.78E-05        !o3 -> o1d
         phot1( 6,lev)=1.63E-3            !hono
         phot1( 7,lev)=7.53e-06        !h2o2
         phot1( 8,lev)=5.42e-06        !ooh1 (cacm)
         phot1( 9,lev)=5.42e-06        !ooh2 (cacm)
         phot1(10,lev)=7.35e-06        !rpr2301 (cacm)
         phot1(11,lev)=4.63e-05        !hcho -> h2 + co
         phot1(12,lev)=3.02e-05        !hcho -> 2 ho2 + co
         phot1(13,lev)=4.87e-06        !ald
         phot1(14,lev)=7.35e-06        !ald2
         phot1(15,lev)=1.58e-06        !ketl
         phot1(16,lev)=6.40e-06        !ooh6101
         phot1(17,lev)=6.40e-06        !ooh6102
         phot1(18,lev)=5.80e-06        !hac
         phot1(19,lev)=5.60e-05        !ap6101
         phot1(20,lev)=3.50e-05        !ap6102
         phot1(21,lev)=5.40e-06        !nald
         phot1(22,lev)=1.10e-05        !hpald
         phot1(23,lev)=7.35e-06        !pacald
         phot1(24,lev)=1.10e-05        !macr
         phot1(25,lev)=1.42e-04        !rp6301
         phot1(26,lev)=6.40e-06        !ooh6301
         phot1(27,lev)=5.92e-05        !gly -> 0.45 hcho + 1.55 co + 0.80 ho2
         phot1(28,lev)=7.35e-06        !rp7102
         phot1(29,lev)=1.42e-04        !mgly (ch3cocho)
         phot1(30,lev)=7.35e-06        !pina
         phot1(31,lev)=7.35e-06        !nrpa
         phot1(32,lev)=7.35e-06        !edlm
         phot1(33,lev)=1.58e-06        !lmkt
         phot1(34,lev)=7.35e-06        !rpr8401
         phot1(35,lev)=7.35e-06        !rpr8501
         phot1(36,lev)=6.40e-06        !ooh8601
         phot1(37,lev)=6.40e-06        !ooh8602
         phot1(38,lev)=3.0E-5        !isop1co2ooh
         phot1(39,lev)=3.0E-5        !isop3ooh4co
         phot1(40,lev)=6.50e-06        !isop3ooh4oh
         phot1(41,lev)=6.50e-06        !isop1oh2ooh
         phot1(42,lev)=2.50e-04        !MVK3OH4CO
         phot1(43,lev)= 2.50e-04       !MVK3CO4OH
         phot1(44,lev)=2.32E-4        !ISOP1CO4OOH
         phot1(45,lev)=2.2E-4        !ISOP1OOH4CO
         phot1(46,lev)=2.5E-4        !MCRENOL ISOP-POW
         phot1(47,lev)=2.50e-04        !MVKENOL ISOP-POW
         phot1(48,lev)=2.5e-04        !MVK3OOH4CO ISOP-POW
         phot1(49,lev)=2.2E-4        !ISOP1N4CO ISOP-POW
         phot1(50,lev)=2.32e-04        !ISOP1CO4N ISOP-POW
         phot1(51,lev)=2.24E-04        !ISOP3CO4N ISOP-POW
      enddo  ! end level loop
c
c Initialise leaf and soil temperature
c 
      do i=1,itotp1
      do j=2,jtot
      templf(i,j)=t(j+levhtr)-273.15+dtleaf !templf here relates to the 
      ! layer and not the level, but here it makes no difference
      enddo
      enddo

c  contot(j) is accumulated condensate on leaf in layer j
c  frwet: fraction of wet leaves
      do ii=1,levcpy-levhtr
      pint1(ii)=0.    ! no intercepted rain
      contot(ii)=0.
      frwet(ii)=0.
      enddo
      twater=t(nlev-1)
c################################################
c Soil 
c################################################ 
      tausum=0.
      vsum=0.
      tsum=0.

c Calculate soil grids
c 
      read(15,*) (zs(i),i=1,6)             ! ######### soil layers
      write(06,'(/,x,"Soil levels (m):")')
      write(06,'(x,6i7)') (ii,ii=1,6)
      write(06,'(x,6f7.3)') (zs(ii),ii=1,6)
      sl=0.005
      dels=1.0
      zs(1)=-1.0e-04
      eps=0.0
      do 30 i=2,nsoil
      eps=eps+dels
   30 zs(i)=-sl*(eps**2)

      read(15,*) bd,sandfc,siltfc,clayfc ! ######### soil properties
      read(15,*) eta                     ! ######### initial soil moisture
      read(15,*) temps                    ! ######### initial soil temperature

      adum=sandfc*alog(1.025)+siltfc*alog(.026)+clayfc*alog(.001)
      bdum=(sandfc*(alog(1.025))**2+siltfc*(alog(.026))**2+clayfc*(alog(
     &.001))**2-adum*adum)**.5
      dgeom=exp(adum)
      sgeom=exp(bdum)
c Compute air entry potential, pe=j/kg
      pes=-.5*dgeom**(-.5)
      bx=-2.*pes+.2*sgeom
      pe=pes*(bd/1.3)**(.67*bx)
      psolid=bd/2.65
      ws=1.-psolid  ! Only possible for previously used psi


c  calc frac of root system extracting water at each node
      froot(1)=0.
      resrot(1)=0.
      rootsm=0.
      do 700 jz=2,nsoil-1    
      zm=0.5*(z(jz-1)+z(jz))
      zmm=0.5*(z(jz)+z(jz+1))
      froot(jz)=exp(-aroot*zm)-exp(-aroot*zmm)
      resrot(jz)=rroot/froot(jz)
 700  rootsm=rootsm+1./resrot(jz)
      froot(nsoil)=froot(nsoil-1)
      resrot(nsoil)=resrot(nsoil-1)


      tsfc=temps(1)-273.15
      tsurf=temps(1)

      close(15)
c
c Read soil properties
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        open(unit=7, file='./fd2dat/lehmst',  ! ######### adapt path if necessary
     &  status='unknown')
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c 
      read(7,4999) ibt
      read(7,5000) isp,isk,ispt,isl
      read(7,5001) (sep(i),i=1,isp)
      read(7,5001) (spp(i),i=1,isp)
      read(7,5001) (sek(i),i=1,isk)
      read(7,5001) (sk(i),i=1,isk)
      read(7,5001) (sept(i),i=1,ispt)
      read(7,5001) (spt(i),i=1,ispt)
      read(7,5001) (sel(i),i=1,isl)
      read(7,5001) (sxl(i),i=1,isl)
      read(7,7510) por,feldka,albdry,albwet,etamin
      read(7,7870) xkr,trefp,trefk,cp3,rho3
        close(7)
 7510 format(5f12.3)
 7870 format(5f10.3)
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        open(unit=8,file='./fd2dat/lehmko',  ! ######### adapt path
     &  status='unknown')
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      read(8,4999) ibt
      do 515 i=1,isp-1
 515  read(8,5222) cpsi(i,1),cpsi(i,2),cpsi(i,3)
      do 516 i=1,isk-1
 516  read(8,5222) ck(i,1),ck(i,2),ck(i,3)
      do 517 i=1,ispt-1
 517  read(8,5222) cpsit(i,1),cpsit(i,2),cpsit(i,3)
      do 518 i=1,isl-1
 518  read(8,5222) clam(i,1),clam(i,2),clam(i,3)
        close(8)
 4999 format(a4)
 5222 format(3e24.16)
 5000 format(6i3)
 5001 format(6e10.3)
c 
      etamax=por-por/30.
      patm=p(1)
      tsurf=temps(1)
      call icsevu(sep,spp,isp,cpsi,isp-1,eta(1),psianf,1,ier)
      psianf=-10.**psianf
      ha=exp(psianf/(r1*tsurf))
      e=aaa+(bbb*tsurf-ccc)/(tsurf-ddd)
      esat=exp(e)
      qgrnd=.622*esat*ha/(p(1)-.378*esat*ha)
      

      return
      end
