c Most of the radiation routine is from CUPID (John M. Norman, Department of Soil Science,
c      University of Wisconsin-Madison, http://www.soils.wisc.edu/~norman/cupid/outline.html)
c**********************************************************
      subroutine declin
      USE parameters, ONLY:pid180
c  calc declin. of sun in rad and eq. of time in fractions of hours
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    DELNU,         SLONG,       XM
      INTEGER             KDAY
      DOUBLE PRECISION    COSDEC,      COSLAT,      DECL,        DECMAX
      DOUBLE PRECISION    DLONG,       EQTM,        SINDEC,      SINLAT
      DOUBLE PRECISION    TANLAT
      REAL    DELTIM,      TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
c 
c  jan 1,1977 at 1 second after midnight is day 28124.0
c 
      common/astron/eqtm,decl,sindec,cosdec,decmax,sinlat,coslat,
     &tanlat,dlong
      common/timpar/timloc,deltim,month,jday,iyear,icumdy

c      pid180=3.1415926537/180.

      kday=(iyear-1977)*365+icumdy+28123
      xm=(-1.+.9856*kday)*pid180
      delnu=2.*.01674*sin(xm)+1.25*.01674*.01674*sin(2.*xm)
      slong=(-79.8280+.9856479*kday)*pid180+delnu
      decl=asin(decmax*sin(slong))
      sindec=sin(decl)
      cosdec=cos(decl)
      eqtm=9.4564*sin(2.*slong)/cosdec-4.*delnu/pid180
      eqtm=eqtm/60.

      return
      end


c**********************************************************
      subroutine zenith 
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    CRTZEN,      PI,          PID2
      DOUBLE PRECISION    COSDEC,      COSLAT,      DECL,        DECMAX
      DOUBLE PRECISION    DLONG,       EQTM,        SINDEC,      SINLAT
      DOUBLE PRECISION    TANLAT
      DOUBLE PRECISION    SUNAZM,      VIEWAZ,      VIEWZN
      DOUBLE PRECISION    XINTV
      INTEGER             NOAZMV,      NOZENV
      DOUBLE PRECISION    HRANG,       TIMSUN
      DOUBLE PRECISION    COSZEN,      FBEAM1,      RADTOP,      RATIO
      DOUBLE PRECISION    RATIOD,      RATION,      ZENANG
      real    DELTIM,      TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"

      common/deg/sunazm,viewzn(99),viewaz(99),nozenv,noazmv
     &,xintv(20,90)
      common/deg2/ hrang,timsun
      common/astron/eqtm,decl,sindec,cosdec,decmax,sinlat,coslat,
     &tanlat,dlong
      common /rad3/radtop(3),fbeam1(3),coszen,zenang
     &,ratiod,ration,ratio
      common/timpar/timloc,deltim,month,jday,iyear,icumdy

      pid2=3.1415926537/2.
      pi=3.1415926537
      timsun=timloc+eqtm+dlong
      hrang=(timsun-12.)*pid2/6.
      zenang=acos(sinlat*sindec+coslat*cosdec*cos(hrang))
      sunazm=asin(cosdec*sin(hrang)/sin(zenang))
c      sunazm=asin(cosdec*sin(hrang)/sin(zenang))+pi
c pi in above stm was added by Chen, 10/31/89.
c  calc crit zenith angle when sun azimuth is 90 deg
      crtzen=acos(sindec/sinlat)
      if(zenang.gt.crtzen)sunazm=(pi-abs(sunazm))*sunazm
     &/abs(sunazm)
      sunazm = sunazm + pi
 600  continue
      coszen=cos(zenang) 
      return
      end


c**********************************************************
      subroutine radiat(kstrt,coszen,radabv,fbeam,thetas,
     & phis,clump,lprin)
c from CUPID file curadia.f
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      use parameters
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER             INT
!      INTEGER             NLEV
      DOUBLE PRECISION    CLUMP,       COSZEN,      FBEAM,       PHIS
      DOUBLE PRECISION    RADABV,      THETAS
      INTEGER             KSTRT
      DOUBLE PRECISION    A,           AD,          ADIF,        ADIF1
      DOUBLE PRECISION    ADIF2,       ADIR,        ADUM,        AUX
      DOUBLE PRECISION    AUX1,        AUX2,        BEAM
      DOUBLE PRECISION    DAVE,        DOWN,        ESOIL,       EXPDIR
      DOUBLE PRECISION    PRAREA,      PROJ,        QITER,       RNLDIV
      DOUBLE PRECISION    SDN,         SHADE,       SOURCE,      SUNLIT
      DOUBLE PRECISION    SUP,         TBEAM,       TLAY2,       UP
      DOUBLE PRECISION    X,           XXX
      INTEGER             I,           ILPRIN,      IOUT,        IREPT
      INTEGER             ITER,        IW,          J,           JJ
      INTEGER             JJJ,         JJM1,        JJP1,        JM1
      INTEGER             JOUT,        JTOTM1,      K,           L
      INTEGER             LOUT
      DOUBLE PRECISION    DSDUM,       DSTNG,       DSTNG2
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
      DOUBLE PRECISION    FDBTM,       XIUPT
      DOUBLE PRECISION    DELT,        PSILF,       TRAN
      DOUBLE PRECISION    ALAM,        EVAP,        GEVAP
      DOUBLE PRECISION    GHEAT,                    HEAT
!      DOUBLE PRECISION    PI,          PID180,      PID2,        SIGMA
      INTEGER             IWRITE,                   KMAX
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    Z,           ZMID
      INTEGER             JZBM1,       JZBOT,       JZCPY,       JZCPY1
      INTEGER             JZCRIT,      JZSFC,       JZSFM1
      DOUBLE PRECISION    XINT,        XINTZ
      INTEGER             ISPHER,      NALPHA,      NXINTZ
      DOUBLE PRECISION    PARABV
      DOUBLE PRECISION    CILEAF,                   CSLEAF
      DOUBLE PRECISION    PSLEAF,                   RDRK,        RGAS
      DOUBLE PRECISION    EAIR,        PHIH,        PHIM,        REFHTE
!ka      DOUBLE PRECISION    REFHTT,      RELH,        TAIR
      DOUBLE PRECISION    REFHTT,        TAIR	!ka - removed relh
      INTEGER             NLABCY,      NLBCPY
      DOUBLE PRECISION    ALEAF,       EMIS,        EMISOL
      DOUBLE PRECISION    EXPDIF,      RLAYR,       RLEAF
      DOUBLE PRECISION    RSOIL,       TLAYR,       TLEAF
      DOUBLE PRECISION    DSTNET,                   DSTRAD
      DOUBLE PRECISION    FRAREA,                   TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    BMFLX,       D,           RNDIV
      DOUBLE PRECISION    RNET,        RNLAM,       TSFC,        U
      DOUBLE PRECISION    SOURDN,                   SOURUP
      DOUBLE PRECISION    EXTF,        FACTE,       POTVIS
      real    DELTIM,      TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
c
c thetas and phis in above stm were added by Chen, 9/7/89.
c 
c     parameter(mh= 1)
      dimension u(3,99), d(3,99), radabv(3), sup(99), sdn(99), adum(99)
     &,fbeam(3)
     &,tbeam(99),beam(3)
     &,dave(3,99),shade(3,99)
      common/balan/dsdum(3,99),dstng(99),dstng2(99)
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common /rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3)
     &,expdif(99),rlayr(3,30),tlayr(3,30)
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common/rad4/d,u,bmflx(3,99),rnet(99),rndiv(99),tsfc,rnlam(3,99)
      common /rad5/ sourup(3,30),sourdn(3,30)
      common /rad7/ facte(99),extf(99),potvis
      common/leaf1/delt(10,99),psilf,tran(10,99)
      common/leaf2/evap(10,99),gevap(10,99),heat(10,99),gheat(10,99)
     &,alam
c     common/resis1/rhleaf(99),rsleaf(10,99),rsnovp(10,99),
c    & rst(10,99),hpsi
      common/misc1/iwrite(9,99),kmax
c pid2,kmax in /misc1/ were added by Chen, 9/4/89.
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/misc6/xint,xintz(99),nalpha,ispher,nxintz
c xintz, extin. coef. at different zenith averaged over azimuth, and
c nxtintz, dimension of xintz, Chen, 05/25/89.
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
c /misc6/ and /indax3/ were added by Chen, 9/7/89.
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
!ka      common/prof1/tair(99),eair(99),phim,phih,refhtt,refhte,relh(99)
!ka     &,nlabcy,nlbcpy
      common/prof1/tair(99),eair(99),phim,phih,refhtt,refhte
     &,nlabcy,nlbcpy	!ka - removed relh(99) from prof1 common block
      common/misc4/z(99),zmid(99),jzcpy,jzcpy1,jzsfc,jzsfm1,jzbot,jzcrit
     &,jzbm1
c     common/prof2/tn(mh,99),akcpy(99),cpcpy(99),u1(99),q(99),et(99),
c    1en(mh,99),qcond(99),econd(99),tcheck(99),esat(99),qwater(99)
      common/photo3/csleaf(10,99),psleaf(10,99),cileaf(10,99),rgas
     &,rdrk(10,99)
!      parameter(nlev=40)
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c this Common is necessary for the boundary layer model
c (Therm. IR at top of Canopy)
      common /irupdn/ xiupt,fdbtm
      common /par/parabv
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      logical lprin
c
c  loop over wavelength, 1=vis; 2=nearir; 3=thermal
c 
      do 5900 k=kstrt, kmax
c  separate solar and thermal
 545  if(k-3) 550, 4000,4000
 550  continue
c  if fbeam or coszen lt .01, don't calculate beam factors
 610  if(coszen-.01)620,635,635
 620  do 630 j=1,jtot
      d(k,j)=0.
      u(k,j)=0.
      do 625 i=1,itotp1
      dstrad(k,i,j)=0.
 625  frarea(i,j)=0.
      frarea(itotp1,j)=1.
 630  continue
      goto 580
 635  if(fbeam(k)-.01)580,565,565
 580  expdir=0.
      do 590 j=1,jtot
      tbeam(j)=0.
      bmflx(k,j)=0.
      sup(j)=0.
      sdn(j)=0.
 590  continue
 585  if(coszen-.01)5000,640,640
c  calculate direct interception for a layer
 565  continue
c  calculate all beam sources
      tbeam(jtot)=fbeam(k)
      sdn(1)=0.
      do  800 j=2,jtot
      jj=jtot-j+1
      jjp1=jj+1
      tbeam(jj)=1.e-20
      expdir=exp(-clump*xint*df(jj)/coszen) !!!!! hopefully correct jj
      if(tbeam(jjp1).gt.1.e-20)tbeam(jj)=tbeam(jjp1)*expdir
      sup(jjp1)=sourup(k,jj)*tbeam(jjp1)
      sdn(jjp1)=sourdn(k,jj)*tbeam(jjp1)
c     xup=(tbeam(jjp1)-tbeam(jj))*rleaf(k)
c     xdn=(tbeam(jjp1)-tbeam(jj))*tleaf(k)

 800  continue
c 
 640  sup(1)=tbeam(1)*rsoil(k)
c   calc up and do wn fluxes with diffuse only equations and add
c    source terms to do wn terms as they are computed.
        aux1=0.
        do 603 j=1,jtot
        aux2=tbeam(j)
        x=sdn(j)+sup(j)
        aux=(aux2-aux1)*(1.-aleaf(k))
        if(j.eq.1) aux=(aux2-aux1)*rsoil(k)
        aux1=aux2
603     continue
      d(k,jtot)=1.-fbeam(k)
      adum(1)=rsoil(k)
      do  600 j=2,jtot
      jm1=j-1
      tlay2=tlayr(k,jm1)*tlayr(k,jm1)
 600     adum(j)=adum(jm1)*tlay2/(1.-adum(jm1)*rlayr(k,jm1))+
     & rlayr(k,jm1)
      do  700 j=2,jtot
      jj=jtot-j+1
      jjp1=jj+1
      d(k,jj)=d(k,jjp1)*tlayr(k,jj)/(1.-adum(jjp1)*
     &  rlayr(k,jj))+sdn(jjp1)
 700  u(k,jjp1)=adum(jjp1)*d(k,jjp1)+sup(jjp1)
      u(k,1)=rsoil(k)*d(k,1)+sup(1)
c  calculate total diffuse up and down considering beam sources.
 1045 iter=0
 900  irept=0
      iter=iter+1
      do  1000 j=2,jtot
      jj=jtot-j+1
      jjp1=jj+1
      down=tlayr(k,jj)*d(k,jjp1)+u(k,jj)*rlayr(k,jj)+sdn(jjp1)
      if(abs(down-d(k,jj))-.0001)1000, 1000, 950
 950  irept=1
 1000 d(k,jj)=down
      u(k,1)=(d(k,1)+tbeam(1))*rsoil(k)
      do  1200 jj=2,jtot
      jjm1=jj-1
      up=rlayr(k,jjm1)*d(k,jj)+u(k,jjm1)*tlayr(k,jjm1)+sup(jj)
      if(abs(up-u(k,jj))-.0001) 1100,1100,1050
 1050 irept=1
 1100 u(k,jj)=up
 1200 continue
        do 1203 j=2,jtot
        a=d(k,j)-u(k,j)+tbeam(j)-(d(k,j-1)-u(k,j-1)+tbeam(j-1))
        adif=(d(k,j)+u(k,j-1))*(1-expdif(j-1))*aleaf(k)
        adir=(tbeam(j)-tbeam(j-1))*aleaf(k)
        ad=adif+adir
1203    continue
 1260 if(irept) 900, 2000, 900
c  calculate fraction of leaf area in class i that is associated with
c  light incident in class i.
c  frarea(i,j) sums to 1. over leaf angl classes including shaded
 2000 jtotm1=jtot-1
      do  3000 j=2,jtot

      jm1=j-1
      if(tbeam(jtot)-.0001)2050,2100,2100
 2050 sunlit=0.
      goto 2200
 2100 sunlit=coszen*(tbeam(j)-tbeam(jm1))/(
     & tbeam(jtot)*df(jm1))/xint   ! correct j?

 2200 do  2500 i=1,itot
       frarea(i,j)=sunlit*distls(i)
 2500 continue
      frarea(itotp1,j)=1-sunlit
 3000 continue

      beam(k)= fbeam(k)*radabv(k)/coszen
      d(k,1)=d(k,1)*radabv(k)
      u(k,1)=u(k,1)*radabv(k)
      bmflx(k,1)=tbeam(1)*radabv(k)
      do 2600 j=2,jtot
      jm1=j-1
c  convert to flux densities
      u(k,j)=u(k,j)*radabv(k)
      d(k,j)=d(k,j)*radabv(k)
 2600 bmflx(k,j)=tbeam(j)*radabv(k)

      do 3100 j=2,jtot

      jm1=j-1
      dsdum(k,j)=0.
c  next line calculates diffuse rad above layer
c      dstrad(k,itotp1,j)=(u(k,jm1)+d(k,j))*(1.-expdif)/df
c above stm was commented because it did not include the multiple
c scattering within the layer j. the absorbed diffuse radiation by
c layer j is (u(k,jm1)+d(k,j))*(1-tlayer-rlayer), and the obsorbed
c direct radiation is (bmflx(k,j)-bmflx(jm1))*aleaf(k). due to the
c multiple scattering within the layer, the total absorption caused
c by direct radiation is (bmflx(k,j)-bmflx(k,jm1))-s(k,j). therefore
c (bmflx(k,j)-bmflx(k,jm1)-s(k,j)-(bmflx(k,j)-bmflx(k,jm1))*aleaf(k)
c =(bmflx(k,j)-bmflx(k,jm1))*(1-aleaf(k))-s(k,j) is the additional
c contribution to absorption, which can be included for simplicity
c in diffuse absorption as follows: (ChEn, 10/3/89)
        adif1=(u(k,jm1)+d(k,j))*(1.-tlayr(k,jm1)-rlayr(k,jm1))
        adif2=(bmflx(k,j)-bmflx(k,jm1))*(1.-aleaf(k))
     &        -(sdn(j)+sup(j))*radabv(k)
        dstrad(k,itotp1,j)=(adif1+adif2)/df(jm1)/aleaf(k)

c  next line calculated average diffuse rad at midpoint of layer.
cc    dstrad(k,itotp1,j)=(d(k,j)+d(k,jm1)+u(k,jm1)+u(k,j))*.5
cc   &*(1.-expdif)/df
      dsdum(k,j)=dsdum(k,j)+dstrad(k,itotp1,j)*frarea(itotp1,j)*df(jm1)
     &*aleaf(k)
        proj=0.
      do  2700 i=1,itot

      dstrad(k,i,j)=dstrad(k,itotp1,j)+ beam(k)*ct(i)
 2700 dsdum(k,j)=dsdum(k,j)+dstrad(k,i,j)*frarea(i,j)*df(jm1)*aleaf(k)
      rnldiv=d(k,j)-u(k,j)+bmflx(k,j)-(d(k,jm1)-u(k,jm1)+bmflx(k,jm1))

 2701 format(' solar ',2i3,2f7.2)
 3100 continue

      goto 5000
c
c  thermal
c 
 4000 continue
c  treat thermal wavelengths.
c  compute thermal source terms
c  check to see if solar has been executed(day or night)
      sdn(1)=0.
      do 4009 j=1,jtot
 4009 bmflx(3,j)=0.
      if(kstrt-3)4020,4010,4010
c  check if night or recycling of leaf energy balance
 4010 if(coszen-.01)4012,4012,4020
 4012 do 4015j=2,jtot
      sdn(j)=emis*sigma*(templf(itotp1,j)+273.)**4
      do 4014 i=1,itotp1
      frarea(i,j)=0.
      dstrad(3,i,j)=sdn(j)
      dstrad(1,i,j)=0.
 4014 dstrad(2,i,j)=0.
      frarea(itotp1,j)=1.
      d(1,j)=0.0
      d(2,j)=0.0
      tbeam(j)=0.0
      bmflx(1,j)=0.0
      bmflx(2,j)=0.0
      sdn(j)=sdn(j)*(1.-expdif(j-1))
 4015 sup(j)=sdn(j)
      d(1,1)=0.0
      d(2,1)=0.0
      tbeam(1)=0.0
      bmflx(1,1)=0.0
      bmflx(2,1)=0.0
      expdir=0.
      goto 4200
 4020 do  4100 j=2,jtot
      source=0.
      do  4050 i=1,itotp1

      dstrad(3,i,j)=emis*sigma*(templf(i,j)+273.)**4
      source=source+frarea(i,j)*dstrad(3,i,j)
 4050 continue
      sdn(j)=source*(1.-expdif(j-1))
 4100 sup(j)=source*(1.-expdif(j-1))
c  thermal boundary conditions
 4200 continue
c  tsfc is soil sfc temp for each iteration.  its calc'd in main
c    program after profl2 is called.
      esoil =emisol*sigma*(tsfc+273.)**4
      d(k,jtot)=radabv(3)
c  compute thermal layer properties
c     tlayer=expdif  !war schon immer kommentiert
c     rlayer=(1.-expdif)*(1.-emis) "
c  assume thermal reflectance =0 and call u+d
      do  4500 j=2,jtot
      jm1=j-1
      jj=jtot-j+1
      jjp1=jj+1
 4500 d(k,jj)=tlayr(k,jj)*d(k,jjp1)+sdn(jjp1)
      sup(1)=d(k,1)*(1.-emisol)
      u(k,1)=esoil+sup(1)
      do 4600j=2,jtot
      jm1=j-1
 4600 u(k,j)=tlayr(k,jm1)*u(k,jm1)+sup(j)
c  iterate 2 times with refl=1-emis
 4580 do  4900 jjj=1,2
      do  4800 j=2,jtot
      jm1=j-1
      jj=jtot-j+1
      jjp1=jj+1
 4800 d(k,jj)=tlayr(k,jj)*d(k,jjp1)+u(k,jj)*rlayr(k,jj)+sdn(jjp1)
      sup(1)=d(k,1)*(1.-emisol)
      u(k,1)=esoil+sup(1)
      do 4810 j=2,jtot
      jm1=j-1
 4810 u(k,j)=rlayr(k,jm1)*d(k,j)+u(k,jm1)*tlayr(k,jm1)+sup(j)
 4900 continue
c
c  sum up thermal contribution to dstnet to be sure it matches u and d
c    thermal. dstnet is calc'd below.
      do 4910 j=2,jtot
      dsdum(k,j)=0.
      jm1=j-1
      do 4905 i=1,itotp1
      xxx=(1.-expdif(jm1))/df(jm1)
      dsdum(k,j)=dsdum(k,j)-2.*dstrad(k,i,j)*xxx*df(jm1)*frarea(i,j)
 4905 continue
      dsdum(k,j)=dsdum(k,j)+(d(k,j)+u(k,jm1))*(1.-tlayr(k,jm1)
     & -rlayr(k,jm1))
      rnldiv=d(k,j)-u(k,j)-(d(k,jm1)-u(k,jm1))

 4915 format(' thermal ',2i3,4f8.2)
 4910 continue
c
c                                                  label file write 5
c 
 5000 continue
      if(k-2)5005,5006,5007
 5005 if(ilprin.eq.5)write(06,5010)k,icumdy
 5010 format(' 150',i3,i4,4x,/14x,'  iter  expdir expdif tlayer',
     &' rlayer  clump ')
      goto 5009
 5006 if(ilprin.eq.5)write(06,5011)k,icumdy
 5011 format(' 150',i3,i3,4x,/14x,'  iter  expdir expdif tlayer',
     &' rlayer  clump ')
      goto 5009
 5007 if(ilprin.eq.5)write(06,5012)k,icumdy
 5012 format(' 150',i3,i4,4x,/14x,'  iter  expdir expdif tlayer',
     &' rlayer  clump ')
 5009 continue
      qiter=iter
      goto 5024
      if(lprin.and.k.eq.1)goto 5014
      if(lprin.and.k.eq.2)goto 5014
      if(lprin.and.k.eq.3)goto 5014
      goto 5024
 5014 write(06,5015)k,icumdy,qiter,expdir,expdif(1),
     & tlayr(k,1),rlayr(k,1)
     &,clump
 5015 format(' 250',i3,i4,4x,9f7.2)
      write(06,5025)
 5025 format(/13x,' clai    adum ',
     &'   d      u    tbeam  bmflx  rn123   rnet ')
 5024 continue
c 
 5030 do 5040 j=1,jtot
      jout=jtot-j+1
      clai(jout)=0.
      do jj=jout,jtot-1
      clai(jout)=clai(jout)+df(jj)
      enddo
c --- altered !!
      if(k.eq.1.or.k.eq.2) rnlam(k,jout)=d(k,jout)-u(k,jout)
     & +bmflx(k,jout)
      if(k.eq.3) rnlam(k,jout)=d(k,jout)-u(k,jout)
      rnet(jout)=rnlam(1,jout)+rnlam(2,jout)+rnlam(3,jout)
      goto 5040
      if(lprin.and.k.eq.1) goto 5038
      if(lprin.and.k.eq.2) goto 5038
      if(lprin.and.k.eq.3) goto 5038
      goto 5040
 5038 continue
      if(k.eq.3) then
      write(06,5026)k,icumdy,jout,clai(jout)
     &,adum(jout),d(k,jout),u(k,jout),tbeam(jout),bmflx(k,jout),rnlam(k,
     2jout),rnet(jout)
      else
      write(06,5026)k,icumdy,jout,clai(jout)
     &,adum(jout),d(k,jout),u(k,jout),tbeam(jout),bmflx(k,jout),rnlam(k,
     2jout)
      endif
 5026 format(i3,i4,2x,i2,  9f7.2)
 5040 continue
      parabv=d(1,jtot)+bmflx(1,jtot)
c 
      if(k.eq.1) then
       do j=1,jtot
        if(radabv(1).le.0.) then
         facte(j)=0.
        else
         facte(j)=rnlam(1,j)/radabv(1)
        endif
       enddo
      endif

c ##########################################################
c Photolysis profile follows PAR profile within the canopy =>
c need to caclulate PAR extinction for this level
c ###############################################################
      do j=1,nlev
       extf(j)=1.
       if(j.le.levhtr) then
       extf(j)=facte(1)
       else
       if((j-levhtr).lt.jtot) extf(j)=facte(j-levhtr+1) 
       endif
      enddo
   
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      xiupt=u(3,jtot)
c  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      do 5970 j=2,jtot
      jm1=j-1

 5970 rndiv(j)=rnet(j)-rnet(jm1)
      rndiv(1)=rnet(1)
 5900 continue 
c 
 5200 do  6500 j=2,jtot
      jm1=j-1
      dstng(j)=0.
c-----------------------------------
      do  6000 i=1,itotp1
c  calc thermal dif at layer midpoint
      dstnet(i,j)=dstrad(1,i,j)*aleaf(1)+dstrad(2,i,j)*aleaf(2)+
     & (d(3,j)+u(3,jm1))*(1.-tlayr(3,jm1)-rlayr(3,jm1))/df(jm1) 
     & -2.*dstrad(3,i,j)*(1.-expdif(j-1))/df(jm1)
      dstng(j)=dstng(j)+dstnet(i,j)*df(jm1)*frarea(i,j)
 6000 continue
c-----------------------------------
      
      dstng2(j)=dsdum(1,j)+dsdum(2,j)+dsdum(3,j)

 6002 format(' total ',2i3,3f8.2)
 6500 continue
c label file write 8
c
c     goto 6550
c 
       iw=0
      jtotm1=jtot-1
      do 5100 iout=1,itotp1
c Class itot+1 is independent of direction and contains the shaded
c Fraction (d.h. nachts 100%)
      do 5100 j=1,jtotm1
      jout=jtot-j+1
      prarea=frarea(iout,jout)*100.
      if(iw.eq.1)write(06,5051)icumdy,iout,jout,dstnet(iout
     &,jout),dstrad(1,iout,jout),dstrad(2,iout,jout),dstrad(3,iout,jout)
     &,prarea,psleaf(iout,jout)
 5051 format(' 2801',i3,i2,i2,i2,  9f7.2)

 5100 continue
 5150 continue
c  calc average flux density downward for shaded leaves. this
c    is for simplified model using only sunlit and shaded leaves.
c    must ave all d(k,j) above lai of interest
      do 5190 k=1,2
      do 5190 j=1,jtot
 5190 dave(k,j)=0.
      do 5350 k=1,2
      jtotm1=jtot-1
      do 5350 j=1,jtotm1
      jout=jtot-j+1
c     shade(k,jout)=d(k,jtot)*exp(-.5*clai(jout)**.7)+radabv(k)*fbeam(k)
c    1*.1*(1.-.1*clai(jout))*exp(-coszen)
c above stm was commented and following three were added by Chen, 9/7/89.
        if(coszen.lt.0.01) goto 5355
c above stm was added by chen, 01/23/90.
      x=0.5
      if(ispher.ne.1) x=xintz(int(thetas/(pi/2.)*90)+1)
      shade(k,jout)=d(k,jtot)*exp(-clump*x*clai(jout)**.7)
     &   +radabv(k)*fbeam(k)*.1*(1.-.1*clai(jout))*exp(-coszen)
 5355   continue
c above stm was added by chen, 01/23/90.
      do 5300l=1,j
      lout=jtot-l+1
 5300 dave(k,jout)=dave(k,jout)+(d(k,lout)+d(k,lout-1))/2.
      dave(k,jout)=dave(k,jout)/j
 5350 continue
      do 5375j=1,jtotm1
      jout=jtot-j+1

 5375 continue
 6550 continue
 7700 continue

      return
      end


c**********************************************************
      subroutine radin4
c based on CUPID file curadin.f
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    A,           AIRMAS,      AXLOG,       EVAL
      DOUBLE PRECISION    FB1,         FB2,         PID180,      PID2
      DOUBLE PRECISION    POTBM1,      POTBM2,      POTDIF,      POTNIR
      DOUBLE PRECISION    POTVIS,      RATIOX,      U,           WATABS
      INTEGER             K
      DOUBLE PRECISION    COSZEN,      FBEAM1,      RADABV,      RATIO
      DOUBLE PRECISION    RATIOD,      RATION,      ZENANG
      real    DELTIM,      TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
c 
c  averad is input as total cal/cm sqd/day
c  compute potential flux for day,pcpt water=1 cm,alpha=1
      common /rad3/radabv(3),fbeam1(3),coszen,zenang
     &,ratiod,ration,ratio
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
c 
      do 10 k=1,2
      radabv(k)=0.0
 10   fbeam1(k)=0.
      pid2=3.1415926537/2.
      pid180=3.1415926737/180.
c  convert ly/hr to w/sq m.
      if(coszen.lt.0.01)goto 900

c  correct for curvature of atmos in airmas
 50   airmas=(sqrt(coszen**2+.0025)-coszen)/.00125
c  correct for refraction(good to 89.5 deg.)
      airmas=airmas-2.8/(90.-zenang/pid180)**2
      potbm1=600.*exp(-.160*airmas)
      potvis=(potbm1+(600.-potbm1)*.4)*coszen
      potdif=(600.-potbm1)*.4*coszen
      u=1/coszen
      axlog=dlog10(u)
      a=10**(-1.195+.4459*axlog-.0345*axlog*axlog)
      watabs=1320.*a
      potbm2=720.*exp(-.05*airmas)-watabs
      if(potbm2.lt.0.)potbm2=0.
 60   eval=(720.-potbm2-watabs)*.54*coszen
 90   potnir=eval+potbm2*coszen
c 
 200  ratio=ratiod
      radabv(1)=potvis*ratio
      radabv(2)=potnir*ratio
 300  fb1=potbm1*coszen/potvis
      fb2=potbm2*coszen/potnir

      ratiox=ratio
      if(ratio.gt..9)ratiox=.9
      fbeam1(1)=fb1*(1.-((.9-ratiox)/.7)**.6667)
      if(ratio.gt.0.88)ratiox=.88
      fbeam1(2)=fb2*(1.-((.88-ratiox)/.68)**.6667)
      if(fbeam1(1).lt.0.)fbeam1(1)=0.
      if(fbeam1(2).lt.0.)fbeam1(2)=0.
      if(fbeam1(1).gt.fb1)fbeam1(1)=fb1
      if(fbeam1(2).gt.fb2)fbeam1(2)=fb2
      goto 1000
 900  ratio=ration
 1000 continue

      return
      end


c**********************************************************
      subroutine skyir(temp,vpaira)
c  calc thermal rad from sky with Brutsaert eq.
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    TEMP,        VPAIRA
      DOUBLE PRECISION    ESKY,        SIGMA
      DOUBLE PRECISION    COSZEN,      FBEAM1,      RADTOP,      RATIO
      DOUBLE PRECISION    RATIOD,      RATION,      ZENANG
      !shc end of adding declarations to allow "implicit none"
      common/rad3/radtop(3),fbeam1(3),coszen,zenang
     1,ratiod,ration,ratio
c
c  use ratio to get weighted average of clear sky and 'clouds'.
      sigma=5.67e-08
      esky=1.24*(vpaira/temp)**(1./7.)
      radtop(3)=sigma*(temp**4)*(esky*ration+1.-ration)
c ratio changed to ration (since an overestimation of the
c two solar parts in radin4 needs not to result in an overestimation
c of the thermal)
      return
      end


c**********************************************************
      subroutine dstlit(thets,phis,kmax,rlfdir,tlfdir)
c  subroutine to calculate the distribution of leaf-normal 
c  to sun angles
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    PHIS,        RLFDIR,      THETS
      DOUBLE PRECISION    TLFDIR
      INTEGER             KMAX
      DOUBLE PRECISION    ALPHA,       ANGL,        BETA,        CA
      DOUBLE PRECISION    CT,          FRADD,       FS,          PI
      DOUBLE PRECISION    PID180,      PID2,        SA,          ST
      DOUBLE PRECISION    SUM,         SUM1,        WBETA
      INTEGER             IALPHA,      IANGL,       IANGLE,      IBETA
      INTEGER             K
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
      DOUBLE PRECISION    CLAI,        COST,        DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    XINT,        XINTZ
      INTEGER             ISPHER,      NALPHA,      NXINTZ
      DOUBLE PRECISION    RLFHEM,      TLFHEM
      !shc end of adding declarations to allow "implicit none"
c
c  for a given direction of the sun or given hour. distls(iangle) is
c  the variable where iangle has the same class intervals as fr(ialpha)
c 
      dimension rlfdir(3),tlfdir(3),ca(9),sa(9)
      common/rad8/rlfhem(3,9),tlfhem(3,9)
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
      common/misc2/fr(99),cost(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
c original ct(99) was changed to cost(99) by Chen, 7/5/89.
      common/misc6/xint,xintz(99),nalpha,ispher,nxintz
c xintz, extin. coef. at different zenith averaged over azimuth, and
c nxtintz, dimension of xintz, Chen, 05/25/89.
c 
      pi=3.1415926537
      pid180=pi/180.
      pid2=pi/2.

      do 40 ialpha=1,itot
 40   distls(ialpha)=0.
        wbeta=2.*pi/nbeta

      do 100 ialpha=1,itot
      alpha=pid2*(ialpha-1)/(itot) + pi/(4.*itot)
      sa(ialpha)=sin(alpha)
      ct=cos(thets)
      st=sin(thets)
      ca(ialpha)=cos(alpha)
c      fradd=fr(ialpha)/180
c      do 100 ibeta=1,20
c      beta=(1.+(ibeta-1)*2.)*pid180 - pi !(5/18/89, chen)
c above three were commented and following three added by Chen,5/24/89
       do 100 ibeta=1,nbeta
        beta=wbeta*(ibeta-.5)
        fradd=fr(ialpha)*fraz(ibeta)
      fs=ct*ca(ialpha) + st*sa(ialpha)*cos(beta-phis)
      angl=acos(abs(fs))
      iangl=angl*itot/pid2 + 1.
c 
 50   format(4f4.1,i3,f5.2,9f6.3) !shc
      if(iangl.lt.1.or.iangl.gt.itot)goto 80
      goto 95
 80   write(6,90)iangl
 90   format(' trouble in sub dstlit iangl= ',i6)
      stop
 95   distls(iangl)=distls(iangl) + fradd

100   continue
        sum=0.
        do 105 iangle=1,itot
          sum=sum+cost(iangle)*distls(iangle)
105     continue
        xint=sum

 51   format(2(f7.2,1x),9(f6.3,1x)) !shc
c
c  compute a mean weighted sunlit leaf reflectance
      do 400 k=1,kmax
      sum1=0.
      sum=0.
c     sum2=0.
c     sum3=0.
c     sum4=0.
c above three statements were commented by Chen, 8/31/89.
      do 300 iangl=1,itot
      angl=pid2*(iangl-1)/itot + pid2/(2.*itot)
      sum=sum+distls(iangl)*rlfhem(k,iangl)
      sum1=sum1+distls(iangl)*tlfhem(k,iangl)
c  calc sums to get weighted leaf refl and trans for diffuse inc rad.
c     sum2=sum2+sa(iangl)*ca(iangl)
c     sum3=sum3+rlfhem(k,iangl)*sa(iangl)*ca(iangl)
c     sum4=sum4+tlfhem(k,iangl)*sa(iangl)*ca(iangl)
c above three statements were commented by Chen, 8/31/89.
 300  continue
      rlfdir(k)=sum
      tlfdir(k)=sum1
c     rlfdif(k)=sum3/sum2
c     tlfdif(k)=sum4/sum2
c above two statements were commented by Chen, 8/31/89, since
c rlfdif(k) & tlfdif(k) were calculated in subroutine inplnt.
 400  continue
      return
      end


c**********************************************************
      subroutine simpsn(theta,phi,gr,nalpha,xint)
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    GR,          PHI,         THETA,       XINT
      INTEGER             NALPHA
      DOUBLE PRECISION    ALPHA,       COSDEL,      CT,          PI
      DOUBLE PRECISION    PID180,      PID2,        ST,          SUM
      INTEGER             K
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
      !shc end of adding declarations to allow "implicit none"
c phi in above and following common were added by Chen, 05/24/89.
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
      dimension gr(99)
c  xint  extinct coeff-nel is in first part of array and nthet is
c          in last part of array and tot size of array is nelpth.
c  nalpha  no of leaf angle classes
c
c   program to calc integrl of ga*cosd (simpson rule)
c     write(11,7)
 7    format(' ','*extinction data')
c
c     nthet=num. of sun angles
c     dt= delta theta
c     nalpha=num. of leaf angle classes
c     da= delta alpha
c
      pi=4.0*atan(1.0e0)
      pid180=pi/180.0
      pid2=pi/2.0
c     nalpha must be odd
      xint=0.
      st=sin(theta)
      ct=cos(theta)
      sum=0.
      do  1001 k=1,nalpha
      alpha=(k-.5)*pid2/nalpha
c      call kernal(cosdel,theta,alpha)
c  above was commented and following added by Chen, 7/9/89.
c ----------------------------------------------------------------------
      call kernal(cosdel,theta,phi,alpha)
c ----------------------------------------------------------------------
      sum=sum+cosdel*gr(k)

 1001 continue
      xint=sum*pid2/nalpha
c  for area fractions weighted by leaf angle distrib,gr(alpha),
c    should not mult by dalpha because integral over gr(alpha)dalpha
c    in denominator will cancel dalpha.
      xint=sum
 1000 continue

      return
      end


c**********************************************************
      subroutine kernal(cosdel,theta,phi,alpha)
c  subroutine to calculate cosine delta
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    ALPHA,       COSDEL,      PHI,         THETA
      DOUBLE PRECISION    BETA,        FS,          PI,          WBETA
      INTEGER             IBETA
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
      !shc end of adding declarations to allow "implicit none"
c  theta=zenith angle of sun or viewer
c  phi =azimuth angle of sun or viewer, south = 0., east = -90.
c  alpha=leaf elevation angle from horizontal
c  beta =azimuthal direction leaf normal points
c  fr =fraction of leaf area proj either in dir of sun or viewer.
c  fraz is the fraction of leaf area as a func of azimuth
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
      pi=4.*atan(1.)
      wbeta=2.*pi/nbeta
      cosdel=0.
      do 1000ibeta=1,nbeta
      beta=(ibeta-.5)*wbeta
      fs=cos(theta)*cos(alpha)+sin(theta)*sin(alpha)*cos(phi-beta)
      cosdel=cosdel+abs(fs)*fraz(ibeta)
 1000 continue
      return
      end


c**********************************************************
      subroutine lad(meu,neu,fr,itot,gamrat,frdeg)
c this subroutine calculates the leaf angle distribution
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER             INT
      INTEGER             ITOT
      DOUBLE PRECISION    AUX,         DANG,        DEG
      DOUBLE PRECISION    FRTEMP
      DOUBLE PRECISION    PI,          PID180,      SUM
      DOUBLE PRECISION    TLD90,       X
      INTEGER             I,           IBETA,       IDEG,        II
      INTEGER             J,           JBETA0,      K,           KK
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
      !shc end of adding declarations to allow "implicit none"
c 
c include distribution both in zenith and in azimuth, Chen, 05/24/89.
c leaf angle distributions for various canopy types:
c   canopy type       f(theta)         theta     meu     neu
c   planophile   2(1+cos(2*theta))/pi  26.76   2.770   1.172
c   erectophile  2(1-cos(2*theta))/pi  63.24   1.172   2.770
c   plagiophile  2(1-cos(4*theta))/pi  45.00   3.326   3.326
c   extremophile 2(1+cos(4*theta))/pi  45.00   0.433   0.433
c   uniform      2/pi                  45.00   1.000   1.000
c   spherical    sin(theta)            57.30   1.101   1.930
c      source: naren goel & do n strebel (1984)  agron j. 76:800-802
c 
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
c above common is added by Chen, 05/23/89, and beta0 is in degrees.
      double precision meu,neu,meuneu,gamrat,lang,fr(99),
     &frdum(91),frdeg(91),frcls(9)
c
c  beta-distribution
c
c      bdistr(gamrat,tld90,neu,meu)=gamrat*(tld90**(neu-1.))*
c     &       ((1.-tld90)**(meu-1.))/(360.*90.)
c above was commented by Chen, 7/5/89.
      pi=4.*atan(1.)
      pid180=pi/180.
      meuneu=meu+neu
c the next statement is used to avoid overflows in the calculation
c of gamma functions for large values of meu and neu.
c      if(meu.gt.35..or.neu.gt.35..or.meuneu.gt.35.) goto 1
c      call gamma(meu,gam1)
c      call gamma(neu,gam2)
c      call gamma(meuneu,gam12)
c      goto 4
c 1    gamrat=((meu*neu)/(2.*pi*meuneu))**0.5
c      gamrat=gamrat*((meuneu/meu)**meu)*(meuneu/neu)**neu
c      goto 5
c 4    gamrat=gam12/(gam1*gam2)
c above 9 lines were commented by Chen, 7/5/89.
      sum = 0.
      dang=90./itot
 5    do  10 i = 1,itot
         lang=(dang*(i-1)+dang/2)/90.
c      fr(i)=(gamrat*(lang**(neu-1))*((1-lang)**(meu-1))/90.)*dang
c above was commented and following added by Chen, 7/5/89.
        fr(i)=lang**(neu-1)*(1-lang)**(meu-1)
      sum = sum + fr(i)
 10   continue
      do  20 i = 1,itot
 20      fr(i)=fr(i)/sum
       write(06,'(x,a8,10f8.3)') 'fr(i): ',(fr(i),i=1,itot)
c ****** 7/7/86 revision  ******
c **  removed section that calculates frspec - it's in a separate
c **  subroutine called "fractn"
c ************************************
c  calculates leaf angle distribution for specular contribution at
c 10 degree increments
c
c     sum = 0
c     do  50 i = 1,9
c        lang=(10.*(i-1)+5)/90.
c        frdum(i) = gamrat*(lang**(neu-1))*((1-lang)**(meu-1))/9.
c     sum=sum+frdum(i)
c 50  continue
c     do  70 i=1,9
c        frdum(i) = frdum(i)/sum
c        do  70 j = 1,36
c  frspec is fraction of leaf area in single azimuth and zenith angle
c    class that contributes to specular reflection.
c           frspec(i,j)=frdum(i)/36.
c 70  continue
c
c  calculate for every 10 deg view zenith angle class
      do  75 i=1,91
         frdum(i)=0.
  75  continue
c for deg classes 1 through 91
      do  200 ideg=1,90
         deg=float(ideg)
         tld90=(deg-0.5)/90.
c         frdum(ideg)=bdistr(gamrat,tld90,neu,meu)
c above was commented and following added by Chen, 7/5/89.
         frdum(ideg)=tld90**(neu-1)*(1-tld90)**(meu-1)
  200 continue
c
c  normalize fraction of leaves
c
      sum=0.
      do  250 i=1,91
         sum=sum+frdum(i)
 250  continue
      do  275 i=1,91
         frdeg(i)=frdum(i)/sum
 275  continue
c  put into 9 classes for screen display
      do  300 i=1,9
         frcls(i)=0.
 300  continue
      do  305 i=2,10
 305     frcls(1)=frcls(1)+frdeg(i)
      do  310 i=2,9
         j=i-1
         frtemp=0.
         do  320 k=1,10
            kk=(k+j*10)
            frtemp=frtemp+frdeg(kk)
 320     continue
         frcls(i)=frtemp
 310  continue
c
c following 12 statments for calculating leaf azimuthal angle distrb.
c were added by Chen, 05/23/89.
c 
      jbeta0=int(beta0/360.*nbeta)
      sum=0.
      do  410 ibeta=1,nbeta
      x=(ibeta-.5)/nbeta
      aux=x**(xneuaz-1)*(1-x)**(xmeuaz-1)
      ii=ibeta+jbeta0
      if (ii.gt.nbeta) ii=ii-nbeta
      fraz(ii)=aux
      sum=sum+aux
 410  continue
      do  420 ii=1,nbeta
      fraz(ii)=fraz(ii)/sum
 420  continue
c 
c following 13 statements for calculating leaf azimuthal angle distrb.
c in each degree class were added by Chen, 07/11/89.
      jbeta0=int(beta0)
      sum=0.
      do  411 ibeta=1,360
      x=(ibeta-.5)/360
      aux=x**(xneuaz-1)*(1-x)**(xmeuaz-1)
      ii=ibeta+jbeta0
      if (ii.gt.360) ii=ii-nbeta
      fradeg(ii)=aux
      sum=sum+aux
 411  continue
      do  421 ii=1,360
      fradeg(ii)=fradeg(ii)/sum
 421  continue
c 
      return
      end


c**********************************************************
      subroutine difint(clump,sunazm)
c  calculate diffuse interception for layer
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      USE parameters
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    CLUMP,       SUNAZM
      DOUBLE PRECISION    A,           ADUM,        CA,          DFDCA
      DOUBLE PRECISION    DIFSUB,      DIRSUB,      DOWN
      DOUBLE PRECISION    RSUBL,       SA,          SDNSUB
      DOUBLE PRECISION    SUPSUB,                   TLAY2
      DOUBLE PRECISION    TSUBL,       UP,          WZ,          X
      INTEGER             II,          IREPT,       ITER,        J
      INTEGER             JJ,          JJP1,        JL,          JM1
      INTEGER             K,           LAYSP1,      LAYSUB
      DOUBLE PRECISION    HCPY,        HTR,         TLAI,        ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    FRDEG,       GAMRAT,      XMEU,        XNEU
      INTEGER             IANGOT,      KPNTSP
      DOUBLE PRECISION    BETA0,       FRADEG,      FRAZ,        XMEUAZ
      DOUBLE PRECISION    XNEUAZ
      INTEGER             NBETA
!      DOUBLE PRECISION    PI,          PID180,      PID2,        SIGMA
      INTEGER             IWRITE,                   KMAX
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    XINT,        XINTZ
      INTEGER             ISPHER,      NALPHA,      NXINTZ
      DOUBLE PRECISION    ALEAF,       EMIS,        EMISOL
      DOUBLE PRECISION    EXPDIF,      RLAYR,       RLEAF
      DOUBLE PRECISION    RSOIL,       TLAYR,       TLEAF
      DOUBLE PRECISION    COSZEN,      FBEAM1,      RADTOP,      RATIO
      DOUBLE PRECISION    RATIOD,      RATION,      ZENANG
      DOUBLE PRECISION    BMFLX,       D,           RNDIV
      DOUBLE PRECISION    RNET,        RNLAM,       TSFC,        U
      DOUBLE PRECISION    SOURDN,                   SOURUP
      DOUBLE PRECISION    RLFDIF,      RLFDIR,      TLFDIF
      DOUBLE PRECISION    TLFDIR
      !shc end of adding declarations to allow "implicit none"
c
c above two statements were added by Chen, 8/31/89.
      dimension sdnsub(9,20),supsub(9,20),rsubl(9),tsubl(9),adum(99)
      common /rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3)
     &,expdif(99),rlayr(3,30),tlayr(3,30)
      common /rad3/radtop(3),fbeam1(3),coszen,zenang
     &,ratiod,ration,ratio
      common /rad4/d(3,99),u(3,99),bmflx(3,99),rnet(99),rndiv(99),tsfc
     &,rnlam(3,99)
      common /rad5/ sourup(3,30),sourdn(3,30)
      common/rad6/rlfdif(3),tlfdif(3),rlfdir(3),tlfdir(3)
      common/indax2/xmeu,xneu,gamrat,frdeg(91),iangot(9),kpntsp
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
      common/misc1/iwrite(9,99),kmax
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/misc6/xint,xintz(99),nalpha,ispher,nxintz
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough
c xintz, extin. coef. at different zenith averaged over azimuth, and
c nxtintz, dimension of xintz, Chen, 05/25/89.
c
c    calc layer trans and refl
c    for thick layers using layer equations and dividing df layer into
c    10 sub-layers and solving for diffuse radiation.
c    should use good precision in diffuse integral so errors do nt
c    compound. thus use 90 intervals because 9 intervals accurate
c    to about .005 and this error becomes .03 in layer non-interception
c    factor because of 10th power from 10 sublayers.
      laysub=10
      laysp1=laysub+1
        wz=pi/2/nxintz
      do 600 jl=1,jtot-1  !neu df hoehenabh.
      x=0.
        do 500 ii=1,nxintz
      a=wz*(ii-.5)
      ca=cos(a)
      sa=sin(a)
      dfdca=df(jl)/ca !!!
        x=x+ca*sa*exp(-clump*xintz(ii)*dfdca)
c parameter clump in above statement was added by Chen, 9/4/89.
 500  continue
        expdif(jl)=2.*x*wz
c  calc diffuse non-interception for a sublayer from non-inter
c    for a df layer because integral takes care of ang dist thru
c    the whole df layer and sublayers assume isotropy of inc rad.
      difsub=expdif(jl)**.1
      do 601 k=1,kmax
      rlayr(k,jl)=(1.-expdif(jl))*rlfdif(k) !!!!!
      tlayr(k,jl)=(1.-expdif(jl))*tlfdif(k)+expdif(jl)
      rsubl(k)=(1.-difsub)*rlfdif(k)
      tsubl(k)=(1.-difsub)*tlfdif(k)+difsub
 601  continue
c  solve sublayer equations
      do 620 k=1,kmax
 603  continue
      d(k,laysp1)=1.
      adum(1)=0.
      tlay2=tsubl(k)*tsubl(k)
      do 605 j=2,laysp1
      jm1=j-1
 605  adum(j)=adum(jm1)*tlay2/(1.-adum(jm1)*rsubl(k)) + rsubl(k)
      do 610 j=2,laysp1
      jj=laysp1-j+1
      jjp1=jj+1
      d(k,jj)=d(k,jjp1)*tsubl(k)/(1.-adum(jj)*rsubl(k))
      u(k,jjp1)=adum(jjp1)*d(k,jjp1)

 610  continue

 615  rlayr(k,jl)=u(k,laysp1)
      tlayr(k,jl)=d(k,1)

 620  continue
 625  continue

c
c  diffuse should be evaluated at layer midpoints
c  calc thick layer sources from laysub sublayers. since source dist
c    depends on sun angle, this effective source is dep on hour and
c    wavelength.
c  zero out all effective sources.
      do 627 k=1,kmax
      sourup(k,jl)=0.
 627  sourdn(k,jl)=0.
c 
      if(coszen.lt.0.01)goto 690
      dirsub=exp(-clump*xint*df(jl)/(coszen*laysub))
c parameter clump in above statement was added by Chen, 9/4/89.
c 
      do 691 k=1,kmax
c original 690 for k loop was changed to 691 since the stm
c if... goto 690 caused trouble. Chen, 9/7/89.
      if(fbeam1(k).lt.0.01)goto 690
      sdnsub(k,laysp1)=(1.-dirsub)*tlfdir(k)
      supsub(k,laysp1)=(1.-dirsub)*rlfdir(k)
      do 630 j=2,laysub
      jj=laysp1-j+1
      jjp1=jj+1
      sdnsub(k,jj)=sdnsub(k,jjp1)*dirsub
      supsub(k,jj)=supsub(k,jjp1)*dirsub

 630  continue
c  zero all u(k,j) and d(k,j).  d(k,laysp1)=u(k,1)=0. are b.c.
      iter=0
      do 640 j=1,laysp1
      d(k,j)=0.
 640  u(k,j)=0.
 645  iter=iter+1
      irept=0
      do 650 j=2,laysp1
      jj=laysp1-j+1
      jjp1=jj+1
      do wn=tsubl(k)*d(k,jjp1)+u(k,jj)*rsubl(k)+sdnsub(k,jjp1)
      if(abs(down-d(k,jj))-.0001)646,646,644
 644  irept=1
 646  d(k,jj)=down
      up=tsubl(k)*u(k,jj)+d(k,jjp1)*rsubl(k)+supsub(k,jjp1)
      if(abs(up-u(k,jjp1))-.0001)649,648,648
 648  irept=1
 649  u(k,jjp1)=up
 650  continue
      if(irept.ne.0) goto 645
c  calc source terms to be used in radiat, these must be multiplied
c  tbeam above each layer to get the sdn(j) and sup(j) needed there.
      sourup(k,jl)=u(k,laysp1)
      sourdn(k,jl)=d(k,1)

 691  continue
c above stm was added by Chen, 9/7/89.
 690  continue
 600  continue
c 
      return
      end
