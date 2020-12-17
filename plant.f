c Leaf energy balance subroutines are from CUPID (John M. Norman, Department of Soil Science,
c University of Wisconsin-Madison, http://www.soils.wisc.edu/~norman/cupid/outline.html)
c**********************************************************
      subroutine cpstruct(imunu,gr)
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      USE parameters
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    GR
      INTEGER             IMUNU
      DOUBLE PRECISION    ALPHA,       ANG2M,       ANGA2M,      ANGAMN
      DOUBLE PRECISION    ANGASD,      ANGL,        ANGMN,       ANGSD
      DOUBLE PRECISION    ANIRR,       ANIRT,       AR,          AT
      DOUBLE PRECISION    AVISR,       AVIST,       BETA,        BNIRR
      DOUBLE PRECISION    BNIRT,       BR,          BT,          BVISR
      DOUBLE PRECISION    BVIST,       CA,          CNIRR,       CNIRT
      DOUBLE PRECISION    CR,          CTR,         CVISR,       CVIST
      DOUBLE PRECISION    DELTHE,      FR1,         FR2
      DOUBLE PRECISION    FRDEG,       GAMRAT,      PHI,         PHIV
      DOUBLE PRECISION    SA,          SUM,         SUM1,        SUM2
      DOUBLE PRECISION    SUM3,        SUM4,        SUMA
      DOUBLE PRECISION    THETA,       THETAS,      THETV,       WBETA
      DOUBLE PRECISION    WZ,          XMEU,        XNEU,        XXINT
      INTEGER             I,           IALPHA,      IANGL,       IM1
      INTEGER             IMUNUA,      IPHI,        IPHIV,       ITHETV
      INTEGER             IXINTZ,      IZ,          K,           L
      DOUBLE PRECISION    COSDEC,      COSLAT,      DECL,        DECMAX
      DOUBLE PRECISION    DLONG,       EQTM,        SINDEC,      SINLAT
      DOUBLE PRECISION    TANLAT
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    SUNAZM,      VIEWAZ,      VIEWZN
      DOUBLE PRECISION    XINTV
      INTEGER             NOAZMV,      NOZENV
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
      DOUBLE PRECISION    DSTNET,                   DSTRAD
      DOUBLE PRECISION    FRAREA,                   TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    BMFLX,       D,           RNDIV
      DOUBLE PRECISION    RNET,        RNLAM,       TSFC
      DOUBLE PRECISION    UI
      DOUBLE PRECISION    RLFDIF,      RLFDIR,      TLFDIF
      DOUBLE PRECISION    TLFDIR
      DOUBLE PRECISION    RLFHEM,      TLFHEM
      !shc end of adding declarations to allow "implicit none"

c leaf angle distribution
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/deg/sunazm,viewzn(99),viewaz(99),nozenv,noazmv
     & ,xintv(20,90)
      common/astron/eqtm,decl,sindec,cosdec,decmax,sinlat,coslat,
     1 tanlat,dlong
      common/indax3/xmeuaz,xneuaz,beta0,fraz(99),fradeg(360),nbeta
      common/misc1/iwrite(9,99),kmax
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/misc6/xint,xintz(99),nalpha,ispher,nxintz
      common /rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3)
     &,expdif(99),rlayr(3,30),tlayr(3,30)
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common/rad4/d(3,99),ui(3,99),bmflx(3,99),rnet(99),rndiv(99),
     & tsfc,rnlam(3,99)
      common/rad6/rlfdif(3),tlfdif(3),rlfdir(3),tlfdir(3)
      common/rad8/rlfhem(3,9),tlfhem(3,9)
      dimension gr(99),frdeg(91),theta(99),ca(9),sa(9)
c
c  input coefficients of parabola describing hemispherical refl and tran
c    as a func of source inc angle(a*x**2 + b*x + c = refl) for vis and
c    nir. other wavebands interpolated from these. c=normal inc value.
      read(15,*)avisr,bvisr,cvisr,anirr,bnirr,cnirr
      read(15,*)avist,bvist,cvist,anirt,bnirt,cnirt
      nxintz=9
      read(15,*)ispher,nalpha,(gr(l),l=1,nalpha)
      read(15,*)imunu,xmeu,xneu
      write(06,'(/,2x,"ispher, imunu, meu, neu: ",2i3,2f6.3)')
     & ispher,imunu,xmeu,xneu

      do 15 k=1,kmax
 15   aleaf(k)=1.-rleaf(k)-tleaf(k)
c
c*****  cupradazm insert begin 2: leaf reflectance/transmittance *******
c  Chen, 8/30/89.
c
c  input coefficients of parabola describing hemispherical refl and tran
c    as a func of source inc angle(a*x**2 + b*x + c = refl) for vis and
c    nir. other wavebands interpolated from these. c=normal inc value.
c  interpolate a,b,and c coeff for leaf refl and trans dep on source inc
      do 26 k=1,kmax
      do 26 i=1,itot
      angl=90./itot*(i-1)+90./(2.*itot)
      ar=anirr-(cnirr-rleaf(k)*100.)*(anirr-avisr)/(cnirr-cvisr)
      br=bnirr-(cnirr-rleaf(k)*100.)*(bnirr-bvisr)/(cnirr-cvisr)
      cr=rleaf(k)*100.
      rlfhem(k,i)=(ar*angl*angl+br*angl+cr)/100.
c  rlfhem is used with angle betw leaf normal and sun from light dist.
      at=anirt-(cnirt-tleaf(k)*100.)*(anirt-avist)/(cnirt-cvist)
      bt=bnirt-(cnirt-tleaf(k)*100.)*(bnirt-bvist)/(cnirt-cvist)
      ctr=tleaf(k)*100.
      tlfhem(k,i)=(at*angl*angl+bt*angl+ctr)/100.
c  if refl +tran of leaf > 1 then reduce trans.
      if(rlfhem(k,i)+tlfhem(k,i).gt.1.)tlfhem(k,i)=1.-rlfhem(k,i)
 26   continue
c  calc sums to get weighted leaf refl and trans for diffuse inc rad.
      do 101 ialpha=1,itot
      alpha=pid2*(ialpha-1)/(itot) + pi/(4.*itot)
      sa(ialpha)=sin(alpha)
      ca(ialpha)=cos(alpha)
 101  continue
      do 400 k=1,kmax
      sum2=0.
      sum3=0.
      sum4=0.
      do 300 iangl=1,itot
      angl=pid2*(iangl-1)/itot + pid2/(2.*itot)
      sum2=sum2+sa(iangl)*ca(iangl)
      sum3=sum3+rlfhem(k,iangl)*sa(iangl)*ca(iangl)
      sum4=sum4+tlfhem(k,iangl)*sa(iangl)*ca(iangl)
 300  continue
      rlfdif(k)=sum3/sum2
      tlfdif(k)=sum4/sum2
 400  continue

c   calc area fraction of leaves in each leaf angle class for
c   spherical leaf angle distrib.
      delthe=pi/(2.*itot)
      theta(1)=delthe*.5
      ct(1)=cos(theta(1))
      fr2=cos(delthe)
      fr(1)=1.-fr2
      fr1=fr2
      do  100 i=2,itot
      im1=i-1
      theta(i)=theta(im1)+delthe
      ct(i)=cos(theta(i))
      fr2=cos(i*delthe)
      fr(i)=fr1-fr2
      fr1=fr2
 100  continue

      if(nalpha.ne.itot) write(6,123)itot,nalpha
 123  format(' itot= ',i3,' nalpha= ',i3)
      if(ispher.ne.0) goto 180
      if(imunu.eq.1) goto 185
      if(imunu.eq.2) angmn=xmeu
c      if(imunu.eq.2)ang2m=xneu
c above was commented and following two added by Chen, 7/6/89.
	if(imunu.eq.2) angsd=xneu
	if(imunu.eq.2) ang2m=angmn*angmn+angsd*angsd
      if(imunu.eq.2) goto 176
      sum=0.
      do 160 i=1,nalpha
 160  sum=sum+gr(i)
      do 170 i=1,nalpha
 170  fr(i)=gr(i)/sum
      do i=1,nalpha

      enddo
      itot=nalpha
      goto 172
c  set nalpha=itot so leaf angle classes are consistent
 180    do 171 ixintz=1,nxintz
 171    xintz(ixintz)=.5
      nalpha=itot
c
c  calc 1st and 2nd moments to fit beta dist to fr(i) for input 
c  distribution [gr(i)] or spherical.
c 
 172  sum1=0.
      sum2=0.
      do 175 i=1,nalpha
      alpha= 90./(2.*itot) + (i-1)*90./itot
      sum1=sum1+fr(i)*alpha
      sum2=sum2+fr(i)*alpha**2
 175  continue
      angmn=sum1
      ang2m=sum2
      angsd=sqrt(ang2m-angmn**2)
 176  xneu=(1.-ang2m/(90.*angmn))/(ang2m/angmn**2 -1.)
      xmeu=xneu*(90./angmn - 1.)
 185  continue
c  set nalpha=itot so leaf angle classes are consistent
      nalpha=itot
      read(15,*) imunua,xmeuaz,xneuaz,beta0,nbeta

	if (ispher.eq.0 .and. imunua.ne.0) goto 194
	  xmeuaz=1.
	  xneuaz=1.
	  beta0=0.
	  angamn=180.
	  angasd=103.92
	  goto 196
 194    if (imunua.ne.2) goto 196
	  angamn=xmeuaz
	  angasd=xneuaz
	  anga2m=angamn*angamn+angasd*angasd
	  xneuaz=(1.-anga2m/(360.*angamn))/(anga2m/angamn**2-1.)
	  xmeuaz=xneuaz*(360./angamn-1.)
 196    continue
       write(06,'(/,x,a)') 'Fraction in each angle class:'
       write(06,'(x,a8,10f8.3)') 'fr(i): ',(fr(i),i=1,itot)

c ----------------------------------------------------------------------
      call lad(xmeu,xneu,fr,itot,gamrat,frdeg)
c ----------------------------------------------------------------------
 
	if (ispher.eq.1) goto 198
	  if (imunu.ne.1) goto 197
	    angmn=0.
	    ang2m=0.
	    do 186 i=1,nalpha
	    alpha=(90./nalpha)*(i-.5)
	    angmn=angmn+alpha*fr(i)
 186        ang2m=ang2m+alpha*alpha*fr(i)
	    angmn=angmn/nalpha
	    ang2m=ang2m/nalpha
	    angsd=sqrt(ang2m-angmn*angmn)
 197      if (imunua.ne.1) goto 198
	    angamn=0.
	    anga2m=0.
	    do 187 i=1,nbeta
	    beta=(360./nbeta)*(i-.5)
	    angamn=angamn+beta*fraz(i)
 187        anga2m=anga2m+beta*beta*fraz(i)
	    angamn=angamn/nbeta
	    anga2m=anga2m/nbeta
	    angasd=sqrt(anga2m-angamn*angamn)
 198    continue

c following 17 statements calculates extinction coefficient at nxintz
c zeniths averaged over azimuth, xintz(iz). Chen, 05/25/89.
	if (ispher.eq.1) goto 290
	wz=pi/2/nxintz
	wbeta=2*pi/nbeta
	do 211 iz=1,nxintz
	thetas=wz*(iz-.5)
	if (abs(xmeuaz-1.+xneuaz-1.).gt.0.0001) goto 263
c ----------------------------------------------------------------------
	call simpsn(thetas,0.,fr,nalpha,xintz(iz))
c ----------------------------------------------------------------------
	goto 211
 263    suma=0.
	do 220 iphi=1,nbeta
	phi=wbeta*(iphi-.5)
c ----------------------------------------------------------------------
	call simpsn(thetas,phi,fr,nalpha,xxint)
c ----------------------------------------------------------------------
 220    suma=suma+xxint
	xintz(iz)=suma/nbeta
 211    continue
 290    continue
c calculation of extinction coef. in viewing directions, xintv(.,.)
c follow. Chen, 9/19/89.
      do 270 ithetv=1,nozenv
      do 280 iphiv=1,noazmv
      xintv(ithetv,iphiv)=0.5
      if(ispher.eq.1) goto 280
      thetv=viewzn(ithetv)*pid180
      phiv=viewaz(iphiv)*pid180
c ----------------------------------------------------------------------
      call simpsn(thetv,phiv,fr,nalpha,xintv(ithetv,iphiv))
c ----------------------------------------------------------------------
 280  continue
 270  continue
      return
      end 

c**********************************************************
      subroutine stoma(lprin,tair,relh)
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NLEV
      DOUBLE PRECISION    RELH,        TAIR
      DOUBLE PRECISION    DVPDS,       G,           RBL,         RCUT
      DOUBLE PRECISION    RSTOM,       TL,          VPD,         VPDS
      DOUBLE PRECISION    VPFAC,       VSAT,        ZNUM
      INTEGER             I,           J
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    DELT,        PSILF,       TRAN
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    BKV,         CAIR,        CONMIN,      CUCOND
      DOUBLE PRECISION    D1,          FACJ,        OX,          RASTOM
      DOUBLE PRECISION    RSFAC,       RXCHAM,      ZNON
      INTEGER             IC3C4
      DOUBLE PRECISION    DSTNET,                   DSTRAD
      DOUBLE PRECISION    FRAREA,                   TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    HPSI,        RHLEAF,      RSLEAF
      DOUBLE PRECISION    RSNOVP,                   RST
      DOUBLE PRECISION    ANSTOM,      PSI1,        PSI2,        RADN
      DOUBLE PRECISION    RCUT20,      RSEXP,       RSM,         RSMIN
      DOUBLE PRECISION    TRSMAX,      TRSMIN,      TRSOPT
      DOUBLE PRECISION    AROOT,       CPYTR,       FROOT,       PSISUM
      DOUBLE PRECISION    PSITOP,      PSIXY,       RESROT,      ROOTSM
      DOUBLE PRECISION    ROOTUP,      RROOT
      !shc end of adding declarations to allow "implicit none"
      parameter (nlev=40)    ! number of levels  
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/resis1/rhleaf(99),rsleaf(10,99),rsnovp(10,99),
     & rst(10,99),hpsi
      common/leaf1/delt(10,99),psilf,tran(10,99)
      common/resis2/radn,anstom,rcut20,rsmin,trsopt,trsmin,trsmax,rsm
     1,rsexp,psi1,psi2
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/root1/froot(99),resrot(99),rootsm,psixy,psitop,rroot
     1,rootup(99),cpytr,psisum,aroot
      common/photo2/ox,cair,rastom,cucond,d1,bkv,facj,rxcham,rsfac
     &             ,znon,conmin,ic3c4

      logical lprin
       dimension relh(nlev)

c  cuticular resistance may be a func of temp
      rcut=rcut20
      do 1000 i=1,itotp1
      do 1000 j=2,jtot
      tl=templf(i,j)
c   calc rs vs radiation
      if(dstrad(1,i,j).gt.0.01) goto 100
      rsleaf(i,j)=rcut
      rst(i,j)=rsleaf(i,j)   ! Includes the cuticular resistance !!!!
      goto 1000
 100  znum=(anstom-1.)*radn/dstrad(1,i,j)
      rsleaf(i,j)=rsmin*(znum+1.)
c  calc temp effect g(templf)
      if(templf(i,j)-trsopt) 400,400,500
 400  g=1.+(((trsopt-templf(i,j))/(trsopt-trsmin))**rsexp)*rsm
      goto 600
 500  g=1.+(((templf(i,j)-trsopt)/(trsmax-trsopt))**rsexp)*rsm
 600  continue
c  calc water pot effect hpsi
      if(psitop.ge.psi1) hpsi=1.
      if(psitop.lt.psi1) hpsi=(psitop-psi2)/(psi1-psi2)
      if(psitop.lt.psi2) hpsi=0.
      if(hpsi.lt.0.001) hpsi=0.001
      rsleaf(i,j)=rsleaf(i,j)*g/hpsi

      vsat=6.108*10**(7.5*tl/(237.3+tl))
      vpd=vsat*(1.-relh(j-1+levhtr))
c     vpd=vsat*(1.-relh(j))
      rstom=rsleaf(i,j)
      rbl=0.0 ! rx*cmo*rsfac Boundary layer resistance is assumed zero
              ! as long as some numbers are still missing
      vpds=vpd*rstom/(rstom+rbl)
      vpfac=rstom/(rstom+rbl)
      dvpds=vpds/vsat
      rsleaf(i,j)=rsleaf(i,j)*(1.+bkv*(dvpds-d1))
      rst(i,j)=rsleaf(i,j)   ! without cuticular resistance

c  calc final resis as parallel comb of rcut and rs
      rsleaf(i,j)=rcut*rsleaf(i,j)/(rcut+rsleaf(i,j)) 
 1000 continue
      return
      end
 
c**********************************************************
      subroutine lfebal(tair,eair,wind,vpaira,lprin)
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      USE parameters
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
!      INTEGER NLEV
      DOUBLE PRECISION    dstradsky, VPAIRA,esky  ! added by Dandan Wei
      DOUBLE PRECISION    EAIR,        TAIR,          WIND
      DOUBLE PRECISION    BAL,         BAL2,        DEN,         DSUM1
      DOUBLE PRECISION    EMEAN,       ESTAIR,      ESTEMP,      ESTLF
      DOUBLE PRECISION    HH,          R1,          RH,          RSFAC2
      DOUBLE PRECISION    RSLST,       S,           SUM,         SUMBAL
      DOUBLE PRECISION    TEMP,        TMEAN,       UMEAN,       ZNUM
      INTEGER             I,           ICOND,       IFLGRS,      IREPT
      INTEGER             ISATE,       J
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    EAVE,        ETOTWT,      EVAPG
      DOUBLE PRECISION    HAVE,        HEATG,       HTOT,        RPLNT
      DOUBLE PRECISION    WATERG
      DOUBLE PRECISION    CONTOT,      CPESTR,      CPHSTR,      ECPYS
      DOUBLE PRECISION    ETOTW,       EVSMIC,      EVTOT,       HCPYS
      DOUBLE PRECISION    HSOIL,       SCOND,                    WCPYS
      INTEGER             IHRWET
      DOUBLE PRECISION    FRWET,       FRWTMX,      PILAST
      DOUBLE PRECISION    PINT,        PINT1,       TWATER,      WTP
      DOUBLE PRECISION    DRIP,        EVIMM,       EVINT,       FRSTEM
      DOUBLE PRECISION    PINTMX,      STEM
      DOUBLE PRECISION    DELT,        PSILF,       TRAN
      DOUBLE PRECISION    ALAM,        EVAP,        GEVAP
      DOUBLE PRECISION    GHEAT,                    HEAT
      DOUBLE PRECISION    GHWATR,                   HWATER
!      DOUBLE PRECISION    PI,          PID180,      PID2,        SIGMA
      INTEGER             IWRITE,                   KMAX
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    CAIR,        CONMIN,      CUCOND,      D1
      DOUBLE PRECISION    D2,          FACJ,        OX,          RASTOM
      DOUBLE PRECISION    RSFAC,       RXCHAM,      ZNON
      INTEGER             IC3C4
      DOUBLE PRECISION    ALEAF,       EMIS,        EMISOL
      DOUBLE PRECISION    EXPDIF,      RLAYR,       RLEAF
      DOUBLE PRECISION    RSOIL,       TLAYR,       TLEAF
      DOUBLE PRECISION    DSTNET,                   DSTRAD
      DOUBLE PRECISION    FRAREA,                   TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    HPSI,        RHLEAF,      RSLEAF
      DOUBLE PRECISION    RSNOVP,      TEMPS,       RST
      REAL    DT,          TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
!      parameter (nlev=40)    ! number of levels  

      logical lprin
      common/timpar/timloc,dt,month,jday,iyear,icumdy
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/resis1/rhleaf(99),rsleaf(10,99),rsnovp(10,99),
     &rst(10,99),hpsi
      common/photo2/ox,cair,rastom,cucond,d1,d2,facj,rxcham,rsfac
     &             ,znon,conmin,ic3c4
      common/leaf1/delt(10,99),psilf,tran(10,99)
      common/leaf2/evap(10,99),gevap(10,99),heat(10,99),gheat(10,99)
     &,alam
      common/leaf3/hwater(10,99),ghwatr(10,99)
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common/misc1/iwrite(9,99),kmax
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/rad1/emis,emisol,rsoil(3),rleaf(3),tleaf(3),aleaf(3),
     1expdif(99),rlayr(3,30),tlayr(3,30)
      common/cpy1/etotwt,htot,rplnt,evapg(99),heatg(99),eave(99)
     &,have(99),waterg(99)
      common/cpy2/hsoil,hcpys,evtot,etotw(99),contot(99),scond(10,99)
     1,ecpys,cphstr,cpestr,wcpys,evsmic,ihrwet(99)
      common/inter2/evint(99),evimm(99),pintmx,frstem,drip(99),stem
      common/inter1/wtp(99),frwet(99),frwtmx,pint(99),pilast(99)
     1,pint1(99),twater
       dimension tair(nlev),eair(nlev),wind(nlev),rh(99)
       
       
c  calc leaf energy bal for each angle class and layer
      do 1100 j=2,jtot
c  tair and eair are at top of layers so ave adjacent values to get
c    layer average
c     tmean=(tair(levhtr+j)+tair(levhtr+j-1))*.5-273.15
c     emean=(eair(levhtr+j)+eair(levhtr+j-1))*.5
c     umean=(wind(levhtr+j)+wind(levhtr+j-1))*.5
      tmean=(tair(levhtr+j-1)+tair(levhtr+j-2))*.5-273.15
      emean=(eair(levhtr+j-1)+eair(levhtr+j-2))*.5
      umean=(wind(levhtr+j-1)+wind(levhtr+j-2))*.5
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      if(umean.le..01) rhleaf(j)=500.
      if(umean.gt..01) rhleaf(j)=180.*sqrt(sizelf/umean) 

      sumbal=0.
      do 1000 i=1,itotp1
      iflgrs=1
      irept=0
      icond=0
      isate=0
 100  temp=templf(i,j) !!!!! (templf(i,j)+tmean)*0.5
c  sat vp calc (in hpa ddw)
c      estemp=6.108*10**(7.5*temps/(temps+237.3))
c      estair=6.108*10**(7.5*tmean/(tmean+237.3)) 
      estemp=6.112*exp(17.67*temp/(temp+243.5))  ! ddw newer version
      estair=6.112*exp(17.67*tmean/(tmean+243.5))
c  latent heat
      alam=597.-0.57*temp
c  slope of sat vp curve
c      s=estemp*alam/(0.1103*(temps+273.15)**2)
       s = estair*17.67*243.5/(tmean+243.5)**2  !ddw dq(Tair)/dt
      if(contot(j).ge.0.)goto 400
 300  isate=1
      iflgrs=0
 400  continue
c  rsleaf is for both sides of leaf  rhleaf is for one side of leaf
c  for amphistom leaves (i.e. stomata on both sides) rsfac=.5  see sub inplnt
c  for hypostom leaves rsfac=1
      if(isate.eq.1) rsfac2=.5
      if(isate.eq.0) rsfac2=rsfac
      r1=iflgrs*rsleaf(i,j)+rhleaf(j)*rsfac2
      rslst=rsleaf(i,j)
      hh=frwet(j)/(rhleaf(j)*0.5)+(1.-frwet(j))/r1 
c     approx correction for heat taken from leaf by cold irrig water
      if(pint1(j).lt..0001) twater=tair(j)
c ddw see Bonan book leaf temperature; changed temp in "den" to tmean
      znum=1825.*hh*(estair-emean)+1.134e-7*emis*(tmean+273.15)**4
     &+pint1(j)*4187.6*(tair(j)-twater)/(dt*df(j-1))
      den=1825.*s*hh+2400./rhleaf(j)+4.536e-7*emis*(tmean+273.15)**3
     &+pint1(j)*4187.6/(dt*df(j-1))
c ddw added absorption emis to dstrad(3, i, j): too much nighttime cooling

      delt(i,j)=(dstnet(i,j)+2.*dstrad(3,i,j)-znum)/den
      if(delt(i,j).lt.-10.) delt(i,j)=-10.
      if(delt(i,j).gt.10.) delt(i,j)=10.

c     write(06,'("delt",2i3,e12.4)') i,j,delt(i,j)
      templf(i,j)=tmean+delt(i,j)
      estlf=6.112*exp(17.67*templf(i,j)/(templf(i,j)+243.5))  ! ddw newer version
c      estlf=6.108*10**(7.5*templf(i,j)/(templf(i,j)+237.3))
      if(isate.eq.1) goto 700
      if(emean-estlf)580,580,550
 550  if(icond.eq.1) goto 700
      icond=1
      iflgrs=0
      goto 100
 580  if(icond-1) 700,595,595
 595  if(irept.eq.1) goto 700
      iflgrs=1
      irept=1
      goto 100
 700  continue
      
      rh(j)=emean/estair

      evap(i,j)=1825.*hh*((estair-emean) +s*delt(i,j))
       
      tran(i,j)=(1.-frwet(j))*1825./r1*(estair-emean+s*delt(i,j))
      gevap(i,j)=evap(i,j)*df(j-1)*frarea(i,j)
      heat(i,j)=2400./rhleaf(j)*delt(i,j)
      gheat(i,j)=heat(i,j)*df(j-1)*frarea(i,j)
      hwater(i,j)=pint1(j)*4187.6*(templf(i,j)-twater)/(dt*df(j-1))
      ghwatr(i,j)=hwater(i,j)*df(j-1)*frarea(i,j)
c ----------------------------------------------------------------

      dstnet(i,j)=dstnet(i,j)+2.0*dstrad(3,i,j)
      dstrad(3,i,j)=emis*sigma*(templf(i,j)+273.)**4
      dstnet(i,j)=dstnet(i,j)-2.*dstrad(3,i,j) 
c
c stms to check leaf energy bal, should be good to .01wm-2.
      bal=dstnet(i,j)*df(j-1)*frarea(i,j)-gevap(i,j)-gheat(i,j)
     &  -ghwatr(i,j)
      bal2=dstnet(i,j)-evap(i,j)-heat(i,j)-hwater(i,j)
      sumbal=sumbal+bal
c     if((iday.eq.3).and.(ihr.eq.9))write(15,999)iday,ihr,i,j,bal,bal2
c    &   ,ghwatr(i,j),hwater(i,j)
 999  format(4i3,4f8.4)
 1000 continue
c     if((iday.eq.3).and.(ihr.eq.9))write(15,1001)ihr,j,sumbal
 1001 format(' leaf ',2i3,f7.2)
 1100 continue

      goto 7777

      do j=jtot,2,-1
       dsum1=0.0
       do i=1,itotp1
        dsum1=dsum1+delt(i,j)*frarea(i,j)
       enddo
      enddo

 7777 continue

      do 2000 j=2,jtot
      sum=0.
      do 1500 i=1,itotp1
      sum=sum+(evap(i,j)-tran(i,j))*frarea(i,j)
 1500 continue
      evint(j)=sum*df(j-1)
 2000 continue
c     if(ihr.eq.21)write(26,405)ihr,isate,iflgrs,irept,icond,eair(jtot),
c    1estlf,tair(jtot)
 405  format(1x,5i3,2(1x,e11.4)) !shc

      return
      end
 
c**********************************************************
      subroutine stress
c**********************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NLEV
      DOUBLE PRECISION    TRANW
      INTEGER             I,           J
      DOUBLE PRECISION    EAVE,        ETOTWT,      EVAPG
      DOUBLE PRECISION    HAVE,        HEATG,       HTOT,        RPLNT
      DOUBLE PRECISION    WATERG
      DOUBLE PRECISION    CONTOT,      CPESTR,      CPHSTR,      ECPYS
      DOUBLE PRECISION    ETOTW,       EVSMIC,      EVTOT,       HCPYS
      DOUBLE PRECISION    HSOIL,       SCOND,                    WCPYS
      INTEGER             IHRWET
      DOUBLE PRECISION    DELT,        PSILF,       TRAN
      DOUBLE PRECISION    ALAM,        EVAP,        GEVAP
      DOUBLE PRECISION    GHEAT,                    HEAT
      DOUBLE PRECISION    GHWATR,                   HWATER
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      !shc end of adding declarations to allow "implicit none"
      parameter (nlev=40)    ! number of levels  
      common/leaf1/delt(10,99),psilf,tran(10,99)
      common/leaf2/evap(10,99),gevap(10,99),heat(10,99),gheat(10,99)
     &,alam
      common/leaf3/hwater(10,99),ghwatr(10,99)
      common/cpy1/etotwt,htot,rplnt,evapg(99),heatg(99),eave(99)
     &,have(99),waterg(99)
      common/cpy2/hsoil,hcpys,evtot,etotw(99),contot(99),scond(10,99)
     1,ecpys,cphstr,cpestr,wcpys,evsmic,ihrwet(99)
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot

c  calc leaf water pot from total evap,soil potent and plant res a


      etotwt=0.
      htot=0.
      tranw=0.
      rplnt=30.
      do 600 j=2,jtot
      evapg(j)=0.
      heatg(j)=0.
      waterg(j)=0.
      etotw(j)=0.


      do 500 i=1,itotp1
c  calc tot flux den per unit of leaf area
      etotwt=etotwt+gevap(i,j)
      htot=htot+gheat(i,j)
      heatg(j)=heatg(j)+gheat(i,j)
      waterg(j)=waterg(j)+ghwatr(i,j)
      evapg(j)=evapg(j)+gevap(i,j)
      etotw(j)=etotw(j)+gevap(i,j)/(alam*4.18e-03)
 500  continue
 600  if(contot(j).ge.0.)tranw=tranw+etotw(j)

      evtot=etotwt
      etotwt=etotwt/(alam*4.18e-3)
c     psilf=psisol(ihr)-tranw/rplnt
      psilf=-.1 ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>

      return
      end
C

c**********************************************************
      subroutine rootex
c**********************************************************
c  subroutine to calculate root extraction from the soil
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NLEV, NSOIL
      INTEGER             I,           J,           JZ
C
C     Common variables
C 
      DOUBLE PRECISION    AXA,         BXB,         CXC
      DOUBLE PRECISION    TEMPS 
      DOUBLE PRECISION    D1ETA,                    D1T
      DOUBLE PRECISION    D2ETA,                    D2T
      DOUBLE PRECISION    DATM,        DXD,         ETA
      DOUBLE PRECISION    H,           P21,         P21T
      DOUBLE PRECISION    PSI,         PSIETA
      DOUBLE PRECISION    PSIT,        PSITET
      DOUBLE PRECISION    RHO0,        RHO1,        RHO22
      DOUBLE PRECISION    TEMP,        XK,          XL11
      DOUBLE PRECISION    XL210,                    XL22
      DOUBLE PRECISION    XLAM
      DOUBLE PRECISION    DELT,        PSILF,       TRAN
      DOUBLE PRECISION    ALAM,        EVAP,        GEVAP
      DOUBLE PRECISION    GHEAT,                    HEAT
      DOUBLE PRECISION    CLAI,        CT,          DF
      DOUBLE PRECISION    DISTLS,      FR,          TOTLAI
      INTEGER             ITOT,        ITOTP1,      JTOT
      DOUBLE PRECISION    Z,           ZMID
      INTEGER             JZBM1,       JZBOT,       JZCPY,       JZCPY1
      INTEGER             JZCRIT,      JZSFC,       JZSFM1
      DOUBLE PRECISION    DSTNET,                   DSTRAD
      DOUBLE PRECISION    FRAREA,                   TEMPLF
      DOUBLE PRECISION    TSOIL
      DOUBLE PRECISION    AROOT,       CPYTR,       FROOT,       PSISUM
      DOUBLE PRECISION    PSITOP,      PSIXY,       RESROT,      ROOTSM
      DOUBLE PRECISION    ROOTUP,      RROOT
      DOUBLE PRECISION    DELS,        SL,          ZS
      DOUBLE PRECISION    AKS,         AN,          ASOIL,       BD
      DOUBLE PRECISION    BSOIL,       BX,          CSOIL,       DSOIL
      DOUBLE PRECISION    ESOIL,       PE,          WS
      DOUBLE PRECISION    TPRECP,      WCOND,       WSTOR
      INTEGER             IPRECP
      !shc end adding declarations to allow "implicit none"

      parameter (nlev=40,nsoil=15)    ! number of levels  
      common/misc2/fr(99),ct(99),totlai,df(99)
     &,clai(99),distls(10),itot,itotp1,jtot
      common/misc4/z(99),zmid(99),jzcpy,jzcpy1,jzsfc,jzsfm1,jzbot,jzcrit
     &,jzbm1
      common/leaf1/delt(10,99),psilf,tran(10,99)
      common/leaf2/evap(10,99),gevap(10,99),heat(10,99),gheat(10,99)
     1,alam
      common/root1/froot(99),resrot(99),rootsm,psixy,psitop,rroot
     1,rootup(99),cpytr,psisum,aroot
      common /rad2/dstrad(3,10,99),dstnet(10,99),frarea(10,99)
     &,templf(10,99),tsoil
      common/soil1/zs(nsoil),sl,dels
      common/feld/ eta(nsoil),psi(nsoil),psit(nsoil),
     &temps(nsoil),xk(nsoil),psieta(nsoil),
     &psitet(nsoil),xlam(nsoil),axa(nsoil),bxb(nsoil),
     &cxc(nsoil),dxd(nsoil),d1t(nsoil),
     &d2t(nsoil),d1eta(nsoil),d2eta(nsoil),
     &xl11(nsoil),xl22(nsoil),rho0(nsoil),
     &rho1(nsoil),rho22(nsoil),xl210(nsoil),
     &p21(nsoil),p21t(nsoil),h(nsoil),datm(nsoil)
      common/soil4/pe,bx,bd,aks,an,ws,asoil,bsoil,csoil,dsoil,esoil
      common/water1/tprecp,wcond(99),wstor(99),iprecp
c
      cpytr=0.
      do 10 j=2,nsoil
      psi(j)=pe*(eta(j)/ws)**(-bx) ! should be replaced by psi from soil model
      if(psi(j).ge.pe) psi(j)=pe
      
      do 10 i=1,itotp1
      cpytr=cpytr+tran(i,j)*frarea(i,j)
 10   continue
c  convert w m-2 to kg m-2 s-1
      cpytr=cpytr*df(j-1)/(alam*4.18e3)
      psisum=0.
 300  do 200 jz=2,nsoil          
      if(psi(jz).ge.-1500) psisum=psisum+psi(jz)/resrot(jz)
c     write(19,*)ihr,jz,psisum,resrot(jz),rootsm,froot(jz),psi(jz)
 200  continue
 
c  rootsm is sum of reciprocals of root resis from main prog
c      if(cpytr.ge.0)psixy=(-cpytr+psisum)/rootsm
c      if(cpytr.lt.0.)psixy=psisum/rootsm
c above 2 stm were commented and following 3 added. chen, 02/19/90.
        psisum=psisum/rootsm
        if(cpytr.ge.0) psixy=-cpytr/rootsm+psisum
        if(cpytr.lt.0.) psixy=psisum
c  assume .6 of plant resis is in root so rleaf=2*rroot/3
      psitop=psixy-cpytr*rroot*0.6667
      if(cpytr.lt.0.) psitop=psixy
      do 400 jz=2,nsoil            
      rootup(jz)=(psi(jz)-psixy)/resrot(jz)
      if(psi(jz).lt.-1500.) rootup(jz)=0.
 400  continue
      rootup(1)=0.
      rootup(nsoil)=rootup(nsoil-1)
      return
      end

