c  Soil model (Forkel et al., 1984)
c 
      subroutine grenzn(fnet,patm,p21,h,rhoair,psi,
     &p21t,psieta,xl111,xl112,xl221,xl222,psit,d1t1,d1t2,d1eta1,
     &d1eta2,d2t1,d2t2,d2eta1,d2eta2,xl210,xlam1,xlam2,ajn1,
     &z1,z2,zb1,zb2,
     &ta1,ta2,t1,vcp1,q2,tb1,tb2,eta1,eta2,qgrnd,qsg,
     &afl1,aflh,bflh,bfl1,bfl2,deltim,tausum,itime,
     &p1a,p1b)


      USE cacm3_Precision, ONLY:dp 
c describes soil-atmosphere boundary

      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    AFL1,        AFLH,        AJN1,        BFL1
      DOUBLE PRECISION    BFL2,        BFLH,        D1ETA1,      D1ETA2
      DOUBLE PRECISION    D1T1,        D1T2,        D2ETA1,      D2ETA2
      DOUBLE PRECISION    D2T1,        D2T2,        ETA1
      DOUBLE PRECISION    ETA2,        FNET,        H,           P1A
      DOUBLE PRECISION    P1B,         P21,         P21T,        PATM
      DOUBLE PRECISION    PSI,         PSIETA,      PSIT,        Q2
      DOUBLE PRECISION    QGRND,       QSG,         RHOAIR,      T1
      DOUBLE PRECISION    TA1,         TA2,         TAUSUM
      DOUBLE PRECISION    VCP1,        XL111,       XL112
      DOUBLE PRECISION    XL210,       XL221,       XL222,       XLAM1
      DOUBLE PRECISION    XLAM2,       ZB1,         ZB2
      REAL(KIND=dp) ::  Z1,  Z2
      DOUBLE PRECISION   TB1, TB2
      REAL   deltim
C 
      INTEGER             ITIME
C
C     Local variables
C 
      DOUBLE PRECISION    AA,          D1ETA,       D1T,         D2ETA
      DOUBLE PRECISION    D2T,         DEST,        DET,         DF1DE
      DOUBLE PRECISION    DF1DT,       DF2DE,       DF2DT,       DQGDE
      DOUBLE PRECISION    DQGDT,       EB,          EBP,         EDEST
      DOUBLE PRECISION    F1,          F2,          GRADE,       GRADTB
      DOUBLE PRECISION    P21N,        S,           TS,          TST
      DOUBLE PRECISION    TSTP,        VERL,        XKE,         XL11
      DOUBLE PRECISION    XL22,        XLAM
      INTEGER             I,           ICP,         IER1,        IER2
      INTEGER             IFALL
      DOUBLE PRECISION    ALBDRY,      ALBWET,      CP3,         ETAMAX
      DOUBLE PRECISION    ETAMIN,      FELDKA,      POR,         RHO3
      DOUBLE PRECISION    TREFK,       TREFP,       XKR
      DOUBLE PRECISION    AAA,         AK,          BBB,         CCC
      DOUBLE PRECISION    CP0,         DDD,         GRAV,        P0
      DOUBLE PRECISION    PII,          R0,          R1,          TCONV
      DOUBLE PRECISION    XL
      DOUBLE PRECISION    CK,          CLAM,        CPSI
      DOUBLE PRECISION    CPSIT,       SEK,         SEL,         SEP
      DOUBLE PRECISION    SEPT,        SK,          SXL,          SPP
      DOUBLE PRECISION    SPT
      INTEGER             ISK,         ISL,         ISP,         ISPT
      !shc end of adding declarations to allow "implicit none"
c
c       interface condition
c       calculate temp(z=0) and eta(z=0) from the flux balance at
c       the soil surface by a newton-iteration
c 
      common/konst/p0,aaa,bbb,ccc,ddd,ak,r0,r1,grav,cp0,xl,tconv
      common/stuetzw/sep(18),spp(18),sek(17),sk(17),sept(10),spt(10)
     &,sel(10),sxl(10),cpsi(17,3),ck(16,3),cpsit(9,3),clam(9,3),
     &isp,isk,ispt,isl
      common/bodk/ etamin,etamax,xkr,por,trefp,trefk,cp3,rho3,
     &feldka,albdry,albwet

      save
      data ifall/0/
      data dest/0./,ts/0./

      icp=isp-1
      s=fnet
      aa=(100000./patm)**.286
      p1a=p1a*100.
c   average transport coefficient of the first level
      xl11=xl111*.95+xl112*.05
      xl22=xl221*.95+xl222*.05
      d1t=d1t1*.95+d1t2*.05
      d1eta=d1eta1*.95+d1eta1*.05
      d2t=d2t1*.95+d2t2*.05
      d2eta=d2eta1*.95+d2eta2*.05
      xlam=xlam1*.95+xlam2*.05
c
c  iterative solution of border area for fluxes of energy f1
c  and water f2
c 
      do 10 i=1,10
      qgrnd=.622*p21*h/(patm-p21*h*.378)
      dqgdt=(h*p21t+p21*h/r1/tb1*(psit-psi/tb1))*.622/patm
      dqgde=.622/patm*p21*h/r1/tb1*psieta
      aflh=-ajn1*rhoair*cp0*(ta2-tb1*aa)/z2
      gradtb=(tb1-tb2)/(zb1-zb2)
      grade=(eta1-eta2)/(zb1-zb2)
      afl1=-ajn1*rhoair/z2*(q2-qgrnd)
      bflh=-xlam*gradtb
      bfl1=-d1t*gradtb-d1eta*grade


      bfl2=-d2t*gradtb-d2eta*grade-xl22*grav

c     case decision
c

      if(i.ne.1.and.ifall.eq.5) goto 50
      if(afl1.gt.0.0) goto 30
      if(p1b.lt.p21) goto 29
c 1.case taubildung und destillation
      ifall=1
      if(p1a.lt.p21) ifall=2
      goto 33
c 2.case wd-fluss nach unten aber keine taubildung
  29  ifall=2
      if(p1a.gt.p21) ifall=1
      goto 33
  30  if(tausum.gt.0.0) goto 34
c  3.case normale evaporation
      ifall=3
      goto 33
c  4.fall verdunstung von tau
  34  verl=(afl1-dest-ts)*deltim
      if(verl.gt.tausum) goto 35
      ifall=4
      goto 33
c  5.fall verbrauch eines taurestes
  35  ifall=5
  33  continue
c
      goto(51,50,50,51,50),ifall
c
c  kein tau
c 
  50  dest=0.0
      ts=0.0
c
c      flussbilanz
c 
      f1=s+aflh-bflh+(xl210-psi)*(afl1-bfl1)
      
      f2=afl1-bfl1-bfl2

c 
      if(abs(f1).gt.0.1) goto 22
      if(abs(f2).le.1.0e-07) goto 11
c vereinfachte ableitungen der grenzflaechenbed. nach t und eta
  22  df1dt=2.+rhoair*cp0*ajn1*aa/z2+xlam/(zb1-zb2)+(xl210-psi)*
     &(ajn1*rhoair/z2*dqgdt+d1t/(zb1-zb2)+afl1/rhoair*.003)
       
      df1de=(xl210-psi)*(ajn1*rhoair/z2*dqgde+d1eta/(zb1-zb2))
     &+7.*gradtb-psieta*(afl1-bfl1)
c 
      df2dt=rhoair*ajn1/z2*dqgdt+(d1t+d2t)/(zb1-zb2)+afl1/
     &rhoair*.003
      df2de=rhoair*ajn1/z2*dqgde+(d1eta+d2eta)/(zb1-zb2)
      goto 55
c
c  tau
c 
  51  edest=.5
      dest=bfl1*edest
      ts=xl22*grav/1000.
      ts=-dmin1(ts,tausum)
      xke=xl22
      f1=s+aflh-bflh+xl210*(afl1-bfl1)-psi*(dest-bfl1-ts)
      f2=dest+ts-bfl1-bfl2
      if(abs(f1).gt.0.1) goto 52
      if(abs(f2).le.1.e-07) goto 11
 52   df1dt=2.+rhoair*cp0*ajn1*aa/z2+xlam/(zb1-zb2)+xl210*
     &(ajn1*rhoair/z2*dqgdt+d1t/(zb1-zb2)+afl1/rhoair*.003)
     &-psi*((1.-edest)*d1eta/(zb1-zb2))
      df1de=7.*gradtb+xl210*d1eta/(zb1-zb2)-
     &psi*((1.-edest)*d1eta/(zb1-zb2)-xke*grav)-
     &psieta*(dest-bfl1-ts)
      df2dt=((1.-edest)*d1t+d2t)/(zb1-zb2)
      df2de=((1.-edest)*d1eta+d2eta)/(zb1-zb2)
c 
  55  det=df1dt*df2de-df1de*df2dt
      eb=eta1-eta1/10.
      ebp=eta1+eta1/10.
      tst=tb1-5.
      tstp=tb1+5.
c
c neue temperatur und wassergehalt an der erdoberflaeche
c 
      tb1=tb1-(df2de*f1-df1de*f2)/det
      
      eta1=eta1-(-df2dt*f1+df1dt*f2)/det
      !shc if(eta1.lt.eb.or.eta1.lt.etamin) eta1=amax1(eb,etamin)
      if(eta1.lt.eb.or.eta1.lt.etamin) eta1=dmax1(eb,etamin) !shc
      if(eta1.lt.eb) eta1=eb
      if(eta1.gt.ebp) eta1=ebp
      if(tb1.lt.tst) tb1=tst 
      if(tb1.gt.tstp) tb1=tstp

      if(eta1.gt.etamax) eta1=etamax

      call icsevu(sep,spp,isp,cpsi,icp,eta1,psi,1,ier1)
      call dcsevu(sep,spp,isp,cpsi,icp,eta1,psieta,1,ier2)
      psi=-10.**psi
      psieta=2.3026*psi*psieta
      p21=610.7*exp(17.1536*(tb1-273.15)/(tb1-38.33))
      p21t=4028./(tb1-38.33)/(tb1-38.33)*p21
      h=exp(psi/r1/tb1)
  10  continue

  11  continue
c
c      calculate amount of dew
c 
      if(ifall.eq.3) then
      tausum=0.0
      goto 60
      endif
      if(ifall.eq.1) then
      tausum=tausum+deltim*(-afl1+dest+ts)
      goto 60
      endif
      if(ifall.eq.4) then
      tausum=tausum+deltim*(-afl1+dest+ts)
      goto 60
      endif
      if(ifall.eq.5) then
      tausum=0.0
      goto 60
      endif
      if(ifall.eq.2) tausum=tausum
  60  continue
c
c  spez. feuchte und temp. im atm. niveau 1 aus constant-flux-bedingung
      vcp1=h*p21/patm
      
      ta1=tb1*aa

      t1=ta1/aa


      p21n=610.7*exp(17.1536*(tb1-273.15)/(tb1-38.33))
      qsg=.622*p21n/(patm-.387*p21n)
      return
      end


      subroutine erdbod
c      calculate new soil temperature (temp) and
c     liquid water content (eta) in the soil

      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NSOIL
      DOUBLE PRECISION    A,           AA,          AB,          ABDC
      DOUBLE PRECISION    ABK,         ABKA,        ABKB,        ABKG
      DOUBLE PRECISION    ABKT,        AD2,         ADQ,         ALPHA
      DOUBLE PRECISION    AT,          B,           BA,          BB
      DOUBLE PRECISION    BETA,        BT,          BXC         
      DOUBLE PRECISION    TEMPS
      DOUBLE PRECISION    C,           CBDA,        CP0,         CP1
      DOUBLE PRECISION    CP2,         D,           DB,          DEDS
      DOUBLE PRECISION    DEDS2,       DELS2,       DELS2I,      DELSQ
      DOUBLE PRECISION    DELSQI,      DL22DS,      DL2PDS,      DT
      DOUBLE PRECISION    DTDS,        DTDS2,       DXA
      DOUBLE PRECISION    E,           F,           F1D1E,       F1D1T
      DOUBLE PRECISION    F1D2E,       F1D2T,       F1DEP,       F1DTP
      DOUBLE PRECISION    F1LAM,       F2D1E,       F2D1T,       F2D2E
      DOUBLE PRECISION    F2D2T,       F2DEP,       F2DTP,       F2LAM
      DOUBLE PRECISION    F3D1E,       F3D1T,       F3D2E,       F3D2T
      DOUBLE PRECISION    F3DEP,       F3DTP,       F3LAM,       FB
      DOUBLE PRECISION    G,           HG1,         HG2,         P21B
      DOUBLE PRECISION    R1T,         RHO2B,       SL2I,        SQSZ
      DOUBLE PRECISION    X,           XL21
      INTEGER             I,           J,           NMIN,        NTYPE
      DOUBLE PRECISION    ALBDRY,      ALBWET,      CP3,         ETAMAX
      DOUBLE PRECISION    ETAMIN,      FELDKA,      POR,         RHO3
      DOUBLE PRECISION    TREFK,       TREFP,       XKR
      DOUBLE PRECISION    AXA,         BXB,         CXC  
      DOUBLE PRECISION    D1ETA,                    D1T       
      DOUBLE PRECISION    D2ETA,                    D2T       
      DOUBLE PRECISION    DATM,        DXD       
      DOUBLE PRECISION    H,           P21,         P21T       
      DOUBLE PRECISION    PSI,         PSIETA       
      DOUBLE PRECISION    PSIT,        PSITET       
      DOUBLE PRECISION    RHO0,        RHO1,        RHO22       
      DOUBLE PRECISION    TEMP,        XK,          XL11       
      DOUBLE PRECISION    XL210,                    XL22       
      DOUBLE PRECISION    XLAM       
      DOUBLE PRECISION    AFL1,        AFLH,        BFL1,        BFL2
      DOUBLE PRECISION    BFLH,        FNET,        PATM,        QGRND
      DOUBLE PRECISION    QSATG,       S,           TAUSUM,      TSUM
      DOUBLE PRECISION    TSURF,       VSUM
      DOUBLE PRECISION    AAA,         AK,          BBB,         CCC
      DOUBLE PRECISION    CP,          DDD,         GRAV,        P0
      DOUBLE PRECISION    PII,          R0,          R1,          TCONV
      DOUBLE PRECISION    XL,          ETA
      DOUBLE PRECISION    DELS,        SL,          ZS       
      DOUBLE PRECISION    CK,          CLAM,        CPSI
      DOUBLE PRECISION    CPSIT,       SEK,         SEL,         SEP
      DOUBLE PRECISION    SEPT,        SK,          SPP,          SPT
      DOUBLE PRECISION    SXL
      INTEGER             ISK,         ISL,         ISP,         ISPT
      REAL    DELTIM,      TIMLOC
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
      parameter (nsoil=15)
      save
      common/soil1/ zs(nsoil),sl,dels
      common/konst/p0,aaa,bbb,ccc,ddd,ak,r0,r1,grav,cp,xl,tconv
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
      common/intf/ qgrnd,qsatg,tausum,tsum,vsum,tsurf,patm,
     &s,fnet,afl1,aflh,bfl1,bfl2,bflh
      common/feld/eta(nsoil),psi(nsoil),psit(nsoil),
     &temps(nsoil),xk(nsoil),psieta(nsoil),
     &psitet(nsoil),xlam(nsoil),axa(nsoil),bxb(nsoil),
     &cxc(nsoil),dxd(nsoil),d1t(nsoil),
     &d2t(nsoil),d1eta(nsoil),d2eta(nsoil),
     &xl11(nsoil),xl22(nsoil),rho0(nsoil),
     &rho1(nsoil),rho22(nsoil),xl210(nsoil),
     &p21(nsoil),p21t(nsoil),h(nsoil),datm(nsoil)
      common/stuetzw/ sep(18),spp(18),sek(17),sk(17),sept(10),
     &spt(10),sel(10),sxl(10),
     &cpsi(17,3),ck(16,3),cpsit(9,3),clam(9,3),
     &isp,isk,ispt,isl
      common/bodk/ etamin,etamax,xkr,por,trefp,trefk,cp3,rho3,
     &feldka,albdry,albwet
      dimension a(nsoil),b(nsoil),c(nsoil),d(nsoil),
     &e(nsoil),f(nsoil),g(nsoil),x(nsoil)
c 
      data cp0/1004.4/,cp1/1845.58/,cp2/4185.0/,rho2b/1000./
      data alpha/2./,beta/0./

      do 5 ntype=1,2
      if(ntype.eq.2) goto 333
      call wert(nsoil,eta,psi,xk,psit,psieta,psitet,xlam)
      do 60 j=1,nsoil
      i=j
      psi(i)=psi(i)+psit(i)*(temps(i)-trefp)
      psieta(i)=psieta(i)+psitet(i)*(temps(i)-trefp)
      xlam(i)=xlam(i)/2.
c 
      r1t=r1*temps(i)
      h(i)=exp(psi(i)/r1t)
      p21(i)=610.7*exp(17.1536*(temps(i)-273.15)/(temps(i)-
     &38.33))
      p21b=h(i)*p21(i)
      xl210(i)=2499700.-2339.42*(temps(i)-273.15)
      xl21=xl210(i)-psi(i)+temps(i)*psit(i)
      p21t(i)=p21(i)*xl210(i)/(r1t*temps(i))
      rho0(i)=(patm-p21b)*(por-eta(i))/(r0*temps(i))
      rho1(i)=p21b*(por-eta(i))/r1t
      rho22(i)=eta(i)*rho2b
      datm(i)=4.42e-8*temps(i)**2.3/patm*133.3
      xl11(i)=p21b*datm(i)*xkr*(por-eta(i))/r1t/r1t
     &*rho0(i)/(rho0(i)+rho1(i))
      xl22(i)=xk(i)*exp(2122.5/trefk)/exp(2122.5/temps(i))
     &*rho2b/grav
      d1t(i)=xl11(i)*xl21/temps(i)
      d1eta(i)=xl11(i)*psieta(i) 
      d1eta(i)=0.1*d1eta(i) ! Fluss wird sonst unrealistisch 
      d2t(i)=xl22(i)*psit(i)
      d2eta(i)=xl22(i)*psieta(i)
      aa=p21b/rho2b/r1t*((por-eta(i))*psieta(i)/r1t-1.)
      abk=p21t(i)-p21(i)/temps(i)+psit(i)*p21(i)/r1t-
     &p21(i)*psi(i)/(r1t*temps(i))
      bb=h(i)*(por-eta(i))/r1t*abk
      axa(i)=rho2b*(1.+aa)
      bxb(i)=bb
      cxc(i)=cp0*rho0(i)+cp1*rho1(i)+cp2*rho22(i)+cp3*rho3+xl210(i)*bb
      dxd(i)=rho2b*(xl210(i)*aa+psi(i)-temps(i)*psit(i))
  60  continue
 333  continue
      dels2=2.*dels
      delsq=dels*dels
      dels2i=1./dels2
      delsqi=1./delsq
      ad2=alpha/dels2
      adq=alpha/delsq
      ba=beta/alpha
c
c-------berechnung von abkuerzungen und koeffizienten------
      do 520 i=2,nsoil-1
c-------berecnung der abkuerzungen-------------------------
      sqsz=sqrt(-sl/zs(i))
      dtds=(temps(i+1)-temps(i-1))*dels2i
      deds=(eta(i+1)-eta(i-1))*dels2i
      dl22ds=(xl22(i+1)-xl22(i-1))*dels2i
      dl2pds=(xl22(i+1)*psi(i+1)-xl22(i-1)*psi(i-1))*dels2i
      dtds2=(temps(i+1)+temps(i-1)-2.*temps(i))*delsqi
      deds2=(eta(i+1)+eta(i-1)-2.*eta(i))*delsqi
c 
      f1d1t=(d1t(i+1)-d1t(i-1))*dels2i-d1t(i)*sqsz
      f1d2t=(d2t(i+1)-d2t(i-1))*dels2i-d2t(i)*sqsz
      f1d1e=(d1eta(i+1)-d1eta(i-1))*dels2i-d1eta(i)*sqsz
      f1d2e=(d2eta(i+1)-d2eta(i-1))*dels2i-d2eta(i)*sqsz
      f1lam=(xlam(i+1)-xlam(i-1))*dels2i-xlam(i)*sqsz
      f1dtp=(d2t(i+1)*psi(i+1)-d2t(i-1)*psi(i-1))*dels2i-d2t(i)*psi(i)*
     &sqsz
      f1dep=(d2eta(i+1)*psi(i+1)-d2eta(i-1)*psi(i-1))*dels2i-d2eta(i)*
     &psi(i)*sqsz
c 
      f2d1t=f1d1t*ad2+d1t(i)*adq
      f2d2t=f1d2t*ad2+d2t(i)*adq
      f2d1e=f1d1e*ad2+d1eta(i)*adq
      f2d2e=f1d2e*ad2+d2eta(i)*adq
      f2lam=f1lam*ad2+xlam(i)*adq
      f2dtp=f1dtp*ad2+d2t(i)*psi(i)*adq
      f2dep=f1dep*ad2+d2eta(i)*psi(i)*adq
c 
      f3d1t=-f1d1t*ad2+d1t(i)*adq
      f3d2t=-f1d2t*ad2+d2t(i)*adq
      f3d1e=-f1d1e*ad2+d1eta(i)*adq
      f3d2e=-f1d2e*ad2+d2eta(i)*adq
      f3lam=-f1lam*ad2+xlam(i)*adq
      f3dtp=-f1dtp*ad2+d2t(i)*psi(i)*adq
      f3dep=-f1dep*ad2+d2eta(i)*psi(i)*adq
c 
      sl2i=1./(2.*sqrt(-sl*zs(i)))
      abkg=deltim*grav*sl2i
      abkt=deltim/(-4.*sl*zs(i))
      abka=deltim*alpha/(2.*dels)*sl2i
      abkb=abka*beta/alpha
      hg1=-cp1*dtds*sl2i+grav
      hg2=-cp2*dtds*sl2i+grav
c----berechnung der koeffizienten fuer die temperatur---------
      if (ntype.eq.2) goto 400
      dxa=dxd(i)/axa(i)
      cbda=cxc(i)-dxa*bxb(i)
      a(i)=abkt*dxa*(f3d1t+f3d2t)+abkt*(-f3lam-xl210(i)*f3d1t-f3dtp)
     &-abka*(d1t(i)*hg1+d2t(i)*hg2)
       
c 
      fb=abkt*2.*adq*(-dxa*d1t(i)-dxa*d2t(i)
     &+xlam(i)+xl210(i)*d1t(i)+d2t(i)*psi(i))
      b(i)=cbda+fb
c 
      c(i)=abkt*dxa*(f2d1t+f2d2t)-abkt*(f2lam+xl210(i)*f2d1t+f2dtp)
     &+abka*(d1t(i)*hg1+d2t(i)*hg2)
c 
      d(i)=cbda-fb*ba
c 
      e(i)=-c(i)*ba
c 
      f(i)=-a(i)*ba
c 
      g(i)=-abkt*dxa*((f1d1e+f1d2e)*deds+(d1eta(i)+d2eta(i))*deds2)
     &+dxa*abkg*dl22ds+
     &abkt*xl210(i)*(f1d1e*deds+d1eta(i)*deds2)
     &+abkt*(f1dep*deds+d2eta(i)*psi(i)*deds2)-abkg*dl2pds
     &-deltim*deds*(d1eta(i)*hg1+d2eta(i)*hg2)+deltim+xl22(i)*grav*hg2
c 
      goto 420
c
c-----berechnung der koeffizienten fuer den wassergehalt-----
 400  bxc=bxb(i)/cxc(i)
      abdc=axa(i)-bxc*dxd(i)
 
      a(i)=abkt*(bxc*xl210(i)*f3d1e+bxc*f3dep-f3d1e-f3d2e)+abka*bxc*
     &(d1eta(i)*hg1+d2eta(i)*hg2)
c 
      fb=-abkt*(bxc*xl210(i)*d1eta(i)+bxc*d2eta(i)*
     &psi(i)-d1eta(i)-d2eta(i))*2.*adq
      b(i)=abdc+fb
c 
      c(i)=abkt*bxc*xl210(i)*f2d1e+abkt*bxc*f2dep-abka*bxc*(d1eta(i)*hg1
     &-d2eta(i)*hg2)-abkt*(f2d1e+f2d2e)
c 
      d(i)=abdc-fb*ba
c 
      e(i)=-c(i)*ba
c 
      f(i)=-a(i)*ba
c 
      g(i)=abkt*bxc*(-xl210(i)*(f1d1t*dtds+d1t(i)*dtds2)-(f1lam*dtds+
     &xlam(i)*dtds2)-(f1dtp*dtds+d2t(i)*psi(i)*dtds2))+bxc*abkg*
     &dl2pds+bxc*deltim*(d1t(i)*dtds/(2.*sqrt(-sl*
     &zs(i))))*hg1+bxc*deltim*(d2t(i)*dtds/(2.*sqrt(-sl
     &*zs(i)))-xl22(i)*grav)*hg2+abkt*((f1d1t+f1d2t)*dtds+(d1t(i)+
     &d2t(i))*dtds2)-abkg*dl22ds
c 
 420  continue
c 
 520  continue
c
c
      nmin=1
      
      if(ntype.eq.1) then
      do 80 i=2,nsoil-1
 80   x(i)=d(i)*temps(i) + e(i)*temps(i+1) +
     & f(i)*temps(i-1)+g(i)
c 
      ab=1.0
      bb=0.0
      db=tsurf
      at=0.0
      bt=1.0
      dt=temps(nsoil)
      call solveb(ab,bb,db,at,bt,dt,a,b,c,x,temps)
      else
      do 81 i=2,nsoil-1
 81   x(i)=d(i)*eta(i) + e(i)*eta(i+1) +
     & f(i)*eta(i-1)+g(i)
      ab=1.0
      bb=0.0
      db=eta(1)
      at=0.0
      bt=1.0
      dt=eta(nsoil-1)
      
      call solveb(ab,bb,db,at,bt,dt,a,b,c,x,eta)
      endif
c 
  5   continue
      return
      end


      subroutine wert(nb,eta,psi,xk,psit,psieta,psitet,xlam)
      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    ETA,         PSI,         PSIETA
      DOUBLE PRECISION    PSIT,        PSITET,      XK
      DOUBLE PRECISION    XLAM
      INTEGER             NB
      INTEGER             I,           ICK,         ICL,         ICP
      INTEGER             ICPT,        IER,         IER1,        IER2
      INTEGER             IER3,        IER4,        IER5
      DOUBLE PRECISION    CK,          CLAM,        CPSI
      DOUBLE PRECISION    CPSIT,       SEK,          SEL,         SEP
      DOUBLE PRECISION    SEPT,        SK,           SPP,          SPT
      DOUBLE PRECISION    SXL
      INTEGER             ISK,         ISL,         ISP,         ISPT
      !shc end of adding declarations to allow "implicit none"
c      calculate actual soil properties for given moisture at
c      each grid
      dimension eta(nb),psi(nb),xk(nb),psit(nb),
     &psieta(nb),psitet(nb),xlam(nb)
      common/stuetzw/ sep(18),spp(18),sek(17),sk(17),sept(10),
     &spt(10),sel(10),sxl(10),
     &cpsi(17,3),ck(16,3),cpsit(9,3),clam(9,3),
     &isp,isk,ispt,isl
c 
      icp=isp-1
      ick=isk-1
      icpt=ispt-1
      icl=isl-1
c      cubic spline interpolation
      call icsevu(sep,spp,isp,cpsi,icp,eta,psi,nb,ier1)
      call icsevu(sek,sk,isk,ck,ick,eta,xk,nb,ier2)
      call icsevu(sept,spt,ispt,cpsit,icpt,eta,psit,nb,ier3)
      call icsevu(sel,sxl,isl,clam,icl,eta,xlam,nb,ier4)
      call dcsevu(sep,spp,isp,cpsi,icp,eta,psieta,nb,ier5)
      call dcsevu(sept,spt,ispt,cpsit,icpt,eta,psitet,nb,ier)
      do 11 i=1,nb
      psi(i)=-10.**psi(i)
      psieta(i)=2.3026*psi(i)*psieta(i)
  11  xk(i)=10.**xk(i)
c 
      return
      end


      subroutine icsevu(xs,ys,ms,c,mc,x,y,m,ier)
      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    C,           X,           XS,          Y
      DOUBLE PRECISION    YS
      INTEGER             IER,         M,           MC,          MS
      DOUBLE PRECISION    D
      INTEGER             I,           J,           JJ
      DOUBLE PRECISION    TIME
      INTEGER             ITIME,       NUMIT,       NUMRAD,      NUMT1
      INTEGER             NUMT2,       NZEIT
      !shc end of adding declarations to allow "implicit none"
      common/iter/time,itime,nzeit,numit,numt1,numt2,numrad
      dimension xs(ms),ys(ms),c(mc,3),x(m),y(m)

      ier=0
      do 10 i=1,m
      do 11 j=1,ms
      if(x(i).lt.xs(j)) goto 12
  11  continue
  12  if(x(i).lt.xs(1)) write(06,1001)i,x(i),xs(1),ys(1),itime
      if(x(i).gt.xs(ms)) then
      write(06,1002)i,x(i),xs(ms),ys(ms),itime
      write(06,1003)
 1003 format(1x,"xxxxxxxxxxxx")
      endif
      jj=j-1
      d=x(i)-xs(jj)
      y(i)=((c(jj,3)*d+c(jj,2))*d+c(jj,1))*d+ys(jj)
  10  continue
 1001 format(1x,' icsevu: x lt xs(1) ',i3,3e12.4,i6)
 1002 format(1x,' icsevu: x gt xs(ms) ',i3,3e12.4,i6)
      return
      end


      subroutine dcsevu(xs,ys,ms,c,mc,x,y,m,ier)
      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      DOUBLE PRECISION    C,           X,           XS,          Y
      DOUBLE PRECISION    YS
      INTEGER             IER,         M,           MC,          MS
      DOUBLE PRECISION    D
      INTEGER             I,           J,           JJ
      !shc end of adding declarations to allow "implicit none"
      dimension xs(ms),ys(ms),c(mc,3),x(m),y(m)
      ier=0
      do 10 i=1,m
      do 11 j=1,ms
      if(x(i).lt.xs(j)) goto 12
  11  continue
  12  if(x(i).lt.xs(1)) write(06,1001)
      if(x(i).gt.xs(ms)) write(06,1002)
      jj=j-1
      d=x(i)-xs(jj)
      y(i)=(3.0*c(jj,3)*d+2.0*c(jj,2))*d+c(jj,1)
  10  continue
 1001 format(1x,' dcsevu: x lt xs(1) ')
 1002 format(1x,' dcsevu: x gt xs(ms) ')
      return
      end


      subroutine solveb(ab,bb,db,at,bt,dt,a,b,c,d,var)
c      solve tridiagonal system to calculate new values (soil)
      !shc implicit double precision(a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER             NSOIL
      parameter(nsoil=15) !shc moved down from higher up
      DOUBLE PRECISION    A,           AB,          AT
      DOUBLE PRECISION    B,           BB,          BT
      DOUBLE PRECISION    C,           D,           DB,          DT
      DOUBLE PRECISION    VAR
      DOUBLE PRECISION    ALPHA,                    BETA
      DOUBLE PRECISION    XX
      INTEGER             J,           JJ,          M1
      !shc end of adding declarations to allow "implicit none"
      dimension a(nsoil),b(nsoil),c(nsoil),d(nsoil),alpha(nsoil),
     &beta(nsoil),xx(nsoil),var(nsoil)


      alpha(1)=-bb/ab
      beta(1)=db/ab
      m1=nsoil-1
      do 4 j=2,m1
      xx(j)=(a(j)*alpha(j-1)+b(j))
      alpha(j)=-c(j)/xx(j)
   4  beta(j)=-(a(j)*beta(j-1)-d(j))/xx(j)
      var(nsoil)=(dt-at*beta(m1))/(bt+alpha(m1)*at)
      do 5 jj=1,m1
      j=nsoil-jj
   5  var(j)=var(j+1)*alpha(j)+beta(j)
      return
      end
