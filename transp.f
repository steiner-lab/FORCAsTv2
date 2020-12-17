c     subroutine wprofil(vg,alphab,wind)
      subroutine wprofil(vg,alphab,wind,ustarb,tstamp)
      USE cacm3_Precision, ONLY:dp 
c Wind profile (roughly based on Baldocchi, 1988)

      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NLEV,NSOIL
      INTEGER lev
      DOUBLE PRECISION    ALPHAB,      TSTAMP,      VG
      DOUBLE PRECISION    WIND, ustarb
      DOUBLE PRECISION    D,           UH
      DOUBLE PRECISION    A, B,  FCT
      DOUBLE PRECISION    Z0M,         ZR
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR

      parameter (nlev=40,nsoil=15)    ! number of levels  
      dimension wind(nlev)
      REAL(KIND=dp) :: Z, ZF, DZ, DZF
      common/height/z(nlev),zf(nlev),dz(nlev),dzf(nlev)
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf

c Get ustarb from measured wind
c          ustarb=vg/30. 
          zr     = 36.4
          d      = 2./3.*hcpy
          z0m    = zrough
          uh = ustarb*dlog((z(levcpy+1)-d)/zrough)/0.4 !estimate    

          a = ustarb*dlog((z(levcpy)-d)/zrough)/0.4
          b = uh*exp(-alphab*(1.-z(levcpy)/(hcpy)))
          fct = a/b   ! fct is the factor to keep continuity

c ######### For test purposes you may modify the parameters of the wind profile
      do lev=1,nlev
       if(z(lev).ge.0.95*hcpy) then
        wind(lev)=ustarb*dlog((z(lev)-d)/zrough)/0.4
       else
        wind(lev)=fct*uh*exp(-alphab*(1.-z(lev)/(hcpy)))  
c        wind(lev)=uh*exp(-alphab*((hcpy-z(lev))/(hcpy-htr)))
       endif
c       wind(levcpy)=0.2*wind(levcpy+1)+0.8*wind(levcpy-1)
      enddo
      wind(1)=0.

c Test for (practically) bare soil
      if(tlai.le.0.05) then
       do lev=1,nlev
        if(lev.ne.1) 
     & wind(lev)=ustarb*dlog((z(lev))/0.01 )/0.4
       enddo
      endif

      return
      end


c***************************************************************
      subroutine atk(vg, rlat,z,wind,akh,ajnn,lprin,tstamp, ustarb)
!ka - added tstamp for diagnostic output
!ddw - added ustarb for akh
c***************************************************************

      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER NLEV
      DOUBLE PRECISION    AJNN,        AKH,         RLAT
      DOUBLE PRECISION    WIND,        Z
      DOUBLE PRECISION    AIRMOL,      AK,          AKHST
      DOUBLE PRECISION    AKM,         AKMST,       ALL,         CORIOL
      DOUBLE PRECISION    DZ,          GRAV,        RI
      DOUBLE PRECISION    STRAT,       THELAM,      TM,          VZBET
      DOUBLE PRECISION    VZBET2,                   XLAM,        XX0
      DOUBLE PRECISION    XXH,         ZM
      DOUBLE PRECISION    zi, zb,  ustarb, vg, acontin, Lmonin
      DOUBLE PRECISION    e, psat,  spechumid 
      INTEGER             IM,          IP,          K
      DOUBLE PRECISION    P,           RELH,        RHOAIR
      DOUBLE PRECISION    T,           THETA,       hvkms, thetav
      DOUBLE PRECISION    HCPY,        HTR,         SIZELF,      TLAI
      DOUBLE PRECISION    ZROUGH
      INTEGER             LEVCPY,      LEVHTR
      DOUBLE PRECISION    ALPHADAY,    ALPHAN,      VGDAY,       VGN
      DOUBLE PRECISION    ZRDAY,       ZRN
      REAL   TIMLOC, DELTIM
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"

c simple formulation of exchange coefficient following Blackadar

      parameter (nlev=40)    ! number of levels  
      common /cpy/ hcpy,htr,tlai,levcpy,levhtr,zrough,sizelf
      common/timpar/ timloc,deltim,month,jday,iyear,icumdy
      dimension z(nlev),wind(nlev),akm(nlev),akh(nlev),hvkms(nlev)
      dimension thetav(nlev)
      common/atmo1/ t(nlev),theta(nlev),p(nlev),rhoair(nlev)
     & ,relh(nlev)
      !shc common/stratpar/ stable,vgday,vgn,zrday,zrn,alphaday,alphan
      common/stratpar/ vgday,vgn,zrday,zrn,alphaday,alphan,stable !shc
      logical stable

      dimension ajnn(nlev),akhst(nlev),akmst(nlev),fct(nlev)
      logical konv,lprin
      data akm/nlev*0./

c ALLOCATE ARRAYS FOR PRINTING TO OUTPUT (AMB, 11/1/10)
      DIMENSION ALL(NLEV), STRAT(NLEV), VZBET2(NLEV), RI(NLEV), 
     +   XLAM(NLEV), KONV(NLEV)
      double precision tstamp,fct 
      
      !dw - added for Stull PBL height
      double precision hsrf,xx0Stull,rhocp  
      DOUBLE PRECISION    EVAPS,       HEATS
      common/cpysrc/ evaps(99),heats(99)
       
c     ========== PARAMETERS =========     
      ak      = 0.4
      grav    = 9.81
      airmol  = 0.1529e-04
      coriol  = 1./3600./6.*3.1415*sin(rlat)
      thelam  = 2.7e-04*wind(nlev-1)/coriol
      
      ! added by ddw for test purpose
      zi      = 1000.0
      xx0     = 0.
      xxh     = 0.
c     ============== SMOOTHING ????  ========= 
      do k=1,nlev
c         e = vcp(k,lh2o)*p(k)
         psat=611.2*exp(17.67*(t(k)-273.15)/(t(k)-29.65))  ! in Pa
         e = relh(k)*psat
         spechumid = 0.622*e/(p(k)-0.378*e)    ! in Kg/Kg
         theta(k)  = t(k)*((100000.0/p(k))** .286)
         thetav(k) = theta(k)*(1.0+0.61*spechumid)
         rhoair(k) = p(k)/287./t(k)
         im=max(k-1,1)
         ip= min(k+1,nlev)
         akhst(k)=0.05*akh(k)+0.90*akh(im)+0.05*akh(ip)+airmol
         akmst(k)=0.05*akm(k)+0.90*akm(im)+0.05*akm(ip)+airmol
      enddo
c     ======================================== 

c     ============== LOOP of xx0 (boundary layer height) at each level ========= 
      do k=1,nlev-1
         im=max(k-1,1)
         ip=min(k+1,nlev)
         dz=z(k+1)-z(k)
         strat(k)=(theta(k+1)-theta(k))/dz
         
         ! Above the canopy
         if(k.ge.levcpy) then
           if(strat(k).le.0.0) then   ! unstable
             xx0=xx0+(z(ip)-z(k))
             konv(k)=(theta(k)+0.1).lt.theta(1)
             ! xxh estimation (don't know what xxh is)
             if(konv(k)) then
               if(xxh.le.0.) xxh=z(k)
             endif
           else                       ! stable 
             konv(k)=.false.
           endif
         endif
      
      enddo

c     =========================    END of the LOOP   =============================      

      ! Stull (1988) boundary layer height based on slab model      
      rhocp = 1200.*p(k)/100000.
      hsrf  = 0.5*(heats(levcpy-1)+heats(levcpy))/rhocp  ! surface heat flux
      xx0Stull = 200.
        
      if(hsrf.gt.0.) then
       xx0Stull=sqrt(2.0*(2.0*0.2+1.0)/0.0065
     & *timloc*3600.*hsrf) 
      !assuming the xx0Stull lasts one more hour after hsrf becomes negative 
      else
         if(timloc.gt.18.0.and.timloc.le.21) xx0Stull = 500.  
      endif

       zi      = xx0Stull 
       zb      = 0.1*zi
       acontin = 1.0/(1.0-zb/zi)**2.0    

c     ====================== LOOP of K at each level =============== 
      do k=1,nlev-1
         dz=z(k+1)-z(k)
         strat(k)=(theta(k+1)-theta(k))/dz
         vzbet=(wind(k+1)-wind(k))/dz
         zm=0.5*(z(k)+z(k+1))
         tm=0.5*(t(k)+t(k+1))
         vzbet2(k)=vzbet*vzbet
         ri(k)=grav*strat(k)/(tm*vzbet2(k))

c        ---------------------- K in boundary layer by ddw ------------ 
         if(zm.lt.zi) then
           ajnn(k) = ak*ustarb*zm
         endif
         if(zm.ge.zb.and.zm.lt.zi) then
           ajnn(k) = acontin*ak*ustarb*zm*(1.0-zm/zi)**2.0
         endif 
                        
c        ---------------------- K in free atmosphere ----------------
         if(zm.ge.zi) then
            xlam(k)=thelam
            all(k)=ak*(zm-0.67*hcpy)/(1.0+ak*(zm-0.67*hcpy)/xlam(k))
            ajnn(k)=(all(k)**2)*sqrt(vzbet2(k)) 
            ajnn(k)= 0.0 
         endif
        
c        ---------------------- K within canopy  ----------------
         if(k.lt.levcpy) then
            xlam(k)=2.
            all(k)=ak*zm/(1.0+ak*zm/xlam(k))
            ajnn(k)=(all(k)**2)*sqrt(vzbet2(k))
         endif
         
cc        ----------------------- ddw Stability correction -----------------      
        ! Calculate Monin-Obukhov length
        hvkms(k)=-akh(k)*(thetav(k+1)-thetav(k))
     &/(z(k+1)-z(k))
          Lmonin  = -(ustarb**3.0*thetav(k))/(ak*grav*hvkms(k)) 
          
          if(Lmonin.gt.0.0) then
            akh(k)=ajnn(k)/(1.0+5.0*zm/Lmonin)
          elseif(Lmonin.lt.0.0) then
            akh(k)=ajnn(k)*(1.0-16.0*zm/Lmonin)**(1.0/2.0)
          endif
c        ---------------------------     smoothing      -----------------      
c         akm(k)=0.1*akmst(k)+0.9*akm(k)
         akh(k)=0.1*akhst(k)+0.9*akh(k)
      
      enddo
c     ======================  END of LOOP  =============== 
      
      if(lprin) then  
        do k  =nlev-1,1,-1
           zm =0.5*(z(k)+z(k+1))
           strat(k)=(theta(k+1)-theta(k))/dz
           vzbet=(wind(k+1)-wind(k))/dz
           vzbet2(k)=vzbet*vzbet
           ri(k)=grav*strat(k)/(tm*vzbet2(k))       
        enddo
      endif
      
      return
      end


      subroutine solve(ab,bb,db,at,bt,dt,a,b,c,d,var,nr)
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER             NLEV
      DOUBLE PRECISION    A,           AB,          AT,          B
      DOUBLE PRECISION    BB,          BT,          C,           D
      DOUBLE PRECISION    DB,          DT,          VAR
      INTEGER             NR
      DOUBLE PRECISION    ALPHA,       BETA,        XX
      INTEGER             I,           J,           JJ,          L
      INTEGER             M1
      !shc end of adding declarations to allow "implicit none"
      save
c
c      solve tridiagonal system to calculate new values (atmosphere)
c 
      parameter(nlev=40)
      dimension a(nlev),b(nlev),c(nlev),d(nlev),
     &alpha(nlev),beta(nlev),xx(nlev),var(nlev)

c      if(var(22).gt.273.)  print*,'ab,bb,at,bt in solve = ', ab,
c     &bb,at,bt

      alpha(nr)=-bb/ab
      beta(nr)=db/ab
      m1=nlev-1

      do 4 j=nr+1,m1
      xx(j)=(a(j)*alpha(j-1)+b(j))
      alpha(j)=-c(j)/xx(j)
   4  beta(j)=-(a(j)*beta(j-1)-d(j))/xx(j)

      var(nlev)=(dt-at*beta(m1))/(bt+alpha(m1)*at)


      do 5 jj=nr,m1
      j=nlev-jj
      var(j)=var(j+1)*alpha(j)+beta(j)
      if(var(j).lt.0.) var(j)=.0001*var(j+1)
   5  continue
      
      

      if(nr.eq.1)return
      l=nr-1
      do 7 i=1,l
   7  var(i)=var(nr)

      return
      end

c **************************************************************
      subroutine newt(akh)                                        
      USE cacm3_Precision, ONLY:dp 
c Subroutine to calculate vertical exchange of energy; updates
c air temperature for each model layer       
c **************************************************************
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
      INTEGER             NLEV
      DOUBLE PRECISION    AKH
      DOUBLE PRECISION    AB,          AT,          BB,          BT
      DOUBLE PRECISION    DB,          DT,          P0,          RHOCP
      DOUBLE PRECISION    SOURCT,                   X
      INTEGER             K,           NR
      DOUBLE PRECISION    P,           RELH,        RHOAIR
      DOUBLE PRECISION    T,           THETA
      DOUBLE PRECISION    EVAPS,       HEATS
      DOUBLE PRECISION    A0,          B0,          C0,          DEL
      DOUBLE PRECISION    ZMIN1,       ZMIN2
      DOUBLE PRECISION    BMFLX,       D,           RNDIV
      DOUBLE PRECISION    RNET,        RNLAM,       TSFC
      DOUBLE PRECISION    UI
      REAL                TIMLOC, deltim
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      REAL(KIND=dp) ::    Z, ZF, DZ, DZF
      !shc end of adding declarations to allow "implicit none"
      save
      parameter(nlev=40)
      common/timpar/timloc,deltim,month,jday,iyear,icumdy
      common/height/z(nlev),zf(nlev),dz(nlev),dzf(nlev)
      common/grid/ del,zmin1,zmin2,a0(nlev),b0(nlev),c0(nlev)
      common/cpysrc/ evaps(99),heats(99)
      common/rad4/d(3,99),ui(3,99),bmflx(3,99),rnet(99),rndiv(99),
     & tsfc,rnlam(3,99)
      common/atmo1/ t(nlev),theta(nlev),p(nlev),rhoair(nlev)
     & ,relh(nlev)
      dimension akh(nlev),
     & sourct(nlev),x(nlev)
      
      
      p0=100000.0
      do k=1,nlev
      theta(k)=t(k)*((p0/p(k))** .286)
      enddo
      
      do k=2,nlev-1
      rhocp=1200.*p(k)/101300.
c     sourct(k)=0.5*(heats(k)+heats(k-1))/rhocp /dz(k)
      sourct(k)=heats(k-1)/rhocp /dz(k)

c a0 coefficient for implicit solver altered by SC/KA (Nov 2014)
c based on Taylor expansion of continuity equation
       a0(k)=-akh(k-1)/dzf(k)/dz(k-1)*deltim
cxc       a0(k)=-akh(k)/dzf(k)/dz(k-1)*deltim
       c0(k)=-akh(k)/dzf(k)/dz(k)*deltim
       b0(k)=1.+(-a0(k)-c0(k))
c
c Crank-Nicholson explicit solution scheme
c but there is likely still an error in here
c  
c       deltak=(akh(k)-akh(k-1))
c     akhm=0.5*akh(k-1)+0.5*akh(k)
c     s=(k-1)*del
c     yy=zmin1*exp(s)
c     xx=yy+zmin2
c     xdsi=1.0/(xx*xx*del*del)
c     xdsih=0.5*xdsi
c New:
c     a0(k)=( deltak*xdsih-akhm*yy/xx*xdsih-akhm*xdsi)*deltim
c     b0(k)=1.+2.*akhm*xdsi*deltim
c     c0(k)=(-deltak*xdsih+akhm*yy/xx*xdsih-akhm*xdsi)*deltim
      x(k)=theta(k)+((p0/p(k))**.286)*sourct(k)*deltim                
      enddo
c      print*, 'sourct(9)=', sourct(9)
c Boundary conditions might still need to be improved   
      nr=1
      at=0.0
      bt=1.0
      dt=theta(nlev)
      ab=1.0
      bb=0.0
      db=(tsfc+273.15)*((p0/p(1))**.286)     !preliminary

!      print*, 'tsfc', tsfc
!      print*, 'dt', dt
!      print*, 'a0', a0
!      print*, 'b0', b0
!      print*, 'c0', c0
!      print*, 'x', sourct(17)
!      print*, 'theta', theta(17)
c **************************************************************
c     solve tridiagonal system to calculate new values of theta    
      call solve(ab,bb,db,at,bt,dt,a0,b0,c0,x,theta,nr)
c **************************************************************
      
!      print*, 'after', theta(17)
      do k=1,nlev
       t(k)=theta(k)/((100000.0/p(k))**.286)
      enddo

      return
      end


c **************************************************************
c subroutine newc(akh,vcnc)
c Subroutine to calculate vertical exchange of matter; updates
c gas- and particle-phase concentrations for each model layer  
c **************************************************************
      subroutine newc(akh,vcp,vcnc,vset) !shc added vset
      USE parameters
!      USE cacm_parameters, ONLY:nspec
      USE cacm3_Precision, ONLY:dp 
      USE module_interf
      !shc implicit double precision (a-h,o-z)
      implicit none !shc
      !shc The following declarations are added so I can use "implicit none"
!      INTEGER             NLEV
      INTEGER             MAXRECT
!      INTEGER             NSPEC, ndspec	!ka - added deposition of meoh and etoh; added ndspec
!      INTEGER             NRECT
      INTEGER             NVARS
      DOUBLE PRECISION    AKH
      DOUBLE PRECISION    VCP
      DOUBLE PRECISION    AB,          AT,          BB,          BT
      DOUBLE PRECISION    CS,          DB,          DTT,          X
      INTEGER             JSPEC,       K,           NR
      DOUBLE PRECISION    P,           RELH,        RHOAIR
      DOUBLE PRECISION    T,           THETA
      DOUBLE PRECISION    A0,          B0,          C0,          DEL
      DOUBLE PRECISION    ZMIN1,       ZMIN2
      real(kind=dp)    DZ,          DZF,         Z,  ZF
      DOUBLE PRECISION    DEPN,                     DEPR
      DOUBLE PRECISION    DEPRAT,                   SOURCE
      REAL                TIMLOC, deltim
      INTEGER             ICUMDY,      IYEAR,       JDAY,        MONTH
      !shc end of adding declarations to allow "implicit none"
      double precision  vset !shc settling velocity for aerosol
!      integer ncacm,nmpmpo,nvarspc,ind_aers12
!      parameter (ncacm=351,nmpmpo=11,nfix=5) !shc ncacm === nspec in cacm_Parameters.f90; updated for cacm2.0 ka
!      parameter (ind_aers12=39) !shc same as cacm_Parameter.f90; updated for cacm2.0 ka
!      parameter (ind_H2O=349) !ka same as cacm_Parameter.f90
      save
!      parameter(nlev=40)
      parameter(maxrect=500)
!      parameter (nspec=84, ndspec=86, nrect=249)!ka - added deposition of meoh and etoh; added ndspec
      parameter (nvars=nspec) !shc ddw 2019
      double precision  a0p(nlev), b0p(nlev), c0p(nlev) !shc Crank-Nicolson coefficients for particles

      common/timpar/timloc,deltim,month,jday,iyear,icumdy
      common/height/z(nlev),zf(nlev),dz(nlev),dzf(nlev)
      common/atmo1/ t(nlev),theta(nlev),p(nlev),rhoair(nlev)
     & ,relh(nlev)
      common/grid/ del,zmin1,zmin2,a0(nlev),b0(nlev),c0(nlev)
      common/quellen/source(nlev,nspec),depn(nlev,nspec),
     & deprat(nlev,nspec),depr(nspec)		!ka - added deposition of meoh and etoh; altered nspec to ndspec

      double precision vcnc
      dimension akh(nlev),vcp(nlev,nspec),vcnc(nlev,nmpmpo)
      dimension x(nlev),cs(nlev)

c Turbulent transport - solved with simple implicit scheme
c Water vapor is included
c---------------------------------------------------------------
c here follows the loop for the gas-phase species
c---------------------------------------------------------------
      do k=2,nlev-1
c a0 coefficient for implicit solver altered by SC/KA (Nov 2014)
c based on Taylor expansion of continuity equation
       a0(k)=-akh(k-1)/dzf(k)/dz(k-1)*deltim
       c0(k)=-akh(k)/dzf(k)/dz(k)*deltim
       b0(k)=1.+(-a0(k)-c0(k))

      enddo

      do 30 jspec=1,nspec                ! species loop for transport
       do k=1,nlev
        cs(k)=vcp(k,jspec)
       enddo

       do k=2,nlev-1
       x(k)= cs(k) + source(k,jspec)*deltim 
       enddo

       at=0.0
       bt=1.0
       dtt=cs(nlev)
       ab=1.0
       bb=0.0
       db=cs(2)  + source(1,jspec)*deltim !preliminary: Lower boundary at soil surface
       if(jspec.eq.lh2o) db=cs(1)  ! Taking into account the boundary conditions!
cxc Boundary conditions might need to be improved  
       nr=1
cxc ----------------------------------------------------------------------
       if(jspec.ne.losd) then ! O3P leads to time-step problems if transported
       call solve(ab,bb,db,at,bt,dtt,a0,b0,c0,x,cs,nr)
       endif
cxc ----------------------------------------------------------------------
  
       do k=1,nlev
        if(cs(k).lt.-1.e-15) write(06,'(" cs < 0. ",2i3,e12.4)')
     &                             jspec, k, cs(k)
        vcp(k,jspec)=max(cs(k),0.D0)
       enddo
   30 continue                     ! End of loop for species
     
cxc----------------------------------------------------------------------
cxc End of the loop for the gases
cxc----------------------------------------------------------------------
c---------------------------------------------------------------
c here follows the loop for the aerosol species
c---------------------------------------------------------------
       do k=2,nlev-1
       a0p(k)= (0.5*vset-akh(k-1)/dzf(k))/dz(k-1)*deltim !shc for particles
       c0p(k)= (-0.5*vset-akh(k)/dzf(k))/dz(k)*deltim !shc for particles
       b0p(k) = b0(k)  !shc for particles, same as for gas
      enddo

      do 31 jspec=1,nmpmpo              ! species loop for transport
       do k=1,nlev
          cs(k)=vcnc(k,jspec)
       enddo

       do k=2,nlev-1
          x(k)= cs(k)  
       enddo

       at=0.0
       bt=1.0
       dtt=cs(nlev)
       ab=1.0
       bb=0.0
       db=cs(2)   !preliminary
       nr=1
cxc ----------------------------------------------------------------------
       call solve(ab,bb,db,at,bt,dtt,a0p,b0p,c0p,x,cs,nr) 
cxc ----------------------------------------------------------------------
  
       do k=1,nlev
        if(cs(k).lt.-1.e-15) write(06,'(" cs < 0. ",2i3,e12.4)')
     & jspec, k, cs(k)
        vcnc(k,jspec)=max(cs(k),0.D0)
       enddo
   31  continue                     ! End of loop for species
     
cxc----------------------------------------------------------------------
cxc End of RACM species loop
cxc----------------------------------------------------------------------

      return
      end
