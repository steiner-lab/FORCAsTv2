

      !This subroutine is used to call MPMPO
      subroutine call_mpmpo(tstamp,deltm,lev,zlev,boxcacm,  &
      tt,parcni,parcno,jjrh0,depr, alwc,is_run_mpmpo) 
      !shc added alwc

      use parameters
      use cacm3_parameters, ONLY:nspec
      use modd_glo, ONLY:naaero !ddw
      !shc use cacm_Parameters !shc
      implicit none !shc

      integer impmpo, icacm, icount, icount2, isoa !shc
      double precision tstamp,  zlev
      double precision depr(nspec)  !ka - added deposition of meoh and etoh 
      real tt, deltm, jjcpt(8), jjcpart(nmpmpo), jjtemp, jjrh, jjlwc, jjhplus
      real jjpar(nmpmpo), jjorganion(nmpmpo), jjgamma_aq_r(1, NAAERO),&
      jjgamma_aq_h(1, NAAERO), diffs(nmpmpo), gas1(nmpmpo), gas2(nmpmpo)
      integer LOGDEV, LAYER,lev       
      integer boxh, option, I 
      real poa
      real conv, conv1
      real totsoa, totdiff
      real totorggas, dummy !shc
      real gas_ratio(nmpmpo) !shc
      double precision parcni(nmpmpo),parcno(nmpmpo),jjrh0
      double precision boxcacm(nspec),boxcacmppm(nspec), cncfct 
      double precision   alwc !shc aerosol liquid water content
      logical is_run_mpmpo !shc
!      logical is_call_mpmpo !shc
      double precision kdep, sum1, sum2, sum3 !shc

! Pass gas-phase concentrations from CACHE to CACM
      cncfct = 1.0E06     
      boxcacmppm = boxcacm*cncfct   ! ppm in CACM
!*********************************************************************
! Deposition for the aerosol phase species 
! At this point, the portion of condensable species that
! are in the gas phase has undergone deposition, but the
! portion that is in the aerosol phase has not; so apply
! deposition to the portion that are in the aerosol phase,
! assuming that they have the same deposition rate as in the gas phase
!*********************************************************************

        if (is_run_mpmpo) then 
           do impmpo = 1, nmpmpo
              sum1=0.0
              sum2=0.0
              sum3=0.0
              icount2 = 0
              do icount = 1, max_nspec_in_surrogate
                 icacm = icacm_in_surrogate(impmpo,icount)
                 if ( icacm > 0 ) then
                       sum1 = sum1 + depr(icacm)*boxcacm(icacm)
                       sum2 = sum2 + boxcacm(icacm)
                       sum3 = sum3 + depr(icacm)
                       icount2 = icount2 + 1
                 end if ! icacm > 0
              end do ! icount
              if ( sum2 > 1.e-18 ) then
                 kdep = sum1/sum2
              else if ( icount2 > 0 ) then
                 kdep = sum3/real(icount2)
              else
                 kdep = 0.
              end if
              ! deposition lost as a first-order loss rate: 
              parcni(impmpo) = parcni(impmpo)*exp(-kdep*deltm)
           end do ! impmpo
        end if

!**********************************
!      Parameters
!**********************************
!      tstamp=tstamp
! Particle phase concentrations passed from CACHE to MPMPO
      do I=1,nmpmpo  !nmpmpo=11
        jjpar(I)=parcni(I)
      enddo

! these are the assumptions
!shc      poa = 2.0 ! assuming that poa conc. is 2 ug/m3
      poa = 0.2 !shc lowered poa conc (~ CMAQ results for July 2006) 
      jjtemp = real(tt) ! use temperature from CACHE !shc uncommented
!shc  if ( option .eq. 1 .and. boxh .ge. 13 .and. boxh .le. 15) then
!shc      jjrh = 1.0
!shc    else 
!shc      jjrh = 0.7
!shc  end if
      jjrh = min(0.99,real(jjrh0)) !shc
!shc  jjlwc = 30. ! assuming jjlwc is 30. ug/m3, liquid water content
      jjlwc = real(alwc) !shc use observed data
      jjhplus = 1.0e-3  ! assuming that the ph value of particle is 3.0
!ddw      jjhplus = 3.91e-2  ! the jjhplus from E-AIM for 0722-0723, 2016 
!conversion factor except the MW from ppm to ug/m3  
      conv = 1.0133e5/(8.314*jjtemp)
      conv1 = 1./conv
      
      jjcpt(1) = poa*0.0047
      jjcpt(2) = poa*0.0010
      jjcpt(3) = poa*0.0053
      jjcpt(4) = poa*0.0080
      jjcpt(5) = poa*0.0012
      jjcpt(6) = poa*0.0
      jjcpt(7) = poa*0.0164
      jjcpt(8) = poa*0.9634
      
! the dimer is treated as POA in the MPMPO ( AERS12)
! dw - the index for aers12 is 30 in CACM-POW, so change 40 to laers12 
      jjcpt(6) = SNGL(boxcacmppm(LAERS12)*conv*186.)

!shc code below has been changed significantly from the original 
      diffs(:) = 0.0
      gas1(:) = 0.0   
      totsoa = 0.0
      do isoa = 1, nmpmpo
         do icount = 1, max_nspec_in_surrogate
            icacm = icacm_in_surrogate(isoa,icount)
            if ( icacm > 0 ) then
               gas1(isoa)   = gas1(isoa) + SNGL(boxcacmppm(icacm)*mw_of_icacm(icacm)*conv)
               diffs(isoa)  = diffs(isoa) + boxcacmppm(icacm)
            end if ! icacm > 0
         end do ! icount
         totsoa  = totsoa + jjpar(isoa)     
         jjcpart(isoa) = gas1(isoa) + jjpar(isoa) ! [=] ug/m3
      end do ! isoa
      totsoa = totsoa + jjcpt(6) ! "AER12 is put in jjcpt(6) but is 2ndary"
      
      if(is_run_mpmpo ) then
         call MAIN_MPMPO(jjcpt, jjcpart, jjtemp, jjrh, jjlwc, jjhplus, &
              jjpar, jjorganion, jjgamma_aq_r,jjgamma_aq_h, LOGDEV, LAYER)

!shc code below has been changed significantly from the original 
         totorggas = 0.0
         do isoa = 1, nmpmpo
            ! SUBROUTINE MPMPO set jjcpart min of 1.0e-12
            if ( jjcpart(isoa) .gt. 1.0e-12 ) then
               gas2(isoa) = jjcpart(isoa) - jjpar(isoa)
               gas_ratio(isoa) = gas2(isoa)/(gas1(isoa)+1.0e-30)
            else
               gas2(isoa) = gas1(isoa)
               gas_ratio(isoa) = 1.0
            endif
            totorggas = gas2(isoa) + totorggas
         end do

         ! pass to main CACHE variable and 
         ! calculate the difference before and after partitioning
         totdiff = totsoa ! total SOA concentration prior to partitioning
         totsoa = 0.0
         do isoa = 1, nmpmpo
            do icount = 1, max_nspec_in_surrogate
               icacm = icacm_in_surrogate(isoa,icount)
               if ( icacm > 0 ) then
                  boxcacmppm(icacm)   = boxcacmppm(icacm)*gas_ratio(isoa)
                  diffs(isoa)  = diffs(isoa) - boxcacmppm(icacm)
               end if ! icacm > 0
            end do ! icount
            totsoa  = totsoa + jjpar(isoa)     
            parcno(isoa) = jjpar(isoa) ! output to main CACHE program; [=] ug/m3
         end do ! isoa
         diffs(:) = -1.0*diffs(:) ! so >0 if condensed; <0 if evaporated
         totsoa = totsoa + jjcpt(6)
         totdiff = totsoa-totdiff
      else ! don't call MPMPO
         gas2(:) = gas1(:)
         parcno(:) = parcni(:)
         jjpar(:) = parcni(:)
         diffs(:) = 0.0
         totdiff=0.0
!ddw         diff(:)=0.0
         totorggas = 0.0
         do isoa = 1, nmpmpo
            totorggas = gas1(isoa) + totorggas
         end do
      endif

! Write out particle phase concentrations for the 11 surrogates
! 85: mpmpo_aer.out
! 86: mpmpo_gas.out
      dummy = 0.0 !shc
      if ((mod(tstamp,dble(1800.0))).eq.0) then 
        write(86,985) tstamp,zlev,(gas2(I) ,I=1,nmpmpo),dummy,totorggas !shc
        if ( is_run_mpmpo ) then
!shc added jjorganion, jjh2oorg; ddw removed jjorganion
           write(85,987) tstamp,zlev,(jjpar(I),I=1,nmpmpo),jjcpt(6),totsoa
           write(87,985) tstamp,zlev,(diffs(I),I=1,nmpmpo),dummy,totdiff   !shc
           write(88,986) tstamp,zlev,(jjgamma_aq_r(1,I),I=1,naaero) !ddw -- output the gamma_raoult
           write(89,986) tstamp,zlev,(jjgamma_aq_h(1,I),I=1,naaero) !ddw -- output the organic phase
           write(90,988) tstamp,zlev,(jjorganion(I),I=1,nmpmpo) !ddw -- output the organic phase
        endif
      endif

!shc 985   format(f10.1,x,f10.4,2x,13(E11.4))
985   format(f10.1,x,f10.4,2x,19(E13.4E3,1x))  !shc
986   format(f10.1,x,f10.4,2x,24(E13.4E3,1x))  !ddw
987   format(f10.1,x,f10.4,2x,24(E13.4E3,1x))  !shc
988   format(f10.1,x,f10.4,2x,19(E13.4E3,1x))  !shc
!shc 986   format(a10,x,a10,2x,13(E11.4))
      do icacm = 1, nspec
         boxcacm(icacm) = boxcacmppm(icacm)/cncfct
      enddo

      end subroutine
