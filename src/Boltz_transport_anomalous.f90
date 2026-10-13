!-------------------------------------------------------------------!
!> This file contains subroutines that calculate anomalous          !
!> transport properties, including anomalous Hall conductivity,     !
!> spin Hall conductivity and anomalous Nernst coefficient.         !
!>                                                                  !
!> References :                                                     !
!> [1] Physical Review B 74, 195118(2006)        (AHC1)             !
!> [2] Physcial Review Letter 95, 156601 (2005)  (SHC1)             !
!> [3] Physcial Review B 98, 214402 (2018)       (SHC2)             !
!> [4] Physical Review Letter 97, 026603 (2006)  (ANC1)             !
!> restructed by Hanqi Pi on May 2023
!-------------------------------------------------------------------!

   subroutine sigma_AHC
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Hall conductivity !
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. sigma_ahc_eta***meV.txt: sigma tensor for several eta values !
   !> 2. sigma_ahc.gnu: gnuplot script to plot sigma tensor           !
   !>                                                                 !
   !> The output sigma is in the unit of S/cm.                        !
   !>                                                                 !
   !> References: eq(35) of AHC1                                      !
   !>                                                                 !
   !> Dec. 05 2017 by Quansheng Wu @ EPFL                             !
   !> modified on June. 26 2018 by Quansheng Wu @ Airplane from       !
   !> Beijing to Zurich                                               !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: i, m, ie, ieta, iT
     integer :: ierr, knv3
     integer :: NumberofEta

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :, :)

     real(dp), allocatable :: eta_array(:)
     real(dp), allocatable :: T_list(:)

     character*40 :: ahcfilename, etaname

     NumberofEta = 9 

     allocate(eta_array(NumberofEta))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum, NumberofEta))
     sigma_tensor_ahc    = 0d0
      
     eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Eta_Arc


     !> energy
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie

     !> temperature
     do iT=1, NumT
         if (NumT>1) then
            T_list(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
         else
            T_list= Tmin
         endif
     enddo ! iT

     call  sigma_ahc_vary_ChemicalPotential(OmegaNum, energy, NumberofEta, eta_array, sigma_tensor_ahc)

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta=1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(ahcfilename, '(7a)')'sigma_ahc_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=ahcfilename)
           write(outfileindex, '("#",10a)')' Anomalous hall conductivity in unit of S/cm,', 'Brodening eta= ',  trim(adjustl(etaname)), ' meV'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 4)
           write(outfileindex, '("#",a13, 20a16)')'Eenergy (eV)', '\sigma_xy', '\sigma_yz', '\sigma_zx'
           do ie=1, OmegaNum
              write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_tensor_ahc(3, ie, ieta), &
                                                                      sigma_tensor_ahc(1, ie, ieta), &
                                                                      sigma_tensor_ahc(2, ie, ieta)
           
           enddo ! ie
           close(outfileindex)
        enddo
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
        write(ahcfilename, '(7a)')'sigma_ahc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_ahc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ahc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "AHC (S/cm)"'
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(ahcfilename)),  "' u 1:2 w l title '\sigma_{xy}' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(ahcfilename)), "' u 1:3 w l title '\sigma_{zx}' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ahcfilename)), "' u 1:4 w l title '\sigma_{xz}' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif

     deallocate(energy, sigma_tensor_ahc)

     return
  end subroutine sigma_AHC

  subroutine alpha_ANE
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Nernst coefficient!
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. alpha_ane_eta***meV.txt: ane tensor for several eta values   !
   !> 2. alpha_ane.gnu: gnuplot script to plot ane tensor             !
   !>                                                                 !
   !> The first type of files are output on a grid of (T, mu)         !
   !> T is the temperature and mu is the chemical potential.          !
   !> The output coefficient is in the unit of A/(mK).                !
   !>                                                                 !
   !> References: eq(8) of ANC1                                       !
   !>                                                                 !
   !> Dec 14 2022  by Hanqi Pi @ iopcas                               !
   !> Modified at Apr 30 2023 by Hanqi Pi @ iopcas                    !
   !------------------------------------------------------------------!
   use para 
   implicit none 
   
   ! energy refers to the variable to be integrated
   ! mu refers to the chemical potential
   real(dp) :: energy,  omega, T, mu, time_start, time_end

   real(dp), parameter :: ANE_int_interval = 1.0d0*eV2Hartree
   real(dp), parameter :: ANE_int_step     = 0.001d0*eV2Hartree
   
   real(dp), allocatable :: T_list(:)
   real(dp), allocatable :: mu_list(:) ! chemical potential alpha_ane
   real(dp), allocatable :: energy_list(:) ! chemical potential for sigma_ahc (num_step+1)
   real(dp), allocatable :: nernst(:,:,:,:) ! (coordinate,Temp, mu, NumberofEta)
   real(dp), allocatable :: ahc(:,:,:) ! (coordinate, energy, NumberofEta) 
   real(dp), allocatable :: minus_dfde(:,:,:) ! derivative of Fermi-Dirac (num_step+1, Temp, mu)
   integer :: i, ie, iT, imu,ieta, NumberofEta

   ! integral region is (OmegaMin-ANE_int_interval,OmegaMax+ANE_int_interval)
   ! num_step = (OmegaMax-OmegaMin+2*ANE_int_interval)/ANE_int_step
   ! Note: make sure integral region includes Efermi
   integer :: num_step 

   real(dp), allocatable :: eta_array(:)

   character*40 :: anefilename, etaname

   num_step = int((OmegaMax-OmegaMin+2*ANE_int_interval)/ANE_int_step)
   NumberofEta = 9 

   allocate(eta_array(NumberofEta))
   allocate( T_list(NumT) )
   allocate( mu_list(OmegaNum) )
   allocate( energy_list(num_step+1) )
   allocate( nernst(3, NumT, OmegaNum, NumberofEta)) 
   allocate( ahc(3, num_step+1, NumberofEta)) 
   allocate( minus_dfde(num_step+1, NumT, OmegaNum) )
   T_list = 0.0d0
   energy_list = 0.0d0
   nernst = 0.0d0
   ahc = 0.0d0
   minus_dfde = 0.0d0

   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Eta_Arc

   !> energy to be integrated
   do ie = 1,  num_step+1
       energy_list(ie) = (OmegaMin-ANE_int_interval)+ANE_int_step*(ie-1)
   enddo

   !> chemical potential
   do imu=1, OmegaNum
       if (OmegaNum>1) then
          mu_list(imu)= OmegaMin+ (OmegaMax-OmegaMin)* (imu-1d0)/dble(OmegaNum-1)
       else
          mu_list= OmegaMin
       endif
   enddo ! imu

   !> temperature
   do iT=1, NumT
       if (NumT>1) then
          T_list(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
       else
          T_list= Tmin
       endif
   enddo ! iT

   !> get the anomalous hall conductivity
   call now(time_start)
   call sigma_ahc_vary_ChemicalPotential(num_step+1, energy_list, NumberofEta, eta_array, ahc)
   call now(time_end)
   call print_time_cost(time_start, time_end, 'calculation of AHC')

   !> get the derivative of Fermi-Dirac
   call now(time_start)
   do imu = 1, OmegaNum
      mu = mu_list(imu)
      do iT = 1, NumT
         T = T_list(iT)
         do ie = 1, num_step+1
            energy = energy_list(ie)
            call minusdfde_calc_single(energy/eV2Hartree, T*8.617333262E-5, &
                                       mu/eV2Hartree, minus_dfde(ie, iT, imu))
         enddo ! imu
      enddo ! iT
   enddo ! ie
   call now(time_end)
   call print_time_cost(time_start, time_end, 'calculation of derivative of Fermi-Dirac')

   ! --------------------------------------------------------------------
   !  nernst = -int dE 1/e* df/dmu * ahc * (E - mu) / T
   !
   !  We use Hartree as the unit of energy and convert them to eV with eV2Hartree
   !  Therefore the unit of nernst is 1/e * eV * 1/eV * S/cm * eV * 1/K = A/cm/K
   !  To convert it to A/m/K, we need to multiply it by 100
   ! --------------------------------------------------------------------

   !> get the anomalous Nernst coefficient
   call now(time_start)
   do ieta = 1, NumberofEta
      do imu = 1, OmegaNum
         mu = mu_list(imu)
         do iT = 1, NumT
            T = T_list(iT)
            do ie = 1, num_step+1
               energy = energy_list(ie)
               omega = energy - mu
               nernst(:, iT, imu,ieta) = nernst(:, iT, imu,ieta)- minus_dfde(ie,iT,imu)*ahc(:,ie,ieta)*&
                                    (omega/eV2Hartree)/T*(ANE_int_step/eV2Hartree)*100
            enddo ! ie
         enddo ! iT
      enddo ! ie
   enddo ! ieta
   call now(time_end)
   call print_time_cost(time_start, time_end, 'integration of AHC')

   if(cpuid .eq. 0) then
      do ieta =1, NumberofEta
         outfileindex = outfileindex+1
         write(etaname,'(f12.2)') eta_array(ieta)*1000d0/eV2Hartree
         write(anefilename,'(7a)') 'alpha_ane_eta',trim(adjustl(etaname)),'meV.txt'
         open(unit=outfileindex, file=anefilename)
         write(outfileindex, '("#",a)')'Anomalous Nernst coefficient  (A(m-1 K-1))'
         write(outfileindex, "('#column', i7, i13,3i16)")(i, i=1, 5)
         write(outfileindex,'("#",2a13, 3a16)')'Mu(eV)', 'T(K)', '\alpha_xy', '\alpha_yz', '\alpha_zx'
         do imu = 1, OmegaNum
            do iT = 1, NumT
                  write(outfileindex,'(F14.6,F13.3,3E16.6)') mu_list(imu)/eV2Hartree,T_list(iT), nernst(3, iT, imu,ieta),&
                                                                                             nernst(1, iT, imu, ieta),&
                                                                                             nernst(2, iT, imu, ieta)
            enddo
         enddo
         close(outfileindex)
      enddo
   endif

   !> write script for gnuplot
   outfileindex = outfileindex+1
   if(cpuid .eq. 0) then
         write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
         write(anefilename, '(7a)')'alpha_ane_eta', trim(adjustl(etaname)), 'meV.txt'
         open(unit=outfileindex, file='alpha_ane.gnu')
         write(outfileindex, '(a)') 'set terminal pdfcairo enhanced color font ",30" size 13, 6'
         write(outfileindex, '(a)') "set output 'ane.pdf'"
         write(outfileindex, '(a)') 'set key outside'
         write(outfileindex, '(a)') "set palette defined (0 'red', 1 'green')"
         write(outfileindex, '(a)') 'unset colorbox'
         write(outfileindex, '(a, f6.2)') 'Tmin = ',Tmin
         write(outfileindex, '(a, f6.2)') 'Tmax = ',Tmax
         write(outfileindex, '(a, I4)') 'NumT = ',NumT
         write(outfileindex, '(a, f6.2)') 'OmegaMin = ',OmegaMin/eV2Hartree
         write(outfileindex, '(a, f6.2)') 'OmegaMax = ',OmegaMax/eV2Hartree
         write(outfileindex, '(a, I4)') 'OmegaNum = ',OmegaNum
         write(outfileindex, '(a)') ''
         write(outfileindex, '(a)') '#plot the temperature-dependent alpha_yx for the first six chemical potentials'
         write(outfileindex, '(a)') 'set xlabel "T (K)"'
         write(outfileindex, '(a)') 'set ylabel "{/Symbol a}_{yx} (A/(mK))"'
         write(outfileindex, '(a)') 'set ylabel offset 0.0,0'
         write(outfileindex, '(5a)')& 
                                    "plot for [i=0:5] '",trim(adjustl(anefilename)),"' every ::i*NumT::i*NumT+(NumT-1)", &
                                    "u 2:(-$3) w l lt palette frac i/5. title sprintf('{/Symbol m}=%.3f eV', ",&
                                    "OmegaMin+(OmegaMax-OmegaMin)/(OmegaNum*1.0-1.0)*i )"
         write(outfileindex, '(a)') ''
         write(outfileindex, '(a)') '#plot the chemical potential dependent alpha_yx'
         write(outfileindex, '(a)') '#set xlabel "E-E_f (eV)"'
         write(outfileindex, '(a)') '#set ylabel "{/Symbol a}_{yx} (A/(mK))"'
         write(outfileindex, '(a)') '#set ylabel offset 0.0,0'
         write(outfileindex, '(5a)')& 
                                    "#plot for [i=0:NumT-1] '",trim(adjustl(anefilename)),"' every NumT::i", &
                                    " u 1:(-$3) w l lt palette frac i/(NumT*1.0-1.0) title sprintf('T=%.3f K',",&
                                    "Tmin+(Tmax-Tmin)/(NumT*1.0-1.0)*i)"
         close(outfileindex)
   endif

   deallocate(mu_list, energy_list, T_list, nernst, ahc)
   return

end subroutine alpha_ANE

subroutine sigma_ahc_vary_ChemicalPotential(NumOfmu, mulist, NumberofEta, eta_array, sigma_tensor_ahc)
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Hall conductivity !
   !> as a function of chemical potential                             !
   !>                                                                 !
   !> This subroutine is called by sigma_AHC and alpha_ANE            !
   !>                                                                 !
   !> Dec. 05 2017 by Quansheng Wu @ EPFL                             !
   !> modified on June. 26 2018 by Quansheng Wu @ Airplane from       !
   !> Beijing to Zurich                                               !
   !> modified on Apr. 29. 2023 by Hanqi Pi @ Beijing                 !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none

     integer, intent(in)    :: NumOfmu
     real(dp), intent(in)    :: mulist(NumOfmu)
     integer, intent(in)    :: NumberofEta
     real(dp), intent(in)    :: eta_array(NumberofEta)
     real(dp), intent(inout) :: sigma_tensor_ahc(3, NumOfmu, NumberofEta)
     
     integer :: iR, ik, ikx, iky, ikz
     integer :: i, m, ie, ieta
     integer :: ierr, knv3

     real(dp) :: mu, Beta_fake, eta_local
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), allocatable :: Dmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_Ham(:, :, :)

     !> Berry curvature vectors for all bands
     real(dp),allocatable :: Omega_BerryCurv(:, :)
     real(dp),allocatable :: Omega_BerryCurv_t(:, :)

     allocate(Dmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Omega_BerryCurv(Num_wann, 3))
     allocate(Omega_BerryCurv_t(Num_wann, 3))

     allocate( W (Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( sigma_tensor_ahc_mpi(3, NumOfmu, NumberofEta))
     sigma_tensor_ahc_mpi= 0d0
     Hamk_bulk=0d0
     UU= 0d0
      
      knv3= Nk1*Nk2*Nk3

     call now(time_start) 
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
           call now(time_end) 
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
       !call ham_bulk_latticegauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Hamiltonian basis
        call dHdk_atomicgauge_Ham(k, UU, Vmn_Ham)
       !call dHdk_latticegauge_Ham(k, W, UU, Vmn_Ham)

        call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)


        !> calculate Berry curvature at a single k point for all bands
        !> \Omega_n^{\gamma}(k)=i\sum_{\alpha\beta}\epsilon_{\gamma\alpha\beta}(D^{\alpha\dag}D^{\beta})_{nn}
        call berry_curvarture_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
 
        do ieta= 1, NumberofEta
           eta_local = eta_array(ieta)
           !> consider the Fermi-distribution according to the broadening Earc_eta
           Beta_fake= 1d0/eta_local
   
           do ie=1, NumOfmu
              mu = mulist(ie)
              do m= 1, Num_wann
                 Omega_BerryCurv_t(m, :)= Omega_BerryCurv(m, :)*fermi(W(m)-mu, Beta_fake)
              enddo
              sigma_tensor_ahc_mpi(:, ie, ieta)= sigma_tensor_ahc_mpi(:, ie, ieta)- &
                 (sum(Omega_BerryCurv_t(:, :), dim=1))
           enddo ! ie
        enddo ! ieta
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif
   ! --------------------------------------------------------------------
   ! At this point ahc contains
   !
   !           sum_k Omega =  -N*V_c*int dk/[8(pi)^3] Omega
   !
   ! (N is the number of kpoints, V_c is the cell volume). We want
   !
   !           ahc = -(e^2/hbar) int dk/[4(pi)^3] Omega
   ! 
   ! in the version older than 2.5.2, we use Angstrom as the unit as length
   ! Hence we need to multiply by  e^2/hbar/N/V_c = 0.00024341/N/V_c
   ! where 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units.
   ! Then it has units of 
   !
   !           C^2 (J*s)^-1 Angstrom^-3  Angstrom^2 = S/Angstrom                         
   ! 
   ! in the latest version, we use the atomic unit
   ! Similary we need to multiply by  e^2/hbar/N/V_c/Bohr_radius 
   ! and obtain the sigma_tensor_ahc in unit of (Ohm*m)^-1
   ! --------------------------------------------------------------------

      !> in the version older than 2.5.2, we use Angstrom as the unit as length
      !sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume  ! in (Omega*cm)^-1
      
      !> in the latest version, we use the atomic unit
      sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*Echarge**2/hbar/&
         Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume
      !> in unit of (Ohm*cm)^-1
      sigma_tensor_ahc= sigma_tensor_ahc/100d0
      
      ! deallocate( W, Hamk_bulk, UU)
      ! deallocate( sigma_tensor_ahc_mpi)
  end subroutine sigma_ahc_vary_ChemicalPotential

  subroutine sigma_SHC
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the spin Hall conductivity      !
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. sigma_shc_eta***meV.txt: sigma tensor for several eta values !
   !> 2. sigma_shc.gnu: gnuplot script to plot sigma tensor           !
   !> 3. shc_berrycurvature_xy_z.dat: BZ map for spin-up/down part    !
   !> 4. shc_berrycurvature_xy_z.gnu: gnuplot script for BZ map        !
   !>                                                                 !
   !> The output sigma is in the unit of (hbar/e)S/cm.                !
   !>                                                                 !
   !> References: SHC1 and SHC2                                       !
   !>                                                                 !
   !> Dec. 05 2022 by Quansheng Wu @ Beijing                          !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none
    
      integer :: ik, ikx, iky, ikz, ieta
      integer :: m, n, i, j, ie, ialpha, ibeta, igamma
      integer :: ierr, knv3, nkmesh2, nwann
      integer :: NumberofEta
 
       real(dp) :: mu, Beta_fake, deno_fac, eta_local, vmin, vmax, fac_nk1, fac_nk2
      real(dp) :: k(3), k_cart(3), kxy_plane(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     !> energy  dim= OmegaNum
     real(dp), allocatable :: energy(:)

     !> sigma^gamma_{alpha, beta}, alpha, beta, gamma=1,2,3 for x, y, z
     !>  sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
     real(dp), allocatable :: sigma_tensor_shc(:, :, :, :, :)
     real(dp), allocatable :: sigma_tensor_shc_mpi(:, :, :, :, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

      complex(dp), allocatable :: Vmn_Ham(:, :, :)
      complex(dp), allocatable :: Vmn_wann(:, :, :)
      complex(dp), allocatable :: spin_sigma(:, :)
      complex(dp), allocatable :: projector_spin_up(:, :)
      complex(dp), allocatable :: projector_spin_dn(:, :)
      complex(dp), allocatable :: j_spin_gamma_alpha(:, :)
      complex(dp), allocatable :: j_spin_up_xy(:, :)
      complex(dp), allocatable :: j_spin_dn_xy(:, :)
      complex(dp), allocatable :: mat_t(:, :)

     !> Berry curvature vectors for all bands
      real(dp),allocatable :: Omega_spin(:)
      real(dp),allocatable :: Omega_spin_t(:)
      real(dp),allocatable :: Omega_spin_up(:)
      real(dp),allocatable :: Omega_spin_dn(:)
      real(dp),allocatable :: Omega_spin_up_t(:)
      real(dp),allocatable :: Omega_spin_dn_t(:)

      real(dp), allocatable :: Omega_k_up(:), Omega_k_dn(:), Omega_k_total(:)
      real(dp), allocatable :: Omega_k_up_mpi(:), Omega_k_dn_mpi(:), Omega_k_total_mpi(:)

     ! spin operator matrix spin_sigma_x,spin_sigma_y in spin_sigma_z representation
     complex(Dp),allocatable :: pauli_matrices(:, :, :) 
    
     real(dp), allocatable :: eta_array(:)
     character(80) :: shcfilename, etaname

     NumberofEta=9

      allocate(Vmn_Ham(Num_wann, Num_wann, 3))
      allocate(Vmn_wann(Num_wann, Num_wann, 3))
      allocate(spin_sigma(Num_wann, Num_wann))
      allocate(projector_spin_up(Num_wann, Num_wann))
      allocate(projector_spin_dn(Num_wann, Num_wann))
      allocate(j_spin_gamma_alpha(Num_wann, Num_wann))
      allocate(j_spin_up_xy(Num_wann, Num_wann))
      allocate(j_spin_dn_xy(Num_wann, Num_wann))
      allocate(Omega_spin(Num_wann))
      allocate(Omega_spin_t(Num_wann))
      allocate(Omega_spin_up(Num_wann))
      allocate(Omega_spin_dn(Num_wann))
      allocate(Omega_spin_up_t(Num_wann))
      allocate(Omega_spin_dn_t(Num_wann))

     allocate( W (Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( mat_t(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate(eta_array(NumberofEta))
     allocate( sigma_tensor_shc    (OmegaNum, 3, 3, 3, NumberofEta))
     allocate( sigma_tensor_shc_mpi(OmegaNum, 3, 3, 3, NumberofEta))

     allocate(pauli_matrices(Num_wann, Num_wann, 3))
     spin_sigma= 0d0
     sigma_tensor_shc    = 0d0
      sigma_tensor_shc_mpi= 0d0
      Hamk_bulk=0d0
      UU= 0d0
      Vmn_wann= 0d0
      pauli_matrices= 0d0
      projector_spin_up=0d0
      projector_spin_dn=0d0
      j_spin_up_xy=0d0
      j_spin_dn_xy=0d0
 
     eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Eta_Arc

     nwann= Num_wann/2
     !> spin operator matrix
     !> this part is package dependent. 
    !if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
    !   .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
        do j=1, nwann
           pauli_matrices(j, nwann+j, 1)=1.0d0
           pauli_matrices(j+nwann, j, 1)=1.0d0
           pauli_matrices(j, nwann+j, 2)=-zi
           pauli_matrices(j+nwann, j, 2)=zi
           pauli_matrices(j, j, 3)= 1d0
           pauli_matrices(j+nwann, j+nwann, 3)=-1d0
        enddo
    !else
    !   if (cpuid.eq.0) write(stdout, *)'Error: please report your software and wannier90.wout to me'
    !   if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
    !   stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
    !endif

   
     !> energy range (chemical potential range)
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie

      knv3= Nk1*Nk2*Nk3
      nkmesh2= Nk1*Nk2
      allocate(Omega_k_up(nkmesh2), Omega_k_dn(nkmesh2), Omega_k_total(nkmesh2))
      allocate(Omega_k_up_mpi(nkmesh2), Omega_k_dn_mpi(nkmesh2), Omega_k_total_mpi(nkmesh2))
      Omega_k_up=0d0
      Omega_k_dn=0d0
      Omega_k_total=0d0
      Omega_k_up_mpi=0d0
      Omega_k_dn_mpi=0d0
      Omega_k_total_mpi=0d0

      call now(time_start) 
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
           call now(time_end) 
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
       !call ham_bulk_latticegauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU= Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Wannier basis
        !> \partial_k H_nm
        call dHdk_atomicgauge(k, Vmn_wann)

         do ialpha= 1, 3
            call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
         enddo

         !> spin axis igamma= x, y, z
         do igamma= 1, 3
           !> set Pauli matrix
           spin_sigma= pauli_matrices(:, :, igamma)
   
           !> calculate spin current operator j_spin_gamma_alpha^l_alpha= 1/2*{Sigma_gamma, v_alpha} 
           ! in order to calculate SHC^z_xy
           do ialpha= 1, 3
              !> in Wannier basis
              j_spin_gamma_alpha= 0d0
              call mat_mul(Num_wann, spin_sigma, Vmn_wann(:, :, ialpha), j_spin_gamma_alpha(:, :))
              call mat_mul(Num_wann, Vmn_wann(:, :, ialpha), spin_sigma, mat_t)
              j_spin_gamma_alpha(:, :)= j_spin_gamma_alpha(:, :)+ mat_t
              j_spin_gamma_alpha= j_spin_gamma_alpha/2d0
         
              !> rotate to Hamiltonian basis
              mat_t= j_spin_gamma_alpha(:, :)
              call rotation_to_Ham_basis(UU, mat_t, j_spin_gamma_alpha(:, :))

              do ieta= 1, NumberofEta
                 eta_local= eta_array(ieta)
                 !> \Omega_spin^l_n^{\gamma}(k)=-2\sum_{m}*aimag(Im({js(\gamma),v(\alpha)}/2)_nm*v_beta_mn))/((w(n)-w(m))^2+eta_arc^2)
                 do ibeta= 1, 3
                    Omega_spin= 0d0
                    do n= 1, Num_wann
                       do m= 1, Num_wann
                          if (abs(W(n)-W(m))<eps9 ) cycle
                          deno_fac= -2d0/((W(n)-W(m))**2+ eta_local**2)
                          Omega_spin(n)= Omega_spin(n)+ &
                             aimag(j_spin_gamma_alpha(n, m)*Vmn_Ham(m, n, ibeta))*deno_fac
                       enddo
                    enddo
        
                    !> consider the Fermi-distribution according to the broadening Earc_eta
                    Beta_fake= 1d0/eta_local
            
                    do ie=1, OmegaNum
                       mu = energy(ie)
                       do n= 1, Num_wann
                          Omega_spin_t(n)= Omega_spin(n)*fermi(W(n)-mu, Beta_fake)
                       enddo
   
                       !> sum over all "spin" Berry curvature below chemical potential mu
                       sigma_tensor_shc_mpi(ie, igamma, ialpha, ibeta, ieta)= &
                          sigma_tensor_shc_mpi(ie, igamma, ialpha, ibeta, ieta)+ &
                          sum(Omega_spin_t(:))
                    enddo ! ie
                 enddo ! ibeta  v
              enddo ! ieta=1, NumberofEta
           enddo ! ialpha  j
        enddo ! igamma  spin
     enddo ! ik

#if defined (MPI)
      call mpi_allreduce(sigma_tensor_shc_mpi,sigma_tensor_shc,size(sigma_tensor_shc),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
      sigma_tensor_shc= sigma_tensor_shc_mpi
#endif

     !> in the latest version, we use the atomic unit
     !> in unit of ((hbar/e)(Ohm*m)^-1
     sigma_tensor_shc= sigma_tensor_shc/dble(knv3)/Origin_cell%CellVolume*&
        Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume/2d0

     !> in unit of ((hbar/e)(Ohm*cm)^-1
     sigma_tensor_shc= sigma_tensor_shc/100d0


     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta=1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(shcfilename, '(7a)')'sigma_shc_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=shcfilename)
           write(outfileindex, '("#",10a)')' Spin hall conductivity in unit of (hbar/e)S/cm,', 'Brodening eta= ',  trim(adjustl(etaname)), ' meV'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 28)
           write(outfileindex, '("#",a13, 27a16)')'Eenergy (eV)', &
             'xx^x', 'xy^x', 'xz^x', 'yx^x', 'yy^x', 'yz^x', 'zx^x', 'zy^x', 'zz^x', &
             'xx^y', 'xy^y', 'xz^y', 'yx^y', 'yy^y', 'yz^y', 'zx^y', 'zy^y', 'zz^y', &
             'xx^z', 'xy^z', 'xz^z', 'yx^z', 'yy^z', 'yz^z', 'zx^z', 'zy^z', 'zz^z'
           do ie=1, OmegaNum
              write(outfileindex, '(E16.8)', advance='no')energy(ie)/eV2Hartree
              do igamma=1, 3
              do ialpha=1, 3
              do ibeta =1, 3
                 if (ialpha*ibeta*igamma/=27)then
                    write(outfileindex, '(200E16.8)', advance='no') sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
                 else
                    write(outfileindex, '(200E16.8)', advance='yes') sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
                 endif
              enddo
              enddo
              enddo
           enddo
        enddo ! ieta
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
      if (cpuid==0) then
         write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
         write(shcfilename, '(7a)')'sigma_shc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_shc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_shc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "SHC (\hbar/e)S/cm"'
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(shcfilename)),  "' u 1:21 w l title '\sigma_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(shcfilename)), "' u 1:17 w l title '\sigma_{zx}^y' lc rgb 'blue' lw 4, \"
         write(outfileindex, '(5a)')"'", trim(adjustl(shcfilename)), "' u 1:13 w l title '\sigma_{xz}^y' lc rgb 'orange' lw 4 "
         close(outfileindex)
      endif

      !> spin Berry curvature map on KPLANE_BULK for sigma_xy^z
      Omega_k_up=0d0
      Omega_k_dn=0d0
      Omega_k_total=0d0
      fac_nk1=0d0
      fac_nk2=0d0
      if (Nk1>1) fac_nk1=1d0/dble(Nk1-1)
      if (Nk2>1) fac_nk2=1d0/dble(Nk2-1)
      eta_local= max(Eta_Arc, eps9)
      Beta_fake= 1d0/eta_local

      do ik= 1+ cpuid, nkmesh2, num_cpu
         ikx= (ik-1)/Nk2+1
         iky= ik-(ikx-1)*Nk2
         k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
          + K3D_vec2*(iky-1)*fac_nk2- (K3D_vec1+K3D_vec2)/2d0

         call ham_bulk_atomicgauge(k, Hamk_bulk)
         UU= Hamk_bulk
         call eigensystem_c( 'V', 'U', Num_wann, UU, W)
         call dHdk_atomicgauge(k, Vmn_wann)
         do ialpha= 1, 3
            call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
         enddo

         spin_sigma= pauli_matrices(:, :, 3)
         projector_spin_up= 0d0
         projector_spin_dn= 0d0
         do n=1, Num_wann
            projector_spin_up(n, n)= 0.5d0
            projector_spin_dn(n, n)= 0.5d0
         enddo
         projector_spin_up= projector_spin_up+ spin_sigma/2d0
         projector_spin_dn= projector_spin_dn- spin_sigma/2d0

         j_spin_up_xy= 0d0
         j_spin_dn_xy= 0d0
         call mat_mul(Num_wann, projector_spin_up, Vmn_wann(:, :, 1), j_spin_up_xy(:, :))
         call mat_mul(Num_wann, Vmn_wann(:, :, 1), projector_spin_up, mat_t)
         j_spin_up_xy= (j_spin_up_xy+ mat_t)/2d0
         call rotation_to_Ham_basis(UU, j_spin_up_xy, j_spin_up_xy)

         call mat_mul(Num_wann, projector_spin_dn, Vmn_wann(:, :, 1), j_spin_dn_xy(:, :))
         call mat_mul(Num_wann, Vmn_wann(:, :, 1), projector_spin_dn, mat_t)
         j_spin_dn_xy= (j_spin_dn_xy+ mat_t)/2d0
         call rotation_to_Ham_basis(UU, j_spin_dn_xy, j_spin_dn_xy)

         Omega_spin_up= 0d0
         Omega_spin_dn= 0d0
         do n= 1, NumOccupied
            do m= NumOccupied+1, Num_wann
               if (abs(W(n)-W(m))<eps9 ) cycle
               deno_fac= -2d0/((W(n)-W(m))**2+ eta_local**2)
               Omega_spin_up(n)= Omega_spin_up(n)+ &
                  aimag(j_spin_up_xy(n, m)*Vmn_Ham(m, n, 2))*deno_fac
               Omega_spin_dn(n)= Omega_spin_dn(n)+ &
                  aimag(j_spin_dn_xy(n, m)*Vmn_Ham(m, n, 2))*deno_fac
            enddo
         enddo
         Omega_k_up(ik)= sum(Omega_spin_up(1:NumOccupied))
         Omega_k_dn(ik)= sum(Omega_spin_dn(1:NumOccupied))
         Omega_k_total(ik)= Omega_k_up(ik)-Omega_k_dn(ik)
      enddo

#if defined (MPI)
      call mpi_allreduce(Omega_k_up,Omega_k_up_mpi,size(Omega_k_up_mpi),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(Omega_k_dn,Omega_k_dn_mpi,size(Omega_k_dn_mpi),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(Omega_k_total,Omega_k_total_mpi,size(Omega_k_total_mpi),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
      Omega_k_up_mpi= Omega_k_up
      Omega_k_dn_mpi= Omega_k_dn
      Omega_k_total_mpi= Omega_k_total
#endif

      !> write BZ distribution of spin Berry curvature for sigma_xy^z
      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='shc_berrycurvature_xy_z.dat')
         write(outfileindex, '(a)')'# Spin Berry curvature for sigma_xy^z, unit: Angstrom^2'
         write(outfileindex, '(a)')'# Columns 7-9: spin-up, spin-down, total(up-down)'
         write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                             "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
                                             'Omega_up', 'Omega_dn', 'Omega_tot'
         write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 9)
         fac_nk1=0d0
         fac_nk2=0d0
         if (Nk1>1) fac_nk1=1d0/dble(Nk1-1)
         if (Nk2>1) fac_nk2=1d0/dble(Nk2-1)
         ik=0
         do ikx= 1, nk1
            do iky= 1, nk2
               ik=ik+1
               k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
                + K3D_vec2*(iky-1)*fac_nk2- (K3D_vec1+K3D_vec2)/2d0
               k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
               call rotate_k3_to_kplane(k_cart, kxy_plane)
               write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, kxy_plane*Angstrom2atomic, &
                  Omega_k_up_mpi(ik)/Angstrom2atomic**2, &
                  Omega_k_dn_mpi(ik)/Angstrom2atomic**2, &
                  Omega_k_total_mpi(ik)/Angstrom2atomic**2
               ! write(outfileindex, *) ' '
            enddo
            write(outfileindex, *) ' '
         enddo
         close(outfileindex)
      endif

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         vmin=minval(Omega_k_total_mpi)/Angstrom2atomic**2
         vmax=maxval(Omega_k_total_mpi)/Angstrom2atomic**2
         open(unit=outfileindex, file='shc_berrycurvature_xy_z.gnu')
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
         write(outfileindex, '(a)')"set output 'shc_berrycurvature_xy_z.png'"
         write(outfileindex, '(a)')'if (!exists("MP_LEFT"))   MP_LEFT = .12'
         write(outfileindex, '(a)')'if (!exists("MP_RIGHT"))  MP_RIGHT = .92'
         write(outfileindex, '(a)')'if (!exists("MP_BOTTOM")) MP_BOTTOM = .12'
         write(outfileindex, '(a)')'if (!exists("MP_TOP"))    MP_TOP = .88'
         write(outfileindex, '(a)')'if (!exists("MP_GAP"))    MP_GAP = 0.08'
         write(outfileindex, '(a)')'set multiplot layout 1,3 rowsfirst \'
         write(outfileindex, '(a)')"              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"
         write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
         write(outfileindex, '(a)')"unset ztics"
         write(outfileindex, '(a)')"unset key"
         write(outfileindex, '(a)')"set pm3d"
         write(outfileindex, '(a, E10.3, a, E10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
         write(outfileindex, '(a)')"set view map"
         write(outfileindex, '(a)')"set size ratio -1"
         write(outfileindex, '(a)')"set border lw 3"
         write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
         write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
         write(outfileindex, '(a)')"set xrange [] noextend"
         write(outfileindex, '(a)')"set yrange [] noextend"
         write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
         write(outfileindex, '(a)')"set pm3d interpolate 2,2"
         write(outfileindex, '(a)')"set title 'Spin-up {/Symbol W}_{xy}^z ({\305}^2)'"
         write(outfileindex, '(a)')"splot 'shc_berrycurvature_xy_z.dat' u 4:5:7 w pm3d"
         write(outfileindex, '(a)')"unset ylabel"
         write(outfileindex, '(a)')"unset ytics"
         write(outfileindex, '(a)')"set title 'Spin-down {/Symbol W}_{xy}^z ({\305}^2)'"
         write(outfileindex, '(a)')"splot 'shc_berrycurvature_xy_z.dat' u 4:5:8 w pm3d"
         write(outfileindex, '(a)')"set colorbox"
         write(outfileindex, '(a)')"set title 'Total {/Symbol W}_{xy}^z ({\305}^2)'"
         write(outfileindex, '(a)')"splot 'shc_berrycurvature_xy_z.dat' u 4:5:9 w pm3d"
         write(outfileindex, '(a)')"unset multiplot"
         close(outfileindex)
      endif



      deallocate( W, Hamk_bulk, UU, energy)
      deallocate( sigma_tensor_shc, sigma_tensor_shc_mpi)
      deallocate(Vmn_Ham, Vmn_wann, spin_sigma, projector_spin_up, projector_spin_dn)
      deallocate(j_spin_gamma_alpha, j_spin_up_xy, j_spin_dn_xy, mat_t)
      deallocate(Omega_spin, Omega_spin_t, Omega_spin_up, Omega_spin_dn, Omega_spin_up_t, Omega_spin_dn_t)
      deallocate(pauli_matrices, eta_array)
      deallocate(Omega_k_up, Omega_k_dn, Omega_k_total, Omega_k_up_mpi, Omega_k_dn_mpi, Omega_k_total_mpi)
 
     return
  end subroutine sigma_SHC

  subroutine Sz_expectation
   !------------------------------------------------------------------!
   !> Calculate <Sz> expectation value summed over occupied bands on  !
   !> a 2D k plane in the BZ.                                         !
   !>                                                                   !
   !> The Pauli sigma_z is built in the Wannier basis with             !
   !> eigenvalues +/-1 for spin up/down. We sum                        !
   !> <u_n| sigma_z |u_n>/2 over occupied bands to get <Sz> in units  !
   !> of hbar.                                                          !
   !>                                                                   !
   !> It produces the following output files:                           !
   !> 1. sz_expectation.dat: <Sz> expectation value map                !
   !> 2. sz_expectation.gnu: gnuplot script for the map                 !
   !------------------------------------------------------------------!
     use wmpi
     use para
     implicit none

     integer :: ik, ikx, iky, n, i
     integer :: ierr, nkmesh2
     real(dp) :: k(3), k_cart(3), kxy_plane(3)
     real(dp) :: fac_nk1, fac_nk2

     real(dp) :: time_start, time_end

     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     complex(dp), allocatable :: sigma_z(:, :)
     complex(dp), allocatable :: Sz_Ham(:, :)

     real(dp), allocatable :: Sz_occ(:)
     real(dp), allocatable :: Sz_occ_mpi(:)
     real(dp), allocatable :: Sz_top(:)
     real(dp), allocatable :: Sz_top_mpi(:)

     !> only valid for SOC systems
     if (SOC==0) then
        if (cpuid==0) write(stdout, *)'ERROR: Sz_expectation_calc requires SOC>0 (set SOC in namelist SYSTEM)'
        stop
     endif

     if (Numoccupied<=0) then
        if (cpuid==0) write(stdout, *)'ERROR: Sz_expectation_calc requires Numoccupied>0, please set it in namelist SYSTEM'
        stop
     endif

     nkmesh2= Nk1*Nk2

     allocate(W(Num_wann))
     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate(sigma_z(Num_wann, Num_wann))
     allocate(Sz_Ham(Num_wann, Num_wann))
     allocate(Sz_occ(nkmesh2))
     allocate(Sz_occ_mpi(nkmesh2))
     allocate(Sz_top(nkmesh2))
     allocate(Sz_top_mpi(nkmesh2))

     W= 0d0
     Hamk_bulk= 0d0
     UU= 0d0
     sigma_z= 0d0
     Sz_Ham= 0d0
     Sz_occ= 0d0
     Sz_occ_mpi= 0d0
     Sz_top= 0d0
     Sz_top_mpi= 0d0

     !> build sigma_z/2 in Wannier basis (so that <sigma_z>/2 = <S_z> in units of hbar)
     !> basis ordering assumed: (spin-up block of size Num_wann/2, spin-down block of size Num_wann/2)
     do n= 1, Num_wann/2
        sigma_z(n, n)= 0.5d0
        sigma_z(n+Num_wann/2, n+Num_wann/2)= -0.5d0
     end do

     if (Nk1>1) fac_nk1= 1d0/dble(Nk1-1)
     if (Nk2>1) fac_nk2= 1d0/dble(Nk2-1)

     time_start= 0d0
     time_end= 0d0

     do ik= 1+cpuid, nkmesh2, num_cpu
        ikx= (ik-1)/Nk2+1
        iky= ik-(ikx-1)*Nk2
        k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1 &
         + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0

        !> calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_atomicgauge(k, Hamk_bulk)

        !> diagonalize
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)

        !> rotate sigma_z to Hamiltonian basis: Sz_Ham = UU_dag * sigma_z * UU
        call rotation_to_Ham_basis(UU, sigma_z, Sz_Ham)

        !> sum <Sz> over occupied bands; also store the top occupied band
        Sz_occ(ik)= 0d0
        do n= 1, Numoccupied
           Sz_occ(ik)= Sz_occ(ik)+ dble(Sz_Ham(n, n))
        end do
        Sz_top(ik)= dble(Sz_Ham(Numoccupied, Numoccupied))

     end do

#if defined (MPI)
     call mpi_allreduce(Sz_occ, Sz_occ_mpi, size(Sz_occ), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
     call mpi_allreduce(Sz_top, Sz_top_mpi, size(Sz_top), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     Sz_occ_mpi= Sz_occ
     Sz_top_mpi= Sz_top
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='sz_expectation.dat')
        write(outfileindex, '(a)')'# <Sz> expectation value summed over occupied bands'
        write(outfileindex, '(a)')'# Units of hbar (i.e. Pauli sigma_z/2, eigenvalues +/-1/2 per electron)'
        write(outfileindex, '(a)')'# Sz_avg = Sz_occ/NumOccupied (per-band average); Sz_top = <Sz> of the highest occupied band'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
                                            'Sz_occ (hbar)', 'Sz_avg (hbar)', 'Sz_top (hbar)'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 9)
        ik= 0
        do ikx= 1, Nk1
           do iky= 1, Nk2
              ik= ik+1
              k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1 &
               + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0
              k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
              call rotate_k3_to_kplane(k_cart, kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, &
                 kxy_plane*Angstrom2atomic, Sz_occ_mpi(ik), &
                 Sz_occ_mpi(ik)/dble(Numoccupied), Sz_top_mpi(ik)
           end do
           write(outfileindex, *) ' '
        end do
        close(outfileindex)
     endif

     !> write gnuplot script
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='sz_expectation.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 2560, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'sz_expectation.png'"
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title '<S_z> expectation value (\hbar)'"
        write(outfileindex, '(a)')"splot 'sz_expectation.dat' u 4:5:7 w pm3d"
        close(outfileindex)
     endif

     deallocate(W, Hamk_bulk, UU, sigma_z, Sz_Ham, Sz_occ, Sz_occ_mpi, Sz_top, Sz_top_mpi)

   return
  end subroutine Sz_expectation


  subroutine orbital_real_basis_coeff(label, l, coeff, supported)
    !------------------------------------------------------------------!
    !> Return coefficients of real cubic harmonics in |l,m> basis.    !
    !> coeff(m) is the coefficient of |l,m>, where m=-l...l.          !
    !------------------------------------------------------------------!
     use para
     implicit none

     character(*), intent(in) :: label
     integer, intent(out) :: l
     complex(dp), intent(out) :: coeff(-3:3)
     logical, intent(out) :: supported

     real(dp) :: invsqrt2
     character(16) :: label_u

     coeff= 0d0
     supported= .true.
     l= -1
     invsqrt2= 1d0/sqrt(2d0)
     label_u= trim(adjustl(label))

     select case(trim(label_u))
     case('S')
        l= 0
        coeff(0)= 1d0
     case('PZ', 'Z')
        l= 1
        coeff(0)= 1d0
     case('PX', 'X')
        l= 1
        coeff(-1)= invsqrt2
        coeff(1)= -invsqrt2
     case('PY', 'Y')
        l= 1
        coeff(-1)= zi*invsqrt2
        coeff(1)= zi*invsqrt2
     case('DZ2', 'Z2', 'D3Z2-R2', 'D3Z2R2')
        l= 2
        coeff(0)= 1d0
     case('DXY', 'XY')
        l= 2
        coeff(-2)= zi*invsqrt2
        coeff(2)= -zi*invsqrt2
     case('DXZ', 'DZX', 'XZ', 'ZX')
        l= 2
        coeff(-1)= invsqrt2
        coeff(1)= -invsqrt2
     case('DYZ', 'YZ')
        l= 2
        coeff(-1)= zi*invsqrt2
        coeff(1)= zi*invsqrt2
     case('DX2-Y2', 'DX2Y2', 'X2-Y2', 'X2Y2')
        l= 2
        coeff(-2)= invsqrt2
        coeff(2)= invsqrt2
     case default
        supported= .false.
     end select

   return
  end subroutine orbital_real_basis_coeff

  subroutine build_orbital_operator_matrices(orbital_matrices)
    !------------------------------------------------------------------!
    !> Build Lx, Ly, Lz matrices in the Wannier basis from PROJECTORS. !
    !------------------------------------------------------------------!
     use wmpi
     use para
     implicit none

     complex(dp), intent(out) :: orbital_matrices(Num_wann, Num_wann, 3)

     integer :: i, j, ibase, jbase
     integer :: isector, jsector, nspinless
     integer :: iatom, jatom, iproj, jproj
     integer :: l_i, l_j, m1, m2
     complex(dp) :: coeff_i(-3:3), coeff_j(-3:3)
     complex(dp) :: lx, ly, lz, lp, lm
     logical :: has_spin, supported_i, supported_j
     character(10) :: label_i, label_j

     orbital_matrices= 0d0

     has_spin= (SOC>0 .or. Add_Zeeman_Field)
     nspinless= Num_wann
     if (has_spin) nspinless= Num_wann/2

     if (has_spin .and. (2*nspinless/=Num_wann)) then
        if (cpuid==0) then
           write(stdout, '(a)')'ERROR in OHC_calc: Num_wann is not compatible with spin duplication.'
        endif
        stop
     endif

     do i= 1, Num_wann
        ibase= i
        isector= 1
        if (has_spin .and. i>nspinless) then
           ibase= i-nspinless
           isector= 2
        endif

        iatom= Origin_cell%spinorbital_to_atom_index(ibase)
        iproj= Origin_cell%spinorbital_to_projector_index(ibase)
        label_i= trim(adjustl(Origin_cell%proj_name(iproj, iatom)))
        call orbital_real_basis_coeff(label_i, l_i, coeff_i, supported_i)
        if (.not.supported_i) then
           if (cpuid==0) then
              write(stdout, '(a, a)')'ERROR in OHC_calc: unsupported projector label: ', trim(label_i)
              write(stdout, '(a)')'Supported labels include s, pz/px/py, dxy/dyz/dxz(dxz)/dx2-y2/dz2.'
           endif
           stop
        endif

        do j= 1, Num_wann
           jbase= j
           jsector= 1
           if (has_spin .and. j>nspinless) then
              jbase= j-nspinless
              jsector= 2
           endif
           if (isector/=jsector) cycle

           jatom= Origin_cell%spinorbital_to_atom_index(jbase)
           if (iatom/=jatom) cycle

           jproj= Origin_cell%spinorbital_to_projector_index(jbase)
           label_j= trim(adjustl(Origin_cell%proj_name(jproj, jatom)))
           call orbital_real_basis_coeff(label_j, l_j, coeff_j, supported_j)
           if (.not.supported_j) then
              if (cpuid==0) then
                 write(stdout, '(a, a)')'ERROR in OHC_calc: unsupported projector label: ', trim(label_j)
                 write(stdout, '(a)')'Supported labels include s, pz/px/py, dxy/dyz/dxz(dxz)/dx2-y2/dz2.'
              endif
              stop
           endif
           if (l_i/=l_j) cycle

           lx= 0d0
           ly= 0d0
           lz= 0d0
           do m1= -l_i, l_i
              do m2= -l_i, l_i
                 lp= 0d0
                 lm= 0d0
                 if (m1==m2+1) lp= sqrt(dble(l_i*(l_i+1)-m2*(m2+1)))
                 if (m1==m2-1) lm= sqrt(dble(l_i*(l_i+1)-m2*(m2-1)))
                 lx= lx+ conjg(coeff_i(m1))*coeff_j(m2)*(lp+lm)/2d0
                 ly= ly+ conjg(coeff_i(m1))*coeff_j(m2)*(lp-lm)/(2d0*zi)
                 if (m1==m2) lz= lz+ conjg(coeff_i(m1))*coeff_j(m2)*dble(m2)
              enddo
           enddo

           orbital_matrices(i, j, 1)= lx
           orbital_matrices(i, j, 2)= ly
           orbital_matrices(i, j, 3)= lz
        enddo
     enddo

   return
  end subroutine build_orbital_operator_matrices

  subroutine atomic_orbital_moment_singlek_allbands(UU, orbital_operator, L_atom)
   !------------------------------------------------------------------!
   !> Calculate <L_atom> for one k point and all bands.               !
   !> orbital_operator is the atom-centered L built from PROJECTORS;  !
   !> UU contains Hamiltonian eigenvectors from eigensystem_c.        !
   !> The output is in units of hbar.                                 !
   !------------------------------------------------------------------!
     use para
     implicit none

     complex(dp), intent(in) :: UU(Num_wann, Num_wann)
     complex(dp), intent(in) :: orbital_operator(Num_wann, Num_wann, 3)
     real(dp), intent(out) :: L_atom(Num_wann, 3)

     integer :: i, n
     complex(dp), allocatable :: L_Ham(:, :)

     allocate(L_Ham(Num_wann, Num_wann))
     L_Ham= 0d0
     L_atom= 0d0

     do i= 1, 3
        call rotation_to_Ham_basis(UU, orbital_operator(:, :, i), L_Ham)
        do n= 1, Num_wann
           L_atom(n, i)= dble(L_Ham(n, n))
        enddo
     enddo

     deallocate(L_Ham)

   return
  end subroutine atomic_orbital_moment_singlek_allbands

  subroutine sigma_OHC
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the orbital Hall conductivity    !
   !> and orbital Berry curvature map (sigma_xy^z) on KPLANE_BULK.    !
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. sigma_ohc_eta***meV.txt: sigma tensor for several eta values !
   !> 2. sigma_ohc.gnu: gnuplot script to plot sigma tensor           !
   !> 3. ohc_berrycurvature_xy_z.dat: BZ map for orbital curvature    !
   !> 4. ohc_berrycurvature_xy_z.gnu: gnuplot script for BZ map       !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none

     integer :: ik, ikx, iky, ikz, ieta
     integer :: m, n, i, ie, ialpha, ibeta, igamma
     integer :: ierr, knv3, nkmesh2
     integer :: NumberofEta

     real(dp) :: mu, Beta_fake, deno_fac, eta_local, vmin, vmax, fac_nk1, fac_nk2
     real(dp) :: k(3), k_cart(3), kxy_plane(3)

     real(dp) :: time_start, time_end

     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)
     real(dp), allocatable :: energy(:)

     real(dp), allocatable :: sigma_tensor_ohc(:, :, :, :, :)
     real(dp), allocatable :: sigma_tensor_ohc_mpi(:, :, :, :, :)

     real(dp), external :: fermi, delta

     complex(dp), allocatable :: Vmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_wann(:, :, :)
     complex(dp), allocatable :: orbital_operator(:, :, :)
     complex(dp), allocatable :: j_orb_gamma_alpha(:, :)
     complex(dp), allocatable :: j_orb_xy(:, :)
     complex(dp), allocatable :: mat_t(:, :)

     real(dp), allocatable :: Omega_orb(:)
     real(dp), allocatable :: Omega_orb_t(:)

     real(dp), allocatable :: Omega_k(:), Omega_k_mpi(:)
     real(dp), allocatable :: L_atom(:, :)
     real(dp), allocatable :: L_atom_allk_Occ(:, :, :), L_atom_allk_Occ_mpi(:, :, :)
     real(dp), allocatable :: L_atom_allk_EF(:, :, :), L_atom_allk_EF_mpi(:, :, :)

     real(dp), allocatable :: eta_array(:)
     character(80) :: ohcfilename, etaname

     NumberofEta= 9

     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_wann(Num_wann, Num_wann, 3))
     allocate(orbital_operator(Num_wann, Num_wann, 3))
     allocate(j_orb_gamma_alpha(Num_wann, Num_wann))
     allocate(j_orb_xy(Num_wann, Num_wann))
     allocate(Omega_orb(Num_wann))
     allocate(Omega_orb_t(Num_wann))
     allocate(L_atom(Num_wann, 3))

     allocate(W(Num_wann))
     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate(mat_t(Num_wann, Num_wann))
     allocate(energy(OmegaNum))
     allocate(eta_array(NumberofEta))
     allocate(sigma_tensor_ohc    (OmegaNum, 3, 3, 3, NumberofEta))
     allocate(sigma_tensor_ohc_mpi(OmegaNum, 3, 3, 3, NumberofEta))

     Vmn_Ham= 0d0
     Vmn_wann= 0d0
     orbital_operator= 0d0
     j_orb_gamma_alpha= 0d0
     j_orb_xy= 0d0
     Omega_orb= 0d0
     Omega_orb_t= 0d0
     L_atom= 0d0
     W= 0d0
     Hamk_bulk= 0d0
     UU= 0d0
     mat_t= 0d0
     sigma_tensor_ohc= 0d0
     sigma_tensor_ohc_mpi= 0d0

     eta_array= (/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Eta_Arc

     call build_orbital_operator_matrices(orbital_operator)

     !> print orbital matrices (Lx, Ly, Lz) in Wannier basis to stdout
     if (cpuid==0) then
        write(stdout, '(a)') ' '
        write(stdout, '(a)') '>> Orbital operator matrices (Lx, Ly, Lz) in Wannier basis'
        write(stdout, '(a)') '   Built from PROJECTORS via build_orbital_operator_matrices'
        write(stdout, '(a)') ' '
        do i= 1, 3
           select case(i)
              case(1); write(stdout, '(a)') '   Lx matrix (real part):'
              case(2); write(stdout, '(a)') '   Ly matrix (real part):'
              case(3); write(stdout, '(a)') '   Lz matrix (real part):'
           end select
           do n= 1, Num_wann
              write(stdout, '(4X, 100F10.5)') (dble(orbital_operator(n, m, i)), m= 1, Num_wann)
           enddo
           write(stdout, '(a)') ' '
        enddo

        do i= 1, 3
           select case(i)
              case(1); write(stdout, '(a)') '   Lx matrix (imag part):'
              case(2); write(stdout, '(a)') '   Ly matrix (imag part):'
              case(3); write(stdout, '(a)') '   Lz matrix (imag part):'
           end select
           do n= 1, Num_wann
              write(stdout, '(4X, 100F10.5)') (aimag(orbital_operator(n, m, i)), m= 1, Num_wann)
           enddo
           write(stdout, '(a)') ' '
        enddo
        write(stdout, '(a)') '>> End of orbital operator matrices'
        write(stdout, '(a)') ' '
     endif

     do ie= 1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)*(ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo

     knv3= Nk1*Nk2*Nk3
     nkmesh2= Nk1*Nk2
     allocate(Omega_k(nkmesh2), Omega_k_mpi(nkmesh2))
     allocate(L_atom_allk_Occ(2, 3, nkmesh2), L_atom_allk_Occ_mpi(2, 3, nkmesh2))
     allocate(L_atom_allk_EF(2, 3, nkmesh2), L_atom_allk_EF_mpi(2, 3, nkmesh2))
     Omega_k= 0d0
     Omega_k_mpi= 0d0
     L_atom_allk_Occ= 0d0
     L_atom_allk_Occ_mpi= 0d0
     L_atom_allk_EF= 0d0
     L_atom_allk_EF_mpi= 0d0

     call now(time_start)
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and.mod(ik/num_cpu, 100).eq.0) then
           call now(time_end)
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= ik-(iky-1)*Nk3-(ikx-1)*Nk2*Nk3
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call ham_bulk_atomicgauge(k, Hamk_bulk)
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)
        call dHdk_atomicgauge(k, Vmn_wann)

        do ialpha= 1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
        enddo

        do igamma= 1, 3
           do ialpha= 1, 3
              j_orb_gamma_alpha= 0d0
              call mat_mul(Num_wann, orbital_operator(:, :, igamma), &
                 Vmn_wann(:, :, ialpha), j_orb_gamma_alpha(:, :))
              call mat_mul(Num_wann, Vmn_wann(:, :, ialpha), &
                 orbital_operator(:, :, igamma), mat_t)
              j_orb_gamma_alpha= (j_orb_gamma_alpha+mat_t)/2d0

              mat_t= j_orb_gamma_alpha
              call rotation_to_Ham_basis(UU, mat_t, j_orb_gamma_alpha)

              do ieta= 1, NumberofEta
                 eta_local= eta_array(ieta)
                 do ibeta= 1, 3
                    Omega_orb= 0d0
                    do n= 1, Num_wann
                       do m= 1, Num_wann
                          if (abs(W(n)-W(m))<eps9 ) cycle
                          deno_fac= -2d0/((W(n)-W(m))**2+eta_local**2)
                          Omega_orb(n)= Omega_orb(n)+ &
                             aimag(j_orb_gamma_alpha(n, m)*Vmn_Ham(m, n, ibeta))*deno_fac
                       enddo
                    enddo

                    Beta_fake= 1d0/eta_local
                    do ie= 1, OmegaNum
                       mu= energy(ie)
                       do n= 1, Num_wann
                          Omega_orb_t(n)= Omega_orb(n)*fermi(W(n)-mu, Beta_fake)
                       enddo
                       sigma_tensor_ohc_mpi(ie, igamma, ialpha, ibeta, ieta)= &
                          sigma_tensor_ohc_mpi(ie, igamma, ialpha, ibeta, ieta)+sum(Omega_orb_t(:))
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ohc_mpi, sigma_tensor_ohc, size(sigma_tensor_ohc), &
                       mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     sigma_tensor_ohc= sigma_tensor_ohc_mpi
#endif

     !> in unit of ((hbar/e)(Ohm*m)^-1)
     sigma_tensor_ohc= sigma_tensor_ohc/dble(knv3)/Origin_cell%CellVolume*&
        Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume

     !> in unit of ((hbar/e)(Ohm*cm)^-1)
     sigma_tensor_ohc= sigma_tensor_ohc/100d0

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta= 1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(ohcfilename, '(7a)')'sigma_ohc_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=ohcfilename)
           write(outfileindex, '("#",10a)')' Orbital hall conductivity in unit of (hbar/e)S/cm,', &
              'Brodening eta= ', trim(adjustl(etaname)), ' meV'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 28)
           write(outfileindex, '("#",a13, 27a16)')'Eenergy (eV)', &
             'xx^x', 'xy^x', 'xz^x', 'yx^x', 'yy^x', 'yz^x', 'zx^x', 'zy^x', 'zz^x', &
             'xx^y', 'xy^y', 'xz^y', 'yx^y', 'yy^y', 'yz^y', 'zx^y', 'zy^y', 'zz^y', &
             'xx^z', 'xy^z', 'xz^z', 'yx^z', 'yy^z', 'yz^z', 'zx^z', 'zy^z', 'zz^z'
           do ie= 1, OmegaNum
              write(outfileindex, '(E16.8)', advance='no')energy(ie)/eV2Hartree
              do igamma= 1, 3
              do ialpha= 1, 3
              do ibeta = 1, 3
                 if (ialpha*ibeta*igamma/=27) then
                    write(outfileindex, '(200E16.8)', advance='no') &
                       sigma_tensor_ohc(ie, igamma, ialpha, ibeta, ieta)
                 else
                    write(outfileindex, '(200E16.8)', advance='yes') &
                       sigma_tensor_ohc(ie, igamma, ialpha, ibeta, ieta)
                 endif
              enddo
              enddo
              enddo
           enddo
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
        write(ohcfilename, '(7a)')'sigma_ohc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_ohc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ohc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', &
           OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "OHC (\hbar/e)S/cm"'
        write(outfileindex, '(5a)')"plot '", trim(adjustl(ohcfilename)), &
           "' u 1:21 w l title '\sigma_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ohcfilename)), &
           "' u 1:17 w l title '\sigma_{zx}^y' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ohcfilename)), &
           "' u 1:13 w l title '\sigma_{xz}^y' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif

     !> orbital Berry curvature map on KPLANE_BULK for sigma_xy^z
     Omega_k= 0d0
     fac_nk1= 0d0
     fac_nk2= 0d0
     if (Nk1>1) fac_nk1= 1d0/dble(Nk1-1)
     if (Nk2>1) fac_nk2= 1d0/dble(Nk2-1)
     eta_local= max(Eta_Arc, eps9)
     Beta_fake= 1d0/eta_local

     do ik= 1+ cpuid, nkmesh2, num_cpu
        ikx= (ik-1)/Nk2+1
        iky= ik-(ikx-1)*Nk2
        k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
         + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0

        call ham_bulk_atomicgauge(k, Hamk_bulk)
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)
        call dHdk_atomicgauge(k, Vmn_wann)
         do ialpha= 1, 3
            call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
         enddo

         call atomic_orbital_moment_singlek_allbands(UU, orbital_operator, L_atom)
         do i= 1, 3
            L_atom_allk_Occ(1, i, ik)= sum(L_atom(1:NumOccupied, i))
            L_atom_allk_Occ(2, i, ik)= L_atom(NumOccupied, i)
            do m= 1, Num_wann
               L_atom_allk_EF(1, i, ik)= L_atom_allk_EF(1, i, ik)+ &
                  L_atom(m, i)*fermi(W(m)-E_arc, Beta_fake)
               L_atom_allk_EF(2, i, ik)= L_atom_allk_EF(2, i, ik)+ &
                  L_atom(m, i)*delta(eta_local, W(m)-E_arc)
            enddo
         enddo

        j_orb_xy= 0d0
        call mat_mul(Num_wann, orbital_operator(:, :, 3), Vmn_wann(:, :, 1), j_orb_xy(:, :))
        call mat_mul(Num_wann, Vmn_wann(:, :, 1), orbital_operator(:, :, 3), mat_t)
        j_orb_xy= (j_orb_xy+mat_t)/2d0
        call rotation_to_Ham_basis(UU, j_orb_xy, j_orb_xy)

        Omega_orb= 0d0
        do n= 1, NumOccupied
           do m= NumOccupied+1, Num_wann
              if (abs(W(n)-W(m))<eps9 ) cycle
              deno_fac= -2d0/((W(n)-W(m))**2+eta_local**2)
              Omega_orb(n)= Omega_orb(n)+ &
                 aimag(j_orb_xy(n, m)*Vmn_Ham(m, n, 2))*deno_fac
           enddo
        enddo
        Omega_k(ik)= sum(Omega_orb(1:NumOccupied))
     enddo

#if defined (MPI)
     call mpi_allreduce(Omega_k, Omega_k_mpi, size(Omega_k_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
     call mpi_allreduce(L_atom_allk_Occ, L_atom_allk_Occ_mpi, size(L_atom_allk_Occ_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
     call mpi_allreduce(L_atom_allk_EF, L_atom_allk_EF_mpi, size(L_atom_allk_EF_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     Omega_k_mpi= Omega_k
     L_atom_allk_Occ_mpi= L_atom_allk_Occ
     L_atom_allk_EF_mpi= L_atom_allk_EF
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='atomic_orbital_moment.dat')
        write(outfileindex, '(a)')'# Atom-centered <L_atom> from PROJECTORS'
        write(outfileindex, '(a)')'# Columns 7-15 are in units of hbar; columns 16-18 are delta-broadened weights in hbar/Ha.'
        write(outfileindex, '(a)')'# Positive values are angular momentum expectation values, not electron magnetic moments.'
        write(outfileindex, '("#", 30X,4a36)')'Sum over 1-NumOccupied bands', &
           "values at NumOccupied'th band", &
           ' Sum below Fermi level','values at Fermi level'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
            'Lx', 'Ly', 'Lz', 'Lx', 'Ly', 'Lz', &
            'Lx', 'Ly', 'Lz', 'Lx', 'Ly', 'Lz'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1,18)
        ik= 0
        do ikx= 1, Nk1
           do iky= 1, Nk2
              ik= ik+1
              k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
               + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0
              k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
              call rotate_k3_to_kplane(k_cart, kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, &
                 kxy_plane*Angstrom2atomic, L_atom_allk_Occ_mpi(1, :, ik), &
                 L_atom_allk_Occ_mpi(2, :, ik), L_atom_allk_EF_mpi(1, :, ik), &
                 L_atom_allk_EF_mpi(2, :, ik)
           enddo
           write(outfileindex, *) ' '
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='ohc_berrycurvature_xy_z.dat')
        write(outfileindex, '(a)')'# Orbital Berry curvature for sigma_xy^z, unit: Angstrom^2'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", 'Omega_orb'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 7)
        ik= 0
        do ikx= 1, Nk1
           do iky= 1, Nk2
              ik= ik+1
              k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
               + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0
              k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
              call rotate_k3_to_kplane(k_cart, kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, &
                 kxy_plane*Angstrom2atomic, Omega_k_mpi(ik)/Angstrom2atomic**2
            !   write(outfileindex, *) ' '
           enddo
           write(outfileindex, *) ' '
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        vmin= minval(Omega_k_mpi)/Angstrom2atomic**2
        vmax= maxval(Omega_k_mpi)/Angstrom2atomic**2
        open(unit=outfileindex, file='ohc_berrycurvature_xy_z.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 2560, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'ohc_berrycurvature_xy_z.png'"
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a, E10.3, a, E10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title 'Orbital {/Symbol W}_{xy}^z ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'ohc_berrycurvature_xy_z.dat' u 4:5:7 w pm3d"
        close(outfileindex)
     endif

     deallocate(W, Hamk_bulk, UU, energy)
     deallocate(sigma_tensor_ohc, sigma_tensor_ohc_mpi)
     deallocate(Vmn_Ham, Vmn_wann, orbital_operator, j_orb_gamma_alpha, j_orb_xy, mat_t)
     deallocate(Omega_orb, Omega_orb_t)
     deallocate(L_atom, L_atom_allk_Occ, L_atom_allk_Occ_mpi, L_atom_allk_EF, L_atom_allk_EF_mpi)
     deallocate(eta_array)
     deallocate(Omega_k, Omega_k_mpi)

   return
  end subroutine sigma_OHC

  subroutine build_L_total_Ham(W, Vmn_Ham, eta, L_Ham_tot)
   !-----------------------------------------------------------------------!
   !> Build the "total" orbital angular momentum operator L_gamma in the   !
   !> Hamiltonian eigenbasis from the modern theory of orbital            !
   !> magnetization. Implements Eq. (5) of                                !
   !> A. Pezo, D. Garcia Ovalle, A. Manchon,                              !
   !> Phys. Rev. B 106, 104414 (2022):                                    !
   !>                                                                    !
   !>   L_gamma^{np}(k) = (1/2) eps_{gamma,alpha,beta} *                  !
   !>      Im sum_{q != n,p} [1/(E_q-E_n)+1/(E_q-E_p)] v_alpha^{nq} v_beta^{qp}!
   !>                                                                    !
   !> Atomic Hartree units (hbar=e=m_e=1, mu_B=1/2) so the formula's     !
   !> prefactor e^2/(4 mu_B) reduces to 1/2 and L comes out in units of  !
   !> hbar, consistent with the intra-atomic L matrices used by sigma_OHC.!
   !-----------------------------------------------------------------------!
     use para, only: dp, Num_wann, eps9, zzero, Eta_Arc
     implicit none

     real(dp),    intent(in)  :: W(Num_wann)
     complex(dp), intent(in)  :: Vmn_Ham(Num_wann, Num_wann, 3)
     real(dp), intent(in)  :: eta
     complex(dp), intent(out) :: L_Ham_tot(Num_wann, Num_wann, 3)

     integer :: n, p, q, ig, ia, ib
     integer :: eps_levi(3, 3, 3)
     real(dp) :: fac, dE
     complex(dp) :: fac2, accum
     complex(dp) :: zi = (0.0_dp, 1.0_dp)
     complex(dp) :: resolvent(Num_wann, Num_wann)

     L_Ham_tot = zzero

     eps_levi = 0
     eps_levi(1, 2, 3) =  1
     eps_levi(2, 3, 1) =  1
     eps_levi(3, 1, 2) =  1
     eps_levi(1, 3, 2) = -1
     eps_levi(3, 2, 1) = -1
     eps_levi(2, 1, 3) = -1

     dE = 1e-9

     !> R(n,q) is reused for every p and Cartesian component.  In the
     !> original implementation the same complex division was repeated
     !> inside the O(Num_wann**3) loop.
     resolvent= zzero
     do n= 1, Num_wann
        do q= 1, Num_wann
           if (q == n) cycle
           if (abs(W(q)-W(n)) < dE) cycle
           resolvent(n, q)= 1d0/(W(q)-W(n)-zi*eta)
        enddo
     enddo

     do ig = 1, 3
        do ia = 1, 3
           do ib = 1, 3
              if (eps_levi(ig, ia, ib) == 0) cycle
              fac = 0.5d0 * dble(eps_levi(ig, ia, ib))
               do n = 1, Num_wann
                  do p = 1, Num_wann
                     accum = zzero
                     do q = 1, Num_wann
                           fac2= resolvent(n, q)+resolvent(p, q)
                           accum= accum+fac2*Vmn_Ham(n, q, ia)*Vmn_Ham(q, p, ib)
                     enddo ! q
                     L_Ham_tot(n, p, ig) = L_Ham_tot(n, p, ig) + &
                                                !  fac * aimag(accum)
                                               (-zi) * fac  * accum
                                             !   cmplx(fac * aimag(accum), 0d0, dp)
                 enddo ! p
              enddo ! n
           enddo
        enddo
     enddo

   return
  end subroutine build_L_total_Ham

  subroutine sigma_OHC_tot
   !------------------------------------------------------------------!
   !> Total orbital Hall conductivity using the modern theory of      !
   !> orbital magnetization. Differs from sigma_OHC by replacing the  !
   !> atom-centered Lz with the Bloch-basis orbital moment operator   !
   !> [Eq. (5), PRB 106, 104414 (2022)], capturing both intra- and    !
   !> inter-atomic contributions.                                     !
   !>                                                                 !
   !> Output files:                                                   !
   !> 1. sigma_ohc_tot_eta***meV.txt: sigma tensor for several eta    !
   !> 2. sigma_ohc_tot.gnu: gnuplot script                            !
   !> 3. ohc_tot_berrycurvature_xy_z.dat / .gnu: BZ map (sigma_xy^z)  !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none

     integer :: ik, ikx, iky, ikz, ieta
     integer :: m, n, i, ie, ialpha, ibeta, igamma
     integer :: ierr, knv3, nkmesh2
     integer :: NumberofEta

     !> Diagnostic Pezo implementation variant: discard diagonal modern-OAM
     !> matrix elements L_nn(k) in the Hamiltonian eigenbasis before forming
     !> J={L,v}/2.  Set false to recover the complete modern-L calculation.
     logical, parameter :: use_offdiagonal_L_only = .false.
     logical, parameter :: run_hermiticity_diagnostics = .false.

     real(dp) :: mu, Beta_fake, deno_fac, eta_local, vmin, vmax, fac_nk1, fac_nk2
     real(dp) :: k(3), k_cart(3), kxy_plane(3)

     real(dp) :: time_start, time_end

     !> [TEMP DIAGNOSTIC] max non-Hermitian deviation over all k-points
     real(dp) :: herm_dev_L1, herm_dev_L2, herm_dev_j, herm_dev_v
     real(dp) :: herm_dev_L1_mpi, herm_dev_L2_mpi, herm_dev_j_mpi, herm_dev_v_mpi

     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)
     real(dp), allocatable :: energy(:)

     real(dp), allocatable :: sigma_tensor_ohc(:, :, :, :, :)
     real(dp), allocatable :: sigma_tensor_ohc_mpi(:, :, :, :, :)

     real(dp), external :: fermi, delta

     complex(dp), allocatable :: Vmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_wann(:, :, :)
     complex(dp), allocatable :: L_Ham_tot(:, :, :)
     complex(dp), allocatable :: j_orb_gamma_alpha(:, :)
     complex(dp), allocatable :: j_orb_xy(:, :)
     complex(dp), allocatable :: mat_t(:, :)

     real(dp), allocatable :: Omega_orb(:)
     real(dp), allocatable :: deno_matrix(:, :)
     real(dp), allocatable :: fermi_weight(:, :)

     real(dp), allocatable :: Omega_k(:), Omega_k_mpi(:)
     real(dp), allocatable :: L_tot(:, :)
     real(dp), allocatable :: L_tot_allk_Occ(:, :, :), L_tot_allk_Occ_mpi(:, :, :)
     real(dp), allocatable :: L_tot_allk_EF(:, :, :), L_tot_allk_EF_mpi(:, :, :)

     real(dp), allocatable :: eta_array(:)
     character(80) :: ohcfilename, etaname

     NumberofEta= 9

     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_wann(Num_wann, Num_wann, 3))
     allocate(L_Ham_tot(Num_wann, Num_wann, 3))
     allocate(j_orb_gamma_alpha(Num_wann, Num_wann))
     allocate(j_orb_xy(Num_wann, Num_wann))
     allocate(Omega_orb(Num_wann))
     allocate(deno_matrix(Num_wann, Num_wann))
     allocate(fermi_weight(Num_wann, OmegaNum))

     allocate(W(Num_wann))
     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate(mat_t(Num_wann, Num_wann))
     allocate(energy(OmegaNum))
     allocate(eta_array(NumberofEta))
     allocate(sigma_tensor_ohc    (OmegaNum, 3, 3, 3, NumberofEta))
     allocate(sigma_tensor_ohc_mpi(OmegaNum, 3, 3, 3, NumberofEta))

     Vmn_Ham= 0d0
     Vmn_wann= 0d0
     L_Ham_tot= 0d0
     j_orb_gamma_alpha= 0d0
     j_orb_xy= 0d0
     Omega_orb= 0d0
     deno_matrix= 0d0
     fermi_weight= 0d0
     W= 0d0
     Hamk_bulk= 0d0
     UU= 0d0
     mat_t= 0d0
     sigma_tensor_ohc= 0d0
     sigma_tensor_ohc_mpi= 0d0

     eta_array= (/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Eta_Arc

     do ie= 1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)*(ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo

     knv3= Nk1*Nk2*Nk3
     nkmesh2= Nk1*Nk2
     allocate(Omega_k(nkmesh2), Omega_k_mpi(nkmesh2))
     allocate(L_tot(Num_wann, 3))
     allocate(L_tot_allk_Occ(2, 3, nkmesh2), L_tot_allk_Occ_mpi(2, 3, nkmesh2))
     allocate(L_tot_allk_EF(2, 3, nkmesh2), L_tot_allk_EF_mpi(2, 3, nkmesh2))
     Omega_k= 0d0
     Omega_k_mpi= 0d0
     L_tot= 0d0
     L_tot_allk_Occ= 0d0
     L_tot_allk_Occ_mpi= 0d0
     L_tot_allk_EF= 0d0
     L_tot_allk_EF_mpi= 0d0

     !> [TEMP DIAGNOSTIC] accumulate max non-Hermitian deviation
     herm_dev_L1= 0d0
     herm_dev_L2= 0d0
     herm_dev_j= 0d0
     herm_dev_v= 0d0

     call now(time_start)
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and.mod(ik/num_cpu, 100).eq.0) then
           call now(time_end)
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= ik-(iky-1)*Nk3-(ikx-1)*Nk2*Nk3
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call ham_bulk_atomicgauge(k, Hamk_bulk)
      !   call ham_bulk_latticegauge(k, Hamk_bulk)          ! H: 纯晶格规范
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)
        call dHdk_atomicgauge(k, Vmn_wann)
        do ialpha= 1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
        enddo

      !   call ham_bulk_latticegauge(k, Hamk_bulk)          ! H: 纯晶格规范
      !   UU= Hamk_bulk
      !   call eigensystem_c('V', 'U', Num_wann, UU, W)
      !   call dHdk_latticegauge_wann(k, Vmn_wann)           ! v: 纯晶格规范
      !   do ialpha= 1, 3
      !      call rotation_to_Ham_basis(UU, Vmn_wann(:,:,ialpha), Vmn_Ham(:,:,ialpha))
      !   enddo
      ! !   call dHdk_latticegauge_Ham(k, W, UU,  Vmn_Ham)           ! v: 纯晶格规范

        if (run_hermiticity_diagnostics) then
           do i= 1, 3
              herm_dev_v= max(herm_dev_v, &
                 maxval(abs(Vmn_Ham(:, :, i)- conjg(transpose(Vmn_Ham(:, :, i))))))
           enddo
        endif

      !   call ham_bulk_latticegauge(k, Hamk_bulk)          ! H: 纯晶格规范
      !   UU= Hamk_bulk
      !   call eigensystem_c('V', 'U', Num_wann, UU, W)
      !   call dHdk_latticegauge_wann(k, Vmn_wann)           ! v: 纯晶格规范
      !   do ialpha= 1, 3
      !      call rotation_to_Ham_basis(UU, Vmn_wann(:,:,ialpha), Vmn_Ham(:,:,ialpha))
      !   enddo
      ! !   call dHdk_latticegauge_Ham(k, W, UU,  Vmn_Ham)           ! v: 纯晶格规范

        !> L_total, the broadened denominator, and the Fermi weights only
        !> depend on eta at fixed k.  Build each of them once per eta.
        do ieta= 1, NumberofEta
           eta_local= eta_array(ieta)
           Beta_fake= 1d0/eta_local

           call build_L_total_Ham(W, Vmn_Ham, eta_local, L_Ham_tot)
           do i= 1, 3
              if (run_hermiticity_diagnostics) then
                 herm_dev_L1= max(herm_dev_L1, &
                    maxval(abs(L_Ham_tot(:, :, i)- conjg(transpose(L_Ham_tot(:, :, i))))))
              endif
              L_Ham_tot(:, :, i)= 0.5d0 * &
                 (L_Ham_tot(:, :, i)+conjg(transpose(L_Ham_tot(:, :, i))))
              if (run_hermiticity_diagnostics) then
                 herm_dev_L2= max(herm_dev_L2, &
                    maxval(abs(L_Ham_tot(:, :, i)- conjg(transpose(L_Ham_tot(:, :, i))))))
              endif
           enddo
           if (use_offdiagonal_L_only) then
              do i= 1, 3
                 do n= 1, Num_wann
                    L_Ham_tot(n, n, i)= zzero
                 enddo
              enddo
           endif

           deno_matrix= 0d0
           do n= 1, Num_wann
              do m= 1, Num_wann
                 if (abs(W(n)-W(m))<eps9) cycle
                 deno_matrix(n, m)= -2d0/((W(n)-W(m))**2+eta_local**2)
              enddo
           enddo

           do ie= 1, OmegaNum
              mu= energy(ie)
              do n= 1, Num_wann
                 fermi_weight(n, ie)= fermi(W(n)-mu, Beta_fake)
              enddo
           enddo

           do igamma= 1, 3
              do ialpha= 1, 3

                  !> J^o_{gamma,alpha} = {v_alpha, L_gamma}/2, all in Ham basis
                  j_orb_gamma_alpha= 0d0
                  call mat_mul(Num_wann, L_Ham_tot(:, :, igamma), &
                     Vmn_Ham(:, :, ialpha), j_orb_gamma_alpha(:, :))
                  call mat_mul(Num_wann, Vmn_Ham(:, :, ialpha), &
                     L_Ham_tot(:, :, igamma), mat_t)
                  j_orb_gamma_alpha= (j_orb_gamma_alpha+mat_t)/2d0

                  if (run_hermiticity_diagnostics) then
                     herm_dev_j= max(herm_dev_j, maxval(abs( &
                        j_orb_gamma_alpha-conjg(transpose(j_orb_gamma_alpha)))))
                  endif

                  do ibeta= 1, 3
                     Omega_orb= 0d0
                     do n= 1, Num_wann
                        do m= 1, Num_wann
                           Omega_orb(n)= Omega_orb(n)+ &
                               aimag(j_orb_gamma_alpha(n, m)*Vmn_Ham(m, n, ibeta))* &
                               deno_matrix(n, m)
                        enddo !m
                     enddo !n

                     do ie= 1, OmegaNum
                        sigma_tensor_ohc_mpi(ie, igamma, ialpha, ibeta, ieta)= &
                           sigma_tensor_ohc_mpi(ie, igamma, ialpha, ibeta, ieta)+ &
                           dot_product(Omega_orb, fermi_weight(:, ie))
                     enddo!ie
                  enddo!ibeta
              enddo!ialpha
           enddo !igamma
        enddo !ieta
     enddo !ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ohc_mpi, sigma_tensor_ohc, size(sigma_tensor_ohc), &
                       mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     sigma_tensor_ohc= sigma_tensor_ohc_mpi
#endif

     !> Optional diagnostic: reduce max non-Hermitian deviation across cpus.
     if (run_hermiticity_diagnostics) then
#if defined (MPI)
        call mpi_allreduce(herm_dev_v, herm_dev_v_mpi, 1, mpi_dp, mpi_max, mpi_cmw, ierr)
        call mpi_allreduce(herm_dev_L1, herm_dev_L1_mpi, 1, mpi_dp, mpi_max, mpi_cmw, ierr)
        call mpi_allreduce(herm_dev_L2, herm_dev_L2_mpi, 1, mpi_dp, mpi_max, mpi_cmw, ierr)
        call mpi_allreduce(herm_dev_j, herm_dev_j_mpi, 1, mpi_dp, mpi_max, mpi_cmw, ierr)
#else
        herm_dev_v_mpi= herm_dev_v
        herm_dev_L1_mpi= herm_dev_L1
        herm_dev_L2_mpi= herm_dev_L2
        herm_dev_j_mpi= herm_dev_j
#endif
        if (cpuid==0) then
           write(stdout, '(a)') ' '
           write(stdout, '(a)') ' >>> [OHC_tot diagnostic] max non-Hermitian deviation over all k-points:'
           write(stdout, '(a, E16.8)') '       Vmn_Ham          |A - A^H|_max = ', herm_dev_v_mpi
           write(stdout, '(a, E16.8)') '       L_Ham_tot_before |A - A^H|_max = ', herm_dev_L1_mpi
           write(stdout, '(a, E16.8)') '       L_Ham_tot_after  |A - A^H|_max = ', herm_dev_L2_mpi
           write(stdout, '(a, E16.8)') '       j_orb {L,v}/2    |A - A^H|_max = ', herm_dev_j_mpi
           write(stdout, '(a)') ' '
        endif
     endif

     !> in unit of ((hbar/e)(Ohm*m)^-1)
     sigma_tensor_ohc= sigma_tensor_ohc/dble(knv3)/Origin_cell%CellVolume*&
        Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume

     !> in unit of ((hbar/e)(Ohm*cm)^-1)
     sigma_tensor_ohc= sigma_tensor_ohc/100d0

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta= 1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(ohcfilename, '(7a)')'sigma_ohc_tot_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=ohcfilename)
           write(outfileindex, '("#",10a)')' Total orbital hall conductivity in unit of (hbar/e)S/cm,', &
              ' Brodening eta= ', trim(adjustl(etaname)), ' meV'
           write(outfileindex, '("#",a)')' Modern theory of orbital magnetization, PRB 106, 104414 (2022) Eq. (5)'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 28)
           write(outfileindex, '("#",a13, 27a16)')'Eenergy (eV)', &
             'xx^x', 'xy^x', 'xz^x', 'yx^x', 'yy^x', 'yz^x', 'zx^x', 'zy^x', 'zz^x', &
             'xx^y', 'xy^y', 'xz^y', 'yx^y', 'yy^y', 'yz^y', 'zx^y', 'zy^y', 'zz^y', &
             'xx^z', 'xy^z', 'xz^z', 'yx^z', 'yy^z', 'yz^z', 'zx^z', 'zy^z', 'zz^z'
           do ie= 1, OmegaNum
              write(outfileindex, '(E16.8)', advance='no')energy(ie)/eV2Hartree
              do igamma= 1, 3
              do ialpha= 1, 3
              do ibeta = 1, 3
                 if (ialpha*ibeta*igamma/=27) then
                    write(outfileindex, '(200E16.8)', advance='no') &
                       sigma_tensor_ohc(ie, igamma, ialpha, ibeta, ieta)
                 else
                    write(outfileindex, '(200E16.8)', advance='yes') &
                       sigma_tensor_ohc(ie, igamma, ialpha, ibeta, ieta)
                 endif
              enddo
              enddo
              enddo
           enddo
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
        write(ohcfilename, '(7a)')'sigma_ohc_tot_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_ohc_tot.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ohc_tot.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', &
           OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "Total OHC (\hbar/e)S/cm"'
        write(outfileindex, '(5a)')"plot '", trim(adjustl(ohcfilename)), &
           "' u 1:21 w l title '\sigma_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ohcfilename)), &
           "' u 1:17 w l title '\sigma_{zx}^y' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ohcfilename)), &
           "' u 1:13 w l title '\sigma_{xz}^y' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif

     !> total orbital Berry curvature map on KPLANE_BULK for sigma_xy^z
     Omega_k= 0d0
     fac_nk1= 0d0
     fac_nk2= 0d0
     if (Nk1>1) fac_nk1= 1d0/dble(Nk1-1)
     if (Nk2>1) fac_nk2= 1d0/dble(Nk2-1)
     eta_local= max(Eta_Arc, eps9)
     Beta_fake= 1d0/eta_local

     do ik= 1+ cpuid, nkmesh2, num_cpu
        ikx= (ik-1)/Nk2+1
        iky= ik-(ikx-1)*Nk2
        k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
         + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0

        call ham_bulk_atomicgauge(k, Hamk_bulk)
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)
        call dHdk_atomicgauge(k, Vmn_wann)
        do ialpha= 1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
        enddo

        call build_L_total_Ham(W, Vmn_Ham, eta_local, L_Ham_tot)
         do i= 1, 3
            L_Ham_tot(:, :, i)= 0.5d0 * (L_Ham_tot(:, :, i) + conjg(transpose(L_Ham_tot(:, :, i))))
            do n= 1, Num_wann
               L_tot(n, i)= real(L_Ham_tot(n, n, i), dp)
            enddo
            L_tot_allk_Occ(1, i, ik)= sum(L_tot(1:NumOccupied, i))
            L_tot_allk_Occ(2, i, ik)= L_tot(NumOccupied, i)
            do m= 1, Num_wann
               L_tot_allk_EF(1, i, ik)= L_tot_allk_EF(1, i, ik)+ &
                  L_tot(m, i)*fermi(W(m)-E_arc, Beta_fake)
               L_tot_allk_EF(2, i, ik)= L_tot_allk_EF(2, i, ik)+ &
                  L_tot(m, i)*delta(eta_local, W(m)-E_arc)
            enddo
         enddo
        if (use_offdiagonal_L_only) then
           do i= 1, 3
              do n= 1, Num_wann
                 L_Ham_tot(n, n, i)= zzero
              enddo
           enddo
        endif

        j_orb_xy= 0d0
        call mat_mul(Num_wann, L_Ham_tot(:, :, 3), Vmn_Ham(:, :, 1), j_orb_xy(:, :))
        call mat_mul(Num_wann, Vmn_Ham(:, :, 1), L_Ham_tot(:, :, 3), mat_t)
        j_orb_xy= (j_orb_xy+mat_t)/2d0

        Omega_orb= 0d0
        do n= 1, NumOccupied
           do m= NumOccupied+1, Num_wann
              if (abs(W(n)-W(m))<eps9 ) cycle
              deno_fac= -2d0/((W(n)-W(m))**2+eta_local**2)
              Omega_orb(n)= Omega_orb(n)+ &
                 aimag(j_orb_xy(n, m)*Vmn_Ham(m, n, 2))*deno_fac
           enddo
        enddo
        Omega_k(ik)= sum(Omega_orb(1:NumOccupied))
     enddo

#if defined (MPI)
     call mpi_allreduce(Omega_k, Omega_k_mpi, size(Omega_k_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
     call mpi_allreduce(L_tot_allk_Occ, L_tot_allk_Occ_mpi, size(L_tot_allk_Occ_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
     call mpi_allreduce(L_tot_allk_EF, L_tot_allk_EF_mpi, size(L_tot_allk_EF_mpi), &
                        mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     Omega_k_mpi= Omega_k
     L_tot_allk_Occ_mpi= L_tot_allk_Occ
     L_tot_allk_EF_mpi= L_tot_allk_EF
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='total_orbital_moment.dat')
        write(outfileindex, '(a)')'# Modern-theory <L_total> from the diagonal of the Pezo orbital-moment operator'
        write(outfileindex, '(a)')'# Columns 7-15 are in units of hbar; columns 16-18 are delta-broadened weights in hbar/Ha.'
        write(outfileindex, '(a)')'# Positive values are angular momentum expectation values, not electron magnetic moments.'
        write(outfileindex, '("#", 30X,4a36)')'Sum over 1-NumOccupied bands', &
           "values at NumOccupied'th band", &
           ' Sum below Fermi level','values at Fermi level'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
           'Lx', 'Ly', 'Lz', 'Lx', 'Ly', 'Lz', &
           'Lx', 'Ly', 'Lz', 'Lx', 'Ly', 'Lz'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1,18)
        ik= 0
        do ikx= 1, Nk1
           do iky= 1, Nk2
              ik= ik+1
              k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
               + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0
              k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
              call rotate_k3_to_kplane(k_cart, kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, &
                 kxy_plane*Angstrom2atomic, L_tot_allk_Occ_mpi(1, :, ik), &
                 L_tot_allk_Occ_mpi(2, :, ik), L_tot_allk_EF_mpi(1, :, ik), &
                 L_tot_allk_EF_mpi(2, :, ik)
           enddo
           write(outfileindex, *) ' '
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='ohc_tot_berrycurvature_xy_z.dat')
        write(outfileindex, '(a)')'# Total orbital Berry curvature for sigma_xy^z, unit: Angstrom^2'
        write(outfileindex, '(a)')'# Modern theory of orbital magnetization, PRB 106, 104414 (2022)'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", 'Omega_orb'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 7)
        ik= 0
        do ikx= 1, Nk1
           do iky= 1, Nk2
              ik= ik+1
              k= K3D_start+ K3D_vec1*(ikx-1)*fac_nk1  &
               + K3D_vec2*(iky-1)*fac_nk2-(K3D_vec1+K3D_vec2)/2d0
              k_cart= k(1)*Origin_cell%Kua+ k(2)*Origin_cell%Kub+ k(3)*Origin_cell%Kuc
              call rotate_k3_to_kplane(k_cart, kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')k_cart*Angstrom2atomic, &
                 kxy_plane*Angstrom2atomic, Omega_k_mpi(ik)/Angstrom2atomic**2
           enddo
           write(outfileindex, *) ' '
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        vmin= minval(Omega_k_mpi)/Angstrom2atomic**2
        vmax= maxval(Omega_k_mpi)/Angstrom2atomic**2
        open(unit=outfileindex, file='ohc_tot_berrycurvature_xy_z.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 2560, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'ohc_tot_berrycurvature_xy_z.png'"
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a, E10.3, a, E10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title 'Total Orbital {/Symbol W}_{xy}^z ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'ohc_tot_berrycurvature_xy_z.dat' u 4:5:7 w pm3d"
        close(outfileindex)
     endif

     deallocate(W, Hamk_bulk, UU, energy)
     deallocate(sigma_tensor_ohc, sigma_tensor_ohc_mpi)
     deallocate(Vmn_Ham, Vmn_wann, L_Ham_tot, j_orb_gamma_alpha, j_orb_xy, mat_t)
     deallocate(Omega_orb, deno_matrix, fermi_weight)
     deallocate(eta_array)
     deallocate(Omega_k, Omega_k_mpi)
     deallocate(L_tot, L_tot_allk_Occ, L_tot_allk_Occ_mpi, L_tot_allk_EF, L_tot_allk_EF_mpi)

   return
  end subroutine sigma_OHC_tot


  subroutine orbital_moment_kpath
   !------------------------------------------------------------------!
   !> Calculate atom-centered and total modern-theory orbital moments !
   !> along KPATH_BULK.  The orbital Berry curvature is Omega_xy^Lz.  !
   !>                                                              !
   !> Output:                                                       !
   !>   orbital_moment_kpath.dat                                    !
   !>   orbital_moment_kpath.gnu                                    !
   !>   orbital_moment_kpath.pdf (after running gnuplot)             !
   !------------------------------------------------------------------!
     use wmpi
     use para
     implicit none

     integer :: ik, i, m, n, ierr
     real(dp) :: k(3), eta_local, deno_fac
     real(dp) :: omega_intra, omega_tot

     real(dp), allocatable :: W(:)
     real(dp), allocatable :: L_atom(:, :), L_tot(:, :)
     real(dp), allocatable :: values(:, :), values_mpi(:, :)

     complex(dp), allocatable :: Hamk_bulk(:, :), UU(:, :)
     complex(dp), allocatable :: Vmn_wann(:, :, :), Vmn_Ham(:, :, :)
     complex(dp), allocatable :: orbital_operator(:, :, :), L_Ham_tot(:, :, :)
     complex(dp), allocatable :: j_intra_wann(:, :), j_intra_ham(:, :)
     complex(dp), allocatable :: j_tot_ham(:, :), mat_t(:, :)

     allocate(W(Num_wann))
     allocate(L_atom(Num_wann, 3), L_tot(Num_wann, 3))
     allocate(values(14, nk3_band), values_mpi(14, nk3_band))
     allocate(Hamk_bulk(Num_wann, Num_wann), UU(Num_wann, Num_wann))
     allocate(Vmn_wann(Num_wann, Num_wann, 3))
     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(orbital_operator(Num_wann, Num_wann, 3))
     allocate(L_Ham_tot(Num_wann, Num_wann, 3))
     allocate(j_intra_wann(Num_wann, Num_wann))
     allocate(j_intra_ham(Num_wann, Num_wann))
     allocate(j_tot_ham(Num_wann, Num_wann))
     allocate(mat_t(Num_wann, Num_wann))

     W= 0d0
     L_atom= 0d0
     L_tot= 0d0
     values= 0d0
     values_mpi= 0d0
     Hamk_bulk= 0d0
     UU= 0d0
     Vmn_wann= 0d0
     Vmn_Ham= 0d0
     orbital_operator= 0d0
     L_Ham_tot= 0d0
     j_intra_wann= 0d0
     j_intra_ham= 0d0
     j_tot_ham= 0d0
     mat_t= 0d0

     eta_local= max(Eta_Arc, eps9)
     call build_orbital_operator_matrices(orbital_operator)

     do ik= 1+cpuid, nk3_band, num_cpu
        k= k3points(:, ik)

        call ham_bulk_atomicgauge(k, Hamk_bulk)
        UU= Hamk_bulk
        call eigensystem_c('V', 'U', Num_wann, UU, W)
        call dHdk_atomicgauge(k, Vmn_wann)
        do i= 1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
        enddo

        !> Atom-centered and total modern-theory orbital moments.
        call atomic_orbital_moment_singlek_allbands(UU, orbital_operator, L_atom)
        call build_L_total_Ham(W, Vmn_Ham, eta_local, L_Ham_tot)
        do i= 1, 3
           L_Ham_tot(:, :, i)= 0.5d0* &
              (L_Ham_tot(:, :, i)+conjg(transpose(L_Ham_tot(:, :, i))))
           do n= 1, Num_wann
              L_tot(n, i)= real(L_Ham_tot(n, n, i), dp)
           enddo
        enddo

        values(1:3, ik)= sum(L_atom(1:NumOccupied, :), dim=1)
        values(4:6, ik)= sum(L_tot(1:NumOccupied, :), dim=1)
        values(7:9, ik)= L_atom(NumOccupied, :)
        values(10:12, ik)= L_tot(NumOccupied, :)

        !> Intra-atomic J_x^Lz in the Wannier basis, then rotate to
        !> the Hamiltonian eigenbasis.
        call mat_mul(Num_wann, orbital_operator(:, :, 3), &
           Vmn_wann(:, :, 1), j_intra_wann)
        call mat_mul(Num_wann, Vmn_wann(:, :, 1), &
           orbital_operator(:, :, 3), mat_t)
        j_intra_wann= 0.5d0*(j_intra_wann+mat_t)
        call rotation_to_Ham_basis(UU, j_intra_wann, j_intra_ham)

        !> Total modern-theory J_x^Lz directly in the Hamiltonian basis.
        call mat_mul(Num_wann, L_Ham_tot(:, :, 3), &
           Vmn_Ham(:, :, 1), j_tot_ham)
        call mat_mul(Num_wann, Vmn_Ham(:, :, 1), &
           L_Ham_tot(:, :, 3), mat_t)
        j_tot_ham= 0.5d0*(j_tot_ham+mat_t)

        omega_intra= 0d0
        omega_tot= 0d0
        do n= 1, NumOccupied
           do m= NumOccupied+1, Num_wann
              if (abs(W(n)-W(m))<eps9) cycle
              deno_fac= -2d0/((W(n)-W(m))**2+eta_local**2)
              omega_intra= omega_intra+aimag( &
                 j_intra_ham(n, m)*Vmn_Ham(m, n, 2))*deno_fac
              omega_tot= omega_tot+aimag( &
                 j_tot_ham(n, m)*Vmn_Ham(m, n, 2))*deno_fac
           enddo
        enddo
        values(13, ik)= omega_intra
        values(14, ik)= omega_tot
     enddo

#if defined (MPI)
     call mpi_allreduce(values, values_mpi, size(values), &
        mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     values_mpi= values
#endif

     outfileindex= outfileindex+1
     if (cpuid==0) then
        open(unit=outfileindex, file='orbital_moment_kpath.dat')
        write(outfileindex, '(a)') &
           '# Orbital moments and orbital Berry curvature along KPATH_BULK'
        write(outfileindex, '(a)') &
           '# L columns are in hbar; Omega_xy^Lz columns are in Angstrom^2.'
        write(outfileindex, '(a)') &
           '# sum: bands 1:NumOccupied; nocc: band NumOccupied.'
        write(outfileindex, '(a)') &
           '# 1:klen(1/A), 2-4:k_frac, 5-7:L_intra_sum(x,y,z),'// &
           ' 8-10:L_tot_sum(x,y,z), 11-13:L_intra_nocc(x,y,z),'// &
           ' 14-16:L_tot_nocc(x,y,z), 17:Omega_intra_xy^Lz,'// &
           ' 18:Omega_tot_xy^Lz'
        write(outfileindex, '("#column", i5, 30i16)')(i, i=1,18)
        do ik= 1, nk3_band
           write(outfileindex, '(f16.8,3f14.8,14E18.8)') &
              k3len(ik)*Angstrom2atomic, k3points(:, ik), &
              values_mpi(1:12, ik), &
              values_mpi(13:14, ik)/Angstrom2atomic**2
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+1
     if (cpuid==0) then
        open(unit=outfileindex, file='orbital_moment_kpath.gnu')
        write(outfileindex, '(a)') &
           'set terminal pdfcairo enhanced color font ",18" size 8,8'
        write(outfileindex, '(a)') &
           "set output 'orbital_moment_kpath.pdf'"
        write(outfileindex, '(a,f14.7,a)') &
           'set xrange [0:', maxval(k3len)*Angstrom2atomic, ']'
        write(outfileindex, '(a)', advance='no')'set xtics ('
        do i= 1, nk3lines+1
           if (i>1) write(outfileindex, '(a)', advance='no')', '
           write(outfileindex, '(a,a,a,f12.6)', advance='no') &
              '"', trim(k3line_name(i)), '" ', &
              k3line_stop(i)*Angstrom2atomic
        enddo
        write(outfileindex, '(a)')')'
        do i= 2, nk3lines
           write(outfileindex, '(a,f12.6,a,f12.6,a)') &
              'set arrow from ', k3line_stop(i)*Angstrom2atomic, &
              ', graph 0 to ', k3line_stop(i)*Angstrom2atomic, &
              ', graph 1 nohead dt 2 lc rgb "gray"'
        enddo
        write(outfileindex, '(a)')'set grid ytics'
        write(outfileindex, '(a)')'set key outside right center'
        write(outfileindex, '(a)')'set multiplot layout 2,1'
        write(outfileindex, '(a)')'set ylabel "<L_z> (hbar)"'
        write(outfileindex, '(a)') &
           "plot 'orbital_moment_kpath.dat' u 1:7 w l lw 2 t 'L^{intra}_{sum}', \"
        write(outfileindex, '(a)') &
           "     '' u 1:10 w l lw 2 t 'L^{tot}_{sum}', \"
        write(outfileindex, '(a)') &
           "     '' u 1:13 w l lw 2 dt 2 t 'L^{intra}_{nocc}', \"
        write(outfileindex, '(a)') &
           "     '' u 1:16 w l lw 2 dt 2 t 'L^{tot}_{nocc}'"
        write(outfileindex, '(a)')'set xlabel "K path"'
        write(outfileindex, '(a)') &
           'set ylabel "{/Symbol W}_{xy}^{L_z} ({\305}^2)"'
        write(outfileindex, '(a)') &
           "plot 'orbital_moment_kpath.dat' u 1:17 w l lw 2 t 'intra', \"
        write(outfileindex, '(a)') &
           "     '' u 1:18 w l lw 2 t 'total'"
        write(outfileindex, '(a)')'unset multiplot'
        close(outfileindex)
     endif

     deallocate(W, L_atom, L_tot, values, values_mpi)
     deallocate(Hamk_bulk, UU, Vmn_wann, Vmn_Ham)
     deallocate(orbital_operator, L_Ham_tot)
     deallocate(j_intra_wann, j_intra_ham, j_tot_ham, mat_t)

   return
  end subroutine orbital_moment_kpath
