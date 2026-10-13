! ====================================================
! ================= Method 2 =========================
! ====================================================
  
  subroutine ek_slab
   !> This subroutine is used for calculating energy 
   !> dispersion with wannier functions for 2D slab system
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  
   use wmpi
   use para
   implicit none 

   ! loop index
   integer :: i, j, l, lwork, ierr, io
   integer :: ig, ia, iq

   real(Dp) :: k(2), emin, emax, maxweight

   ! time measurement
   real(dp) :: time_start, time_end, time_start0

   ! parameters for zheev
   real(Dp), allocatable ::  rwork(:)
   complex(Dp), allocatable :: work(:)

   ! eigenvalue 
   real(Dp), allocatable :: eigenvalue(:)

   ! energy dispersion
   real(Dp),allocatable :: ekslab(:,:), ekslab_mpi(:,:)

   !> color for plot, surface state weight
   real(dp), allocatable :: surf_l_weight(:, :), surf_l_weight_mpi(:, :)
   real(dp), allocatable :: surf_r_weight(:, :), surf_r_weight_mpi(:, :)

   ! hamiltonian slab
   complex(Dp),allocatable ::CHamk(:,:)
   
   ! ---------------- NEW ----------------
   real(dp) ::  norm, minsx, maxsx, minsy, maxsy, minsz, maxsz
   real(dp), allocatable :: &
   Proj_Spin_Slab_x(:,:), Proj_Spin_Slab_x_mpi(:,:),&
   Proj_Spin_Slab_y(:,:), Proj_Spin_Slab_y_mpi(:,:), &
   Proj_Spin_Slab_z(:,:), Proj_Spin_Slab_z_mpi(:,:)
   real(dp), allocatable :: Proj_Atom_Slab(:, :, :), Proj_Atom_Slab_mpi(:, :, :)
   real(dp), allocatable :: Proj_Orbital_Slab(:, :, :), Proj_Orbital_Slab_mpi(:, :, :)
   complex(Dp),allocatable :: spin_x(:,:), spin_y(:,:), spin_z(:,:)
   type(int_array1D), allocatable :: AtomProjOrbitals(:)
   integer  :: num_wann_half, ib
   ! ---------------- END ----------------

   lwork= 16*Nslab*Num_wann
   ierr = 0


   allocate(eigenvalue(nslab*Num_wann))
   allocate( surf_l_weight (Nslab* Num_wann, knv2))
   allocate( surf_l_weight_mpi (Nslab* Num_wann, knv2))
   allocate( surf_r_weight (Nslab* Num_wann, knv2))
   allocate( surf_r_weight_mpi (Nslab* Num_wann, knv2))
   allocate(ekslab(Nslab*Num_wann,knv2))
   allocate(ekslab_mpi(Nslab*Num_wann,knv2))
   allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
   allocate(work(lwork))
   allocate(rwork(lwork))

   ! ---------------- NEW ----------------
   if (SlabBand_spin_proj) then
       allocate(spin_x(Num_wann*Nslab, Num_wann*Nslab))
       allocate(spin_y(Num_wann*Nslab, Num_wann*Nslab))
       allocate(spin_z(Num_wann*Nslab, Num_wann*Nslab))
       spin_x = 0.0d0; spin_y = 0.0d0; 
       spin_z = 0.0d0
       num_wann_half = Num_wann / 2
      !  do j = 1, num_wann_half
      !     spin_x(j, num_wann_half+j)=1.0d0
      !     spin_x(j+num_wann_half, j)=1.0d0
      !     spin_y(j, num_wann_half+j)=-zi
      !     spin_y(j+num_wann_half, j)=zi
      !     spin_z(j, j) = 1.0d0        ! spin up
      !     spin_z(j+num_wann_half, j+num_wann_half) = -1.0d0 ! spin down
      !  enddo
       do i=1, Nslab
         do j=1, num_wann_half
            spin_x(Num_wann*(i-1)+j,               Num_wann*(i-1)+j+num_wann_half)=1.0d0
            spin_x(Num_wann*(i-1)+j+num_wann_half, Num_wann*(i-1)+j)=1.0d0
            spin_y(Num_wann*(i-1)+j,               Num_wann*(i-1)+j+num_wann_half)=-zi
            spin_y(Num_wann*(i-1)+j+num_wann_half, Num_wann*(i-1)+j)=zi
            spin_z(Num_wann*(i-1)+j,               Num_wann*(i-1)+j)= 1d0
            spin_z(Num_wann*(i-1)+j+num_wann_half, Num_wann*(i-1)+j+num_wann_half)= -1d0
         enddo
         ! do j=1, num_wann_half
         !    spin_x(Num_wann*(i-1)+2*j-1, Num_wann*(i-1)+2*j  )=1.0d0
         !    spin_x(Num_wann*(i-1)+2*j  , Num_wann*(i-1)+2*j-1)=1.0d0
         !    spin_y(Num_wann*(i-1)+2*j-1, Num_wann*(i-1)+2*j  )=-zi
         !    spin_y(Num_wann*(i-1)+2*j  , Num_wann*(i-1)+2*j-1)=zi
         !    spin_z(Num_wann*(i-1)+2*j-1, Num_wann*(i-1)+2*j-1)= 1d0
         !    spin_z(Num_wann*(i-1)+2*j,   Num_wann*(i-1)+2*j  )= -1d0
         ! enddo
       enddo
       allocate(Proj_Spin_Slab_x(Nslab*Num_wann, knv2))
       allocate(Proj_Spin_Slab_x_mpi(Nslab*Num_wann, knv2))
       allocate(Proj_Spin_Slab_y(Nslab*Num_wann, knv2))
       allocate(Proj_Spin_Slab_y_mpi(Nslab*Num_wann, knv2))
       allocate(Proj_Spin_Slab_z(Nslab*Num_wann, knv2))
       allocate(Proj_Spin_Slab_z_mpi(Nslab*Num_wann, knv2))
       Proj_Spin_Slab_x = 0.d0; Proj_Spin_Slab_y = 0.d0; 
       Proj_Spin_Slab_z = 0.d0
   endif

   if (SlabBand_atom_proj) then
      allocate(AtomProjOrbitals(NumberofSelectedAtoms_groups))
      do ig=1, NumberofSelectedAtoms_groups
         io= 0
         do i=1, NumberofSelectedAtoms(ig)
            ia= Selected_Atoms(ig)%iarray(i)
            io= io+ Origin_cell%nprojs(ia)
         enddo
         if (SOC>0) io= io*2
         AtomProjOrbitals(ig)%length= io
         allocate(AtomProjOrbitals(ig)%iarray(io))
         AtomProjOrbitals(ig)%iarray= 0
         io= 0
         do i=1, NumberofSelectedAtoms(ig)
            ia= Selected_Atoms(ig)%iarray(i)
            do j=1, Origin_cell%nprojs(ia)
               io= io+ 1
               AtomProjOrbitals(ig)%iarray(io)= index_start(ia)+ j- 1
            enddo
         enddo
         if (SOC>0) then
            do i=1, NumberofSelectedAtoms(ig)
               ia= Selected_Atoms(ig)%iarray(i)
               do j=1, Origin_cell%nprojs(ia)
                  io= io+ 1
                  AtomProjOrbitals(ig)%iarray(io)= index_start(ia)+ j- 1+ Num_wann/2
               enddo
            enddo
         endif
      enddo

      if (cpuid==0) then
         write(stdout, '(a)') '>> AtomProjOrbitals groups and their iarray indices'
         do ig=1, NumberofSelectedAtoms_groups
            write(stdout, '(a, i4, a, i6)') '   group ', ig, ', length =', AtomProjOrbitals(ig)%length
            write(stdout, '(a, 500i6)') '      iarray =', &
               (AtomProjOrbitals(ig)%iarray(i), i=1, AtomProjOrbitals(ig)%length)
         enddo
      endif

      allocate(Proj_Atom_Slab(Nslab*Num_wann, knv2, NumberofSelectedAtoms_groups))
      allocate(Proj_Atom_Slab_mpi(Nslab*Num_wann, knv2, NumberofSelectedAtoms_groups))
      Proj_Atom_Slab= 0d0
      Proj_Atom_Slab_mpi= 0d0
   endif

   if (SlabBand_orbital_proj) then
      allocate(Proj_Orbital_Slab(Nslab*Num_wann, knv2, NumberofSelectedOrbitals_groups))
      allocate(Proj_Orbital_Slab_mpi(Nslab*Num_wann, knv2, NumberofSelectedOrbitals_groups))
      Proj_Orbital_Slab= 0d0
      Proj_Orbital_Slab_mpi= 0d0
      if (cpuid==0) then
         write(stdout, '(a)') 'WARNING: SlabBand_orbital_proj is currently a placeholder; output columns are zero.'
      endif
   endif
   ! ---------------- END ----------------

   surf_l_weight= 0d0
   surf_l_weight_mpi= 0d0
   surf_r_weight= 0d0
   surf_r_weight_mpi= 0d0

   ! sweep k
   ekslab=0.0d0
   ekslab_mpi=0.0d0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do i= 1+cpuid, knv2, num_cpu
      if (cpuid==0.and. mod(i/num_cpu, 4)==0) &
         write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         ' Slabek: ik', i, knv2, ' time left', &
         (knv2-i)*(time_end- time_start)/num_cpu, &
         ' time elapsed: ', time_end-time_start0 

      call now(time_start)

      k= k2_path(i, :)
      chamk=0.0d0 

      call ham_slab(k,Chamk)


      eigenvalue=0.0d0

      ! diagonal Chamk
      call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
     
      ekslab(:,i)=eigenvalue   ! (Nband,Nk)

      ! ------------ NEW ----------------
      if (SlabBand_spin_proj) then
          do ib = 1, Num_wann*Nslab
             !  norm = sum(abs(CHamk(:,ib))**2)
             ! CHamk(行, 列), 列是本征矢的索引
             Proj_Spin_Slab_x(ib, i) = &
                real(dot_product(CHamk(:, ib), matmul(spin_x, CHamk(:, ib))), dp) !/ norm
             Proj_Spin_Slab_y(ib, i) = &
                real(dot_product(CHamk(:, ib), matmul(spin_y, CHamk(:, ib))), dp) !/ norm
             Proj_Spin_Slab_z(ib, i) = &
                real(dot_product(CHamk(:, ib), matmul(spin_z, CHamk(:, ib))), dp) !/ norm
         enddo
      endif
      ! ------------ END ----------------

      ! H*chamk(:,n)=E(n)*chamk(:,n)
      !> Nslab*Num_wann
      !> rho(:)=abs(chamk(:,n))**2
      !> (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), 
      do j=1, Nslab* Num_wann
         !> left is the bottom surface
         do l= 1, NBottomOrbitals
            io= BottomOrbitals(l) 
            surf_l_weight(j, i)= surf_l_weight(j, i) &
               + abs(CHamk(io, j))**2  ! first slab
         enddo ! l sweeps the selected orbitals

         !> right is the top surface
         do l= 1, NTopOrbitals
            io= Num_wann*(Nslab-1)+ TopOrbitals(l)   
            surf_r_weight(j, i)= surf_r_weight(j, i) &
               + abs(CHamk(io, j))**2  ! first slab
         enddo ! l sweeps the selected orbitals

        !do l=1, Num_wann
        !   surf_l_weight(j, i)= surf_l_weight(j, i) &
        !      + abs(CHamk(l, j))**2  ! first slab
        !     !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
        !   surf_r_weight(j, i)= surf_r_weight(j, i) &
        !      + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
        !     !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
        !enddo ! l
      enddo ! j 

      if (SlabBand_atom_proj) then
         do ib=1, Num_wann*Nslab
            do ig=1, NumberofSelectedAtoms_groups
               do iq=1, Nslab
                  do l=1, AtomProjOrbitals(ig)%length
                     io= Num_wann*(iq-1)+ AtomProjOrbitals(ig)%iarray(l)
                     Proj_Atom_Slab(ib, i, ig)= Proj_Atom_Slab(ib, i, ig)+ abs(CHamk(io, ib))**2
                  enddo
               enddo
            enddo
         enddo
      endif
      call now(time_end)
   enddo ! i

#if defined (MPI)
   call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(surf_l_weight, surf_l_weight_mpi,size(surf_l_weight),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(surf_r_weight, surf_r_weight_mpi,size(surf_r_weight),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
         
   ! ------------ NEW ----------------
   if (SlabBand_spin_proj) then
      call mpi_allreduce(Proj_Spin_Slab_x, Proj_Spin_Slab_x_mpi,size(Proj_Spin_Slab_x),&
                         mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(Proj_Spin_Slab_y, Proj_Spin_Slab_y_mpi,size(Proj_Spin_Slab_y),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(Proj_Spin_Slab_z, Proj_Spin_Slab_z_mpi,size(Proj_Spin_Slab_z),&
                         mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif
   if (SlabBand_atom_proj) then
      call mpi_allreduce(Proj_Atom_Slab, Proj_Atom_Slab_mpi,size(Proj_Atom_Slab),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif
   if (SlabBand_orbital_proj) then
      call mpi_allreduce(Proj_Orbital_Slab, Proj_Orbital_Slab_mpi,size(Proj_Orbital_Slab),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif
   ! ------------ END ----------------
#else
   ekslab_mpi= ekslab
   surf_l_weight_mpi= surf_l_weight
   surf_r_weight_mpi= surf_r_weight

   ! ------------ NEW ----------------
   if (SlabBand_spin_proj) then
      Proj_Spin_Slab_x_mpi= Proj_Spin_Slab_x
      Proj_Spin_Slab_y_mpi= Proj_Spin_Slab_y
      Proj_Spin_Slab_z_mpi= Proj_Spin_Slab_z
   endif
   if (SlabBand_atom_proj) then
      Proj_Atom_Slab_mpi= Proj_Atom_Slab
   endif
   if (SlabBand_orbital_proj) then
      Proj_Orbital_Slab_mpi= Proj_Orbital_Slab
   endif
   ! ------------ END ----------------
#endif

   

   !> deal with phonon system
   if (index(Particle,'phonon')/=0) then
      do i=1, knv2
         do j=1, Num_wann*Nslab
            ekslab_mpi(j, i)= sqrt(abs(ekslab_mpi(j, i)))*sign(1d0, ekslab_mpi(j, i))
         enddo
      enddo
   endif
   ekslab_mpi= ekslab_mpi/eV2Hartree

   ekslab=ekslab_mpi

   maxweight=maxval(surf_r_weight_mpi+ surf_l_weight_mpi)
   surf_l_weight= surf_l_weight_mpi/ maxweight
   surf_r_weight= surf_r_weight_mpi/ maxweight
   
   ! ---------- NEW ---------------
   ! if (SlabBand_spin_proj) then
   !    minsx = minval(Proj_Spin_Slab_x_mpi)
   !    maxsx = maxval(Proj_Spin_Slab_x_mpi)
   !    Proj_Spin_Slab_x= 2.0_dp * (Proj_Spin_Slab_x_mpi - minsx) / (maxsx - minsx) - 1.0_dp
   !    minsy = minval(Proj_Spin_Slab_y_mpi)
   !    maxsy = maxval(Proj_Spin_Slab_y_mpi)
   !    Proj_Spin_Slab_y= 2.0_dp * (Proj_Spin_Slab_y_mpi - minsy) / (maxsy - minsy) - 1.0_dp
   !    minsz = minval(Proj_Spin_Slab_z_mpi)
   !    maxsz = maxval(Proj_Spin_Slab_z_mpi)
   !    Proj_Spin_Slab_z= 2.0_dp * (Proj_Spin_Slab_z_mpi - minsz) / (maxsz - minsz) - 1.0_dp
   ! endif
   if (SlabBand_spin_proj) then
      Proj_Spin_Slab_x= Proj_Spin_Slab_x_mpi
      Proj_Spin_Slab_y= Proj_Spin_Slab_y_mpi
      Proj_Spin_Slab_z= Proj_Spin_Slab_z_mpi
   endif
   if (SlabBand_atom_proj) Proj_Atom_Slab= Proj_Atom_Slab_mpi
   if (SlabBand_orbital_proj) Proj_Orbital_Slab= Proj_Orbital_Slab_mpi
   ! ---------- END ---------------
   
   outfileindex= outfileindex+ 1
   if(cpuid==0)then
      open(unit=outfileindex, file='slabek.dat')
      write(outfileindex, '(a)', advance='no') '# k  E  BS weight  TS weight'
      if (SlabBand_spin_proj) write(outfileindex, '(a)', advance='no') '  s_x  s_y  s_z'
      if (SlabBand_atom_proj) then
         do ig=1, NumberofSelectedAtoms_groups
            write(outfileindex, '(a, i3)', advance='no') '  atom_g', ig
         enddo
      endif
      if (SlabBand_orbital_proj) then
         do ig=1, NumberofSelectedOrbitals_groups
            write(outfileindex, '(a, i3)', advance='no') '  orbital_g', ig
         enddo
      endif
      write(outfileindex, *)

      do j=1, Num_wann*Nslab
         do i=1, knv2
            write(outfileindex,'(2f15.7, 2f16.7)', advance='no') &
               k2len(i)*Angstrom2atomic, ekslab(j,i), surf_l_weight(j,i), surf_r_weight(j,i)
            if (SlabBand_spin_proj) then
               write(outfileindex,'(3f16.7)', advance='no') &
                  Proj_Spin_Slab_x(j,i), Proj_Spin_Slab_y(j,i), Proj_Spin_Slab_z(j,i)
            endif
            if (SlabBand_atom_proj) then
               do ig=1, NumberofSelectedAtoms_groups
                  write(outfileindex,'(f16.7)', advance='no') Proj_Atom_Slab(j,i,ig)
               enddo
            endif
            if (SlabBand_orbital_proj) then
               do ig=1, NumberofSelectedOrbitals_groups
                  write(outfileindex,'(f16.7)', advance='no') Proj_Orbital_Slab(j,i,ig)
               enddo
            endif
            write(outfileindex,*)
         enddo
         write(outfileindex , *)''
      enddo

      close(outfileindex)
      write(stdout,*) 'calculate energy band  done'
   endif

   emin= minval(ekslab)-0.5d0
   emax= maxval(ekslab)+0.5d0
   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='slabek.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'slabek.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'slabek.png'"
      write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set ylabel font ",36"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'set ylabel offset -1, 0 '
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
      endif
      ! ! NEW ---------
      ! if (SlabBand_spin_proj) then
      !     write(outfileindex, '(a)') 'set palette defined(0 "blue", 1 "white", 2 "red")'
      !     write(outfileindex, '(a)') 'set cblabel "<{/Symbol s}_z>"'
      !     write(outfileindex, '(a)') 'set cbrange [-1:1]'
      !     write(outfileindex, '(a)') '# To plot with Sx or Sy, change the column index in "using"'
      !     write(outfileindex, '(a)') '# e.g., using 1:2:3 for Sx, 1:2:4 for Sy, 1:2:5 for Sz'
      !     write(outfileindex, '(a)') "plot 'slab_band.dat' using 1:2:5 with lines palette"
      ! else
      !     write(outfileindex, '(a)') "plot 'slab_band.dat' using 1:2 with lines"
      ! endif
      ! ! END -------------
      write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
      write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)*Angstrom2atomic

      do i=1, nk2lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, 0.0, k2line_stop(i+1)*Angstrom2atomic, emax
         else
            write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
         endif
      enddo
      write(outfileindex, '(a)')'#rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
          "w lp lw 2 pt 7  ps 1 lc rgb variable"
      write(outfileindex, '(2a)')"# (a) "
      write(outfileindex, '(2a)')"# plot the top and bottom surface's weight together"
      write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:($3+$4) ",  &
          "w lp lw 2 pt 7  ps 1 lc palette"
      write(outfileindex, '(2a)')"# (b) "
      write(outfileindex, '(2a)') &
         "# plot top and bottom surface's weight with red and blue respectively"
      write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
         '0 "grey", 1 "red" )'
      write(outfileindex, '(2a)')"plot 'slabek.dat' u 1:2:($4-$3) ",  &
          "w lp lw 2 pt 7  ps 1 lc palette"

      !write(outfileindex, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
      !   "w lp lw 2 pt 13 palette"
      close(outfileindex)
   endif

   202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
   203 format(A3,'" ',F10.5,')')
   204 format('set arrow from ',F10.5,',',F10.5, &
      ' to ',F10.5,',',F10.5, ' nohead')
 
   deallocate(eigenvalue)
   deallocate( surf_l_weight )
   deallocate( surf_l_weight_mpi )
   deallocate( surf_r_weight )
   deallocate( surf_r_weight_mpi )
   deallocate(ekslab)
   deallocate(ekslab_mpi)
   deallocate(CHamk)
   deallocate(work)
   deallocate(rwork)

   ! ---------- NEW ---------------
   if (allocated(Proj_Spin_Slab_x)) deallocate( Proj_Spin_Slab_x )
   if (allocated(Proj_Spin_Slab_x_mpi)) deallocate( Proj_Spin_Slab_x_mpi )
   if (allocated(Proj_Spin_Slab_y)) deallocate( Proj_Spin_Slab_y )
   if (allocated(Proj_Spin_Slab_y_mpi)) deallocate( Proj_Spin_Slab_y_mpi )
   if (allocated(spin_x)) deallocate( spin_x )
   if (allocated(spin_y)) deallocate( spin_y )
   if (allocated(Proj_Spin_Slab_z)) deallocate( Proj_Spin_Slab_z )
   if (allocated(Proj_Spin_Slab_z_mpi)) deallocate( Proj_Spin_Slab_z_mpi )
   if (allocated(spin_z)) deallocate( spin_z )
   if (allocated(Proj_Atom_Slab)) deallocate(Proj_Atom_Slab)
   if (allocated(Proj_Atom_Slab_mpi)) deallocate(Proj_Atom_Slab_mpi)
   if (allocated(Proj_Orbital_Slab)) deallocate(Proj_Orbital_Slab)
   if (allocated(Proj_Orbital_Slab_mpi)) deallocate(Proj_Orbital_Slab_mpi)
   if (allocated(AtomProjOrbitals)) then
      do ig=1, size(AtomProjOrbitals)
         if (allocated(AtomProjOrbitals(ig)%iarray)) deallocate(AtomProjOrbitals(ig)%iarray)
      enddo
      deallocate(AtomProjOrbitals)
   endif
   ! ---------- END ---------------
 
return
end subroutine ek_slab

subroutine ek_slab_sparseHR
 use para
 use sparse
 implicit none

 !> some temporary integers
 integer :: ik, ia1, ia2, i, j, ierr, ib, iq, ig

 ! wave vector
 real(dp) :: k(2)

 !> dim= Ndimq, knv3
 integer :: Ndimq
 real(dp), allocatable :: W(:)
 real(dp), allocatable :: eigv(:, :)
 real(dp), allocatable :: eigv_mpi(:, :)

 real(dp) :: emin, emax
 real(dp) :: time_start, time_end, time_start0

 integer :: nnzmax, nnz
 complex(dp), allocatable :: acoo(:)
 integer, allocatable :: jcoo(:)
 integer, allocatable :: icoo(:)

 !> eigenvector of the sparse matrix acoo. Dim=(Ndimq, neval)
 complex(dp), allocatable :: psi(:)
 complex(dp), allocatable :: zeigv(:, :)

 !> print the weight for the Selected_WannierOrbitals
 real(dp), allocatable :: dos_selected(:, :, :)
 real(dp), allocatable :: dos_selected_mpi(:, :, :)

 real(dp), allocatable :: dos_l_selected(:, :, :)
 real(dp), allocatable :: dos_l_selected_mpi(:, :, :)
 real(dp), allocatable :: dos_r_selected(:, :, :)
 real(dp), allocatable :: dos_r_selected_mpi(:, :, :)

 !number of ARPACK eigenvalues
 integer :: neval

 ! number of Arnoldi vectors
 integer :: nvecs

 !shift-invert sigma
 complex(dp) :: sigma

 !> time measurement
 real(dp) :: time1, time2, time3

 logical :: ritzvec
 
 Ndimq= Num_wann* Nslab
 nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
 if(Is_Sparse_Hr) nnzmax=splen*Nslab+Ndimq
 if (NumSelectedEigenVals==0) NumSelectedEigenVals=Ndimq
 neval=NumSelectedEigenVals
 if (neval>=Ndimq) neval= Ndimq- 2

 !> ncv
 nvecs=int(2*neval)

 ! if (nvecs<50) nvecs= 50
 if (nvecs>Ndimq) nvecs= Ndimq


 sigma=(1d0,0d0)*E_arc
 
 allocate( acoo(nnzmax), stat= ierr)
 call printallocationinfo('acoo', ierr)
 allocate( jcoo(nnzmax), stat= ierr)
 call printallocationinfo('jcoo', ierr)
 allocate( icoo(nnzmax), stat= ierr)
 call printallocationinfo('icoo', ierr)
 allocate( W( neval), stat= ierr)
 allocate( eigv( neval, knv2))
 call printallocationinfo('eigv', ierr)
 allocate( eigv_mpi( neval, knv2), stat= ierr)
 call printallocationinfo('eigv_mpi', ierr)
 allocate( psi(ndimq))
 allocate( zeigv(ndimq,nvecs), stat= ierr)
 call printallocationinfo('zeigv', ierr)
 allocate( dos_selected     (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_selected', ierr)
 allocate( dos_selected_mpi (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_selected_mpi', ierr)
 allocate( dos_l_selected     (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_l_selected', ierr)
 allocate( dos_l_selected_mpi (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_l_selected_mpi', ierr)
 allocate( dos_r_selected     (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_r_selected', ierr)
 allocate( dos_r_selected_mpi (neval,   knv2, NumberofSelectedOrbitals_groups), stat= ierr)
 call printallocationinfo('dos_r_selected_mpi', ierr)
 dos_l_selected= 0d0
 dos_l_selected_mpi= 0d0
 dos_r_selected= 0d0
 dos_r_selected_mpi= 0d0
 dos_selected= 0d0
 dos_selected_mpi= 0d0

 eigv_mpi= 0d0
 eigv    = 0d0
 acoo= 0d0

 !> calculate the along special k line
 time_start= 0d0
 time_start0= 0d0
 call now(time_start0)
 time_start= time_start0
 time_end  = time_start0
 ritzvec= .true.
 do ik=1+ cpuid, knv2, num_cpu
    if (cpuid==0.and. mod(ik/num_cpu, 4)==0) &
       write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
       ' Slabek: ik', ik, knv2, ' time left', &
       (knv2-ik)*(time_end- time_start)/num_cpu, &
       ' time elapsed: ', time_end-time_start0 

    call now(time_start)
    
    k= k2_path(ik, :)
    nnz= nnzmax
    call now(time1)
    call ham_slab_sparseHR(nnz, k, acoo,jcoo,icoo)
    call now(time2)

    !> diagonalization by call zheev in lapack
    W= 0d0

#if defined (INTELMKL)
    call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
#endif

    call now(time3)
    eigv(1:neval, ik)= W(1:neval)

    !> calculate the weight on the selected orbitals
    do ib= 1, neval
       psi(:)= zeigv(:, ib)  !> the eigenvector of ib'th band
       do ig=1, NumberofSelectedOrbitals_groups
          do iq=1, Nslab
             do i= 1, NumberofSelectedOrbitals(ig)
                j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2
                
             enddo ! sweep the selected orbitals
          enddo ! iq sweep the magnetic supercell

          do iq=1, 2  ! edge states
             if (iq>Nslab) cycle
             do i= 1, NumberofSelectedOrbitals(ig)
                j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                dos_l_selected(ib, ik, ig)= dos_l_selected(ib, ik, ig)+ abs(psi(j))**2
             enddo ! sweep the selected orbitals
          enddo ! iq sweep the magnetic supercell
          do iq=Nslab-1, Nslab ! edge states
             if (iq<1) cycle
             do i= 1, NumberofSelectedOrbitals(ig)
                j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                dos_r_selected(ib, ik, ig)= dos_r_selected(ib, ik, ig)+ abs(psi(j))**2
             enddo ! sweep the selected orbitals
          enddo ! iq sweep the magnetic supercell
       enddo ! ig
    enddo ! ib sweep the eigenvalue

    if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time1, ' s'
    if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time3-time2, ' s'
    call now(time_end)
 enddo !ik

#if defined (MPI)
 call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
    mpi_dp,mpi_sum,mpi_cmw,ierr)
 call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
    mpi_dp,mpi_sum,mpi_cmw,ierr)
    
 call mpi_allreduce(dos_l_selected, dos_l_selected_mpi,size(dos_l_selected),&
    mpi_dp,mpi_sum,mpi_cmw,ierr)
 call mpi_allreduce(dos_r_selected, dos_r_selected_mpi,size(dos_r_selected),&
    mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
 eigv_mpi= eigv
 dos_selected_mpi= dos_selected
 dos_l_selected_mpi= dos_l_selected
 dos_r_selected_mpi= dos_r_selected
#endif

 !> minimum and maximum value of energy bands
 emin= minval(eigv_mpi)-0.5d0
 emax= maxval(eigv_mpi)+0.5d0

 eigv_mpi= eigv_mpi/eV2Hartree

 outfileindex= outfileindex+ 1
 if (cpuid.eq.0) then
    open(unit=outfileindex, file='slabek.dat')
    write(outfileindex, '("#", a14, a15, a)')'k ', ' E', ' Weight on the selected orbitals'
    do j=1, neval
       do i=1,knv2
          write(outfileindex,'(200f16.8)')k2len(i)*Angstrom2atomic, eigv_mpi(j, i), &
             (dos_l_selected_mpi(j, i, ig), dos_r_selected_mpi(j, i, ig), &
             ig=1, NumberofSelectedOrbitals_groups)
       enddo
       write(outfileindex , *)''
    enddo
    close(outfileindex)
    write(stdout,*) 'calculate landau level done'
 endif

 outfileindex= outfileindex+ 1

 emin= minval(eigv_mpi)-0.5d0
 emax= maxval(eigv_mpi)+0.5d0
 !> write script for gnuplot
 outfileindex= outfileindex+ 1
 if (cpuid==0) then
    open(unit=outfileindex, file='slabek.gnu')
    write(outfileindex, '(a)')"set encoding iso_8859_1"
    write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
    write(outfileindex, '(a)')"#set output 'slabek.eps'"
    write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
       '  font ",60" size 1920, 1680'
    write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
       '  font ",60" size 1920, 1680'
    write(outfileindex, '(a)')"set output 'slabek.png'"
    write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
       '5 "yellow", 10 "red" )'
    write(outfileindex, '(a)')'set style data linespoints'
    write(outfileindex, '(a)')'unset ztics'
    write(outfileindex, '(a)')'unset key'
    write(outfileindex, '(a)')'set pointsize 0.8'
    write(outfileindex, '(a)')'set border lw 3 '
    write(outfileindex, '(a)')'set view 0,0'
    write(outfileindex, '(a)')'#set xtics font ",36"'
    write(outfileindex, '(a)')'#set ytics font ",36"'
    write(outfileindex, '(a)')'#set ylabel font ",36"'
    write(outfileindex, '(a)')'#set xtics offset 0, -1'
    write(outfileindex, '(a)')'set ylabel offset -1, 0 '
    write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len*Angstrom2atomic), ']'
    if (index(Particle,'phonon')/=0) then
       write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
       write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
    else
       write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
       write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
    endif
    write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
    write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)*Angstrom2atomic

    do i=1, nk2lines-1
       if (index(Particle,'phonon')/=0) then
          write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, 0.0, k2line_stop(i+1)*Angstrom2atomic, emax
       else
          write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
       endif
    enddo
    write(outfileindex, '(a)')'#rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
    write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
        "w lp lw 2 pt 7  ps 1 lc rgb variable"
    write(outfileindex, '(2a)')"# plot the top and bottom surface's weight together"
    write(outfileindex, '(2a)')"plot 'slabek.dat' u 1:2:($3+$4) ",  &
        "w lp lw 2 pt 7  ps 1 lc palette"
    write(outfileindex, '(2a)')"#"
    write(outfileindex, '(2a)')"# plot top and bottom surface's weight with different color"
    write(outfileindex,'(2a)') '#set palette defined ( -1  "blue", ', &
       '0 "grey", 1 "red" )'
    write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:($4-$3) ",  &
        "w lp lw 2 pt 7  ps 1 lc palette"

    !write(outfileindex, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
    !   "w lp lw 2 pt 13 palette"
    close(outfileindex)
 endif

 202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
 203 format(A3,'" ',F10.5,')')
 204 format('set arrow from ',F10.5,',',F10.5, &
      ' to ',F10.5,',',F10.5, ' nohead')
 


#if defined (MPI)
 call mpi_barrier(mpi_cmw, ierr)
#endif

 deallocate( acoo)
 deallocate( jcoo)
 deallocate( icoo)
 deallocate( W)
 deallocate( eigv)
 deallocate( eigv_mpi)
 deallocate( zeigv)
 deallocate( dos_selected)
 deallocate( dos_selected_mpi)

 return
end subroutine ek_slab_sparseHR


subroutine ek_slab_kplane
   !> This subroutine is used for calculating energy 
   !> dispersion with wannier functions for 2D slab system
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  
   use wmpi
   use para
   implicit none 


   ! loop index
   integer :: i, j, l, lwork, ierr, kn12, ik, istart, iend

   ! wave vector 
   real(Dp) :: k(2)

   real(Dp) :: time_start, time_end

   real(Dp), allocatable :: eigenvalue(:)
 
   ! energy dispersion
   real(Dp),allocatable :: ekslab(:,:)
   real(Dp),allocatable :: ekslab_mpi(:,:)

   real(dp), allocatable :: k12(:,:)
   real(dp), allocatable :: k12_shape(:,:)

   !> color for plot, surface state weight
   real(dp), allocatable :: surf_weight(:, :)
   real(dp), allocatable :: surf_weight_mpi(:, :)

   complex(Dp),allocatable ::CHamk(:,:)

   kn12= nk1*nk2
   lwork= 16*Nslab*Num_wann
   ierr = 0

   allocate(eigenvalue(nslab*Num_wann))
   allocate( surf_weight (Nslab* Num_wann, kn12))
   allocate( surf_weight_mpi (Nslab* Num_wann, kn12))
   allocate(ekslab(Nslab*Num_wann,kn12))
   allocate(ekslab_mpi(Nslab*Num_wann,kn12))
   allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
   allocate(k12(2,kn12))
   allocate(k12_shape(2,kn12))

   surf_weight= 0d0
   surf_weight_mpi= 0d0

   !> set up k slice
   ik =0
   do i= 1, nk1
      do j= 1, nk2
         ik =ik +1
         k12(:, ik)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                    + (j-1)*K2D_vec2/dble(nk2-1)
         k12_shape(:, ik)= k12(1, ik)* Ka2+ k12(2, ik)* Kb2
      enddo
   enddo

   ! sweep k
   ekslab=0.0d0
   ekslab_mpi=0.0d0
   time_start= 0d0
   time_end= 0d0
   do i=1+cpuid, kn12, num_cpu
      if (cpuid==0.and. mod(i/num_cpu, 100)==0) &
         write(stdout, *) 'SlabBand_plane, ik ', i, 'Nk',nk1*nk2, 'time left', &
         (nk1*nk2-i)*(time_end- time_start)/num_cpu, ' s'
      call now(time_start)
      
      k= k12(:, i)
      chamk=0.0d0 

      call ham_slab(k,Chamk)

      eigenvalue=0.0d0

      ! diagonal Chamk
      call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
     
      ekslab(:,i)=eigenvalue

      do j=1, Nslab* Num_wann
         do l=1, Num_wann
            surf_weight(j, i)= surf_weight(j, i) &
               + abs(CHamk(l, j))**2 & ! first slab
               + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
              !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
              !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
         enddo ! l
        !surf_weight(j, i)= (surf_weight(j, i))
      enddo ! j 
      call now(time_end)
   enddo ! i

#if defined (MPI)
   call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   ekslab_mpi= ekslab
   surf_weight_mpi= surf_weight
#endif

   !> deal with phonon system
   if (index(Particle,'phonon')/=0) then
      do i=1, kn12
         do j=1, Num_wann*Nslab
            ekslab_mpi(j, i)= sqrt(abs(ekslab_mpi(j, i)))*sign(1d0, ekslab_mpi(j, i))
         enddo
      enddo
   endif
   ekslab_mpi= ekslab_mpi/eV2Hartree

   ekslab=ekslab_mpi
   if (maxval(surf_weight_mpi)<0.00001d0)surf_weight_mpi=1d0
   surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)
   

   istart= Numoccupied*Nslab-1
   iend= Numoccupied*Nslab+2
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='slabek_plane.dat')
      write(outfileindex, '(4a16, a)')'# kx', ' ky', ' k1', ' k2', ' (E(ib), dos(ib)), ib=1, NumberofSelectedOrbitals'
      write(outfileindex, '(a, 2i10)')'# Nk1, Nk2=', Nk1, Nk2
      do i=1, kn12
         write(outfileindex,'(2000f16.7)')k12_shape(:,i), k12(:,i), &
            (ekslab(j,i), (255-surf_weight(j, i)*255d0), j=istart, iend)
         if (mod(i, nk1)==0) write (outfileindex, *)' '
      enddo
      close(outfileindex)
      write(stdout,*) 'calculate energy band  done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='slabek_plane-matlab.dat')
      write(outfileindex, '(4a16, a)')'% kx', ' ky', ' k1', ' k2', ' (E(ib), dos(ib)), ib=1, NumberofSelectedOrbitals'
      write(outfileindex, '(a, 2i10)')'% Nk1, Nk2=', Nk1, Nk2
      do i=1, kn12
         write(outfileindex,'(2000f16.7)')k12_shape(:,i), k12(:,i), &
            (ekslab(j,i), (255-surf_weight(j, i)*255d0), j=istart, iend)
      enddo
      close(outfileindex)
   endif


   !> write out a script that can be used for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='slabek_plane.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'slabek_plane.eps'"
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         ' size 1920, 1680 font ",36"'
      write(outfileindex, '(a)')"set output 'slabek_plane.png'"
      write(outfileindex, '(a)')'set palette rgbformulae 33,13,10'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set origin 0.2, 0'
      write(outfileindex, '(a)')'set size 0.8, 1'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'set size ratio -1'
      write(outfileindex, '(a)')'set xtics'
      write(outfileindex, '(a)')'set ytics'
      write(outfileindex, '(a)')'set view 80,60'
      write(outfileindex, '(a)')'set xlabel "k_1"'
      write(outfileindex, '(a)')'set ylabel "k_2"'
      write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
      write(outfileindex, '(a)')'unset colorbox'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a)')'set pm3d interpolate 4,4'
      write(outfileindex, '(2a)')"splot 'slabek_plane.dat' u 1:2:7 w pm3d, \"
      write(outfileindex, '(2a)')"      'slabek_plane.dat' u 1:2:9 w pm3d"

      close(outfileindex)

   endif ! cpuid

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate(eigenvalue)
   deallocate( surf_weight )
   deallocate( surf_weight_mpi )
   deallocate(ekslab)
   deallocate(ekslab_mpi)
   deallocate(CHamk)
 
   return
end subroutine ek_slab_kplane
