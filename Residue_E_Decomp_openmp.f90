  program Residue_E_Decomp
!
!  written by Lalith Perera for energy analysis
!  June 19 2009
!  modified by G. Andres Cisneros to f90
!  June 30 2009
!  attempted to mpi-ify by EML
!  ?? ?? 2019
!  OMP parallel GAC
!  Aug 25 2020
!compile: ifort Residue_E_Decomp_openmp.f90 -o Residue_E_Decomp_openmp.x -qopenmp
!
!
!  constances....
!

!! Commence original
      implicit none

! nat_max-total number of atoms
! nres_max - total number of residues

      integer nat_max,nres_max,nprotat,ntype_max,nres,npairs,natpairs
      !parameter (nat_max=62607,nprotat=4623,nres=273)
      !parameter (nres_max=19601,ntype_max=1000)
      !parameter (npairs=(nres*(nres-1))/2)
      !parameter (natpairs=(nprotat*(nprotat+1))/2)
        real*8 coulomb_const_kcal_per_mole
        !real *8 at_chg(nat_max)
  double precision,allocatable::at_chg(:)
        !real *8 xlj_12(ntype_max),xlj_6(ntype_max)
  double precision,allocatable::xlj_12(:),xlj_6(:)
        !real *8 charge_ij(natpairs)
  double precision,allocatable::charge_ij(:)
        !real *8 vdw_12(natpairs)
        !real *8  vdw_6(natpairs)
  double precision,allocatable::vdw_12(:),vdw_6(:)
        !real *8 LJ_12(ntype_max,ntype_max)
        !real *8 LJ_6(ntype_max,ntype_max)
  double precision,allocatable::LJ_12(:,:),LJ_6(:,:)
        real *8 xd,yd,zd,rd,rd6,rd12
        real *8 E_coul,E_vdW
        !real *8 E_coul_bin(npairs),E_vdW_bin(npairs)
        !real *8 SE_coul_bin(npairs),SE_vdW_bin(npairs)
  double precision,allocatable::E_coul_bin(:),E_vdW_bin(:),SE_coul_bin(:),SE_vdW_bin(:)
        real *8 sigma,eps
        real *8 sig1,eps1
        real *8 sig2,eps2
        real *8 re,re1,re2
        !real *8 pos(3,nat_max)
  double precision,allocatable::pos(:,:)
        real *8 bx,by,bz


        integer n_res,n_files,allochk
        integer natom
        !integer i_res_point(nres_max)
        !integer at_type_index(nat_max)
  integer,allocatable::i_res_point(:),at_type_index(:)
        integer nlast_at
        integer k1,k2,k3,k4,k5,k6,k7,k8,k9
        integer m1,m2,i1,i2
        integer k_bin
        !integer nat_in_res(nres_max)
  integer,allocatable::nat_in_res(:)
        integer at_id_i,at_type_i
        integer at_id_j,at_type_j
        !integer iflag(natpairs)
        !integer at_index_i(natpairs),at_index_j(natpairs)
  integer,allocatable::iflag(:),at_index_i(:),at_index_j(:)
        integer icc
        integer k_bin_tot,kbin
        integer ntypes,nttyp
        integer nsam
        integer rdstat, n_file_sam

!! Now the ones for OpenMP
  integer proc_num, thread_num, my_thread_num, nomptasks !GAC
  integer omp_get_num_procs, omp_get_thread_num, omp_get_max_threads !GAC
          ! define omp functions that return ints

        !character *3 res_name(nres_max)
  !character*3,allocatable::res_name(:)
  character*4,allocatable::res_name(:)
        !character *4 am_at_type(nat_max),std_am_at_type(ntype_max)
  character*4,allocatable::am_at_type(:),std_am_at_type(:)
        character *80 filename,inpfile,prmfile,junk


      coulomb_const_kcal_per_mole = 332.0538200

!!! Create an OpenMP outfile
!  open(unit=872,file='EDA_openmp.log')
!  write(872, '(a)' ) ' '
!  write(872, '(a)' ) 'EDA Program: F90 openmp'
!  write(872, '(a)' ) ' '
!  write(872, '(a,i8,a)' ) ' Program found ', proc_num, ' available processors.'
!  write(872, '(a,i8,a)' ) ' Program found ', thread_num, ' available threads.'
!  close(872)
!!!


!
        !open(801,file='NEW_RES_ENERGY_analysis.inp')
        write(6,*)'Name of input file?'
        read(5,*),inpfile
        inpfile=trim(inpfile)
        !inpfile='Residue_E_Decomp_AA_check_10ns.inp'
        write(6,*)'Name of prmtop file?'
        read(5,*),prmfile
        prmfile=trim(prmfile)
        open(801,file=inpfile)
        !prmfile='AA_check_neutral.prmtop'
        open(801,file=inpfile)
!    n_res = total number of non-water residues
!    n_files = number of mdcrd files....
        read(801,*) n_res,junk
        print *,n_res
        read(801,*) n_files,junk
        print *,n_files
        read(801,*) nat_max,junk
        print *,nat_max
        read(801,*) nprotat,junk
        print *,nprotat
        !read(801,*) nres,junk
        read(801,*) nres_max,junk
        print *,nres_max
        read(801,*) ntype_max,junk
        print *,ntype_max

        nres = n_res
        npairs=(nres*(nres-1))/2
        natpairs=(nprotat*(nprotat+1))/2

! GAC: allocate arrays

  allocate(at_chg(nat_max), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate at_chg, exiting'
     stop
  endif
  allocate(xlj_12(ntype_max),xlj_6(ntype_max), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate xlj_12, xlj_6, exiting'
     stop
  endif
  allocate(charge_ij(natpairs),vdw_12(natpairs),vdw_6(natpairs), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate charge_ij,vdw_12,vdw6, exiting'
     stop
  endif
  allocate(LJ_12(ntype_max,ntype_max),LJ_6(ntype_max,ntype_max), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate LJ_12, LJ_6, exiting'
     stop
  endif
  allocate(E_coul_bin(npairs),E_vdW_bin(npairs),SE_coul_bin(npairs),&
           SE_vdW_bin(npairs), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate E_coul_bin,E_vdW_bin,SE_coul_bin,&
     &SE_vdW_bin, exiting'
     stop
  endif
  allocate(pos(3,nat_max), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate pos, exiting'
     stop
  endif
  allocate(i_res_point(nres_max),at_type_index(nat_max),nat_in_res(nres_max),&
           stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate i_res_point,at_type_index,nat_in_res,&
                & exiting'
     stop
  endif
  allocate(iflag(natpairs),at_index_i(natpairs),at_index_j(natpairs), stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate iflag,at_index_i,at_index_j, exiting'
     stop
  endif
  allocate(res_name(nres_max),am_at_type(nat_max),std_am_at_type(ntype_max),&
           stat = allochk)
  if (allochk .gt. 0 ) then
     write(6,*)'E_decomp:could not allocate res_name,am_at_type,std_atm_at_type,&
                & exiting'
     stop
  endif

!!  if ( proc_id == 0 ) then
    write(6,*)' Finished allocating arrays'
!!  end if

! GAC finished allocation

!  get information from "prmtop"
!  natom = total number of atoms
! ntypes = total atom types (for van der Waals)
! nttyp = (ntypes*(ntypes+1))/2
! i_res_point = residue pointer (starting atoms)
! at_type_index = atom type index for vdW
! at_chg = charges
! xlj_12 = LJ -A12 terms (from prmtop after combining)
! xlj_6 = LJ -6 terms (from prmtop after combining)
! res_name = residue names
! am_at_type = amber atom types

        call rdparm(natom,nat_max,ntypes,i_res_point,at_type_index,&
             at_chg,xlj_12,xlj_6,&
              res_name,am_at_type,prmfile)

! nlast_at = The last atom to be considered
        nlast_at = i_res_point(n_res+1) -1
        nttyp = ((ntypes+1)*ntypes)/2
        print *, 'nttyp = ',nttyp

        std_am_at_type(1) = am_at_type(1)
        k3 = 1
        do 10 k1=2,natom
        icc   = 0
        do 20 k2=1,k3
        if(am_at_type(k1).eq.std_am_at_type(k2)) icc =  1
 20     continue
        if(icc.eq.0) then
        k3 = k3 + 1
        std_am_at_type(k3) = am_at_type(k1)
        print *, k3,std_am_at_type(k3)
        endif
 10     continue

!
!! Change the output name
        open(unit=804,file='fort_sanity_check.txt')
!!
        do 30 k1=1,k3
        write(804,901) k1,std_am_at_type(k1)
 30     continue
        write(804,*) ' done writing std_am_at_types'
 901    format(i10,a4)

        do 40 k1=1,ntypes
        do 50 k2=1,natom
        if(k1.eq.at_type_index(k2)) then
        write(804,*) k1,am_at_type(k2)
        goto 40
        endif
 50     continue
 40     continue
        write(804,*) ' done writing atom types'

!
!   van der Waals parameters in terms of pairs....
!
        k1 = 0
        print *, 'ntypes = ',ntypes
        do 60 k2=1,ntypes
        do 60 k3=1,k2
        k1 = k1 + 1
        LJ_12(k2,k3) = xlj_12(k1)
        LJ_12(k3,k2) = xlj_12(k1)
        LJ_6(k2,k3) = xlj_6(k1)
        LJ_6(k3,k2) = xlj_6(k1)
        if(abs(xlj_6(k1)).gt.0.000000001.or. &
        abs(xlj_12(k1)).gt.0.000000001) then
        sigma = (xlj_12(k1)/xlj_6(k1))**(1.0/6.0)
        eps = (xlj_6(k1)**2)/(4.0*xlj_12(k1))
        else
        sigma = 0.0
        eps = 0.0
        endif
        re = sigma * 2.0**(1.0/6.0)
        if(k2.eq.k3) then
        write(804,902) k2,k3,re/2.0,eps
        endif
 60     continue
 902    format(2i10,4f10.4)
        write(804,*) ' done writing LJ parameters ',k1
!
! total number of atoms in a residue .......
!
        do 100 k1=1,nres
        nat_in_res(k1) = i_res_point(k1+1) - i_res_point(k1)
 100    continue
        write(804,910) (nat_in_res(k1),k1=1,nres)
 910    format(20i4)
         write(804,*) ' done writing number of atoms in a residue'
!
!   set up bins...
!
        do 110 k1=1,npairs
        E_vdW_bin(k1) = 0.0
        E_coul_bin(k1) = 0.0
        SE_vdW_bin(k1) = 0.0
        SE_coul_bin(k1) = 0.0
 110    continue
!
!   extract charges and assign them in a new array
!
        k1 = 0
        k_bin = 0

!    residue k2 ...
        do 120 k2=1,nres-1

!    residue k4 ...
!   ********WARNING************************
!   K4 should start two residues away from k2. Otherwise, 1-2,1-3 and 1-4
!   interactions will result in erroneous energies.....
!
        do 130 k4=k2+1,nres
         k_bin = k_bin + 1

!     atom k3 of residue k2
          do 140 k3=1,nat_in_res(k2)
        at_id_i = i_res_point(k2) + k3 - 1
        at_type_i = at_type_index(at_id_i)

!     atom k5 of residue k4
           do 150 k5=1,nat_in_res(k4)
        at_id_j = i_res_point(k4) + k5 - 1
        at_type_j = at_type_index(at_id_j)

        k1 = k1 + 1

!     assign the charge
        charge_ij(k1) =  at_chg(at_id_i) * at_chg(at_id_j)
!     assign LJ parameters for each atom
        vdw_12(k1) =  LJ_12(at_type_i,at_type_j)
        vdw_6(k1)  =   LJ_6(at_type_i,at_type_j)
        iflag(k1) = k_bin
        at_index_i(k1) = at_id_i
        at_index_j(k1) = at_id_j
 150    continue
 140    continue
 130    continue
 120    continue
 911    format(5i10,4i5)
        k_bin_tot = k_bin
        print *, 'k_bin_tot = ',k_bin_tot

!
! the real calculation starts here
!

       nsam = 0
       do 200 i1 =1,n_files
       read(801,*) filename
       print *, filename
       open(802,file=filename)
       read(802,*)
       n_file_sam = 0
       do while (rdstat == 0)
          read(802,903,iostat=rdstat) pos
          read(802,903,iostat=rdstat) bx,by,bz
          n_file_sam = n_file_sam + 1
       enddo
       rdstat = 0
       rewind(802)
       read(802,*)
       print *,'number of samples in file: n_file_sam = ', n_file_sam
       !do 210 i2=1,250
       do 210 i2=1,n_file_sam
       read(802,903,end=220,err=220) pos
       read(802,903,end=220,err=220) bx,by,bz
 903   format(10f8.3)
       nsam = nsam + 1

! GAC, this is what we want to paralellize
!$omp parallel do

    do k1 = 1, natpairs

        m1 = at_index_i(k1)
        m2 = at_index_j(k1)

        xd = pos(1,m1) - pos(1,m2)
        yd = pos(2,m1) - pos(2,m2)
        zd = pos(3,m1) - pos(3,m2)
        rd = xd**2 + yd**2 + zd**2
        rd6 = rd * rd * rd
        rd12 = rd6 * rd6

! Coulomb energy
        E_coul = charge_ij(k1) / sqrt(rd)
! van der Waals Energy
        E_vdW  = vdw_12(k1)/rd12 - vdw_6(k1)/rd6
!
!  binning.....
!
        k_bin = iflag(k1)
        E_coul_bin(k_bin) = E_coul_bin(k_bin) + E_coul
        E_vdW_bin(k_bin) = E_vdW_bin(k_bin) + E_vdW
        SE_coul_bin(k_bin) = SE_coul_bin(k_bin) + E_coul**2
        SE_vdW_bin(k_bin) = SE_vdW_bin(k_bin) + E_vdW**2

   enddo
!$omp end parallel do
! GAC finished parallel part

 210    continue

 220    continue
        close(802)
 200    continue

!! Change the names of 803 and 806
        open(unit=803,file='fort_coulomb_interaction.dat')
        open(unit=806,file='fort_vdw_interaction.dat')
!!

        write(803,904) nsam
!! Add descriptive column names
        write(803,543) 'Index','ResA','ResB','CoulEnergy','StdErr'
        write(806,543) 'Index','ResA','ResB','VdWEnergy','StdErr'
!! Format them accordingly, assign that format a unit number
!! Ex: 3a10 means first 3 columns, character, 10 max positions
 543    format(3a10,2a19)
!! The 2e20.12 format below means 20 total numbers, 12 b4 decimal
!!
        k1 = 0
        do 500 k2=1,nres-1
        do 500 k4=k2+1,nres
        k1 = k1 + 1
        !write(805,906) k1,k2,k4,res_name(k2),res_name(k4)
        write(803,905) k1,k2,k4,E_coul_bin(k1)/nsam,SE_coul_bin(k1)/nsam
        write(806,906) k1,k2,k4,E_vdW_bin(k1)/nsam,SE_vdW_bin(k1)/nsam
 500    continue
 904    format(i10)
 905    format(3i10,2e20.12)
 906    format(3i10,2e20.12)
! 906    format(3i10,2a6)
        stop
        end

      subroutine rdparm(natom,nat,ntypes,i_res_point,at_type_index,&
             at_chg,xlj_12,xlj_6,&
              res_name,am_at_type,prmfile)

!
!  subroutine reading information from prmtop
!
      integer natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia,nat
      integer nhparm,nparm,nnb,nres,nbona,ntheta,nphia
      integer numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper
      integer mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
      integer nf,nttyp,ntype,allochk
      integer ix_2(nat),ix_3(nat),ix_4(nat),ix_5(nat)
      integer i_res_point(nat), at_type_index(nat)
      integer, allocatable::ix_1(:)
      real *8 ax_1(nat*3)
      real *8 at_chg(nat),at_mass(nat)
      real *8 xlj_12(1000),xlj_6(1000)
      character *4 at_name(nat),res_name(nat),am_at_type(nat)
      character *4 tree_ch_cl(nat)
      character *80 ixxxxxx,prmfile


      nf = 22

      !open(nf,file='prmtop')
      open(nf,file=prmfile)


         read(nf,*)
         read(nf,*)
         read(nf,*)
         read(nf,*)
         read(nf,*)
         read(nf,*)
      read(nf,101) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia,&
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,&
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper,&
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
      print  201, natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia,&
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd,&
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
 101   format(10I8)

 201   format(t2,&
         'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7,&
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7,&
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7,&
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7,&
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7,&
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7&
         ,/' NCOPY  = ',i7/)

         allocate(ix_1(nnb), stat = allochk)
         if (allochk .gt. 0 ) then
            write(6,*)'rdparm:could not allocate ax_1, exiting'
            close(nf)
            stop
         endif
         if (natom .ne. nat) then
            print *,'natom in prmtop != nat in input, exiting'
            close (nf)
            stop
         endif
!  type = 'ATOM_NAME'
         read(nf,*)
         read(nf,*)
         read(nf,102) (at_name(i),i = 1,natom)
          print *, 'ATOM_NAME'
         print 102, at_name(1),at_name(natom)
 102     format(20a4)

!  type = 'CHARGE'
          read(nf,*)
          read(nf,*)
          read(nf,103) (at_chg(i),i = 1,natom)
          print *, 'CHARGE',natom,nat
          print 103, at_chg(1),at_chg(natom)
 103      format(5E16.8)

!  type = 'ATOMIC_NUMBER'
          read(nf,*)
          read(nf,*)
          print *, 'ATOMIC_NUMBER'
          read(nf,104) (at_type_index(i),i = 1,natom)
          print 104, at_type_index(1),at_type_index(natom)

!  type = 'MASS'
          read(nf,*)
          read(nf,*)
          read(nf,103) (at_mass(i),i = 1,natom)
          print *, 'MASS'
          print 103, at_mass(1),at_mass(natom)

!  type = 'ATOM_TYPE_INDEX'
          read(nf,*)
          read(nf,*)
          print *, 'ATOM_TYPE_INDEX'
          read(nf,104) (at_type_index(i),i = 1,natom)
          print 104, at_type_index(1),at_type_index(natom)
 104      format(10I8)

!  type = 'NUMBER_EXCLUDED_ATOMS'
          read(nf,*)
          read(nf,*)
          read(nf,104) (ix_1(i),i = 1,natom)
          print *, 'NUMBER_EXCLUDED_ATOMS'
          print 104, ix_1(1),ix_1(natom)

!  type = 'NONBONDED_PARM_INDEX'
           nttyp = ntypes*(ntypes+1)/2
           ntype = ntypes*ntypes

          read(nf,*)
          read(nf,*)
          read(nf,104) (ix_1(i),i = 1,ntype)
          print *, 'NONBONDED_PARM_INDEX'
          print 104, ntype,ix_1(1),ix_1(ntype),ix_1(ntype-1)
          print *, 'ntype = ',ntype

!  type = 'RESIDUE_LABEL'
          read(nf,*)
          read(nf,*)
          print *, 'RESIDUE_LABEL'
          print *, 'nres = ',nres
          read(nf,102) (res_name(i),i=1,nres)
          !print *,nres,res_name(1),res_name(15000)
          print 102, res_name(1),res_name(nres)

!  type = 'RESIDUE_POINTER'
          read(nf,*)
          read(nf,*)
          print *, 'RESIDUE_POINTER'
          print *, 'nres = ',nres
          read(nf,104) (i_res_point(i),i=1,nres)
          print 104, nres,i_res_point(1),i_res_point(nres),i_res_point(nres-1)
          print *, 'nres = ',nres

!    ----- READ THE PARAMETERS -----

!  type = 'BOND_FORCE_CONSTANT'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),i=1,numbnd)
         print *, 'BOND_FORCE_CONSTANT'
         print *, 'test1'
         print 103, ax_1(1),ax_1(numbnd)

!  type = 'BOND_EQUIL_VALUE'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),   i = 1,numbnd)
          print *, 'BOND_EQUIL_VALUE'
          print *, 'test2'
          print 103, ax_1(1),ax_1(numbnd)

!  type = 'ANGLE_FORCE_CONSTANT'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),    i = 1,numang)
          print *, 'ANGLE_FORCE_CONSTANT'
          print *, 'test3'
          print 103, ax_1(1),ax_1(numang)

!  type = 'ANGLE_EQUIL_VALUE'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),   i = 1,numang)
          print *, 'ANGLE_EQUIL_VALUE'
          print *, 'test4'
          print 103, ax_1(1),ax_1(numang)

!  type = 'DIHEDRAL_FORCE_CONSTANT'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),    i = 1,nptra)
          print *, 'DIHEDRAL_FORCE_CONSTANT'
          print *, 'test5'
          print 103, ax_1(1),ax_1(nptra )

!  type = 'DIHEDRAL_PERIODICITY'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),    i = 1,nptra)
          print *, 'DIHEDRAL_PERIODICITY'
          print *, 'test6'
          print 103, ax_1(1),ax_1(nptra )

!  type = 'DIHEDRAL_PHASE'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i), i = 1,nptra)
          print *, 'DIHEDRAL_PHASE'
          print *, 'test6'
          print 103, ax_1(1),ax_1(nptra )

!  type = 'SCEE_SCALE_FACTOR'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i), i = 1,nptra)
          print *, 'SCEE_SCALE_FACTOR'
          print *, 'test6a'
          print 103, ax_1(1),ax_1(nptra )

!  type = 'SCNB_SCALE_FACTOR'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i), i = 1,nptra)
          print *, 'SCNB_SCALE_FACTOR'
          print *, 'test6b'
          print 103, ax_1(1),ax_1(nptra )

!  type = 'SOLTY'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i), i = 1,natyp)
          print *, 'SOLTY'
          print *, 'test7'
          print 103, ax_1(1),ax_1(natyp)

!  type = 'LENNARD_JONES_ACOEF'
         read(nf,*)
         read(nf,*)
         print *,'LJ nttyp = ',nttyp
         read(nf,103) (xlj_12(i),   i = 1,nttyp)
          print *, 'LENNARD_JONES_ACOEF'
          print *, 'test8'
          print 103, xlj_12(1),xlj_12(nttyp)

!  type = 'LENNARD_JONES_BCOEF'
         read(nf,*)
         read(nf,*)
         read(nf,103) (xlj_6(i),   i = 1,nttyp)
          print *, 'LENNARD_JONES_BCOEF'
          print *, 'test9'
          print 103, xlj_6(1),xlj_6(nttyp)

!     ----- READ THE BONDING INFORMATIONS -----

!  type = 'BONDS_INC_HYDROGEN'
         read(nf,*)
         read(nf,*)
         read(nf,104) (ix_1(i),ix_2(i),ix_3(i),i=1,nbonh)
          print *, 'BONDS_INC_HYDROGEN'
          print *, 'test10'
          print 104, ix_1(1),ix_1(nbonh)

!  type = 'BONDS_WITHOUT_HYDROGEN'
          read(nf,*)
          read(nf,*)
          read(nf,104)(ix_1(i),ix_2(i),ix_3(i),i=1,nbona)
          print *, 'BONDS_WITHOUT_HYDROGEN'
          print *, 'test11'
          print 104, ix_1(1),ix_1(nbona)

!  type = 'ANGLES_INC_HYDROGEN'
          read(nf,*)
          read(nf,*)
          read(nf,104) (ix_1(i),ix_2(i),ix_3(i),ix_4(i),i=1,ntheth)
          print *, 'ANGLES_INC_HYDROGEN'
          print *, 'test12'
          print 104, ix_1(1),ix_1(ntheth)

!  type = 'ANGLES_WITHOUT_HYDROGEN'
         read(nf,*)
         read(nf,*)
         read(nf,104) (ix_1(i),ix_2(i),ix_3(i),ix_4(i),i=1,ntheta)
          print *, 'ANGLES_WITHOUT_HYDROGEN'
          print *, 'test13'
          print 104, ix_1(1),ix_1(ntheta)

!  type = 'DIHEDRALS_INC_HYDROGEN'
         read(nf,*)
         read(nf,*)
         read(nf,104) &
          (ix_1(i),ix_2(i),ix_3(i),ix_4(i),ix_5(i),i=1,nphih)
          print *, 'DIHEDRALS_INC_HYDROGEN'
          print *, 'test14'
          print 104, ix_1(1),ix_1(nphih)

!  type = 'DIHEDRALS_WITHOUT_HYDROGEN'
         read(nf,*)
         read(nf,*)
         read(nf,104)&
           (ix_1(i),ix_2(i),ix_3(i),ix_4(i),ix_5(i),i=1,nphia)
          print *, 'DIHEDRALS_WITHOUT_HYDROGEN'
          print *, 'test15', nnb
          print 104, ix_1(1),ix_1(nphia)

!  type = 'EXCLUDED_ATOMS_LIST'
         read(nf,*)
         read(nf,*)
         read(nf,104) (ix_1(i),i=1,nnb)
          print *, 'EXCLUDED_ATOMS_LIST'
          print *, 'test16'
          print 104, ix_1(1),ix_1(nnb)

!     ----- READ THE H-BOND PARAMETERS -----

!  type = 'HBOND_ACOEF'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),i=1,nphb)
          print *, 'HBOND_ACOEF'
          print *, 'test17'
          print 103, ax_1(1),ax_1(nphb)

!  type = 'HBOND_BCOEF'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),i=1,nphb)
          print *, 'HBOND_BCOEF'
          print *, 'test18'
          print 103, ax_1(1),ax_1(nphb)

!  type = 'HBCUT'
         read(nf,*)
         read(nf,*)
         read(nf,103) (ax_1(i),i=1,nphb)
          print *, 'HBCUT'
          print *, 'test19'
          print 103, ax_1(1),ax_1(nphb)

!     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----

!  type = 'AMBER_ATOM_TYPE'
         read(nf,*)
         read(nf,*)
         read(nf,102) (am_at_type(i),i=1,natom)
          print *, 'AMBER_ATOM_TYPE'
          print *, 'test20'
          print 102, am_at_type(1),am_at_type(natom)

!  type = 'TREE_CHAIN_CLASSIFICATION'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,102) (tree_ch_cl(i),i=1,natom)

!  type = 'JOIN_ARRAY'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,104) (ix_1(i),i=1,natom)

!  type = 'IROTAT'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,104) (ix_1(i), i=1,natom)

!    ----- READ THE BOUNDARY CONDITION STUFF -----

!     type = 'SOLVENT_POINTERS'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,104) iptres,nspm,nspsol

!     type = 'ATOMS_PER_MOLECULE'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,104) (ix_1(i),i=1,nspm)

!     type = 'BOX_DIMENSIONS'
!C       read(nf,*)
!C       read(nf,*)
!C       read(nf,103) oldbeta,duma,dumb,dumc

        return
        end
