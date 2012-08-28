! Copyright 2012 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##################################################################################################
!* $Id$
!##################################################################################################
! Material subroutine for BVP solution using spectral method
!
! Run 'DAMASK_spectral.exe --help' to get usage hints
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts,
!            D.D. Tjahjanto,
!            C. Kords,
!            M. Diehl,
!            R. Lebensohn
!
! MPI fuer Eisenforschung, Duesseldorf

#include "spectral_quit.f90"

program DAMASK_spectral
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   geometryFile, &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   appendToOutFile
   
 use prec, only: &
   pInt, &
   pReal, &
   DAMASK_NaN
   
 use IO, only: &
   IO_isBlank, &
   IO_open_file, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_error, &
   IO_lc, &
   IO_read_jobBinaryFile, &
   IO_write_jobBinaryFile
   
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_levelBasic, &
   debug_spectralDivergence, &
   debug_spectralRestart, &
   debug_spectralFFTW, &
   debug_reset, &
   debug_info
   
 use math
 
 use mesh,  only : &
   homog, &
   res, &
   res1_red, &
   mesh_NcpElems, &
   wgt, &
   geomdim, &
   virt_dim, &
   deformed_FFT
   
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_initAll
   
 use FEsolving, only: &
   restartWrite, &
   restartInc
   
 use numerics, only: &
   err_div_tol, &
   err_stress_tolrel, &
   err_stress_tolabs, &
   rotation_tol, &
   itmax,&
   itmin, &
   memory_efficient, &                            
   DAMASK_NumThreadsInt, &
   fftw_planner_flag, &
   fftw_timelimit
   
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results

 implicit none
!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: & 
   temp_valueVector                                                                                 !> temporarily from loadcase file when reading in tensors
 logical,     dimension(9) :: &
   temp_maskVector                                                                                  !> temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: maxNchunksLoadcase = (1_pInt + 9_pInt)*3_pInt +&                      ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt +&                      ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                    1_pInt, &                                       ! dropguessing
                              maxNchunksGeom     = 7_pInt, &                                        ! 4 identifiers, 3 values
                              myUnit             = 234_pInt
 integer(pInt), dimension(1_pInt + maxNchunksLoadcase*2_pInt) :: positions                          ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_l    = 0_pInt, &
   N_t    = 0_pInt, &
   N_n    = 0_pInt, &
   N_Fdot = 0_pInt
 
 character(len=1024) :: &
   line

 type bc_type
   real(pReal), dimension (3,3) :: deformation            = 0.0_pReal, &                            ! applied velocity gradient or time derivative of deformation gradient
                                   stress                 = 0.0_pReal, &                            ! stress BC (if applicable)
                                   rotation               = math_I3                                 ! rotation of BC (if applicable)
   real(pReal) ::                  time                   = 0.0_pReal, &                            ! length of increment
                                   temperature            = 300.0_pReal                             ! isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               ! number of increments
                                   outputfrequency        = 1_pInt, &                               ! frequency of result writes
                                   restartfrequency       = 0_pInt, &                               ! frequency of restart writes
                                   logscale               = 0_pInt                                  ! linear/logaritmic time inc flag
   logical ::                      followFormerTrajectory = .true., &                               ! follow trajectory of former loadcase
                                   velGradApplied         = .false.                                 ! decide wether velocity gradient or fdot is given 
   logical, dimension(3,3) ::      maskDeformation        = .false., &                              ! mask of deformation boundary conditions
                                   maskStress             = .false.                                 ! mask of stress boundary conditions
   logical, dimension(9) ::        maskStressVector       = .false.                                 ! linear mask of boundary conditions    
 end type
 
 type(bc_type), allocatable, dimension(:) ::  bc

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) :: &
   P_av, &
   F_aim = math_I3, &
   F_aim_lastInc = math_I3, &
   Favg = 0.0_pReal, &
   mask_stress, &
   mask_defgrad, &
   deltaF_aim, &
   F_aim_lab, &
   F_aim_lab_lastIter, &
   P_av_lab
 
 real(pReal), dimension(3,3,3,3) :: &
   dPdF, &
   C_ref = 0.0_pReal, &
   C = 0.0_pReal, &
   S_lastInc, &
   C_lastInc                                                                                        ! stiffness and compliance
 
 real(pReal), dimension(6) ::             sigma                                                     ! cauchy stress
 real(pReal), dimension(6,6) ::           dsde
 real(pReal), dimension(9,9) ::           temp99_Real                                               ! compliance and stiffness in matrix notation 
 real(pReal), dimension(:,:), allocatable ::  s_reduced, c_reduced                                  ! reduced compliance and stiffness (only for stress BC)
 integer(pInt) ::                         size_reduced = 0_pInt                                     ! number of stress BCs

!--------------------------------------------------------------------------------------------------
! pointwise data 
 type(C_PTR) :: tensorField                                                                         ! field in real an fourier space
 real(pReal),    dimension(:,:,:,:,:), pointer :: P_real, deltaF_real                               ! field in real space (pointer)

 complex(pReal), dimension(:,:,:,:,:), pointer :: P_fourier,deltaF_fourier                          ! field in fourier space (pointer)

 real(pReal),    dimension(:,:,:,:,:), allocatable ::  F, F_lastInc
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 type(C_PTR) ::  plan_stress, plan_correction                                                       ! plans for fftw
 real(pReal), dimension(3,3) ::                         xiDyad                                      ! product of wave vectors
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                   ! gamma operator (field) for spectral method
 real(pReal), dimension(:,:,:,:), allocatable ::        xi                                          ! wave vector field for divergence and for gamma operator
 integer(pInt), dimension(3) ::                         k_s

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval 
 real(pReal) :: guessmode, err_div, err_stress, err_stress_tol              
 real(pReal), dimension(3,3), parameter ::  ones = 1.0_pReal, zeroes = 0.0_pReal
 complex(pReal), dimension(3)   ::          temp3_Complex
 complex(pReal), dimension(3,3) ::          temp33_Complex
 real(pReal),    dimension(3,3) ::          temp33_Real
 integer(pInt) :: i, j, k, l, m, n, p, errorID
 integer(pInt) :: N_Loadcases, loadcase = 0_pInt, inc, iter, ielem, CPFEM_mode=1_pInt, &
                  ierr, totalIncsCounter = 0_pInt,&
                  notConvergedCounter = 0_pInt, convergedCounter = 0_pInt
 logical :: errmatinv
 real(pReal) :: defgradDet
 character(len=6)   :: loadcase_string

!--------------------------------------------------------------------------------------------------
!variables controlling debugging
 logical :: debugGeneral, debugDivergence, debugRestart, debugFFTW

!--------------------------------------------------------------------------------------------------
!variables for additional output due to general debugging
 real(pReal) :: defgradDetMax, defgradDetMin, maxCorrectionSym, maxCorrectionSkew

!--------------------------------------------------------------------------------------------------
! variables for additional output of divergence calculations
 type(C_PTR) :: divergence, plan_divergence
 real(pReal),    dimension(:,:,:,:), pointer :: divergence_real
 complex(pReal), dimension(:,:,:,:), pointer :: divergence_fourier
 real(pReal),    dimension(:,:,:,:), allocatable :: divergence_post
 real(pReal) :: pstress_av_L2, err_div_RMS, err_real_div_RMS, err_post_div_RMS,&
                err_div_max, err_real_div_max

!--------------------------------------------------------------------------------------------------
! variables for debugging fft using a scalar field
 type(C_PTR) :: scalarField_realC, scalarField_fourierC,&
                plan_scalarField_forth, plan_scalarField_back
 complex(pReal), dimension(:,:,:), pointer :: scalarField_real
 complex(pReal), dimension(:,:,:), pointer :: scalarField_fourier
 integer(pInt) :: row, column

!##################################################################################################
! reading of information from load case file and geometry file
!##################################################################################################
#ifdef PETSC
 integer :: ierr_psc
 call PetscInitialize(PETSC_NULL_CHARACTER, ierr_psc)
#endif

!-------------------------------------------------------------------------------------------------- 
! initialization of all related DAMASK modules (e.g. mesh.f90 reads in geometry)
 call CPFEM_initAll(temperature = 300.0_pReal, element = 1_pInt, IP= 1_pInt)

 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectral init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ' Working Directory:    ',trim(getSolverWorkingDirectoryName())
 write(6,'(a)') ' Solver Job Name:      ',trim(getSolverJobName())
 write(6,'(a)') ''

!--------------------------------------------------------------------------------------------------
! reading the load case file and allocate data structure containing load cases
 call IO_open_file(myUnit,trim(loadCaseFile))
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt, maxNchunksLoadcase, 1_pInt                                                        ! reading compulsory parameters for loadcase
       select case (IO_lc(IO_stringValue(line,positions,i)))
            case('l','velocitygrad','velgrad','velocitygradient')
                 N_l = N_l + 1_pInt
            case('fdot','dotf')
                 N_Fdot = N_Fdot + 1_pInt
            case('t','time','delta')
                 N_t = N_t + 1_pInt
            case('n','incs','increments','steps','logincs','logincrements','logsteps')
                 N_n = N_n + 1_pInt
        end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

100 N_Loadcases = N_n
 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                     ! sanity check
   call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                               ! error message for incomplete loadcase
 allocate (bc(N_Loadcases))

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)

 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   loadcase = loadcase + 1_pInt
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do j = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,positions,j)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient')                          ! assign values for the deformation BC matrix
         bc(loadcase)%velGradApplied = &
                     (IO_lc(IO_stringValue(line,positions,j)) == 'l'.or. &                          ! in case of given L, set flag to true
                      IO_lc(IO_stringValue(line,positions,j)) == 'velocitygrad'.or.&
                      IO_lc(IO_stringValue(line,positions,j)) == 'velgrad'.or.&
                      IO_lc(IO_stringValue(line,positions,j)) == 'velocitygradient')
         temp_valueVector = 0.0_pReal
         temp_maskVector = .false.
         forall (k = 1_pInt:9_pInt) temp_maskVector(k) = IO_stringValue(line,positions,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (temp_maskVector(k)) temp_valueVector(k) = IO_floatValue(line,positions,j+k)
         enddo
         bc(loadcase)%maskDeformation = transpose(reshape(temp_maskVector,[ 3,3]))
         bc(loadcase)%deformation = math_plain9to33(temp_valueVector)
       case('p','pk1','piolakirchhoff','stress')
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) bc(loadcase)%maskStressVector(k) =&
                                                          IO_stringValue(line,positions,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (bc(loadcase)%maskStressVector(k)) temp_valueVector(k) =&
                                                          IO_floatValue(line,positions,j+k)         ! assign values for the bc(loadcase)%stress matrix
         enddo
         bc(loadcase)%maskStress = transpose(reshape(bc(loadcase)%maskStressVector,[ 3,3]))
         bc(loadcase)%stress = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         bc(loadcase)%time = IO_floatValue(line,positions,j+1_pInt)
       case('temp','temperature')                                                                   ! starting temperature
         bc(loadcase)%temperature = IO_floatValue(line,positions,j+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         bc(loadcase)%incs = IO_intValue(line,positions,j+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         bc(loadcase)%incs = IO_intValue(line,positions,j+1_pInt)
         bc(loadcase)%logscale = 1_pInt
       case('f','freq','frequency','outputfreq')                                                    ! frequency of result writings
         bc(loadcase)%outputfrequency = IO_intValue(line,positions,j+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         bc(loadcase)%restartfrequency = max(0_pInt,IO_intValue(line,positions,j+1_pInt))                
       case('guessreset','dropguessing')
         bc(loadcase)%followFormerTrajectory = .false.                                              ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of loadcase given in euler angles
         p = 1_pInt                                                                                 ! assuming values given in radians
         l = 1_pInt                                                                                 ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,positions,j+1_pInt)))
           case('deg','degree')                                                                     ! for conversion from degree to radian
           case('rad','radian') 
             p = 0_pInt
           case default           
             l = 0_pInt                                                                             ! immediately reading in angles, assuming radians
         end select
         forall(k = 1_pInt:3_pInt) temp33_Real(k,1) = IO_floatValue(line,positions,j+l+k) 
         if (p==1_pInt) temp33_Real = temp33_Real * inRad
         bc(loadcase)%rotation = math_EulerToR(temp33_Real(:,1))
       case('rotation','rot')                                                                       ! assign values for the rotation of loadcase matrix
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) temp_valueVector(k) = IO_floatValue(line,positions,j+k)
         bc(loadcase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
101 close(myUnit)

!--------------------------------------------------------------------------------------------------
! output of geometry
 write(6,'(a)')          ''
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'DAMASK spectral:'
 write(6,'(a)')          'The spectral method boundary value problem solver for'
 write(6,'(a)')          'the Duesseldorf Advanced Material Simulation Kit'
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'geometry file:        ',trim(geometryFile)
 write(6,'(a)')          '============================================================='
 write(6,'(a,3(i12  ))') 'resolution a b c:', res
 write(6,'(a,3(f12.5))') 'dimension  x y z:', geomdim
 write(6,'(a,i5)')       'homogenization:       ',homog
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'loadcase file:        ',trim(loadCaseFile)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 bc(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first loadcase
 errorID = 0_pInt
 do loadcase = 1_pInt, N_Loadcases
   write (loadcase_string, '(i6)' ) loadcase

   write(6,'(a)') '============================================================='
   write(6,'(a,i6)') 'loadcase:            ', loadcase

   if (.not. bc(loadcase)%followFormerTrajectory) write(6,'(a)') 'drop guessing along trajectory'
   if (bc(loadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .true.) .and. &
           any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .false.)) errorID = 832_pInt               ! each row should be either fully or not at all defined
     enddo
     write(6,'(a)')'velocity gradient:'
   else
     write(6,'(a)')'deformation gradient rate:'
   endif
   write (6,'(3(3(f12.7,1x)/))',advance='no') merge(math_transpose33(bc(loadcase)%deformation),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskDeformation))
   write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' stress / GPa:',&
        1e-9_pReal*merge(math_transpose33(bc(loadcase)%stress),&
                         reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskStress))
   if (any(bc(loadcase)%rotation /= math_I3)) &
     write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' rotation of loadframe:',&
                                                          math_transpose33(bc(loadcase)%rotation)
   write(6,'(a,f12.6)') 'temperature:', bc(loadcase)%temperature
   write(6,'(a,f12.6)') 'time:       ', bc(loadcase)%time
   write(6,'(a,i5)')    'increments: ', bc(loadcase)%incs
   write(6,'(a,i5)')    'output  frequency:  ', bc(loadcase)%outputfrequency
   write(6,'(a,i5)')    'restart frequency:  ', bc(loadcase)%restartfrequency

   if (any(bc(loadcase)%maskStress .eqv. bc(loadcase)%maskDeformation)) errorID = 831_pInt          ! exclusive or masking only
   if (any(bc(loadcase)%maskStress .and. transpose(bc(loadcase)%maskStress) .and. &
     reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
                                               errorID = 838_pInt                                   ! no rotation is allowed by stress BC
   if (any(abs(math_mul33x33(bc(loadcase)%rotation,math_transpose33(bc(loadcase)%rotation))&
                                      -math_I3) > reshape(spread(rotation_tol,1,9),[ 3,3]))&
                    .or. abs(math_det33(bc(loadcase)%rotation)) > 1.0_pReal + rotation_tol)&
                                               errorID = 846_pInt                                   ! given rotation matrix contains strain
   if (bc(loadcase)%time < 0.0_pReal)          errorID = 834_pInt                                   ! negative time increment
   if (bc(loadcase)%incs < 1_pInt)             errorID = 835_pInt                                   ! non-positive incs count
   if (bc(loadcase)%outputfrequency < 1_pInt)  errorID = 836_pInt                                   ! non-positive result frequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo

!--------------------------------------------------------------------------------------------------
! debugging parameters
 debugGeneral    = iand(debug_level(debug_spectral),debug_levelBasic)         /= 0
 debugDivergence = iand(debug_level(debug_spectral),debug_spectralDivergence) /= 0
 debugRestart    = iand(debug_level(debug_spectral),debug_spectralRestart)    /= 0
 debugFFTW       = iand(debug_level(debug_spectral),debug_spectralFFTW)       /= 0
 
!##################################################################################################
! initialization 
!##################################################################################################

 allocate (F          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (F_lastInc  (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (xi         (3,res1_red,res(2),res(3)),      source = 0.0_pReal)
 allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
 allocate (temperature(  res(1),  res(2),res(3)),      source = bc(1)%temperature)                  ! start out isothermally
 tensorField = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))                      ! allocate continous data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, P_real,         [ res(1)+2_pInt,res(2),res(3),3,3])                  ! place a pointer for a real representation on tensorField
 call c_f_pointer(tensorField, deltaF_real,    [ res(1)+2_pInt,res(2),res(3),3,3])                  ! place a pointer for a real representation on tensorField
 call c_f_pointer(tensorField, P_fourier,      [ res1_red,     res(2),res(3),3,3])                  ! place a pointer for a complex representation on tensorField
 call c_f_pointer(tensorField, deltaF_fourier, [ res1_red,     res(2),res(3),3,3])                  ! place a pointer for a complex representation on tensorField

!--------------------------------------------------------------------------------------------------
! general initialization of fftw (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=808_pInt)                         ! check for correct precision in C
!$ if(DAMASK_NumThreadsInt > 0_pInt) then
!$   ierr = fftw_init_threads()
!$   if (ierr == 0_pInt) call IO_error(error_ID = 809_pInt)
!$   call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
!$ endif
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans
 plan_stress =    fftw_plan_many_dft_r2c(3,[ res(3),res(2) ,res(1)],9,&                             ! dimensions , length in each dimension in reversed order
                                    P_real,[ res(3),res(2) ,res(1)+2_pInt],&                        ! input data , physical length in each dimension in reversed order
                                         1,  res(3)*res(2)*(res(1)+2_pInt),&                        ! striding   , product of physical lenght in the 3 dimensions
                                    P_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 plan_correction =fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],9,&
                            deltaF_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,&
                               deltaF_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                         1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)

!--------------------------------------------------------------------------------------------------
! depending on (debug) options, allocate more memory and create additional plans 
 if (debugDivergence) then
   divergence = fftw_alloc_complex(int(res1_red*res(2)*res(3)*3_pInt,C_SIZE_T))
   call c_f_pointer(divergence, divergence_real,    [ res(1)+2_pInt,res(2),res(3),3])
   call c_f_pointer(divergence, divergence_fourier, [ res1_red,     res(2),res(3),3])
   allocate (divergence_post(res(1),res(2),res(3),3));  divergence_post = 0.0_pReal
   plan_divergence = fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],3,&
                           divergence_fourier,[ res(3),res(2) ,res1_red],&
                                            1,  res(3)*res(2)* res1_red,&
                              divergence_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                            1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 endif

 if (debugFFTW) then
   scalarField_realC    = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))                    ! do not do an inplace transform  
   scalarField_fourierC = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))
   call c_f_pointer(scalarField_realC,    scalarField_real,    [res(1),res(2),res(3)])
   call c_f_pointer(scalarField_fourierC, scalarField_fourier, [res(1),res(2),res(3)])
   plan_scalarField_forth = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_real,scalarField_fourier,-1,fftw_planner_flag)
   plan_scalarField_back  = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_fourier,scalarField_real,+1,fftw_planner_flag)
 endif 

 if (debugGeneral) write(6,'(a)') 'FFTW initialized'
  
!--------------------------------------------------------------------------------------------------
! init fields                 
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     F(i,j,k,1:3,1:3) = math_I3
     F_lastInc(i,j,k,1:3,1:3) = math_I3
     coordinates(i,j,k,1:3) = geomdim/real(res,pReal)*real([i,j,k],pReal) &
                            - geomdim/real(2_pInt*res,pReal)
   enddo; enddo; enddo
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (debugRestart) write(6,'(a,i6,a)') 'Reading values of increment ',&
                                             restartInc - 1_pInt,' from file' 
   call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',&
                                                trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   call IO_read_jobBinaryFile(777,'F_aim',trim(getSolverJobName()),size(F_aim))
   read (777,rec=1) F_aim
   close (777)
   call IO_read_jobBinaryFile(777,'F_aim_lastInc',trim(getSolverJobName()),size(F_aim_lastInc))
   read (777,rec=1) F_aim_lastInc
   close (777)
   coordinates = 0.0 ! change it later!!!
   CPFEM_mode = 2_pInt
   if (debugRestart) write(6,'(a)') 'Data read in'
 endif
 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt 
   call CPFEM_general(3_pInt,coordinates(i,j,k,1:3),F(i,j,k,1:3,1:3),F(i,j,k,1:3,1:3),temperature(i,j,k),&                         
                        0.0_pReal,ielem,1_pInt,sigma,dsde,P_real(i,j,k,1:3,1:3),dPdF)
 enddo; enddo; enddo
 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt 
   call CPFEM_general(2_pInt,coordinates(i,j,k,1:3),F(i,j,k,1:3,1:3),F(i,j,k,1:3,1:3),temperature(i,j,k),&
                        0.0_pReal,ielem,1_pInt,sigma,dsde,P_real(i,j,k,1:3,1:3),dPdF)
   C = C + dPdF 
 enddo; enddo; enddo
 if (debugGeneral) write(6,'(a)') 'First call to CPFEM finished'
 C = C * wgt
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)

 do k = 1_pInt, res(3)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1_pInt, res1_red
           k_s(1) = i - 1_pInt
           xi(1:3,i,j,k) = real(k_s, pReal)/virt_dim
 enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! calculate the gamma operator
 if (restartInc == 1_pInt) then
   C_ref = C 
   call IO_write_jobBinaryFile(777,'C_ref',size(C_ref))
   write (777,rec=1) C_ref
   close(777)
 elseif (restartInc > 1_pInt) then
   call IO_read_jobBinaryFile(777,'C_ref',trim(getSolverJobName()),size(C_ref))
   read (777,rec=1) C_ref
   close (777)
 endif
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(res1_red ,res(2),res(3),3,3,3,3), source =0.0_pReal)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, p=1_pInt:3_pInt)&
         gamma_hat(i,j,k, l,m,n,p) =  temp33_Real(l,n)*xiDyad(m,p)
     endif  
   enddo; enddo; enddo
   gamma_hat(1,1,1, 1:3,1:3,1:3,1:3) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif

!--------------------------------------------------------------------------------------------------
! write header of output file
 if (appendToOutFile) then
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED', position='APPEND', status='OLD')
   if (debugRestart) write(6,'(a)') 'Result File opened for appending'
 else
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED',status='REPLACE')
   write(538) 'load',       trim(loadCaseFile)
   write(538) 'workingdir', trim(getSolverWorkingDirectoryName())
   write(538) 'geometry',   trim(geometryFile)
   write(538) 'resolution', res
   write(538) 'dimension',  geomdim
   write(538) 'materialpoint_sizeResults', materialpoint_sizeResults
   write(538) 'loadcases',        N_Loadcases
   write(538) 'frequencies', bc(1:N_Loadcases)%outputfrequency                                      ! one entry per loadcase
   write(538) 'times', bc(1:N_Loadcases)%time                                                       ! one entry per loadcase
   write(538) 'logscales',  bc(1:N_Loadcases)%logscale         
   write(538) 'increments', bc(1:N_Loadcases)%incs                                                  ! one entry per loadcase
   write(538) 'startingIncrement', restartInc - 1_pInt                                              ! start with writing out the previous inc
   write(538) 'eoh'                                                                                 ! end of header
   write(538) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:mesh_NcpElems)              ! initial (non-deformed or read-in) results
   if (debugGeneral) write(6,'(a)') 'Header of result file written out'
 endif


!##################################################################################################
! Loop over loadcases defined in the loadcase file
!##################################################################################################
 do loadcase = 1_pInt,  N_Loadcases
   time0 = time
   if (bc(loadcase)%followFormerTrajectory) then 
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                                          ! change of load case, homogeneous guess for the first inc
   endif

!--------------------------------------------------------------------------------------------------
! arrays for mixed boundary conditions
   mask_defgrad = merge(ones,zeroes,bc(loadcase)%maskDeformation)                                   
   mask_stress  = merge(ones,zeroes,bc(loadcase)%maskStress)
   size_reduced = int(count(bc(loadcase)%maskStressVector), pInt)
   allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)

!##################################################################################################
! loop oper incs defined in input file for current loadcase
!##################################################################################################
   do inc = 1_pInt,  bc(loadcase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 
     
!--------------------------------------------------------------------------------------------------
! forwarding time
     timeinc_old = timeinc
     if (bc(loadcase)%logscale == 0_pInt) then                                                      ! linear scale
       timeinc = bc(loadcase)%time/bc(loadcase)%incs                                                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (loadcase == 1_pInt) then                                                                 ! 1st loadcase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st loadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(    1_pInt-bc(1)%incs ,pReal))                     ! assume 1st inc is equal to 2nd 
         else                                                                                       ! not-1st inc of 1st loadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(inc-1_pInt-bc(1)%incs ,pReal))
         endif
       else                                                                                         ! not-1st loadcase of logarithmic scale
           timeinc = time0 *( (1.0_pReal + bc(loadcase)%time/time0 )**(real(          inc,pReal)/&
                                                                  real(bc(loadcase)%incs ,pReal))&
                             -(1.0_pReal + bc(loadcase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                                   real(bc(loadcase)%incs ,pReal)) )
       endif
     endif
     time = time + timeinc

     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 
       if (bc(loadcase)%velGradApplied) then                                                        ! calculate deltaF_aim from given L and current F
         deltaF_aim = timeinc * mask_defgrad * math_mul33x33(bc(loadcase)%deformation, F_aim)
       else                                                                                         ! deltaF_aim = fDot *timeinc where applicable
         deltaF_aim = timeinc * mask_defgrad * bc(loadcase)%deformation
       endif

!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
       temp33_Real = F_aim                                            
       F_aim = F_aim &                                                                         
                  + guessmode * mask_stress * (F_aim - F_aim_lastInc)*timeinc/timeinc_old &      
                  + deltaF_aim
       F_aim_lastInc = temp33_Real

!--------------------------------------------------------------------------------------------------
! update local deformation gradient and coordinates
       deltaF_aim = math_rotate_backward33(deltaF_aim,bc(loadcase)%rotation)
       Favg = 0.0_pReal                                            
       do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
         temp33_Real = F(i,j,k,1:3,1:3)
         F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) &                                                      ! decide if guessing along former trajectory or apply homogeneous addon
                                + guessmode * (F(i,j,k,1:3,1:3) - F_lastInc(i,j,k,1:3,1:3))&        ! guessing... 
                                            *timeinc/timeinc_old &
                                + (1.0_pReal-guessmode) * deltaF_aim                                ! if not guessing, use prescribed average deformation where applicable
         F_lastInc(i,j,k,1:3,1:3) = temp33_Real 
         Favg = Favg + F(i,j,k,1:3,1:3)
       enddo; enddo; enddo
       deltaF_aim =  guessmode *(Favg*wgt -math_rotate_backward33(F_aim,bc(loadcase)%rotation))                                                   ! average correction in case of guessing to 
       do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
         F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) - deltaF_aim                                           ! correct in case avg of F is not F_aim
       enddo; enddo; enddo
       call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,bc(loadcase)%rotation),&          ! calculate current coordinates
                                                          1.0_pReal,F_lastInc,coordinates)

!--------------------------------------------------------------------------------------------------
! calculate reduced compliance
       if(size_reduced > 0_pInt) then                                                               ! calculate compliance in case stress BC is applied
         C_lastInc = math_rotate_forward3333(C,bc(loadcase)%rotation)                               ! calculate stiffness from former inc
         temp99_Real = math_Plain3333to99(C_lastInc)
         k = 0_pInt                                                                                 ! build reduced stiffness
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
               if(bc(loadcase)%maskStressVector(m)) then
                 j = j + 1_pInt
                 c_reduced(k,j) = temp99_Real(n,m)
         endif; enddo; endif; enddo
         call math_invert(size_reduced,c_reduced, s_reduced, errmatinv)                             ! invert reduced stiffness
         if(errmatinv) call IO_error(error_ID=400_pInt)
         temp99_Real = 0.0_pReal                                                                    ! build full compliance
         k = 0_pInt
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
             if(bc(loadcase)%maskStressVector(m)) then
                   j = j + 1_pInt
                   temp99_Real(n,m) = s_reduced(k,j)
         endif; enddo; endif; enddo
         S_lastInc = (math_Plain99to3333(temp99_Real))
       endif

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       
       write(6,'(a)') '##################################################################'
       write(6,'(A,I5.5,A,es12.5)') 'Increment ', totalIncsCounter, ' Time ',time
       
       iter = 0_pInt
       err_div = huge(err_div_tol)                                                                  ! go into loop 

!##################################################################################################
! convergence loop (looping over iterations)
!##################################################################################################
       do while((iter < itmax .and. (err_div > err_div_tol .or. err_stress > err_stress_tol))&
                 .or. iter < itmin)
         iter = iter + 1_pInt

!--------------------------------------------------------------------------------------------------
! report begin of new iteration
         write(6,'(a)') ''
         write(6,'(a)') '=================================================================='
         write(6,'(6(a,i6.6))') 'Loadcase ',loadcase,' Inc. ',inc,'/',bc(loadcase)%incs,&
                                                     ' @ Iter. ',itmin,' < ',iter,' < ',itmax
         write(6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient aim =',&
                                                             math_transpose33(F_aim)
         write(6,'(a)') ''
         write(6,'(a)') '... update stress field P(F) .....................................'
         if (restartWrite) write(6,'(a)') 'writing CP restart information for last increment'
         
         F_aim_lab_lastIter = math_rotate_backward33(F_aim,bc(loadcase)%rotation)
!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
         ielem = 0_pInt
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(3_pInt,&                                                              ! collect cycle
                              coordinates(i,j,k,1:3), F_lastInc(i,j,k,1:3,1:3),F(i,j,k,1:3,1:3), &
                              temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde,&
                              P_real(i,j,k,1:3,1:3),dPdF)
         enddo; enddo; enddo

         P_real = 0.0_pReal                                                                         ! needed because of the padding for FFTW
         C = 0.0_pReal
         ielem = 0_pInt 
         call debug_reset()
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(CPFEM_mode,&                                                          ! first element in first iteration retains CPFEM_mode 1, 
                              coordinates(i,j,k,1:3),F_lastInc(i,j,k,1:3,1:3), F(i,j,k,1:3,1:3), &  ! others get 2 (saves winding forward effort)
                              temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde, &
                              P_real(i,j,k,1:3,1:3),dPdF)
           CPFEM_mode = 2_pInt
           C = C + dPdF
         enddo; enddo; enddo
         call debug_info()
         restartWrite = .false.

! for test of regridding
        ! if( bc(loadcase)%restartFrequency > 0_pInt .and. &
         !    mod(inc-1,bc(loadcase)%restartFrequency) == 0_pInt .and. &
          !   restartInc/=inc) call quit(-1*(restartInc+1))                                          ! trigger exit to regrid

!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
         if (debugFFTW) then
           row =    (mod(totalIncsCounter+iter-2_pInt,9_pInt))/3_pInt + 1_pInt                      ! go through the elements of the tensors, controlled by totalIncsCounter and iter, starting at 1
           column = (mod(totalIncsCounter+iter-2_pInt,3_pInt))        + 1_pInt
           scalarField_real(1:res(1),1:res(2),1:res(3)) =&                                          ! store the selected component
                  cmplx(P_real(1:res(1),1:res(2),1:res(3),row,column),0.0_pReal,pReal)
         endif

!--------------------------------------------------------------------------------------------------
! call function to calculate divergence from math (for post processing) to check results
         if (debugDivergence) &
           divergence_post = math_divergenceFFT(virt_dim,P_real(1:res(1),1:res(2),1:res(3),1:3,1:3))      ! padding
              
!--------------------------------------------------------------------------------------------------
! doing the FT because it simplifies calculation of average stress in real space also
         call fftw_execute_dft_r2c(plan_stress,P_real,P_fourier)

         P_av_lab = real(P_fourier(1,1,1,1:3,1:3),pReal)*wgt
         P_av = math_rotate_forward33(P_av_lab,bc(loadcase)%rotation)
         write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'Piola-Kirchhoff stress / MPa =',&
                                                          math_transpose33(P_av)/1.e6_pReal

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
         if (debugFFTW) then
           call fftw_execute_dft(plan_scalarField_forth,scalarField_real,scalarField_fourier)
           write(6,'(a,i1,1x,i1)') 'checking FT results of compontent ', row, column
           write(6,'(a,2(es11.4,1x))')  'max FT relative error = ',&
             maxval( real((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-& 
                                     P_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                           scalarField_fourier(1:res1_red,1:res(2),1:res(3)))), &
             maxval(aimag((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-&
                                     P_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                           scalarField_fourier(1:res1_red,1:res(2),1:res(3))))
         endif

!--------------------------------------------------------------------------------------------------
! removing highest frequencies
         P_fourier  (  res1_red,1:res(2) ,             1:res(3)              ,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         P_fourier  (1:res1_red,  res(2)/2_pInt+1_pInt,1:res(3)              ,1:3,1:3)& 
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         if(res(3)>1_pInt) &
          P_fourier (1:res1_red,1:res(2),                res(3)/2_pInt+1_pInt,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! stress BC handling
         if(size_reduced > 0_pInt) then                                                             ! calculate stress BC if applied
           err_stress = maxval(abs(mask_stress * (P_av - bc(loadcase)%stress)))                     ! maximum deviaton (tensor norm not applicable)
           err_stress_tol = min(maxval(abs(P_av)) * err_stress_tolrel,err_stress_tolabs)            ! don't use any tensor norm for the relative criterion because the comparison should be coherent
           write(6,'(a)') '' 
           write(6,'(a)') '... correcting deformation gradient to fulfill BCs ...............'
           write(6,'(a,f6.2,a,es11.4,a)') 'error stress = ', err_stress/err_stress_tol, &
                                                                            ' (',err_stress,' Pa)'  
           F_aim = F_aim - math_mul3333xx33(S_lastInc, ((P_av - bc(loadcase)%stress)))              ! residual on given stress components
           write(6,'(a,1x,es11.4)')'determinant of new deformation = ',math_det33(F_aim)
         else
           err_stress_tol = +huge(1.0_pReal)
         endif
                                            
         F_aim_lab = math_rotate_backward33(F_aim,bc(loadcase)%rotation)                            ! boundary conditions from load frame into lab (Fourier) frame

!--------------------------------------------------------------------------------------------------
! actual spectral method         
         write(6,'(a)') ''
         write(6,'(a)') '... calculating equilibrium with spectral method .................'

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
         pstress_av_L2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(P_av_lab,&                    ! L_2 norm of average stress (http://mathworld.wolfram.com/SpectralNorm.html)
                                                     math_transpose33(P_av_lab)))))
         err_div_RMS = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2)
           do i = 2_pInt, res1_red -1_pInt                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
             err_div_RMS = err_div_RMS &
                   + 2.0_pReal*(sum (real(math_mul33x3_complex(P_fourier(i,j,k,1:3,1:3),&           ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                   xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&            ! --> sum squared L_2 norm of vector 
                               +sum(aimag(math_mul33x3_complex(P_fourier(i,j,k,1:3,1:3),& 
                                                                  xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
           enddo
           err_div_RMS = err_div_RMS &                                                              ! Those two layers (DC and Nyquist) do not have a conjugate complex counterpart
                         + sum( real(math_mul33x3_complex(P_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(P_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum( real(math_mul33x3_complex(P_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(P_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)
         enddo; enddo

         err_div_RMS = sqrt(err_div_RMS)*wgt                                                        ! RMS in real space calculated with Parsevals theorem from Fourier space

         if (err_div_RMS/pstress_av_L2 >  err_div &
                     .and.  err_stress <  err_stress_tol &
                     .and.        iter >= itmin ) then
           write(6,'(a)') 'Increasing divergence, stopping iterations'
           iter = itmax
         endif
         err_div = err_div_RMS/pstress_av_L2                                                        ! criterion to stop iterations

!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
         if (debugDivergence) then                                                                  ! calculate divergence again
           err_div_max = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
             temp3_Complex = math_mul33x3_complex(P_fourier(i,j,k,1:3,1:3)*wgt,&                    ! weighting P_fourier
                                                     xi(1:3,i,j,k))*TWOPIIMG
             err_div_max = max(err_div_max,sum(abs(temp3_Complex)**2.0_pReal))
             divergence_fourier(i,j,k,1:3) = temp3_Complex                                          ! need divergence NOT squared
           enddo; enddo; enddo
           
           call fftw_execute_dft_c2r(plan_divergence,divergence_fourier,divergence_real)            ! already weighted

           err_real_div_RMS = 0.0_pReal
           err_post_div_RMS = 0.0_pReal
           err_real_div_max = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             err_real_div_RMS = err_real_div_RMS +   sum(divergence_real(i,j,k,1:3)**2.0_pReal)     ! avg of squared L_2 norm of div(stress) in real space
             err_post_div_RMS = err_post_div_RMS +   sum(divergence_post(i,j,k,1:3)**2.0_pReal)     ! avg of squared L_2 norm of div(stress) in real space
             err_real_div_max = max(err_real_div_max,sum(divergence_real(i,j,k,1:3)**2.0_pReal))    ! max of squared L_2 norm of div(stress) in real space
           enddo; enddo; enddo

           err_real_div_RMS = sqrt(wgt*err_real_div_RMS)                                            ! RMS in real space
           err_post_div_RMS = sqrt(wgt*err_post_div_RMS)                                            ! RMS in real space
           err_real_div_max = sqrt(    err_real_div_max)                                            ! max in real space
           err_div_max      = sqrt(    err_div_max)                                                 ! max in Fourier space
           
           write(6,'(a,es11.4)')        'error divergence  FT  RMS = ',err_div_RMS
           write(6,'(a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
           write(6,'(a,es11.4)')        'error divergence post RMS = ',err_post_div_RMS
           write(6,'(a,es11.4)')        'error divergence  FT  max = ',err_div_max
           write(6,'(a,es11.4)')        'error divergence Real max = ',err_real_div_max
         endif
         write(6,'(a,f6.2,a,es11.4,a)') 'error divergence = ', err_div/err_div_tol,&
                                                           ' (',err_div,' N/m³)'

!--------------------------------------------------------------------------------------------------
! to the actual spectral method calculation (mechanical equilibrium)
         if(memory_efficient) then                                                                  ! memory saving version, on-the-fly calculation of gamma_hat
           
           do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res1_red
               if(any([i,j,k] /= 1_pInt)) then                                                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
                 forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                   xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
                 forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                   temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
                 temp33_Real = math_inv33(temp33_Real)
                 forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, p=1_pInt:3_pInt)&
                   gamma_hat(1,1,1, l,m,n,p) =  temp33_Real(l,n)*xiDyad(m,p)
                 forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                   temp33_Complex(l,m) = sum(gamma_hat(1,1,1, l,m, 1:3,1:3) *&
                                                                P_fourier(i,j,k,1:3,1:3))
                 deltaF_fourier(i,j,k,1:3,1:3) = temp33_Complex 
             endif             
           enddo; enddo; enddo
   
         else                                                                                       ! use precalculated gamma-operator
           
           do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
             forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt) &
               temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n, 1:3,1:3) *&
                                                                P_fourier(i,j,k,1:3,1:3))
             deltaF_fourier(i,j,k, 1:3,1:3) = temp33_Complex
           enddo; enddo; enddo

         endif
         deltaF_fourier(1,1,1,1:3,1:3) = cmplx((F_aim_lab_lastIter - F_aim_lab) &                   ! assign (negative) average deformation gradient change to zero frequency (real part)
                                                * real(mesh_NcpElems,pReal),0.0_pReal,pReal)              ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
         if (debugFFTW) then
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
              scalarField_fourier(i,j,k) = deltaF_fourier(i,j,k,row,column)
           enddo; enddo; enddo
           do i = 0_pInt, res(1)/2_pInt-2_pInt                                                      ! unpack fft data for conj complex symmetric part
            m = 1_pInt
            do k = 1_pInt, res(3)
              n = 1_pInt
              do j = 1_pInt, res(2)
                scalarField_fourier(res(1)-i,j,k) = conjg(scalarField_fourier(2+i,n,m))
                if(n == 1_pInt) n = res(2) + 1_pInt
                n = n-1_pInt
             enddo
             if(m == 1_pInt) m = res(3) + 1_pInt
             m = m -1_pInt
           enddo; enddo
         endif
!--------------------------------------------------------------------------------------------------
! doing the inverse FT
         call fftw_execute_dft_c2r(plan_correction,deltaF_fourier,deltaF_real)                      ! back transform of fluct deformation gradient

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
         if (debugFFTW) then
           write(6,'(a,i1,1x,i1)') 'checking iFT results of compontent ', row, column
           call fftw_execute_dft(plan_scalarField_back,scalarField_fourier,scalarField_real)
           write(6,'(a,es11.4)') 'max iFT relative error = ',&
               maxval((real(scalarField_real(1:res(1),1:res(2),1:res(3)))-&
                       deltaF_real(1:res(1),1:res(2),1:res(3),row,column))/&
                       real(scalarField_real(1:res(1),1:res(2),1:res(3))))
         endif

!--------------------------------------------------------------------------------------------------
! calculate some additional output
         if(debugGeneral) then
           maxCorrectionSkew = 0.0_pReal
           maxCorrectionSym  = 0.0_pReal
           temp33_Real = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             maxCorrectionSym  = max(maxCorrectionSym,&
                                     maxval(math_symmetric33(deltaF_real(i,j,k,1:3,1:3))))
             maxCorrectionSkew = max(maxCorrectionSkew,&
                                     maxval(math_skew33(deltaF_real(i,j,k,1:3,1:3))))
             temp33_Real = temp33_Real + deltaF_real(i,j,k,1:3,1:3)
           enddo; enddo; enddo
           write(6,'(a,1x,es11.4)') 'max symmetric correction of deformation =',&
                                         maxCorrectionSym*wgt
           write(6,'(a,1x,es11.4)') 'max skew      correction of deformation =',&
                                         maxCorrectionSkew*wgt
           write(6,'(a,1x,es11.4)') 'max sym/skew of avg correction =         ',&
                                         maxval(math_symmetric33(temp33_real))/&
                                         maxval(math_skew33(temp33_real))
         endif

!--------------------------------------------------------------------------------------------------
! updated deformation gradient
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           F(i,j,k,1:3,1:3) = F(i,j,k,1:3,1:3) - deltaF_real(i,j,k,1:3,1:3)*wgt                       ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization
         enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
         if(debugGeneral) then
           defgradDetMax = -huge(1.0_pReal)
           defgradDetMin = +huge(1.0_pReal)
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             defgradDet = math_det33(F(i,j,k,1:3,1:3))
             defgradDetMax = max(defgradDetMax,defgradDet)
             defgradDetMin = min(defgradDetMin,defgradDet) 
           enddo; enddo; enddo

           write(6,'(a,1x,es11.4)') 'max determinant of deformation =', defgradDetMax
           write(6,'(a,1x,es11.4)') 'min determinant of deformation =', defgradDetMin
         endif
         
       enddo    ! end looping when convergency is achieved 
           
       CPFEM_mode = 1_pInt                                                                          ! winding forward
       C = C * wgt
       write(6,'(a)') ''
       write(6,'(a)') '=================================================================='
       if(err_div > err_div_tol .or. err_stress > err_stress_tol) then
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       else
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' converged'
       endif

       if (mod(inc,bc(loadcase)%outputFrequency) == 0_pInt) then                                    ! at output frequency
         write(6,'(a)') ''
         write(6,'(a)') '... writing results to file ......................................'
         write(538)  materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:mesh_NcpElems)       ! write result to file
         flush(538)
       endif
       
       if( bc(loadcase)%restartFrequency > 0_pInt .and. &
                      mod(inc,bc(loadcase)%restartFrequency) == 0_pInt) then                        ! at frequency of writing restart information set restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
         restartInc=totalIncsCounter
         restartWrite = .true.
         write(6,'(a)') 'writing converged results for restart'
         call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F))                        ! writing deformation gradient field to file
         write (777,rec=1) F
         close (777)
         call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',size(F_lastInc))        ! writing F_lastInc field to file
         write (777,rec=1) F_lastInc
         close (777)
         call IO_write_jobBinaryFile(777,'F_aim',size(F_aim))
         write (777,rec=1) F_aim
         close(777)
         call IO_write_jobBinaryFile(777,'F_aim_lastInc',size(F_aim_lastInc))
         write (777,rec=1) F_aim_lastInc
         close(777)
       endif 
       
     endif ! end calculation/forwarding
     guessmode = 1.0_pReal                                                                        ! keep guessing along former trajectory during same loadcase

   enddo  ! end looping over incs in current loadcase
   deallocate(c_reduced)
   deallocate(s_reduced)
   enddo    ! end looping over loadcases
   write(6,'(a)') ''
   write(6,'(a)') '##################################################################'
   write(6,'(i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                     notConvergedCounter + convergedCounter, ' (', &
                                     real(convergedCounter, pReal)/&
                                     real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                     ' %) increments converged!'
 close(538)
 call fftw_destroy_plan(plan_stress); call fftw_destroy_plan(plan_correction)
 if (debugDivergence) call fftw_destroy_plan(plan_divergence)
 if (debugFFTW) then
   call fftw_destroy_plan(plan_scalarField_forth)
   call fftw_destroy_plan(plan_scalarField_back)
 endif
 if (notConvergedCounter > 0_pInt) call quit(3_pInt)
 call quit(0_pInt)
end program DAMASK_spectral
