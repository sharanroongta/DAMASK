!* $Id$
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) :: cycleCounter = 0_pInt, theInc = -1_pInt
 real(pReal)   :: theTime = 0.0_pReal, theDelta = 0.0_pReal
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false.,outdatedFFN1 = .false.,terminallyIll = .false.
 logical :: symmetricSolver = .false. 
 logical :: parallelExecution = .true. 
 logical :: restartWrite = .false.
 logical :: restartRead  = .false.
 logical :: lastMode = .true., cutBack = .false.
 logical, dimension(:,:), allocatable :: calcMode
 integer(pInt), dimension(:,:), allocatable :: FEsolving_execIP
 integer(pInt), dimension(2) :: FEsolving_execElem
 character(len=1024) restartJob

 CONTAINS

!***********************************************************
! determine wether a symmetric solver is used
! and whether restart is requested
!***********************************************************
 subroutine FE_init()
 
 use prec, only: pInt
 use IO
 implicit none
 
 integer(pInt), parameter :: fileunit = 222
 integer(pInt), parameter :: maxNchunks = 6
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=64) tag
 character(len=1024) line


 if (IO_open_inputFile(fileunit)) then
 
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=100) line
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('solver')
         read (fileunit,'(a1024)',END=100) line  ! next line
         positions = IO_stringPos(line,maxNchunks)
         symmetricSolver = (IO_intValue(line,positions,2) /= 1_pInt)
       case ('restart')
         read (fileunit,'(a1024)',END=100) line  ! next line
         positions = IO_stringPos(line,maxNchunks)
         restartWrite = iand(IO_intValue(line,positions,1),1_pInt) > 0_pInt
         restartRead  = iand(IO_intValue(line,positions,1),2_pInt) > 0_pInt
     end select
   enddo
 else
   call IO_error(101) ! cannot open input file
 endif

100 close(fileunit)
 
 if (restartRead .and. IO_open_logFile(fileunit)) then
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=200) line
     positions = IO_stringPos(line,maxNchunks)
     if ( IO_lc(IO_stringValue(line,positions,1)) == 'restart' .and. &
          IO_lc(IO_stringValue(line,positions,2)) == 'file' .and. &
          IO_lc(IO_stringValue(line,positions,3)) == 'job' .and. &
          IO_lc(IO_stringValue(line,positions,4)) == 'id' ) &
       restartJob = IO_StringValue(line,positions,6)
   enddo
 endif

200 close(fileunit)

!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  FEsolving init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 write(6,*) 'restart writing:    ', restartWrite
 write(6,*) 'restart reading:    ', restartRead
 if (restartRead) write(6,*) 'restart Job:        ', trim(restartJob)
 write(6,*)
!$OMP END CRITICAL (write2out)
 

 return

 end subroutine

 END MODULE FEsolving
