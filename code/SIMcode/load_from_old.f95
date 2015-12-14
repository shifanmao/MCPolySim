!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain within
!     the capsid
!     
!     Andrew Spakowitz
!     Written 4-16-04
      
      SUBROUTINE load_from_old(R,U,AB,CHI,DCHI,NT,N,NP,IDUM,FRMFILE,PARA,LBOX,INDEND)

      use mt19937, only : grnd, init_genrand, rnorm, mt, mti
      
      PARAMETER (PI=3.141593)
      INTEGER NR                 ! Total number of replicas      
      DOUBLE PRECISION R(NT,3)
      DOUBLE PRECISION U(NT,3)
      INTEGER AB(NT)            ! Chemical identity of beads
      DOUBLE PRECISION CHI      ! Chi interaction parameter
      DOUBLE PRECISION DCHI     ! Chi interaction parameter

      INTEGER N,NP,NT           ! Number of beads
      DOUBLE PRECISION GAM       ! Equil bead separation
      DOUBLE PRECISION LBOX     ! Box edge length
      INTEGER I,J,IB            ! Index Holders
      INTEGER FRMFILE           ! Is conformation in file?
      INTEGER INPUT             ! Is input file set?
      DOUBLE PRECISION RMIN
      DOUBLE PRECISION R0(3)
      DOUBLE PRECISION PARA(10)
      
      character*4 fileind       ! Index of output
      character*4 repind        ! Replica index for output
      INTEGER INDEND
      INTEGER INDTEST
      INTEGER INDNOW
      character*16 snapnm       ! File for output      
      INTEGER TENS              ! Decimal of index
      INTEGER TENR              ! Decimal of index
      INTEGER stat     ! Status output for erase current output
      
!     Variables in the simulation
      
      DOUBLE PRECISION KAP,EPS  ! Elastic properties
      DOUBLE PRECISION XI       ! Drag coefficients

!     Random number generator initiation

      integer IDUM
      character*8 datedum
      character*10 timedum
      character*5 zonedum
      integer seedvalues(8)
      
!     Setup the choice parameters
      
      INPUT=1

!     Seed the random number generator off the computer clock

      call date_and_time(datedum,timedum,zonedum,seedvalues)
	
! concatenate filename, time within mins, secs, millisecs to seed random number generator	

      IDUM=-seedvalues(5)*1E7-seedvalues(6)*1E5-seedvalues(7)*1E3-seedvalues(8)
      call init_genrand(IDUM)
      
!     Input the conformation if FRMFILE=1
         call read_to_end(1, 'data/out2',INDTEST)
         PRINT*, 'restart from replica # ', ' conformation # ', INDTEST
         
         ! erase current unfinished index
         TENS=nint(log10(1.*INDTEST)-0.4999)+1
         write (fileind,'(I4)'), INDTEST

         snapnm= 'data/r'//fileind((4-TENS+1):4)
         OPEN (UNIT = 1, IOSTAT = stat, FILE = snapnm, STATUS = 'OLD')
         IF (stat.EQ.0) CLOSE(1, STATUS = 'delete')

         snapnm= 'data/u'//fileind((4-TENS+1):4)
         OPEN (UNIT = 2, IOSTAT = stat, FILE = snapnm, STATUS = 'OLD')
         IF (stat.EQ.0) CLOSE(2, STATUS = 'delete')

         snapnm= 'data/p'//fileind((4-TENS+1):4)
         OPEN (UNIT = 3, IOSTAT = stat, FILE = snapnm, STATUS = 'OLD')
         IF (stat.EQ.0) CLOSE(3, STATUS = 'delete')

         ! read in one index before
         INDEND = INDTEST-1

         TENS=nint(log10(1.*INDEND)-0.4999)+1
         write (fileind,'(I4)'), INDEND

         snapnm= 'data/r'//fileind((4-TENS+1):4)
         OPEN (UNIT = 1, FILE = snapnm, STATUS = 'OLD')
         DO 10 I=1,NT
            READ(1,"(3f7.2,I2)") R(I,1),R(I,2),R(I,3),AB(I)
10          CONTINUE
         CLOSE(1)
         
         snapnm= 'data/u'//fileind((4-TENS+1):4)
         OPEN (UNIT = 2, FILE = snapnm, STATUS = 'OLD')
         DO 20 I=1,NT
            READ(2,"(3f7.2,I2)") U(I,1),U(I,2),U(I,3),AB(I)
20          CONTINUE 
         CLOSE(2)

         OPEN (UNIT = 3, FILE = 'data/out2', STATUS = 'OLD')
         DO I = 1,INDEND-1
            READ(3,*)
         ENDDO
         READ(3,*) J, CHI, DCHI
         CLOSE(3)

RETURN     
END
      
!---------------------------------------------------------------*
