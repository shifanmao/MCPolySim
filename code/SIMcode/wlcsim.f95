!---------------------------------------------------------------*
      
PROGRAM wlcsim
      
!     
!     This simulation tracks the dynamics of a single polymer
!     chain modeled as a discrete wormlike chain with bending
!     and stretching energy.
!     
!     Andrew Spakowitz
!     Written 9-2-13
!     
      
!     Variables within the simulation

  use mt19937, only : grnd, sgrnd, rnorm, mt, mti

  PARAMETER (PI=3.141592654) ! Value of pi
  INTEGER PTON                  ! Parallel Tempering on
  INTEGER NRABOVE                 ! Total number of replicas  
  INTEGER NRBELOW                 ! Total number of replicas  
  INTEGER NRNOW                 ! Total number of replicas  
  INTEGER NRMAX                 ! Total number of replicas  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: R0	 ! Conformation of polymer chains
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: U0	 ! Conformation of polymer chains
  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIA ! Volume fraction of A
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: PHIB ! Volume fraction of B
  INTEGER, ALLOCATABLE, DIMENSION(:):: AB            ! Chemical identity of beads

  INTEGER NT                ! Number of beads in simulation
  INTEGER N                 ! Number of beads in simulation
  INTEGER NB                 ! Number of beads in simulation
  INTEGER NP                ! Number of polymers in simulation
  DOUBLE PRECISION ENERGY   ! Total energy
  DOUBLE PRECISION DT       ! Time step size
  INTEGER I,J,IB            ! Index
  INTEGER INDMAX            ! Maximum index in series
  INTEGER IND               ! Ind in series
  INTEGER TENS              ! Decimal of index
  character*4 fileind       ! Index of output
  character*4 repind        ! Replica index for output
  character*16 snapnm       ! File for output
  INTEGER INDEND            ! Restart index
  logical restart           ! Restart from previous?
  INTEGER WRTON

!     Simulation input variables
  
  INTEGER FRMFILE           ! Initial condition
  INTEGER BROWN             ! Include Brownian forces
  INTEGER INTON             ! Include polymer interactions
  INTEGER NSTEP				! Number of MC steps between save
  INTEGER NINIT				! Number of initialization MC steps
  INTEGER NNOINT			! Number of initialization MC steps without interactions
  INTEGER NCHI			        ! Number of savepoints between chi change
  INTEGER NPT			        ! Number of savepoints between tempering
  INTEGER FRMCHEM           ! Initial chemical sequence

!     Monte Carlo variables

  DOUBLE PRECISION MCAMP(6) ! Amplitude of random change
  INTEGER MOVEON(6)			! Is the move active
  INTEGER WINDOW(6)			! Size of window for bead selection
  INTEGER SUCCESS(6)        ! Number of successes
  DOUBLE PRECISION PHIT(6)     ! % hits per total steps
      
!     Energy variables
      
  DOUBLE PRECISION EELAS(3) ! Elastic force
  DOUBLE PRECISION ECHI   ! CHI energy
  DOUBLE PRECISION EKAP   ! KAP energy
      
!     Structure analysis
      
  DOUBLE PRECISION RCOM(3)  ! Center of mass
  DOUBLE PRECISION DELR(3)  ! Mag of gyration tensor
  DOUBLE PRECISION RCOM0(3) ! Init val RCOM
  DOUBLE PRECISION DELR0(3) ! Init val DELR
  DOUBLE PRECISION DRCOM    ! Change in RCOM
  DOUBLE PRECISION SIG(3,3)
  DOUBLE PRECISION COR
  
  INTEGER  SON               !calculate Structure Factors
  INTEGER, PARAMETER:: XNUM = 11
  INTEGER, PARAMETER:: KNUM = (XNUM)**3
  INTEGER NVEC              !number of times calculating SVEC
  DOUBLE PRECISION KVEC(KNUM)
  DOUBLE PRECISION SVEC(KNUM)
      
!     Variables in the simulation
      
  DOUBLE PRECISION PARA(10)
  INTEGER NBIN              ! Number of bins
  INTEGER NBINX              ! Number of bins

!     Variables for the random number generators

  INTEGER IDUM              ! Seed for the generator
  DOUBLE PRECISION MOM(6)

!     Simulation parameters
  
  INTEGER G					! Beads per monomer
  DOUBLE PRECISION LBOX		! Box length (approximate)
  DOUBLE PRECISION DEL      ! Discretization size (approximate)
  DOUBLE PRECISION V		! Bead volume
  DOUBLE PRECISION FA		! Fraction of A beads
  DOUBLE PRECISION LAM		! Chemical correlation parameter
  DOUBLE PRECISION EPS		! Elasticity l0/(2lp)
  DOUBLE PRECISION CHI       ! Chi parameter value
  DOUBLE PRECISION DCHI       ! Chi parameter value
  DOUBLE PRECISION CHI0         ! Initial Chi before ramping
  DOUBLE PRECISION DELCHI       ! Chi ramping rate
  DOUBLE PRECISION KAP		! Incompressibility parameter
  DOUBLE PRECISION L0       ! Equilibrium segment length
  INTEGER PTID              ! ID to pair up replicas for PT
  INTEGER ACCBELOW

!     Load in the parameters for the simulation

  open (unit=5, file='input/input')
  read (unit=5, fmt='(4(/))')
  read (unit=5, fmt=*) PTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) N
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) G
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) LBOX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) DEL
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) L0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) V
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FA
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) LAM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMCHEM
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) EPS
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) CHI0
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) DELCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) KAP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INDMAX
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) FRMFILE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) BROWN
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) INTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NCHI
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NPT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NNOINT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NINIT
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NSTEP
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) WRTON
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) PTID
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRABOVE
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRBELOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRNOW
  read (unit=5, fmt='(2(/))')
  read (unit=5, fmt=*) NRMAX
  close(5)

  NP=nint(LBOX**3./(N*G*V))
  LBOX=(V*N*G*NP)**(1./3.)
  call getpara(PARA,EPS,L0,LBOX)
  
  NT=N*NP*G
  ALLOCATE(R(NT,3))
  ALLOCATE(U(NT,3))
  ALLOCATE(R0(NT,3))
  ALLOCATE(U0(NT,3))
  ALLOCATE(AB(NT))

  NBINX=nint(LBOX/DEL)
  NBIN=NBINX**3.
  DEL=LBOX/NBINX
  
  ALLOCATE(PHIA(NBIN))
  ALLOCATE(PHIB(NBIN))

!  Monte-Carlo simulation parameters
  MCAMP(1)=0.5*PI
  MCAMP(2)=0.3*L0
  MCAMP(3)=0.5*PI
  MCAMP(4)=0.5*PI
  MCAMP(5)=0.5*PI
  MCAMP(6)=5.0*L0
  MOVEON(1)=1
  MOVEON(2)=1
  MOVEON(3)=1
  MOVEON(4)=1
  if (INTON.EQ.1) then
     MOVEON(5)=1
     MOVEON(6)=1
  else
     MOVEON(5)=1
     MOVEON(6)=1
  endif

!     Initial segment window for MC moves
  WINDOW(1)=N*G
  WINDOW(2)=N*G
  WINDOW(3)=N*G
  WINDOW(4)=N*G
  WINDOW(5)=N*G
  WINDOW(6)=N*G

  INQUIRE (FILE = 'data/out1', exist = restart)
  if (.NOT.restart) then

    PRINT*, '-----new simulation-----'

!     Setup the initial condition

    NB=N*G
    call initcond(R,U,AB,NT,NB,NP,IDUM,FRMFILE,PARA,LBOX)

!     Load in AB sequence
    IF (FRMCHEM.EQ.1) THEN
       OPEN (UNIT = 1, FILE = 'input/ab', STATUS = 'OLD')      
       IB=1
       DO 1 I=1,NP
          DO 2 J=1,NB
             READ(1,"(I2)") AB(IB)
             IB=IB+1
2         CONTINUE
1      CONTINUE 
       CLOSE(1)
    ELSE
       call initchem(AB,NT,N,G,NP,FA,LAM)
    ENDIF

!     Perform an initialization MC simulation
    SON=0
    call MCsim(R,U,PHIA,PHIB,AB,NT,NB,NP,NBIN,NNOINT,BROWN,0,PARA,V,CHI0,KAP,LBOX,L0,DEL,MCAMP,SUCCESS,MOVEON,WINDOW,PHIT,&
         KVEC,SVEC,SON,0,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID)
    call MCsim(R,U,PHIA,PHIB,AB,NT,NB,NP,NBIN,NINIT,BROWN,INTON,PARA,V,CHI0,KAP,LBOX,L0,DEL,MCAMP,SUCCESS,MOVEON,WINDOW,PHIT,&
         KVEC,SVEC,SON,0,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID)

!     Save initial conformation and PSI angles 
    OPEN (UNIT = 1, FILE = 'data/r0', STATUS = 'NEW')      
    IB=1
    DO 10 I=1,NP
       DO 20 J=1,NB
          R0(IB,1)=R(IB,1)
          R0(IB,2)=R(IB,2)
          R0(IB,3)=R(IB,3)
          U0(IB,1)=U(IB,1)
          U0(IB,2)=U(IB,2)
          U0(IB,3)=U(IB,3)
          WRITE(1,"(3f8.3,I2)") R(IB,1),R(IB,2),R(IB,3),AB(IB)
          IB=IB+1
20     CONTINUE
10  CONTINUE 
    CLOSE(1)
      
    OPEN (UNIT = 1, FILE = 'data/u0', STATUS = 'NEW')
    IB=1
    DO 30 I=1,NP
       DO 40 J=1,NB
          WRITE(1,"(3f8.3,I2)") U(IB,1),U(IB,2),U(IB,3),AB(IB)
          IB=IB+1
40     CONTINUE
30  CONTINUE 
    CLOSE(1)

!     Open the output files
    INDEND = 0
    OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'NEW')
    OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'NEW')
    OPEN (UNIT = 3, FILE = 'data/out3', STATUS = 'NEW')

 else

    PRINT*, '-----load simulation-----'
    NB=N*G
    call load_from_old(R,U,AB,CHI,DCHI,NT,NB,NP,IDUM,FRMFILE,PARA,LBOX,INDEND)

 endif

!     Begin simulation
  IND=1
  TIME=0.

  !part 7 - PT
  OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/calcpnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)
  OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='new')
  CLOSE(unit=1,IOSTAT=IOStatus)

  DO WHILE ((IND+INDEND).LE.INDMAX)
!     Parallel tempering + chi annealing
     call chisched(IND,NR,CHI,DCHI,CHI0,DELCHI,NCHI,NPT,restart)

!     Perform a MC simulation
     SON=0

     call MCsim(R,U,PHIA,PHIB,AB,NT,NB,NP,NBIN,NSTEP,BROWN,INTON,PARA,V,CHI,KAP,LBOX,L0,DEL,MCAMP,SUCCESS,MOVEON,WINDOW,PHIT,&
          KVEC,SVEC,SON,PTON,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,ECHI,NPT,PTID)

!     Save the conformation and the metrics
     TENS=nint(log10(1.*INDEND+IND)-0.4999)+1
     write (fileind,'(I4)'), INDEND+IND

      !part 1 - energy
     ECHI=0.
     EKAP=0
     EELAS(1)=0.
     EELAS(2)=0.
     EELAS(3)=0.
     call energy_elas(EELAS(1:3),R,U,NT,NB,NP,PARA)
     if (INTON.EQ.1) then
        call energy_int(R,AB,NP,NB,NT,NBIN,V,CHI,KAP,LBOX,DEL,ECHI,EKAP)
     endif
     
     OPEN (UNIT = 1, FILE = 'data/out1', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(1,*)
     ENDDO
     WRITE(1,"(I4,5f20.3)") INDEND+IND,EELAS(1:3),ECHI,EKAP
     CLOSE(1)
      
     !part 2 - CHI
     OPEN (UNIT = 2, FILE = 'data/out2', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(2,*)
     ENDDO
     !  WRITE(2,"(I15,1f10.2)") INDPT,CHI
     WRITE(2,"(I4,1f20.8)") INDEND+IND,CHI*G
     CLOSE(2)

     !part 2.5 - adaptations
     OPEN (UNIT = 3, FILE = 'data/out3', STATUS = 'OLD')
     DO K = 1,IND+INDEND-1
        READ(3,*)
     ENDDO
     WRITE(3,"(I4,18f8.2)") INDEND+IND,REAL(WINDOW(1)),MCAMP(1),PHIT(1), &
          REAL(WINDOW(2)),MCAMP(2),PHIT(2), &
          REAL(WINDOW(3)),MCAMP(3),PHIT(3), &
          REAL(MOVEON(4)),MCAMP(4),PHIT(4), &
          REAL(MOVEON(5)),MCAMP(5),PHIT(5), &
          REAL(MOVEON(6)),MCAMP(6),PHIT(6)
     CLOSE(3)
      
     IF (WRTON.EQ.1) THEN

     !part 3 - R
     snapnm= 'data/r'//fileind((4-TENS+1):4)
     OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
     IB=1
     DO 50 I=1,NP
        DO 60 J=1,NB
           WRITE(1,"(3f7.2,I2)") , &
                R(IB,1)-0.*nint(R(IB,1)/LBOX-0.5)*LBOX, &
                R(IB,2)-0.*nint(R(IB,2)/LBOX-0.5)*LBOX, &
                R(IB,3)-0.*nint(R(IB,3)/LBOX-0.5)*LBOX,AB(IB)
           IB=IB+1
60         CONTINUE
50         CONTINUE
     CLOSE(1)

     !part 4 - U
!     snapnm= 'data/u'//fileind((4-TENS+1):4)
!     OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
!     IB=1
!     DO 70 I=1,NP
!        DO 80 J=1,NB
!           WRITE(1,"(3f7.2,I2)") U(IB,1),U(IB,2),U(IB,3),AB(IB)
!           IB=IB+1
!80         CONTINUE
!70         CONTINUE 
!     CLOSE(1)

     !part 5 - PHI
!     snapnm= 'data/p'//fileind((4-TENS+1):4)
!     OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
!     DO K=1,NBIN
!        WRITE(1,"(2f7.2)") PHIA(K),PHIB(K)
!     ENDDO
!     CLOSE(1)

     ENDIF

     !part 6 - S(k)
     IF (SON.EQ.1) THEN
        snapnm= 'data/s'//fileind((4-TENS+1):4)
        OPEN (UNIT = 1, FILE = snapnm, STATUS = 'NEW')
        DO K=1,KNUM
           WRITE(1,"(2f20.3)") KVEC(K),SVEC(K)
        ENDDO
        CLOSE(1)
     ENDIF

     
     PRINT*, '________________________________________'
     PRINT*, 'Time point ',IND+INDEND, ' out of', INDMAX
     PRINT*, 'Bending energy ', EELAS(1)
     PRINT*, 'Par compression energy ', EELAS(2)
     PRINT*, 'Perp compression energy ', EELAS(3)
     PRINT*, 'Chi interaction energy ', ECHI
     PRINT*, 'Kap compression energy ', EKAP
     !   PRINT*, 'Current number of beads ', N*G
     !   PRINT*, 'Number of polymers ', NP
         
     IND=IND+1
         
  ENDDO
  
END
      
!---------------------------------------------------------------*
