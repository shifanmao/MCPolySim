!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain within
!     the capsid
!     
!     Andrew Spakowitz
!     Written 4-16-04
      
SUBROUTINE chisched(IND,NR,CHI,DCHI,CHI0,DELCHI,NCHI,NPT,restart)

  PARAMETER (PI=3.141593)
  INTEGER NR                ! Total number of replicas
  INTEGER IR                ! Replica index
  INTEGER IND               ! Ind in series
  DOUBLE PRECISION CHI      ! Chi interaction parameter
  DOUBLE PRECISION DCHI     ! Chi interaction parameter
  DOUBLE PRECISION CHI0     ! Initial Chi before ramping
  DOUBLE PRECISION DELCHI   ! Chi ramping rate
  INTEGER NPT               ! Number of savepoints between tempering
  logical restart           ! Restart from previous?

  IF ((IND.LE.1).AND.(.NOT.restart)) THEN
     PRINT*, 'restart chi sched'
     DCHI = DELCHI
     CHI = CHI0
  ELSEIF (MOD(IND,NCHI).EQ.0) THEN
     CHI = CHI+DCHI
  ENDIF

  RETURN
END
