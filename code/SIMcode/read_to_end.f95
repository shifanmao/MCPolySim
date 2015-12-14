!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain within
!     the capsid
!     
!     Andrew Spakowitz
!     Written 4-16-04
      
SUBROUTINE read_to_end(UNITIND,snapnm,INDEND)
  INTEGER UNITIND
  character*11 snapnm       ! File for output
  INTEGER INDEND
  integer testread          ! End of line from previous simulation
  
  OPEN (UNIT = UNITIND,FILE = snapnm,status = 'old')
  testread = 0
  INDEND = -1
  DO while(testread==0)
     READ(UNITIND,*,IOSTAT = testread)
     INDEND = INDEND+1
  ENDDO

  RETURN
END
