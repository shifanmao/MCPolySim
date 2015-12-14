!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain within
!     the capsid
!     
!     Andrew Spakowitz
!     Written 4-16-04
      
SUBROUTINE initchem(AB,NT,N,G,NP,FA,LAM)

  use mt19937, only : grnd, init_genrand, rnorm, mt, mti
  
  PARAMETER (PI=3.141593)
  
  INTEGER AB(NT)            ! Chemical identity of beads
  INTEGER N					! Number of monomers per polymer
  INTEGER G					! Number of beads per monomer
  INTEGER NP				! Number of polymer chains
  INTEGER NT				! Total number of beads
  
  INTEGER I,J,K,IB
  DOUBLE PRECISION TEST
  INTEGER ABVAL

  DOUBLE PRECISION FA		! Fraction of A beads
  DOUBLE PRECISION LAM		! Chemical correlation parameter
  DOUBLE PRECISION PAA,PBB,PAB,PBA	! Chemical identity statistics

  !		Translate LAM and FA to probabilities
	  
  PAA=FA*(1.-LAM)+LAM
  PBB=FA*(LAM-1.)+1.
  PBA=1.-PAA
  PAB=1.-PBB

  !		Determine the bead identities
  
  IB=1
  DO 10 I=1,NP
     TEST=grnd()
     if (TEST.LE.FA) then
        AB(IB)=1
     else
        AB(IB)=0
     endif
     IB=IB+1
     DO 15 K=2,G
        AB(IB)=AB(IB-1)
        IB=IB+1
15   CONTINUE
        
     DO 20 J=2,N
        TEST=grnd()		 
        if (AB(IB-1).EQ.1) then
           if (TEST.LE.PAA) then
              AB(IB)=1
           else
              AB(IB)=0
           endif
        else
           if (TEST.LE.PAB) then
              AB(IB)=1
           else
              AB(IB)=0
           endif
        endif
        IB=IB+1
           
     DO 30 K=2,G
        AB(IB)=AB(IB-1)
        IB=IB+1
30      CONTINUE
20      CONTINUE
10      CONTINUE
      
RETURN     
END
      
!---------------------------------------------------------------*
