
!---------------------------------------------------------------*
      
!     Find the volume fractions from R
      
!     Revised 6-22-04
      
      SUBROUTINE r_to_phi(R,AB,NT,N,NP,NTOT,NBIN, &
      V,CHI,KAP,LBOX,DEL,PHIA,PHIB)
      
      PARAMETER (PI=3.141592654) ! Value of pi
      
      DOUBLE PRECISION R(NT,3)	! Conformation of polymer chains
      INTEGER AB(NT)            ! Chemical identity of beads
      INTEGER IND(NT,2)         ! Start and end integers      
      INTEGER NT                ! Number of beads
      INTEGER NP              ! Total number of polymers
      INTEGER N                ! Number of beads per polymer
      
!     Variables for density calculation
      
      DOUBLE PRECISION PHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION PHIB(NBIN) ! Volume fraction of B	  
      
!     Simulation input variables
      
      DOUBLE PRECISION V       ! Volume of monomer A
      DOUBLE PRECISION CHI      ! Chi parameter value
      DOUBLE PRECISION KAP      ! Compressibility value
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)
      INTEGER POLY              ! Polydisperse (step-growth stat)
      DOUBLE PRECISION P        ! Degree of polymerization (polydisperse)
      DOUBLE PRECISION F        ! Fraction of A monomers
      INTEGER FRMFILE           ! Initial condition      
      INTEGER I,J
      
      INTEGER IX(2),IY(2),IZ(2)
      INTEGER IB
      INTEGER NBINX
      
      DOUBLE PRECISION WX(2),WY(2),WZ(2)
      DOUBLE PRECISION WTOT
      DOUBLE PRECISION RBIN(3)
      INTEGER INDBIN
      INTEGER ISX,ISY,ISZ
      
!     Initialize the bins
      
      do 10 I=1,NBIN
         PHIA(I)=0.
         PHIB(I)=0.
 10   continue
      NBINX=nint(LBOX/DEL)
      
!     Cycle through the beads
      
      IB=1
      do 20 I=1,NP
         do 30 J=1,N
            RBIN(1)=R(IB,1)-nint(R(IB,1)/LBOX-0.5)*LBOX
            RBIN(2)=R(IB,2)-nint(R(IB,2)/LBOX-0.5)*LBOX
            RBIN(3)=R(IB,3)-nint(R(IB,3)/LBOX-0.5)*LBOX
            
            IX(1)=nint(RBIN(1)/DEL+0.5)
            IY(1)=nint(RBIN(2)/DEL+0.5)
            IZ(1)=nint(RBIN(3)/DEL+0.5)
            
            IX(2)=IX(1)-1
            IY(2)=IY(1)-1
            IZ(2)=IZ(1)-1
            
!     Calculate the bin weighting
            
            WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL)
            WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL)
            WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL)
            WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL)
            WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL)
            WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL)
            
            IX(1)=IX(1)-floor(REAL((IX(1)-1))/REAL(NBINX)) * NBINX
            IX(2)=IX(2)-floor(REAL((IX(2)-1))/REAL(NBINX)) * NBINX
            IY(1)=IY(1)-floor(REAL((IY(1)-1))/REAL(NBINX)) * NBINX
            IY(2)=IY(2)-floor(REAL((IY(2)-1))/REAL(NBINX)) * NBINX
            IZ(1)=IZ(1)-floor(REAL((IZ(1)-1))/REAL(NBINX)) * NBINX
            IZ(2)=IZ(2)-floor(REAL((IZ(2)-1))/REAL(NBINX)) * NBINX
            
!     Add volume fraction with weighting to each bin
            
            do 40 ISX=1,2
               do 50 ISY=1,2
                  do 60 ISZ=1,2
                     WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                     INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                     if (AB(IB).EQ.1) then
                        PHIA(INDBIN)=PHIA(INDBIN)+WTOT*V/DEL**3.
                     else
                        PHIB(INDBIN)=PHIB(INDBIN)+WTOT*V/DEL**3.
                     endif
 60               continue
 50            continue
 40         continue
            IB=IB+1
            
 30      continue
 20   continue
      
      RETURN            
      END
      
!---------------------------------------------------------------*
