!---------------------------------------------------------------!
      
!     
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!     
!     Corrections to force magnitude made 6-3-04.
!     
!     Andrew Spakowitz
!     Written 6-29-04
      
      SUBROUTINE MC_int(DEINT,R,AB,NT,NBIN, &
     	 	V,CHI,KAP,LBOX,DEL,PHIA,PHIB,DPHIA,DPHIB, &
     		INDPHI,NPHI,RP,I1,I2)
      
      DOUBLE PRECISION R(NT,3)	! Conformation of polymer chains
      DOUBLE PRECISION RP(NT,3)	! Conformation of polymer chains
      INTEGER AB(NT)            ! Chemical identity of beads    
      INTEGER NT                ! Number of beads
      
!     Simulation input variables
      

      DOUBLE PRECISION V        ! Monomer volume
      DOUBLE PRECISION CHI      ! Chi parameter value
      DOUBLE PRECISION KAP      ! Compressibility value
      DOUBLE PRECISION LBOX     ! Simulation box size (approximate)
      DOUBLE PRECISION DEL      ! Discretization size (approximate)
      INTEGER POLY              ! Polydisperse (step-growth stat)
      DOUBLE PRECISION P        ! Degree of polymerization (polydisperse)
      INTEGER N                 ! Number of beads (monodisperse)
      DOUBLE PRECISION F        ! Fraction of A monomers
      INTEGER FRMFILE           ! Initial condition
      
!     Variables for density calculation
      
      DOUBLE PRECISION PHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION PHIB(NBIN) ! Volume fraction of B
      DOUBLE PRECISION DPHIA(NBIN) ! Volume fraction of A
      DOUBLE PRECISION DPHIB(NBIN) ! Volume fraction of B
      INTEGER INDPHI(NBIN)      ! Indices of the phi
      INTEGER NPHI		! Number of phi values that change 
      DOUBLE PRECISION DR(NT,3) ! Change in bead position
      DOUBLE PRECISION DEINT    ! Change in Self-interaction energy
      INTEGER I1                ! Test bead position 1
      INTEGER I2                ! Test bead position 2
      DOUBLE PRECISION RT(NT,3)
      INTEGER I,J
      INTEGER IB
      INTEGER STOP
      INTEGER NBINX

      INTEGER IX(2),IY(2),IZ(2)      
      DOUBLE PRECISION WX(2),WY(2),WZ(2)
      DOUBLE PRECISION WTOT
      DOUBLE PRECISION RBIN(3)
      INTEGER INDBIN
      INTEGER ISX,ISY,ISZ
      
      do 10 I=1,NBIN
         DPHIA(I)=0.
         DPHIB(I)=0.
	 INDPHI(I)=0
 10   continue
      NBINX=nint(LBOX/DEL)

      DO 15 I=I1,I2
         DR(I,1)=RP(I,1)-R(I,1)
         DR(I,2)=RP(I,2)-R(I,2)
         DR(I,3)=RP(I,3)-R(I,3)
 15   CONTINUE
      
      NPHI=0
      
      do 20 IB=I1,I2
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
         
         do 30 ISX=1,2
            do 40 ISY=1,2
               do 50 ISZ=1,2
                  WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                  INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                  
                  I=1
                  STOP=0
                  do while (STOP.EQ.0)
                     if (INDBIN.EQ.INDPHI(I).AND.I.LE.NPHI) then
                        STOP=1
                        if (AB(IB).EQ.1) then
                           DPHIA(I)=DPHIA(I)-WTOT*V/DEL**3.
                        elseif (AB(IB).EQ.0) then
                           DPHIB(I)=DPHIB(I)-WTOT*V/DEL**3.
                        endif
                     elseif (I.GT.NPHI) then
                        STOP=1
                        NPHI=NPHI+1
                        INDPHI(I)=INDBIN

                        if (AB(IB).EQ.1) then
                           DPHIA(I)=DPHIA(I)-WTOT*V/DEL**3.
                        elseif (AB(IB).EQ.0) then
                           DPHIB(I)=DPHIB(I)-WTOT*V/DEL**3.
                        endif
                     else
                        I=I+1
                     endif                     
                  enddo

 50            continue
 40         continue
 30      continue
 20   continue

      do 60 IB=I1,I2
         RBIN(1)=R(IB,1)+DR(IB,1)
         RBIN(2)=R(IB,2)+DR(IB,2)
         RBIN(3)=R(IB,3)+DR(IB,3)

         RBIN(1)=RBIN(1)-nint(RBIN(1)/LBOX-0.5)*LBOX
         RBIN(2)=RBIN(2)-nint(RBIN(2)/LBOX-0.5)*LBOX
         RBIN(3)=RBIN(3)-nint(RBIN(3)/LBOX-0.5)*LBOX
         
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
         
         do 70 ISX=1,2
            do 80 ISY=1,2
               do 90 ISZ=1,2
                  WTOT=WX(ISX)*WY(ISY)*WZ(ISZ)
                  INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX**2
                  
                  I=1
                  STOP=0
                  do while (STOP.EQ.0)
                     if (INDBIN.EQ.INDPHI(I).AND.I.LE.NPHI) then
                        STOP=1
                        if (AB(IB).EQ.1) then
                           DPHIA(I)=DPHIA(I)+WTOT*V/DEL**3.
                        elseif (AB(IB).EQ.0) then
                           DPHIB(I)=DPHIB(I)+WTOT*V/DEL**3.
                        endif
                     elseif (I.GT.NPHI) then
                        STOP=1
                        NPHI=NPHI+1
                        INDPHI(I)=INDBIN
                        if (AB(IB).EQ.1) then
                           DPHIA(I)=DPHIA(I)+WTOT*V/DEL**3.
                        elseif (AB(IB).EQ.0) then
                           DPHIB(I)=DPHIB(I)+WTOT*V/DEL**3.
                        endif
                     else
                        I=I+1
                     endif
                     
                  enddo
 90            continue
 80         continue
 70      continue
 60   continue
      
      DEINT=0.
      do 100 I=1,NPHI
         J=INDPHI(I)
         DEINT=DEINT+(DEL**3.)*(CHI/V)*((PHIA(J)+DPHIA(I))*(PHIB(J)+DPHIB(I))-PHIA(J)*PHIB(J)) &
     		    +(DEL**3.)*(KAP/V)*((PHIA(J)+DPHIA(I)+PHIB(J)+DPHIB(I)-1.)**2.-(PHIA(J)+PHIB(J)-1.)**2.)
 100  continue
      
      RETURN
      END
      
!---------------------------------------------------------------!
