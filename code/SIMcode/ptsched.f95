!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain within
!     the capsid
!     
!     Andrew Spakowitz
!     Written 4-16-04
      
SUBROUTINE ptsched(INDPT,IND,NRABOVE,NRBELOW,NRMAX,NRNOW,CHI,ECHI,NPT,PTID)
  
  use mt19937, only : grnd, sgrnd, rnorm, mt, mti

  PARAMETER (PI=3.141593)
  INTEGER NR                 ! Total number of replicas  
  INTEGER NRABOVE                 ! Total number of replicas  
  INTEGER NRBELOW                 ! Total number of replicas  
  INTEGER NRNOW                 ! Total number of replicas  
  INTEGER NRMAX                 ! Total number of replicas  
  INTEGER NSORT(NRMAX+1)
  DOUBLE PRECISION CHISORT(NRMAX+1)

  INTEGER IR                ! Replica index
  INTEGER K,J
  INTEGER IND               ! Ind in series
  INTEGER INDPT               ! Ind in series
  INTEGER INDWRT
  DOUBLE PRECISION CHI  ! Chi interaction parameter
  DOUBLE PRECISION CHIOLD
  DOUBLE PRECISION ECHI ! Chi interaction parameter
  DOUBLE PRECISION EPT      ! Exchange energy difference
  DOUBLE PRECISION PPT
  DOUBLE PRECISION PPTABOVE1, PPTABOVE2, PPTABOVE3      ! Exchange probability
  DOUBLE PRECISION PPTBELOW1, PPTBELOW2, PPTBELOW3      ! Exchange probability

  DOUBLE PRECISION TPT      ! Random number
  DOUBLE PRECISION CHITEMP  ! Temporary CHI for exchange
  INTEGER NPT	            ! Number of savepoints between tempering

  INTEGER PTALL(NRMAX+1)
  INTEGER PTCHECK
  INTEGER PTNOW, PTABOVE, PTBELOW
  DOUBLE PRECISION ECHIABOVE, ECHIBELOW
  DOUBLE PRECISION CHIABOVE, CHIBELOW
  INTEGER SWAPNOW, SWAPABOVE, SWAPBELOW
  INTEGER CALCPNOW, CALCPABOVE, CALCPBELOW
  INTEGER ACCBELOW, ACCABOVE
  DOUBLE PRECISION DUMB1, DUMB2
  INTEGER PTID              ! ID to pair up replicas for PT
  INTEGER PTIDM             ! Modify signs of PTID
  INTEGER TENSABOVE              ! Decimal of index
  INTEGER TENSBELOW              ! Decimal of index
  INTEGER TENSTEST              ! Decimal of index
  character*4 numabove       ! Index of output
  character*4 numbelow       ! Index of output
  character*4 numtest       ! Index of output
  character*30 fileabove       ! Index of output
  character*30 filebelow       ! Index of output
  character*30 filetest       ! Index of output

  PTNOW=INDPT
  SWAPNOW=INDPT
  CALCPNOW=INDPT
  INDWRT=1

  PTABOVE=0
  PTBELOW=0
  CALCPABOVE=0
  CALCPBELOW=0
  SWAPABOVE=0
  SWAPBELOW=0
  ACCBELOW=0
  ACCABOVE=0
  PTIDM=1

!  PRINT*, 'PTSTEP NOW ', INDPT, NPT

  IF (MOD(INDPT,NPT).EQ.0) THEN
     !Step 1: Synchronize before calculating probability
!     PRINT*, ' Parallel Tried'
     OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='old')
     WRITE(unit=1,fmt=*,IOSTAT=IOStatus) PTNOW, ECHI, CHI, PTALL(1:(NRMAX+1))
     CLOSE(unit=1,IOSTAT=IOStatus)

     ! Sort CHIs to update NRABOVE and NRBELOW
     DO K=0,NRMAX
        PTALL(K+1)=0
        CHISORT(K+1)=0
     ENDDO

     DO WHILE (PTNOW.GT.(MINVAL(PTALL)))
        DO K=0,NRMAX
           TENSTEST=nint(log10(1.*K)-0.4999)+1
           write (numtest,'(I4)'), K
           filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/ptnow'
           
           OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
           READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1), ECHIABOVE, CHISORT(K+1)
           CLOSE(unit=1,IOSTAT=IOStatus)
        ENDDO

        OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) PTNOW, ECHI, CHI, PTALL(1:(NRMAX+1))
        CLOSE(unit=1,IOSTAT=IOStatus)
     ENDDO

     DO K=1,NRMAX+1
        NSORT(K)=0
     ENDDO

     DO K=1,NRMAX+1
        DO J=1,NRMAX+1
           IF (CHISORT(J).LT.CHISORT(K)) THEN
              NSORT(K)=NSORT(K)+1
           ENDIF
        ENDDO
     ENDDO

     NRABOVE=-1
     NRBELOW=NRMAX+1

     DO J=1,NRMAX+1
        IF (NSORT(J).EQ.(NSORT(NRNOW+1)-1)) THEN
           NRABOVE=J-1
        ENDIF
        
        IF (NSORT(J).EQ.(NSORT(NRNOW+1)+1)) THEN
           NRBELOW=J-1
        ENDIF
     ENDDO

     IF (MOD(INDPT/NPT,10).EQ.0) THEN
        PTIDM=1
     ELSE
        PTIDM=-1
     ENDIF

     IF (MOD(NSORT(NRNOW+1),2).EQ.0) THEN
        PTID=-1*PTIDM
     ELSE
        PTID=1*PTIDM
     ENDIF
     ! End of sorting

     ! Re-read neighbor CHIs
     TENSABOVE=nint(log10(1.*NRABOVE)-0.4999)+1
     TENSBELOW=nint(log10(1.*NRBELOW)-0.4999)+1
     write (numabove,'(I4)'), NRABOVE
     write (numbelow,'(I4)'), NRBELOW
     fileabove='../rand-wlc-'//numabove((4-TENSABOVE+1):4)//'/data/ptnow'
     filebelow='../rand-wlc-'//numbelow((4-TENSBELOW+1):4)//'/data/ptnow'

     ECHIABOVE=ECHI
     ECHIBELOW=ECHI
     CHIABOVE=CHI
     CHIBELOW=CHI
     IF (NRABOVE.GT.0) THEN
        OPEN(unit=1,file=fileabove,IOSTAT=IOStatus,status='old')
        READ(unit=1,fmt=*,IOSTAT=IOStatus) PTABOVE, ECHIABOVE, CHIABOVE
        CLOSE(unit=1,IOSTAT=IOStatus)
     ENDIF

     IF (NRBELOW.LT.(NRMAX+1)) THEN
        OPEN(unit=1,file=filebelow,IOSTAT=IOStatus,status='old')
        READ(unit=1,fmt=*,IOSTAT=IOStatus) PTBELOW, ECHIBELOW, CHIBELOW
        CLOSE(unit=1,IOSTAT=IOStatus)
     ENDIF

     !Step 2: Calculate swap probability
!     PRINT*, 'PTSTEP 2'
     if (CHI.EQ.0) then
        EPT = 0
     else
        EPT = (CHI)*(ECHI/CHI-ECHIABOVE/CHIABOVE)
     endif
     PPTABOVE1 = EPT

     if (CHI.EQ.0) then
        EPT = 0
     else
        EPT = (CHI)*(ECHI/CHI-ECHIBELOW/CHIBELOW)
     endif
     PPTBELOW1 = EPT
     
     !Step 3: Synchronize before calculating probability
!     PRINT*, 'PTSTEP 3'
     OPEN(unit=1,file='data/calcpnow',IOSTAT=IOStatus,status='old')
     WRITE(unit=1,fmt=*,IOSTAT=IOStatus) CALCPNOW, PPTABOVE1, PPTBELOW1, PTALL(1:(NRMAX+1))
     CLOSE(unit=1,IOSTAT=IOStatus)

     DO K=0,NRMAX
        PTALL(K+1)=0
     ENDDO

     DO WHILE (CALCPNOW.GT.(MINVAL(PTALL)))
        DO K=0,NRMAX
           TENSTEST=nint(log10(1.*K)-0.4999)+1
           write (numtest,'(I4)'), K
           filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/calcpnow'

           OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
           READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1)
           CLOSE(unit=1,IOSTAT=IOStatus)
        ENDDO

        OPEN(unit=1,file='data/calcpnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) CALCPNOW, PPTABOVE1, PPTBELOW1, PTALL(1:(NRMAX+1))
        CLOSE(unit=1,IOSTAT=IOStatus)
     ENDDO

     TENSABOVE=nint(log10(1.*NRABOVE)-0.4999)+1
     TENSBELOW=nint(log10(1.*NRBELOW)-0.4999)+1
     write (numabove,'(I4)'), NRABOVE
     write (numbelow,'(I4)'), NRBELOW
     fileabove='../rand-wlc-'//numabove((4-TENSABOVE+1):4)//'/data/calcpnow'
     filebelow='../rand-wlc-'//numbelow((4-TENSBELOW+1):4)//'/data/calcpnow'

     IF (NRBELOW.LT.(NRMAX+1)) THEN
        OPEN(unit=1,file=filebelow,IOSTAT=IOStatus,status='old')
        READ(unit=1,fmt=*,IOSTAT=IOStatus) CALCPBELOW, PPTABOVE2, PPTBELOW2
        CLOSE(unit=1,IOSTAT=IOStatus)
     ENDIF

     !Step 4: Decide whether or not to swap
!     PRINT*, 'PTSTEP 4'
     IF (PTID.EQ.1) THEN
        TTP = grnd()
        PPT = PPTBELOW1+PPTABOVE2

        if (TTP.le.exp(PPT)) then
           ACCBELOW=1
        endif

        OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW, PTALL(1:(NRMAX+1))
        CLOSE(unit=1,IOSTAT=IOStatus)

        DO K=0,NRMAX
           PTALL(K+1)=0
        ENDDO

        DO WHILE (SWAPNOW.GT.(MINVAL(PTALL)))
           DO K=0,NRMAX
              TENSTEST=nint(log10(1.*K)-0.4999)+1
              write (numtest,'(I4)'), K
              filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/swapnow'

              OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
              READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1)
              CLOSE(unit=1,IOSTAT=IOStatus)
           ENDDO

           OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
           WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW, PTALL(1:(NRMAX+1))
           CLOSE(unit=1,IOSTAT=IOStatus)
        ENDDO
        OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW, PTALL(1:(NRMAX+1))
        CLOSE(unit=1,IOSTAT=IOStatus)

        TENSBELOW=nint(log10(1.*NRBELOW)-0.4999)+1
        write (numbelow,'(I4)'), NRBELOW
        filebelow='../rand-wlc-'//numbelow((4-TENSBELOW+1):4)//'/data/swapnow'

        IF (NRBELOW.LE.NRMAX) THEN
           OPEN(unit=1,file=filebelow,IOSTAT=IOStatus,status='old')
           READ(unit=1,fmt=*,IOSTAT=IOStatus) SWAPBELOW, CHIBELOW
           CLOSE(unit=1,IOSTAT=IOStatus)

           if (ACCBELOW.EQ.1) then
              CHI=CHIBELOW
           endif
        ENDIF

        !another check pt. before moving on to MC sim
!        OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='old')
!        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW
!        CLOSE(unit=1,IOSTAT=IOStatus)

!        DO K=0,NRMAX
!           PTALL(K+1)=0
!        ENDDO

!        DO WHILE (SWAPNOW.GT.(MINVAL(PTALL)))
!           DO K=0,NRMAX
!              TENSTEST=nint(log10(1.*K)-0.4999)+1
!              write (numtest,'(I4)'), K
!              filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/swapend'

!              OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
!              READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1)
!              CLOSE(unit=1,IOSTAT=IOStatus)
!           ENDDO

!           OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='old')
!           WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW
!           CLOSE(unit=1,IOSTAT=IOStatus)
!        ENDDO

     ELSEIF (PTID.EQ.-1) THEN

        OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW
        CLOSE(unit=1,IOSTAT=IOStatus)
        
        DO K=0,NRMAX
           PTALL(K+1)=0
        ENDDO

        DO WHILE (SWAPNOW.GT.(MINVAL(PTALL)))
           DO K=0,NRMAX
              TENSTEST=nint(log10(1.*K)-0.4999)+1
              write (numtest,'(I4)'), K
              filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/swapnow'

              OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
              READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1)
              CLOSE(unit=1,IOSTAT=IOStatus)
           ENDDO

           OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
           WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW, PTALL(1:(NRMAX+1))
           CLOSE(unit=1,IOSTAT=IOStatus)
        ENDDO
        OPEN(unit=1,file='data/swapnow',IOSTAT=IOStatus,status='old')
        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW, CHI, ACCBELOW
        CLOSE(unit=1,IOSTAT=IOStatus)

        TENSABOVE=nint(log10(1.*NRABOVE)-0.4999)+1
        write (numabove,'(I4)'), NRABOVE
        fileabove='../rand-wlc-'//numabove((4-TENSABOVE+1):4)//'/data/swapnow'

        IF (NRABOVE.GE.0) THEN
           OPEN(unit=1,file=fileabove,IOSTAT=IOStatus,status='old')
           READ(unit=1,fmt=*,IOSTAT=IOStatus) SWAPABOVE, CHIABOVE, ACCABOVE
           CLOSE(unit=1,IOSTAT=IOStatus)
           if (ACCABOVE.EQ.1) then
              CHI=CHIABOVE
           endif
        ENDIF

        !another check pt. before moving on to MC sim
!        OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='old')
!        WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW
!        CLOSE(unit=1,IOSTAT=IOStatus)

!        DO K=0,NRMAX
!           PTALL(K+1)=0
!        ENDDO

!        DO WHILE (SWAPNOW.GT.(MINVAL(PTALL)))
!           DO K=0,NRMAX
!              TENSTEST=nint(log10(1.*K)-0.4999)+1
!              write (numtest,'(I4)'), K
!              filetest='../rand-wlc-'//numtest((4-TENSTEST+1):4)//'/data/swapend'

!              OPEN(unit=1,file=filetest,IOSTAT=IOStatus,status='old')
!              READ(unit=1,fmt=*,IOSTAT=IOStatus) PTALL(K+1)
!              CLOSE(unit=1,IOSTAT=IOStatus)
!           ENDDO

!           OPEN(unit=1,file='data/swapend',IOSTAT=IOStatus,status='old')
!           WRITE(unit=1,fmt=*,IOSTAT=IOStatus) SWAPNOW
!           CLOSE(unit=1,IOSTAT=IOStatus)
!        ENDDO

     ENDIF

     OPEN(unit=1,file='data/ptnow',IOSTAT=IOStatus,status='old')
     WRITE(unit=1,fmt=*,IOSTAT=IOStatus) PTNOW, ECHI, CHI, PTALL(1:(NRMAX+1))
     CLOSE(unit=1,IOSTAT=IOStatus)
  ENDIF

  RETURN
END
