      PROGRAM FLUX7
C PROGAM FLUX      VERSION 7.9.0 April 2013
C
C 7.53: fixed bug in call CROSSECT (type of first argument)
C 7.54: eliminated questionable handling of 'assigned' labels.
C 7.55: support tetragonal lattice.
C 7.6:  new option STARTCOORD
C 7.61: optional: Moliere potential instead of ZBL
C 7.62: fixed a few floating point exception problems
C 7.7:  new option MIXED: mixed lattice, such as Ge(1-x)Si(x)
C       option CROSS SECTION made more general
C       optional Hartree Fock potential
C 7.71: fixed long standing bug in subroutine FFIELD
C       update HF potential
C 7.72: optional energy dependent stopping.
C 7.8 : major revision  -  merging of flux7 with Luis Rebouta's flux4
C 7.8.5: option MIXED now also suitable for more then 2 kinds of atoms
C        on the same type of lattice site.
C 7.8.6: fixed bug in overestimating the continuum string force
C       in cases from COMBI4.FIG.
C 7.8.7: added support for GaN (Wurtzite) <11-20> and <1-100>
C 7.9.0: Adaptation to gfortran, some fortran77 constructs eliminated
C        Multi-processor implementation: OpenMP
C        Flux files now written in decimal format. Programs can read
C        the old hexadecimal format as well as decimal format.
C
C This version is for a multilayer structure of two alternating layers.
C At each interface a transformation is applied to the beam coordinates:
C 1. a rotation of the velocity vector, specified by a 3 by 3
C    rotation matrix
C 2. a translation, specified by delta (x,y,z) in angstrom units
C The close encounter probability etc. is calculated for each layer
C The flux is only calculated for one of the layers
C
C ===---===---===---===---===---===---===---===--- MAIN PROGRAM FLUX
C
      IMPLICIT NONE

#include "FLUX7.COM"

      integer ier, K, KREP

C SET THE VERSION STRING:
      VERSION(1:5)='V7.8'
      VERSION(16:28)='Rev8 2012'

C DO SOME INPUT AND INITIALIZING :
      CALL BATINIT
      CALL INP(INLU,ier)
      if(ier.ne.0)STOP 'error(s) in INP'
      CALL PRELIMI
C CALCULATE THE FIELD OF THE SURROUNDING STRINGS
      CALL FFIELD
C OUTPUT THINGS FIT TO PRINT :
      CALL OUTP1
C NOW REALLY GO TO WORK:
      DO 10 K=   NSTART,NANG
      DO 10 KREP= NRSTART,NRUN
C BINNENSTE BINNENLOOP:
      CALL BBLOOP(K)
C OUTPUT RESULTS FOR CURRENT ANGLE:
      CALL OUTP2(K)
      CALL UPDATE(K,KREP)
  10  CONTINUE
C FINISHED
      END
C ------------------------------------------------- SUBROUTINE BBLOOP
      SUBROUTINE BBLOOP(K)
      IMPLICIT NONE
      INTEGER K
C In this subroutine the ion trajectories are followed

#include "FLUX7.COM"

C
C PFACTOR and TFACTOR are scaling factors for output of PX, PY, T
C Should have the same values in FLUXLIB/getcoord.f !
      REAL PFACTOR, TFACTOR
      PARAMETER(PFACTOR=1000.0, TFACTOR=1000.0)

      INTEGER IZ2
      INTEGER I,ier,IMPS
      INTEGER J,KL,LX
      INTEGER LOGLU
      INTEGER NTRACK

      LOGICAL RANDOM

      REAL COSTHETA
      REAL PHI
      REAL RLX,RNX,RNY
      REAL SINTHETA
      REAL THETA
      REAL Xoff
      REAL Yoff
      REAL PXA,PYA,PZA,XAV,YAV,T0AV,PXoff,PYoff,Toff

      DOUBLE PRECISION  DAVDPR(0:NBIN-1,ML)
      DOUBLE PRECISION DDAVDPR(0:NBIN-1,ML)
      DOUBLE PRECISION  CEYDPR(0:NBIN-1,ML)
      COMMON /DOUBLEDAV/DAVDPR, DDAVDPR, CEYDPR

      REAL     ARCCOS,CROSSECT,FF,URAND
      EXTERNAL ARCCOS,CROSSECT,FF,URAND
c
C FFX,FFY: components of the force field of the surrounding atom rows.
      REAL FFX,FFY
      COMMON /COMFIELD/ FFX(NF*NF,NLAYER),FFY(NF*NF,NLAYER)
c==============================================================

C Inner loop of the ion trajectories simulation


      CHARACTER*12 DATUM

      INTEGER LIND
      INTEGER IB,IDELTAZ,IL,IL1,IL2,INDIC(6),INDX,INDXT
      INTEGER INVERTX,INVERTY,ISWITCH,IX,IXIY,IY,IZ
      INTEGER JIND,JXYOUT,JY,JZ,JZN
      INTEGER KLTRUE, LAYER, LY, LAYER2
      INTEGER NEF,NPPX,NPPXY,NPPY,NZ,NZINTERF
      INTEGER NINDKL
      INTEGER MODZKL
      INTEGER IS
      INTEGER IXUM,IYUM

      LOGICAL BENT
      LOGICAL QZZKL
      LOGICAL QUIT
      LOGICAL DPERR
      LOGICAL DPTEST
      EXTERNAL DPTEST

      REAL ANORM(ML,NLAYER)
      REAL CHANCE
      REAL COSTET
      REAL DECO,DELTAX,DELTAY,DELX,DELY,DENU
      REAL DNZINTERF,RMAXNZ
      REAL FDEVL,FRAC
      REAL PHI0,PX,PX0,PXN,PXX,PY,PY0,PYN,PYY,PYZ,PZ,PZ0,PZN,PZSQ
      REAL QTH,QX,QY,QZ
      REAL RLY,RP,RWBEAM
      REAL S,S0,SINTET,SPARE,SX,SY,SZ
      REAL TEE1,TEE2,TET,TEX,TEY,THICKNESS,TRUEX,TRUEY,TX,TY
      REAL TRUEPX, TRUEPY
      REAL X0,XT,XTH,XX,XXN
      REAL Y0,YT,YTH,YY,YYN
      REAL ZINTFSIG(NLAYER),ZTH,ZW,ZZ
      REAL TT
      REAL DEVLKL, DEVPKL
      REAL UDKL, XUNKL, YUNKL, ZUNKL, XUMKL, YUMKL, TWOXUNKL, TWOYUNKL
      REAL STARTX, STARTY, STARTPX, STARTPY, STARTT
      REAL C1,C2,FRAC2
      REAL DELZ1,DELZ2,DELZZ,DZZ,DIFSUP,DIFINF,DELZZ1,DELZZ2,DELZZ3
      INTEGER JZZ,IJJ
      INTEGER IRAND, NRAND
      REAL RANBUF(NBINF*NDELZ*4+4), UNIRAN(2)
**TEST save old coordinates, see array bound test
      REAL XSAVE,YSAVE,PXSAVE,PYSAVE
      INTEGER JZSAV
**

      DOUBLE PRECISION CLENCY,T,TINV

      INTEGER IXIYZ(0:NBINF*NDELZ-1)
      REAL DELFLUX(0:NBINF*NDELZ-1)
      DOUBLE PRECISION  DELDAVDPR(0:NBINF-1,MLA)
      DOUBLE PRECISION DELDDAVDPR(0:NBINF-1,MLA)
      DOUBLE PRECISION  DELCEYDPR(0:NBINF-1,MLA)
      DOUBLE PRECISION DELRAAK(0:NBINF-1,MLA,NLAY)
      INTEGER DELLEFTOVER(0:NBINF-1), DELNCHANNEL(0:NBINF-1)
c
C
Cleanup:
      DO  2 I =  1,NF*NF
      FLUX(I) = 0.0
  2   CONTINUE
      DO  21 I =  1,NV*NV
      VFINAL(I) = 0.0
  21  CONTINUE
      DO 32 IL = 1,ML
      DO 31 J = 0,NBIN-1
      DDAVDPR(J,IL) = 0.0
      DAVDPR(J,IL)  = 0.0
      CEYDPR(J,IL)  = 0.0
      DO 30 KL = 1,NLAY
      RAAK(J,IL,KL) = 0.0
  30  CONTINUE
  31  CONTINUE
  32  CONTINUE
      DO 33 J = 0,NBIN-1
      NCHANNEL(J) = 0
      LEFTOVER(J) = 0
      SNUCLEAR(J) = 0
      SCORE(J) = 0
      SVALENCE(J) = 0
  33  CONTINUE
      DO 34 J = 0, NEBIN
      EFINAL(J) = 0
  34  CONTINUE
      DO 35 KL=1,NLAY
      DO 35 IL=1,NLS(KL)
      DO 35 J=1,NL
      HITLIST(J,IL,KL)=0
  35  CONTINUE
      TFINAL = 0.0
      DTFINAL = 0.0
      PXFINAL = 0.0
      DPXFINAL = 0.0
      PYFINAL = 0.0
      DPYFINAL = 0.0
      PTFINAL = 0.0
      DPTFINAL = 0.0
      NFINAL = 0
      NRAAK=0
      RENRGY= 0.0

      NAWBX = 0

      RNX=0.999995*NFX
      RNY=0.999995*NFY
C Note. According to U.Wahl the above precision causes problems on
C certain "no name" 66 MHz 486 PC's. The remedy was to use:
c      RNX=0.9995*NFX
c      RNY=0.9995*NFY

C Set switch IMPS according to impurity cross section type, i.e. for flux
      IF (.NOT. JOPTION(LFFLUX)) THEN
C  no flux calculation
        IMPS = 509
      ELSE IF (NCROSS(0,1).LT.0) THEN
C Rutherford cross section
                       IMPS= 500
          IF(ZSIG.EQ.0)IMPS= 503
        ELSE IF (NCROSS(0,1).EQ.0) THEN
C No cross section in flux
                       IMPS= 501
          IF(ZSIG.EQ.0)IMPS= 504
        ELSE
C Interpolated cross section
                       IMPS= 502
          IF(ZSIG.EQ.0)IMPS= 505
      ENDIF

      NTRACK= 0
      T0AV= T0
      CALL THETAPHI(THETA, PHI, K)
      IF ( PHI .GT. 360.01 ) PHI = 0
      RANDOM = .FALSE.
C IF PHI>360 DEGREES WAS GIVEN, USE RANDOM PHI.
C IF THETA>360 DEGREES WAS GIVEN, USE RANDOM MEDIUM (randomize x,y).
      IF ( THETA .GT. 360.01 ) THEN
         write(6,*)'theta=',THETA,': set theta=0, but randomize x,y'
         write(LOGLU,*)'theta=',THETA,': set theta=0, but randomize x,y'
         THETA = 0
         RANDOM = .TRUE.
      ENDIF
      THETA =  THETA*DEGREE
      PHI   =    PHI*DEGREE
      SINTHETA = SIN(THETA)
      COSTHETA = COS(THETA)
      PXA = SINTHETA * COS ( PHI )
      PYA = SINTHETA * SIN ( PHI )
      PZA = SQRT(1.0 - PXA**2 - PYA**2)
      XAV =  (ZMAX - ZZ0)*PXA/PZA
      YAV =  (ZMAX - ZZ0)*PYA/PZA
C Note: PXA,PYA,PZA are constants of a simulation.
C       PX0,PY0,PZ0 are recalculated for each trajectory:
C        they follow the curvature of the crystal and are
C        possibly subject to rotation at a interface.

      IF(LUOPEN(2).eq.1)THEN
        LOGLU=LUN(2)
      ELSE
        LOGLU=6
      ENDIF

      IF (JOPTION(LFXYOUT)) THEN
       write(LOGLU,897)
  897  FORMAT(/'Output during each trajectory:',
     + /'''*xyz'' x, y, z (in angstroms), px,py',
     + /'each new track starts with 5 zeroes.'/)
      ENDIF
      IF (JOPTION(LFCFILE)) THEN
C skip to next record starting with offsets line
C Note that these offsets are also saved in 'getcoord' and applied
C automatically at subsequent calls.
  204   call getcoord(LUN(4),1, Xoff,Yoff,PXoff, PYoff, Toff, ier)
        if(ier .eq. 2) goto 204
        if(ier .gt. 1) then
          write(6,'(''error'',I2,'' in coordinates file'')') ier
          stop 'error in coordinates file'
        endif
        XAV= XAV+Xoff
        YAV= YAV+Yoff
        PXA= PXoff
        PYA= PYoff
        PZA = SQRT(1.0 - PXA**2 - PYA**2)
        T0AV= Toff
      ENDIF
      if(K .gt. 1)WRITE(LOGLU,317)
  317 FORMAT('*new angle')
C

      if(JOPTION(LEXITCO))then
       write(LOGLU,890)
  890   FORMAT('OUTPUT AT EXIT FOR EACH TRAJECTORY:',/ 12X,
     +  5X,'X-X0',
     +  5X,'Y-Y0',
     +  3X,'PX-PX0',
     +  3X,'PY-PY0',
     +  3X,'T0-T')
       call wrt5(LOGLU,'offsets:    ',XAV,YAV,
     +            (PXA)*PFACTOR,(PYA)*PFACTOR,
     +            REAL(T0AV)*TFACTOR)
       call wrt5(LOGLU,'scalefactor:',1.0,1.0,PFACTOR,PFACTOR,TFACTOR)
      endif

C if "#ifdef" or "omp_set_num_threads" are not recognized modify or
C remove the following.
C All we try to achieve is:
C For option XYOUT there may be no parallel threads.
#ifdef OpenMP
      IF (JOPTION(LFXYOUT)) THEN
        call omp_set_num_threads(1)
      ENDIF
#endif

c==============================================================
C loop through the tracks
!$OMP PARALLEL DO
!$OMP& DEFAULT(shared)
!$OMP& PRIVATE(LX,RLX,LY,RLY)
!$OMP& PRIVATE(IZ2)
!$OMP& PRIVATE(ier)
!$OMP& PRIVATE(CLENCY,T,TINV)
!$OMP& PRIVATE(DELCEYDPR)
!$OMP& PRIVATE(DELDAVDPR)
!$OMP& PRIVATE(DELDDAVDPR)
!$OMP& PRIVATE(DELRAAK)
!$OMP& PRIVATE(DELLEFTOVER, DELNCHANNEL)
!$OMP& PRIVATE(I,IB,IDELTAZ,IL,IL1,IL2,INDIC,INDX,INDXT)
!$OMP& PRIVATE(INVERTX,INVERTY,ISWITCH,IX,IXIY,IY,IZ)
!$OMP& PRIVATE(IRAND, NRAND)
!$OMP& PRIVATE(IS)
!$OMP& PRIVATE(IXIYZ)
!$OMP& PRIVATE(IXUM,IYUM)
!$OMP& PRIVATE(J,JIND,JXYOUT,JY,JZ,JZN)
!$OMP& PRIVATE(JZSAV)
!$OMP& PRIVATE(JZZ,IJJ)
!$OMP& PRIVATE(KL, KLTRUE, LAYER, LAYER2)
!$OMP& PRIVATE(LIND)
!$OMP& PRIVATE(MODZKL)
!$OMP& PRIVATE(NEF,NPPX,NPPXY,NPPY,NZ,NZINTERF)
!$OMP& PRIVATE(NINDKL)
!$OMP& PRIVATE(BENT)
!$OMP& PRIVATE(DPERR)
!$OMP& PRIVATE(QUIT)
!$OMP& PRIVATE(QZZKL)
!$OMP& PRIVATE(ANORM)
!$OMP& PRIVATE(C1,C2,FRAC2)
!$OMP& PRIVATE(CHANCE)
!$OMP& PRIVATE(COSTET,COSTHETA)
!$OMP& PRIVATE(DECO,DELTAX,DELTAY,DELX,DELY,DENU)
!$OMP& PRIVATE(DELFLUX)
!$OMP& PRIVATE(DELZ1,DELZ2,DELZZ,DZZ,DIFSUP,DIFINF,DELZZ1,DELZZ2,DELZZ3)
!$OMP& PRIVATE(DEVLKL, DEVPKL)
!$OMP& PRIVATE(DNZINTERF,RMAXNZ)
!$OMP& PRIVATE(FDEVL,FRAC,PHI)
!$OMP& PRIVATE(PHI0,PX,PX0,PXN,PXX,PY,PY0,PYN,PYY,PYZ,PZ,PZ0,PZN,PZSQ)
!$OMP& PRIVATE(QTH,QX,QY,QZ)
!$OMP& PRIVATE(RANBUF, UNIRAN)
!$OMP& PRIVATE(RP,RWBEAM)
!$OMP& PRIVATE(S,S0,SINTET,SINTHETA,SPARE,SX,SY,SZ)
!$OMP& PRIVATE(STARTX, STARTY, STARTPX, STARTPY, STARTT)
!$OMP& PRIVATE(TEE1,TEE2,TET,TEX,TEY,THETA,THICKNESS,TRUEX,TRUEY,TX,TY)
!$OMP& PRIVATE(TRUEPX, TRUEPY)
!$OMP& PRIVATE(TT)
!$OMP& PRIVATE(UDKL,XUNKL,YUNKL,ZUNKL,XUMKL,YUMKL,TWOXUNKL,TWOYUNKL)
!$OMP& PRIVATE(X0,XT,XTH,XX,XXN)
!$OMP& PRIVATE(XSAVE,YSAVE,PXSAVE,PYSAVE)
!$OMP& PRIVATE(Y0,YT,YTH,YY,YYN)
!$OMP& PRIVATE(ZINTFSIG,ZTH,ZW,ZZ)
!$OMP& SCHEDULE(DYNAMIC)

c==============================================================
      DO 8 LX=1,NXI
c==============================================================
      RLX=LX

      ier=0

      DO 1120 IZ=0,NBINF-1
      DO 1119 IL=1,MLA
        DELDAVDPR(IZ,IL) = 0.0
       DELDDAVDPR(IZ,IL) = 0.0
        DELCEYDPR(IZ,IL) = 0.0
      DO 1118 KL=1,NLAY
        DELRAAK(IZ,IL,KL)= 0.0
 1118 CONTINUE
 1119 CONTINUE
      DELLEFTOVER(IZ) = 0
 1120 CONTINUE

c the following is to prevent complaints from compiler:
      ZUNKL=1.0e38
      YUNKL=1.0e38
      YUMKL=1.0e38
      XUNKL=1.0e38
      XUMKL=1.0e38
      UDKL=1.0e38
      TWOYUNKL=1.0e38
      TWOXUNKL=1.0e38
      DEVPKL=1.0e38
      DEVLKL=1.0e38
      MODZKL=2140000000
      LAYER=2140000000
      JZZ=2140000000
      ISWITCH =2140000000
      JXYOUT=2140000000
      NINDKL=2140000000
      QZZKL=.FALSE.
      BENT=.FALSE.

c==============================================================
      DO 88 LY=1,NYI
c==============================================================
      RLY=LY
C Initialize a new track:
      IF(JOPTION(LFEFILE))THEN
C get T0 from the energy spectrum, recalculate dE/dx_valence
        CALL TSPEC(T0,E0,Es,ICRATE)
        RENRGY=RENRGY+T0
        DO 205 KL=1,NLAY
        CALL SDEDX(DEDXVL(KL,1),DEDXVP(KL,1),T0,Z1,A1,
     +   AVZVAL(KL),CELL(KL)/TAPC(KL))
        DEVL(KL,1) = DEDXVL(KL,1) * 1.0E-6 * UD(KL)
        DEVP(KL,1) = DEDXVP(KL,1) * 1.0E-6 * UD(KL)
  205   CONTINUE
      ENDIF

C Initialize RANBUF
      NRAND=NBINF*NDELZ*4+4
      CALL RANMAR(RANBUF,NRAND)
      IRAND=1

      IF(JOPTION(LFCFILE))THEN
C caution T is double precision , use TT instead
  248   call getcoord(LUN(4),0, XX, YY, PX, PY, TT, ier)
        if(ier .eq. 1) then
          write(6,*)'coordinates file: short line!?'
          goto 248
        else if(ier.ne.0) then
          write(6,*)'coordinates file: premature end of record/file'
          QUIT= .TRUE.
          GOTO 8889
        endif
        PX0 = PX
        PY0 = PY
        PZ = SQRT(1.0 - PX*PX - PY*PY)
        PZ0 = PZ
        COSTHETA = PZ
        T= TT
      ELSE

        XX=(RLX-RANBUF(IRAND))*DXI
        YY=(RLY-RANBUF(IRAND+1))*DYI
        IRAND=IRAND+2
        T=T0
C get THETA, PHI from angles in ANG(3,K) of FLUX7.COM
        CALL THETAPHI(THETA, PHI, K)
        IF ( PHI .GT. 360.01 ) THEN
          PHI = RANBUF(IRAND) * 360.0
          IRAND=IRAND+1
        ENDIF
C IF PHI>360 DEGREES WAS GIVEN, USE RANDOM PHI.
        IF ( THETA .GT. 360.01 ) THETA = 0.0
C IF THETA>360 DEGREES WAS GIVEN, randomize x,y, and set theta=0
        THETA =    THETA*DEGREE
        PHI   =    PHI*DEGREE
C (PX, PY, PZ) is unit vector in the direction of the current velocity:
        SINTHETA = SIN(THETA)
        COSTHETA = COS(THETA)
        PX = SINTHETA * COS ( PHI )
        PY = SINTHETA * SIN ( PHI )
        PZ = COSTHETA

        PX0 = PX
        PY0 = PY
        PZ0 = PZ
C BEAM DIVERGENCE
**
        UNIRAN(1)=RANBUF(IRAND)
        UNIRAN(2)=RANBUF(IRAND+1)
        PXX = UNIRAN(1)-0.4
        PYY = UNIRAN(2)-0.5    
        IRAND=IRAND+2
        PXX = PXX * XSDR 
        PYY = PYY * YSDR
        PHI0= PHI + ANGR
        PX  = PX + PXX * COS(PHI0) - PYY * SIN(PHI0)
        PY  = PY + PXX * SIN(PHI0) + PYY * COS(PHI0)

      ENDIF
      ZZ=ZZ0
      INVERTX = +1
      INVERTY = +1
      TRUEX=XX
      TRUEY=YY
      IF (JOPTION(LFXYOUT)) THEN
        JXYOUT=IXYOUT
        CALL EXPORTXYZ(LOGLU,ZERO, ZERO,ZERO, ZERO,ZERO)
        CALL EXPORTXYZ(LOGLU,TRUEX,TRUEY,ZERO, PX,PY)
      ENDIF
      PZ  = SQRT(1.0 - PX*PX - PY*PY)
**TEST check if PZ is a positive number
*     CALL TEST(PZ,'PZ1')
**
C
C COSTHETA in ANORM is angle of axis with average beam direction,
C do not 'correct' for beam divergence.
      DO 257 KL=1,NLAY
      DO 256 IL = 1,NLA(KL)
      ANORM(IL,KL) = ANORM0(IL,KL) * COSTHETA
 256  CONTINUE
      ZINTFSIG(KL)=WINTERF(KL)/SQRT(ALOG(256.0))
 257  CONTINUE
      JZ = 0
**TEST output PX,PY,PZ
*      WRITE(6,691)'------',PX,PY,PZ, XX,YY,JZ
**

      STARTX=TRUEX
      STARTY=TRUEY
      STARTPX=(PX-PXA)*PFACTOR
      STARTPY=(PY-PYA)*PFACTOR
      STARTT=REAL(T-T0AV)*TFACTOR

      RMAXNZ= NBINF*NDELZ

      QUIT= .FALSE.
      DO 120 NZ=0,NBINF*NDELZ-1
        IXIYZ(NZ)=0
  120 CONTINUE
C
C From here on, KL is the kind of layer (1 or 2)
C KLTRUE numbers the layers from 1 upwards
      KL=0
      KLTRUE=0
      NZINTERF = 0
C
c==============================================================
      DO 123 NZ=0,NBINF*NDELZ-1
c==============================================================
C
C Here starts the loop where a track is followed from one collision to the next.
C
C select the interface depth
C set KL and NZINTERF, first time through, and when crossing interface.
C Note that KL=0 first time through, so KL now set to 1.
      IF(NZ .EQ. NZINTERF) THEN
        LAYER= KL+1
        IF (LAYER.gt.NLAY)LAYER=1

        BENT= (CURRADIUS(LAYER) .NE. 0.0)
C curved crystal: rotate (PX0,PY0,PZ0) over the total bend angle
        IF(BENT)THEN
         PXN=CURMATT(1,1,LAYER)*PX0+ CURMATT(1,2,LAYER)*PY0+
     +      CURMATT(1,3,LAYER)*PZ0
         PYN=CURMATT(2,1,LAYER)*PX0+ CURMATT(2,2,LAYER)*PY0+
     +      CURMATT(2,3,LAYER)*PZ0
         PZN=CURMATT(3,1,LAYER)*PX0+ CURMATT(3,2,LAYER)*PY0+
     +      CURMATT(3,3,LAYER)*PZ0
         PX0 = PXN
         PY0 = PYN
         PZ0 = PZN
        ENDIF

        if(KL.ne.0)then

C crossing of layer boundary. Apply translation and rotation
          DELTAX= TRANSL(1,KL)
          DELTAY= TRANSL(2,KL)
          IDELTAZ= NINT(TRANSL(3,KL)/UD(KL))
C set PX to  true PX for a while
          IF(INVERTX .LT. 0) THEN
             DELTAX= -DELTAX
             PX= -PX
          ENDIF
          IF(INVERTY .LT. 0) THEN
             DELTAY= -DELTAY
             PY= -PY
          ENDIF
          PXN = ROTM(1,1,KL)*PX + ROTM(1,2,KL)*PY + ROTM(1,3,KL)*PZ
          PYN = ROTM(2,1,KL)*PX + ROTM(2,2,KL)*PY + ROTM(2,3,KL)*PZ
          PZN = ROTM(3,1,KL)*PX + ROTM(3,2,KL)*PY + ROTM(3,3,KL)*PZ
          XXN = XX + DELTAX
          YYN = YY + DELTAY
          JZN = JZ + IDELTAZ
          TRUEX = TRUEX + DELTAX*INVERTX
          TRUEY = TRUEY + DELTAY*INVERTY
          IF (JOPTION(LFXYOUT)) THEN
            CALL EXPORTXYZ(LOGLU,TRUEX,TRUEY,ZZ-ZZ0, PX,PY)
            CALL EXPORTXYZ(LOGLU,TRUEX,TRUEY,ZZ-ZZ0, PXN,PYN)
            JXYOUT=IXYOUT
          ENDIF
**TEST test transformation
*      WRITE(6,691)'before',PX,PY,PZ,XX,YY,JZ
*      WRITE(6,691)'after ',PXN,PYN,PZN,XXN,YYN,JZN
*691   FORMAT(A6,'  PX,PY,PZ=',3(F9.4,1X),' XX,YY=',2F9.3,' JZ=',I2)
**
          IF(INVERTX .LT. 0)PXN= -PXN
          IF(INVERTY .LT. 0)PYN= -PYN
          PX  = PXN
          PY  = PYN
          PZ  = PZN
C now PX is again INVERTX * (true PX)
          XX= XXN
          YY= YYN
          JZ= JZN

C Also rotate (PX0,PY0,PZ0) over the rotation angle
          PXN=ROTM(1,1,KL)*PX0+ ROTM(1,2,KL)*PY0+ ROTM(1,3,KL)*PZ0
          PYN=ROTM(2,1,KL)*PX0+ ROTM(2,2,KL)*PY0+ ROTM(2,3,KL)*PZ0
          PZN=ROTM(3,1,KL)*PX0+ ROTM(3,2,KL)*PY0+ ROTM(3,3,KL)*PZ0
          PX0 = PXN
          PY0 = PYN
          PZ0 = PZN
        endif  !!(KL.ne.0)

        KL=KL+1
        IF(KL.GT.NLAY)KL=1
        KLTRUE= KLTRUE+1
        IF(ZINTFSIG(KL) .GT. 0)THEN
          IF(IRAND .gt. NRAND-1) THEN
            CALL RANMAR(RANBUF,NRAND)
            IRAND=1
          ENDIF
          UNIRAN(1)=RANBUF(IRAND)
          UNIRAN(2)=RANBUF(IRAND+1)
          CALL PNORM(ZW, SPARE, ZINTFSIG(KL), UNIRAN)
          IRAND=IRAND+2
        ELSE
          ZW=0
        ENDIF
C get the thickness in angstrom of this layer
        if(NTHICK .gt. 0)then
          if(KLTRUE.le.NTHICK)then
            THICKNESS= THICK(KLTRUE)
          else
            write(LOGLU,*)'error: list of thicknesses exhausted'
            stop 'error: list of thicknesses exhausted'
          endif
        else
          THICKNESS= ZINTERF(KL)
        endif
        DNZINTERF=(THICKNESS+ZW)/UD(KL)
        IF(DNZINTERF.LT.1.0)DNZINTERF=1.0
        IF(DNZINTERF .GT. RMAXNZ)THEN
          NZINTERF=INT(RMAXNZ)
        ELSE
          NZINTERF = NZINTERF + INT(DNZINTERF)
        ENDIF

**TEST output of current layer selection
*     WRITE(6,'(A,2I9)')'New layer. KL and NZINTERF are',
*    + KL,NZINTERF
**
C Set switch ISWITCH according to lattice type
        ISWITCH= LATTICE(KL)

        NINDKL=NIND(KL)
        MODZKL=MODZ(KL)
        XUNKL=XUN(KL)
        YUNKL=YUN(KL)
        UDKL=UD(KL)
        ZUNKL=UD(KL)*MODZ(KL)
        XUMKL=XUM(KL)
        YUMKL=YUM(KL)
        TWOXUNKL=TWOXUN(KL)
        TWOYUNKL=TWOYUN(KL)
        QZZKL=(NLZ(KL).ne.0)
C       true if flux4 case, ZZ steps of varying length
C END IF(NZ .EQ. NZINTERF) THEN
      ENDIF

      IF(RANDOM) THEN
        IF(IRAND .gt. NRAND-1) THEN
          CALL RANMAR(RANBUF,NRAND)
          IRAND=1
        ENDIF
        XX = (RANBUF(IRAND  ) - 0.5)*XUNKL + XX
        YY = (RANBUF(IRAND+1) - 0.5)*YUNKL + YY
        IRAND = IRAND+2
      ENDIF

**TEST save old coordinates, see array bound test
       XSAVE=XX
       YSAVE=YY
       JZSAV=JZ
       PYSAVE=PY
       PXSAVE=PX
**
      IZ = NZ / NDELZ
C  MAP XX,YY,JZ,PX,PY TO UNIT CELL
C BRING XX INTO THE RANGE [0..MODX*XUN] BY TRANSLATION
      IXUM = INT(XX/ ( XUMKL ))
      IF ( XX.LT.0 ) IXUM = IXUM - 1
      XX = XX - IXUM*XUMKL
C SAME FOR YY
      IYUM = INT(YY/ ( YUMKL ))
      IF ( YY.LT.0 ) IYUM = IYUM - 1
      YY = YY - IYUM*YUMKL
C if IXUM or IYUM is non-zero the ion has left the cell.
C Used in SYMMETRY.NEW
C preset DELZZ with some big number to flag use without being set
      DELZZ=1.0E38

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SELLATTICE: SELECT CASE (ISWITCH)
C
C INTERNAL ROUTINES FOR MAPPING COORDINATES ACCORDING TO LATTICE TYPE

#include "SYMMETRY.NEW"

      CASE DEFAULT
        STOP 'Unknown lattice type'
      END SELECT SELLATTICE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      GOTO 7

77    CONTINUE
      IF(DELZZ.EQ.0.0) GOTO 7
      IF(DELZZ.gt. 1.0E10) STOP 'at label 77: DELZZ undefined'
      J=INT(DELZZ/(ZUNKL))
      IF(DELZZ.LT.0.0) J=J-1
      DELZZ=DELZZ-J*ZUNKL
      DZZ=0.0
C      Choose the nearest plane from the new coordinate
      DO 119 I=1,MODZKL
      JZZ=JZ+I
      IF(JZZ.GT.MODZKL) JZZ=JZZ-MODZKL
      DZZ=DUZ(JZZ-1,KL) + DZZ
      IF(DZZ.GE.DELZZ) GOTO 122
  119 CONTINUE
  122 CONTINUE
C
      DIFSUP=DZZ-DELZZ
      DIFINF=DELZZ-(DZZ-DUZ(JZZ-1,KL))
      IF(DIFSUP.LT.DIFINF) THEN
      JZ=JZZ
      ELSE
      JZ=JZZ-1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  7      CONTINUE
C BRING JZ IN THE RANGE [0..MODZKL-1]
      JZ = JZ - JZ/MODZKL*MODZKL
      IF ( JZ.LT.0 ) JZ = JZ + MODZKL
**TEST
*     IF(XX.ne.XSAVE .or. YY.ne.YSAVE)
*    +  WRITE(6,452) XSAVE/XUNKL,YSAVE/YUNKL,JZSAV,XX/XUNKL,YY/YUNKL,
*    + (UZ(JZ,KL)-UZ(JZSAV,KL))/ZUNKL,JZ
* 452 FORMAT(' BEFORE MAPPING:',2F8.3,I4,'  AFTER:',3F8.3,I4)
**

C SET ARRAY INDICES IX,IY ACCORDING TO CURRENT POSITION
      IX = INT(XX/XUNKL*RNX) + 1
      IY = INT(YY/YUNKL*RNY) + 1

**TEST: uncomment the following to get range checking.
*      IF(IX.LE.0 .OR. IX.GT.NFX.OR. IY.LE.0 .OR. IY.GT.NFY.OR.
*     +   JZ.LT.0 .OR. JZ.GE.MODZKL)THEN
*      WRITE(6,1006) IX,IY,JZ,JZSAV,DELX
*     +  XSAVE/XUNKL,YSAVE/YUNKL,
*     +  XX/XUNKL,YY/YUNKL
* 1006  FORMAT(
*     +  'ILLEGAL INDEX! IX, IY, JZ, JZS,DELX, XS1, XS, YS, X, Y :'/
*     +  4(1X,I2), 6(1X,E14.6) )
*         IF(IX.LE.0)write(6,*)'ix <= 0'
*         IF(IX.GT.NFX)write(6,*)'ix > ',NFX
*         IF(IY.LE.0)write(6,*)'iy <= 0'
*         IF(IY.GT.NFY)write(6,*)'iy > ',NFY
*         IF(JZ.LT.0)write(6,*)'jz < 0'
*         IF(JZ.GE.MODZKL)write(6,*)'jz >=',MODZKL
*         IF(IX.LE.0)IX=1
*         IF(IX.GT.NFX)IX=NFX
*         IF(IY.LE.0)IY=1
*         IF(IY.GT.NFY)IY=NFY
*         IF(JZ.LT.0)JZ=0
*         IF(JZ.GE.MODZKL)JZ=0
*         IXIY=1
*      ENDIF
**

      IXIY=IX + (IY-1)*NFX
      TINV=1.0/T
*      DPERR= DPTEST('TINV',TINV)
C
C Update FLUX array:
      IF(KL.EQ.LFLX) THEN
        IXIYZ(NZ)=IXIY
      SELCROSS: SELECT CASE (IMPS)
C One of the following is chosen for option FLUX:
        CASE (500)
          DELFLUX(NZ) = PROFIEL(IZ) * (REAL(TINV)*T0)**2
        CASE (501)
          DELFLUX(NZ) = PROFIEL(IZ)
        CASE (502)
          DELFLUX(NZ) = PROFIEL(IZ) *
     +    CROSSECT( REAL(T), NCROSS(0,1), ECROSS(1,0,1), CROSS(1,0,1) )
        CASE (503)
          DELFLUX(NZ) = (REAL(TINV)*T0)**2
        CASE (504)
          DELFLUX(NZ) = 1.0
        CASE (505)
          DELFLUX(NZ) =
     +    CROSSECT( REAL(T), NCROSS(0,1), ECROSS(1,0,1), CROSS(1,0,1) )
      END SELECT SELCROSS
      ENDIF
C
C LOOP THROUGH THE POSSIBLE ATOMS TO COLLIDE WITH
      DECO=0.0
      DENU=0.0
c==============================================================
      DO 6 IL=1,NLS(KL)
c==============================================================
        INDX=IATOM(JZ,IL,KL)
        IF(JOPTION(LMIXED))THEN
C set IL1 to the index of the randomly selected atom type
C with weight ATOMFRAC(IL1,KL)
          IL1= ILMAP(IL,KL)
          IF(IRAND .gt. NRAND) THEN
            CALL RANMAR(RANBUF,NRAND)
            IRAND=1
          ENDIF
          CHANCE= RANBUF(IRAND)
          IRAND=IRAND+1
          do 27 IL2= 1,NLA1(IL,KL)
            if (CHANCE .LT. ATOMFRAC(IL1,KL) ) then
              goto 28
            else
              CHANCE = CHANCE - ATOMFRAC(IL1,KL)
              IL1=1+IL1
            endif
  27      continue
  28      continue
          HITLIST(IL2,IL,KL)= 1+HITLIST(IL2,IL,KL)
C vacancy when the random nr was larger than the single fraction:
          if(IL2 .gt. NLA1(IL,KL)) goto 6
        ELSE
          IL1= ILMAP(IL,KL)
        ENDIF
        LIND=0
  726   CONTINUE
        IF( INDX .GT. 0) THEN
          LIND= LIND+1
          INDXT= (INDX-1)/NINDKL
          INDIC(LIND)= INDX - (INDXT*NINDKL)
          INDX= INDXT
          GOTO 726
        ENDIF
c==============================================================
        DO 599 JIND=1, LIND
c==============================================================
        INDX= INDIC(JIND)
C UPDATE CLOSE ENCOUNTER YIELD:
        XT= XATOM(INDX,KL) - XX
        YT= YATOM(INDX,KL) - YY
**TEST output
*       write(6,7501) JZ, IL1, JIND, INDX,
*    +   A2(IL1,KL), XATOM(INDX,KL)/XUNKL, YATOM(INDX,KL)/YUNKL
*7501   format(1X, 'jz,il1,jind,indx',4I7, 3F13.4)
** end of test output
        RP= XT*PX + YT*PY
C sqrt(S0) is the collision parameter for an atom in equilibrium position
C now calculated as perpendicular distance, not, as before, in XY plane:
        S0= ( XT*XT + YT*YT -RP*RP ) / U2SQ(IL1,KL)
        IF (S0.GT.60.0) THEN
            CLENCY=0.0
        ELSE
            CLENCY=ANORM(IL1,KL) * EXP ( -S0)
        ENDIF
        IF (NCROSS(IL1,KL).NE.0) THEN
          IF (NCROSS(IL1,KL).LT.0) THEN
C Rutherford cross section
            CLENCY=CLENCY * (TINV*T0)**2
          ELSE
C Interpolated cross section
            CLENCY=CLENCY *
     +      CROSSECT(REAL(T), NCROSS(IL1,KL),
     +               ECROSS(1,IL1,KL), CROSS(1,IL1,KL))
          ENDIF
        ENDIF
*        DPERR= DPTEST('CLENCY', CLENCY)
        DELRAAK(IZ,IL1,KL) = DELRAAK(IZ,IL1,KL) + CLENCY
*        DPERR= DPTEST('DELRAAK(IZ,IL1,KL)', DELRAAK(IZ,IL1,KL))
C sample random displacement (X0,Y0,ZTH) of thermally vibrating atom.
        IF(IRAND .gt. NRAND-3) THEN
          CALL RANMAR(RANBUF,NRAND)
          IRAND=1
        ENDIF
        UNIRAN(1)=RANBUF(IRAND)
        UNIRAN(2)=RANBUF(IRAND+1)
        CALL PNORM(X0,Y0,U1(IL1,KL),UNIRAN)
        UNIRAN(1)=RANBUF(IRAND+2)
        UNIRAN(2)=RANBUF(IRAND+3)
        CALL PNORM(ZTH,SPARE,U1(IL1,KL),UNIRAN)
        IRAND= IRAND+4
        XTH= XATOM(INDX,KL) + X0 - XX
        YTH= YATOM(INDX,KL) + Y0 - YY
C (XTH,YTH,ZTH) = (radius vector of atom) - (XX,YY,0)
C (QX,QY,QZ)    = (radius vector of scattering center) - (XX,YY,0)
C with length = inner poduct of (XTH,YTH,ZTH) with (PX,PY,PZ)
        QTH= XTH*PX + YTH*PY + ZTH*PZ
        QX= QTH*PX
        QY= QTH*PY
        QZ= QTH*PZ
Collision parameter:
        SX= QX - XTH
        SY= QY - YTH
        SZ= QZ - ZTH
        S = SQRT ( SX*SX + SY*SY + SZ*SZ)
C
C Scattering angle in lab system (impulse approximation):
        IZ2= NINT(Z2(IL1,KL))
        TET= 0.5 * ZZA2(IL1,KL) * FF(S/RSCREEN(IL1,KL),IZ2) * REAL(TINV)
        SINTET= SIN(TET)
        COSTET= COS(TET)
C small angle approximation only if angle is small! aug. 2018        
        IF (ABS(TET) .LT. 0.1)THEN
          DENU= DENU + RM(IL1,KL)*TET*TET*REAL(T)
        ELSE
          DENU= DENU + RM(IL1,KL)*(1.0-COSTET)*REAL(T)
        ENDIF
        PX= COSTET*PX + SINTET*SX/S
        PY= COSTET*PY + SINTET*SY/S
        PZ= COSTET*PZ + SINTET*SZ/S
**TEST
*      WRITE(6,691)'scattr',PX,PY,PZ, XX,YY,JZ
**
        IF(PZ .LT. 0)THEN
!$OMP CRITICAL
          NRAAK = NRAAK + 1
          IF ( NRAAK .LE. MAXRAAK  ) THEN
              KRAAK ( NRAAK ) = NZ
          ENDIF
!$OMP END CRITICAL
          WRITE(6,630)NTRACK,'Backscattering!',
     +      T*1.0E3,ZZ-ZZ0,ARCCOS(PZ)/DEGREE
C  LEAVE THIS TRACK:
          QUIT= .TRUE.
          GOTO 8889
        ENDIF
**TEST test if PZ is in the range of floating point numbers
*      CALL TEST(PZ,'PZ2')
**
C z component of vector (QX, QY, QZ) does not change:
        QTH = QZ / PZ
C intersection with z=0 plane of updated ion path:
        DELX = QX - QTH*PX
        DELY = QY - QTH*PY
        XX = XX + DELX
        YY = YY + DELY
        TRUEX = TRUEX + DELX*INVERTX
        TRUEY = TRUEY + DELY*INVERTY
C Energy loss due to CORE electrons
        if(NSTOP.gt.1)then
          IS=1
  982     IF(ESTOP(IS+1).LT.T)goto 983
          if(IS.lt.NSTOP-1)then
            IS=IS+1
            goto 982
          endif
  983     continue
          FRAC2= ESTOP(IS)-REAL(T)
          IB=INT(S*25.0)
          IF(IB.LT.50)THEN
              IF(IB.EQ.0)IB=1
              FRAC=S*25.0-FLOAT(IB)
              C1 = EXP((1.0-FRAC)*ELCORE(IB,IL1,KL,IS) +
     +                 FRAC*ELCORE(IB+1,IL1,KL,IS))
              C2 = EXP((1.0-FRAC)*ELCORE(IB,IL1,KL,IS+1) +
     +                 FRAC*ELCORE(IB+1,IL1,KL,IS+1))
              DECO = DECO + (1.0-FRAC2)*C1 + FRAC2*C2
          ENDIF
          DEVLKL = (1.0-FRAC2)*DEVL(KL,IS) + FRAC2*DEVL(KL,IS+1)
          DEVPKL = (1.0-FRAC2)*DEVP(KL,IS) + FRAC2*DEVP(KL,IS+1)
        else
          IB=INT(S*25.0)
          IF(IB.LT.50)THEN
              IF(IB.EQ.0)IB=1
              FRAC=S*25.0-FLOAT(IB)
              C1 = EXP((1.0-FRAC)*ELCORE(IB,IL1,KL,1) +
     +                 FRAC*ELCORE(IB+1,IL1,KL,1))
              DECO = DECO + C1
          ENDIF
          DEVLKL=DEVL(KL,1)
          DEVPKL=DEVP(KL,1)
        endif

Collect current energy loss:

C for WEIGHTED energy loss (Backscattering), use the following:
        IF(WEIGHTDE .GE. 0 .and. CLENCY.gt. 0.0)THEN

          DELCEYDPR(IZ,IL1) = DELCEYDPR(IZ,IL1)  + CLENCY
          DELDAVDPR(IZ,IL1) = DELDAVDPR(IZ,IL1)  + (T0-T)*CLENCY
          DELDDAVDPR(IZ,IL1)= DELDDAVDPR(IZ,IL1) + ((T0-T)**2)*CLENCY
*          DPERR= DPTEST('DELCEYDPR(IZ,IL1) ',DELCEYDPR(IZ,IL1) )
*          DPERR= DPTEST('DELDAVDPR(IZ,IL1)',DELDAVDPR(IZ,IL1))
*          DPERR= DPTEST('DELDDAVDPR(IZ,IL1)',DELDDAVDPR(IZ,IL1))

        ENDIF

C The UNWEIGHTED version follows outside the current loop (loop 6)
c==============================================================
  599 CONTINUE
c==============================================================
  6   CONTINUE
c==============================================================

      RWBEAM= PX*PX0*REAL(INVERTX) + PY*PY0*REAL(INVERTY) + PZ*PZ0
C for UNWEIGHTED energy loss (Transmission), use the following:
      IF(.not.JOPTION(LFCOLLI) .or. RWBEAM.gt.CMAWBEAM)THEN
        DELLEFTOVER(IZ)= 1+DELLEFTOVER(IZ)
        IF(WEIGHTDE .LT. 0)THEN
           DELDAVDPR(IZ,1)= (T0-T)    + DELDAVDPR(IZ,1)
          DELDDAVDPR(IZ,1)= (T0-T)**2 +DELDDAVDPR(IZ,1)
        ENDIF
      ENDIF
C
      IF(ABS(PX) .LT. CRITX .AND. ABS(PY)  .LT. CRITY)THEN
        DELNCHANNEL(IZ)= 1+DELNCHANNEL(IZ)
      ENDIF

C angular deviation due to rest field of surrounding atom rows.
C (This becomes nonsense for large deviations from axial direction,
C  Let's hope this whole effect is negligible)
      if(QZZKL)then
        TEX = FFX(IXIY,KL) * DUZ(JZ,KL) * REAL(TINV) * 0.5
        TEY = FFY(IXIY,KL) * DUZ(JZ,KL) * REAL(TINV) * 0.5
      else
        TEX = FFX(IXIY,KL) * UDKL * REAL(TINV) * 0.5
        TEY = FFY(IXIY,KL) * UDKL * REAL(TINV) * 0.5
      endif

C Energy loss due to local valence electrons, DEVL
C Energy loss due to distant valence electrons, DEVP
Choose random scattering angle from gaussian:
C Note: make sure arg 3 of ENORM is single precision!
      FDEVL= (0.5*AME/(A1*AMU))*DEVLKL
      IF(IRAND .gt. NRAND-1) THEN
        CALL RANMAR(RANBUF,NRAND)
        IRAND=1
      ENDIF
      UNIRAN(1)=RANBUF(IRAND)
      UNIRAN(2)=RANBUF(IRAND+1)
      CALL PNORM(TEE1,TEE2,SQRT(REAL(FDEVL/PZ*TINV)),UNIRAN)
      IRAND=IRAND+2

C Update direction vector (PX,PY,PZ), energy and position
C The following is not free from small angle approximations
C Construct orthogonal unit vectors, S and T, perpendicular to P
      PYZ = SQRT(PY*PY + PZ*PZ)
C     SX = 0
      SY = PZ/PYZ
      SZ = -PY/PYZ
      TX = PY*SZ - PZ*SY
      TY = - PX*SZ
C     TY = PZ*SX - PX*SZ
C     TZ = PX*SY - PY*SX
C
      PX = PX + TEX*PZ           + TEE2*TX
      PY = PY + TEY*PZ + TEE1*SY + TEE2*TY
C
C curved crystal correction:
      LAYER2=LAYER
      IF(BENT.and. LAYER.le.NLAY)THEN
        IF(INVERTX .LT. 0) PX= -PX
        IF(INVERTY .LT. 0) PY= -PY
        PXN=CURMAT(1,1,LAYER)*PX + CURMAT(1,2,LAYER)*PY +
     +      CURMAT(1,3,LAYER)*PZ
        PYN=CURMAT(2,1,LAYER)*PX + CURMAT(2,2,LAYER)*PY +
     +      CURMAT(2,3,LAYER)*PZ
        PZN=CURMAT(3,1,LAYER)*PX + CURMAT(3,2,LAYER)*PY +
     +      CURMAT(3,3,LAYER)*PZ
        IF(INVERTX .LT. 0)PXN= -PXN
        IF(INVERTY .LT. 0)PYN= -PYN
        PX  = PXN
        PY  = PYN
        PZ  = PZN
      ELSE IF(BENT)THEN
C the following should be impossible but happened to me anyways
C at this point LAYER was 1 (OK) but LAYER2 was 3 (!?)
C these variables were equated a few statements ago!
C This mystery has to do with OMP
          write(6,*) 'error: LAYER was too big!??', LAYER, LAYER2, NLAY
      ENDIF
      PZSQ = 1.0 - PX*PX -PY*PY
C
**TEST
*     if(PZSQ.ge.0)then
*       PZ = SQRT(PZSQ)
*     else
*       PZ = SQRT(-PZSQ)
*     endif
*     write(6,'(A,4F12.7)')'tex,tey,tee1,tee2:',TEX,TEY,TEE1,TEE2
*      WRITE(6,691)'elos',PX,PY,PZ, XX,YY,JZ
**
C Test if angle between velocity and string became too large
C Save z coordinate where this occurred for later treatment.
      IF ( PZSQ .LT. TESQMAX ) THEN
!$OMP CRITICAL
        NRAAK = NRAAK + 1
        IF ( NRAAK .LE. MAXRAAK  ) THEN
            KRAAK ( NRAAK ) = NZ
        ENDIF
!$OMP END CRITICAL
        WRITE(6,630)NTRACK,'theta too large',
     +   T*1.0E3,ZZ-ZZ0,ARCCOS(PZ)/DEGREE
  630     FORMAT(I7,A18,' T=',F9.1,'keV  z=',F10.0,
     +     ' theta=',F9.1)
C         LEAVE THIS TRACK:
          QUIT=.TRUE.
        GOTO 8889
      ENDIF
      IF( JOPTION(LFMAXAN) ) THEN
        RWBEAM= PX*PX0*FLOAT(INVERTX) + PY*PY0*FLOAT(INVERTY) + PZ*PZ0
        IF(RWBEAM.lt.CMAWBEAM)THEN
!$OMP CRITICAL
          NRAAK = NRAAK + 1
          IF ( NRAAK .LE. MAXRAAK  ) THEN
              KRAAK ( NRAAK ) = NZ
          ENDIF
!$OMP END CRITICAL
          WRITE(6,630)NTRACK,'angle with beam too large',
     +      T*1.0E3,ZZ-ZZ0,ARCCOS(RWBEAM)/DEGREE
C         LEAVE THIS TRACK:
          QUIT=.TRUE.
          GOTO 8889
        ENDIF
      ENDIF

      PZ = SQRT(PZSQ)
**TEST test if PZ is in the range of floating point numbers
*     CALL TEST(PZ,'PZ3')
**
      if(QZZKL)then
        DELX = PX * DUZ(JZ,KL) / PZ
        DELY = PY * DUZ(JZ,KL) / PZ
      else
        DELX = PX * UDKL / PZ
        DELY = PY * UDKL / PZ
      endif
      XX = XX + DELX
      YY = YY + DELY
      TRUEX = TRUEX + DELX*INVERTX
      TRUEY = TRUEY + DELY*INVERTY
      IF (JOPTION(LFXYOUT)) THEN
        JXYOUT=JXYOUT-1
        IF (JXYOUT .LE. 0) THEN
          TRUEPX= PX*INVERTX
          TRUEPY= PY*INVERTY
          CALL EXPORTXYZ(LOGLU,TRUEX,TRUEY,ZZ-ZZ0, TRUEPX,TRUEPY)
          JXYOUT=IXYOUT
        ENDIF
      ENDIF

C LOSS.GT.0 IS THE "ION BEAM VERSION" WHERE AN INCOMING PARTICLE LOSES
C  ENERGY BY COLLISIONS WITH ELECTRONS AND NUCLEI
C LOSS.LT.0 IS THE "EMISSION CHANNELING VERSION" WHERE AN INCOMING
C  PARTICLE GAINS ENERGY BY COLLISIONS WITH ELECTRONS AND NUCLEI
C IF LOSS.EQ.0 THE PARTICLE ONLY LOSES ENERGY BY COLLISIONS WITH NUCLEI
C  COLLISIONS WITH ELECTRONS STILL DEFLECT THE PARTICLE
      IF (LOSS.GT.0) THEN
        T = T - DENU - DECO - (DEVLKL + DEVPKL)/PZ
      ELSE IF (LOSS.LT.0) THEN
        T = T + DENU + DECO + (DEVLKL + DEVPKL)/PZ
      ELSE
        T = T - DENU
      ENDIF
      if(JOPTION(LFSCORE))then
       SNUCLEAR(IZ)= SNUCLEAR(IZ) + DENU
       SCORE(IZ)= SCORE(IZ) + DECO
       SVALENCE(IZ) = SVALENCE(IZ) + (DEVLKL + DEVPKL) / PZ
      endif

C TEST IF ENERGY TOO SMALL
      IF ( T .LT. TMIN ) THEN
!$OMP CRITICAL
        NRAAK = NRAAK + 1
        IF ( NRAAK .LE. MAXRAAK  ) THEN
            KRAAK ( NRAAK ) = NZ
        ENDIF
!$OMP END CRITICAL
        WRITE(6,630)NTRACK,'energy too small',
     +    T*1.0E3,ZZ-ZZ0,ARCCOS(PZ)/DEGREE
C         LEAVE THIS TRACK:
        QUIT=.TRUE.
        GOTO 8889
      ENDIF
**TEST test if T is in the range of floating point numbers
*      DPERR= DPTEST('T',T)
**
C
C
      if(QZZKL)then
        ZZ    =  ZZ  +  DUZ(JZ,KL)
      else
        ZZ    =  ZZ  +  UDKL
      endif
      JZ    =  JZ  +  1
C
      IF(QUIT)GOTO 8889
C =====================================================================
  123 CONTINUE
C end of this track
C =====================================================================
C determine angle with the beam direction:
      IF(INVERTX .LT. 0)PX= -PX
      IF(INVERTY .LT. 0)PY= -PY

      NTRACK= NTRACK+1
      if(JOPTION(LSTARTCO))then
       write(DATUM,891)'start',NTRACK
       call wrt5(LOGLU,DATUM,
     +           STARTX,STARTY,STARTPX,STARTPY,STARTT)
      endif
      if(JOPTION(LEXITCO))then
       write(DATUM,891)'exit ',NTRACK
       call wrt5(LOGLU,DATUM,TRUEX-XAV,TRUEY-YAV,
     +            (PX-PXA)*PFACTOR,(PY-PYA)*PFACTOR,
     +            REAL(T0AV-T)*TFACTOR)
  891   FORMAT(A5,I6,':')
      endif

      IF(JOPTION(LFCOLLI))THEN
        RWBEAM = PX*PX0 + PY*PY0 +PZ*PZ0
**TEST test output
*       WRITE(6,778) ARCCOS(RWBEAM)*DEGREE
* 778 FORMAT(' Final angle with beam:',F9.2)
**
        IF (RWBEAM .LT. CMAWBEAM) THEN
!$OMP CRITICAL
          NAWBX = NAWBX+1
!$OMP END CRITICAL
          WRITE(6,631)NTRACK,'maximum exit angle exceeded',
     +    T*1.0E3,ZZ-ZZ0,ARCCOS(PZ)/DEGREE
  631     FORMAT(I7,A28,' T=',F7.1,'keV  z=',F10.0,' theta=',F7.1)
          GOTO 8889
        ENDIF
      ENDIF
      IF (JOPTION(LFEFINAL))THEN
          NEF = INT((T - TFMIN) / EBIN)
          IF(NEF .LE. 0) NEF=0
          IF(NEF .GT. NEBINA) NEF=NEBINA
!$OMP CRITICAL
          EFINAL(NEF)= 1 + EFINAL(NEF)
!$OMP END CRITICAL
      ENDIF
!$OMP CRITICAL
      NFINAL = 1 + NFINAL
      TFINAL = TFINAL + T
      DTFINAL = DTFINAL + T*T
      PXFINAL = PXFINAL + PX
      DPXFINAL = DPXFINAL + PX*PX
      PYFINAL = PYFINAL + PY
      DPYFINAL = DPYFINAL + PY*PY
      PTFINAL = PTFINAL + SQRT(PX*PX +PY*PY)
      DPTFINAL = DPTFINAL + PX*PX + PY*PY
*      DPERR= DPTEST('TFINAL',TFINAL)
*      DPERR= DPTEST('DTFINAL',DTFINAL)
*      DPERR= DPTEST('PXFINAL',PXFINAL)
*      DPERR= DPTEST('DPXFINAL',DPXFINAL)
*      DPERR= DPTEST('PYFINAL',PYFINAL)
*      DPERR= DPTEST('DPYFINAL',DPYFINAL)
*      DPERR= DPTEST('PTFINAL',PTFINAL)
*      DPERR= DPTEST('DPTFINAL',DPTFINAL)
!$OMP END CRITICAL

      if(JOPTION(LFVELO))then
        NPPX=(PX - PX0)/PBIN +NP0
        NPPY=(PY - PY0)/PBIN +NP0
        IF(NPPX .LT. 1) NPPX=1
        IF(NPPX .GT. NVA)NPPX=NVA
        IF(NPPY .LT. 1) NPPY=1
        IF(NPPY .GT. NVA)NPPY=NVA
        NPPXY=NPPX + (NPPY-1)*NVA
!$OMP CRITICAL
        VFINAL(NPPXY)= 1.0+VFINAL(NPPXY)
!$OMP END CRITICAL
      endif

 8889 CONTINUE
**TEST output PX,PY,PZ
C     WRITE(6,691)'at end',PX,PY,PZ, XX,YY,JZ
**
!$OMP CRITICAL
      DO 124 NZ=0,NBINF*NDELZ-1
      IXIY=IXIYZ(NZ)
      IF(IXIY.gt.0)THEN
        FLUX(IXIY)= FLUX(IXIY) + DELFLUX(NZ)
      ENDIF
  124 CONTINUE
!$OMP END CRITICAL
C =====================================================================
  88  CONTINUE
C end loop LY
C =====================================================================
!$OMP CRITICAL
      DO 127 IZ=0,NBINF-1
       LEFTOVER(IZ) = DELLEFTOVER(IZ) + LEFTOVER(IZ)
       NCHANNEL(IZ) = DELNCHANNEL(IZ) + NCHANNEL(IZ)
      DO 126 IL=1,MLA
*      DPERR= DPTEST('DAVDPR(IZ,IL)',DAVDPR(IZ,IL))
*      DPERR= DPTEST('DDAVDPR(IZ,IL)',DDAVDPR(IZ,IL))
*      DPERR= DPTEST('CEYDPR(IZ,IL)',CEYDPR(IZ,IL))
        CEYDPR(IZ,IL) = DELCEYDPR(IZ,IL)  + CEYDPR(IZ,IL)
        DAVDPR(IZ,IL) = DELDAVDPR(IZ,IL)  + DAVDPR(IZ,IL)
       DDAVDPR(IZ,IL) = DELDDAVDPR(IZ,IL) + DDAVDPR(IZ,IL)
*      DPERR= DPTEST('DAVDPR(IZ,IL)',DAVDPR(IZ,IL))
*      DPERR= DPTEST('DDAVDPR(IZ,IL)',DDAVDPR(IZ,IL))
*      DPERR= DPTEST('CEYDPR(IZ,IL)',CEYDPR(IZ,IL))
      DO 125  KL=1,NLAY
        RAAK(IZ,IL,KL)= RAAK(IZ,IL,KL) + REAL(DELRAAK(IZ,IL,KL))
  125 CONTINUE
  126 CONTINUE
  127 CONTINUE
!$OMP END CRITICAL
C =====================================================================
  8   CONTINUE
C end loop LX
C =====================================================================
C     All threads join master thread and disband
!$OMP END PARALLEL DO
!$OMP FLUSH
C
C the end of this simulation
      WRITE(6,830) NFINAL
      IF(NAWBX .GT. 0) THEN
        WRITE(6,831) NAWBX, AMAWBEAM
      ENDIF
  830 FORMAT(I9,' ions accepted')
  831 FORMAT(I9,' ions exited with angle with beam >', F7.3)
      END
C -------------------------------------------------- SUBROUTINE OUTP2
      SUBROUTINE OUTP2(K)
      IMPLICIT NONE
      INTEGER K

C OUTPUT RESULTS FOR ANGLE K.
#include "FLUX7.COM"

      INTEGER I,J,L, KL,IL, MRAAK, KR, IZ, NFX1, NFY1
      REAL FNORM, AVTH, DAVTH, SUM, SAMPLES, AV
      REAL AVRAAK(ML,NLAYER)
      REAL THETA,PHI,THIK,WID
      DOUBLE PRECISION DAVS, DDAVSQ, DPTEMP

      DOUBLE PRECISION  DAVDPR(0:NBIN-1,ML)
      DOUBLE PRECISION DDAVDPR(0:NBIN-1,ML)
      DOUBLE PRECISION  CEYDPR(0:NBIN-1,ML)
      COMMON /DOUBLEDAV/DAVDPR, DDAVDPR, CEYDPR

      LOGICAL DPERR
      LOGICAL DPTEST
      EXTERNAL DPTEST

      CHARACTER*160 FLKOP
      CHARACTER*20 TIME
      DATA FLKOP/'
     +
     +                                        '/
      DATA TIME/'                    '/

C GET CURRENT DATE AND TIME:
      CALL SYSTIM(TIME)

      AVTH=0
      DAVTH=0

C NORMALIZE FLUX TO AN AVERAGE OF UNITY
      IF(JOPTION(LFFLUX)) THEN
        SUM = 0.0
        DO 11 I= 1,NFX*NFY
  11    SUM = SUM + FLUX (I)
        SUM = SUM / FLOAT(NFX*NFY)
        DO 12 I= 1,NFX*NFY
  12    FLUX(I) = FLUX(I) / SUM
      ENDIF

C AVERAGE YIELD:
      DO 101 KL=1,NLAY
      DO 10 IL=1,NLA(KL)
      AV=0.0
      DO 9 J=0,NBINF-1
      AV=AV + RAAK(J,IL,KL)
      DPTEMP= RAAK(J,IL,KL)
      DPERR= DPTEST('RAAK',DPTEMP)
      IF(DPERR)write(6,*)'J,IL,KL, RAAK',J,IL,KL,RAAK(J,IL,KL)
  9   CONTINUE
      AVRAAK(IL,KL) = AV / NBINF
   10 CONTINUE

  101 CONTINUE

C NORMALIZE AVERAGE ENERGY LOSS AS A FUNCTION OF DISTANCE:
C Make it a non-negative quantity, even when LOSS < 0
C FOR UNWEIGHTED ENERGY LOSS USE THE FOLLOWING (TRANSMITTED IONS)
C
      IF(WEIGHTDE .LT. 0) THEN
        DO 17 J=0,NBINF-1
          SAMPLES=LEFTOVER(J)
          IF (SAMPLES .LE. 0) GOTO 17
          DAVS=1000.0*DAVDPR(J,1)/SAMPLES
          IF(DAVS.LT.0.0)DAVS= -DAVS
          DAV(J,1)=REAL(DAVS)
          DDAVSQ= 1.0E6* DDAVDPR(J,1)/SAMPLES - DAVS**2
          IF(DDAVSQ .LT. 0) DDAVSQ=0
          DDAV(J,1)=REAL(SQRT(DDAVSQ))
  17    CONTINUE
      ELSE
C FOR WEIGHTED ENERGY LOSS USE THE FOLLOWING (BACKSCATTERED IONS)
        DO 181 IL=1,MLA
        DO 183 J=0,NBINF-1
          DPERR= DPTEST('CEYDPR(J,IL) ',CEYDPR(J,IL) )
         IF(CEYDPR(J,IL) .gt. 1.0E-20)THEN
          DAVS=1000.0*DAVDPR(J,IL)/CEYDPR(J,IL)
          IF(DAVS.LT.0.0)DAVS= -DAVS
          DDAVSQ= 1.0E6* DDAVDPR(J,IL)/CEYDPR(J,IL) - DAVS**2
          IF(DDAVSQ .LT. 0) DDAVSQ=0
          DAV(J,IL)=REAL(DAVS)
          DDAV(J,IL)=REAL(SQRT(DDAVSQ))
          DPERR= DPTEST('CEYDPR(J,IL) ',CEYDPR(J,IL) )
          DPERR= DPTEST('DAVDPR(J,IL)',DAVDPR(J,IL))
          DPERR= DPTEST('DDAVDPR(J,IL)',DDAVDPR(J,IL))
          DPERR= DPTEST('DDAVSQ',DDAVSQ)
         ELSE
          DAV(J,IL)=0.0
          DDAV(J,IL)=0.0
         ENDIF
 183    CONTINUE
 181    CONTINUE

c [ version 7.7: questionable averaging procedures removed ]
      ENDIF

      IF(NFINAL .le.0)THEN
        do 1443 L=1,3
        IF(LUOPEN(L).ne.1)goto 1443
        write(LUN(L),*)'No ions left ...? NFINAL=',NFINAL
 1443   continue
      ELSE
      TFINAL=TFINAL/NFINAL
      DTFINAL= SQRT(DTFINAL/NFINAL - TFINAL*TFINAL)
      PXFINAL=PXFINAL/NFINAL
      DPXFINAL= SQRT(DPXFINAL/NFINAL - PXFINAL*PXFINAL)
      PYFINAL=PYFINAL/NFINAL
      DPYFINAL= SQRT(DPYFINAL/NFINAL - PYFINAL*PYFINAL)
      PTFINAL=PTFINAL/NFINAL
      AVTH = REAL(ASIN(PTFINAL))
      DPTFINAL= SQRT(DPTFINAL/NFINAL - PTFINAL*PTFINAL)
      DAVTH = REAL(DPTFINAL)/COS(AVTH)
      ENDIF
C OUTPUT TO UNIT LUN(1)..LUN(3)
      DO 30 L=1,3
      IF(LUOPEN(L).ne.1)goto 30
      CALL THETAPHI(THETA,PHI,K)
      WRITE(LUN(L),1)TIME,K,THETA,PHI,NRAAK
  1   FORMAT(/1X,A20,'ANGLE',I4,': THETA=',F7.3,' PHI=',F7.3,
     +I4,' LOST IONS')
      IF(ANG(3,K).NE.0.0) WRITE(LUN(L),1288)(ANG(J,K),J=1,3)
 1288 FORMAT('input angles:',3F7.3)
      MRAAK=NRAAK
      IF(MRAAK.GT.MAXRAAK)MRAAK=MAXRAAK
      IF(NRAAK.NE.0)WRITE(LUN(L),100)(KRAAK(KR),KR=1,MRAAK)
  100 FORMAT(' IONS ''LOST'' AT PLANE NUMBERS:'/,(7I6))
      WRITE(LUN(L),2)((IL,KL,AVRAAK(IL,KL),IL=1,NLA(KL)),KL=1,NLAY)
  2   FORMAT(' AVERAGE CLOSE ENCOUNTER PROBABILITY OF ATOMS',2I2,':',
     +         F8.4)
      WRITE(LUN(L),21)'ENERGY', TFINAL, DTFINAL
  21  FORMAT(' AVERAGE FINAL ',A7,F9.4,'  STANDARD DEV:', F9.4)
      WRITE(LUN(L),21)'ANGP(X)', PXFINAL/DEGREE, DPXFINAL/DEGREE
      WRITE(LUN(L),21)'ANGP(Y)', PYFINAL/DEGREE, DPYFINAL/DEGREE
      WRITE(LUN(L),21)'THETA  ', AVTH/DEGREE, DAVTH/DEGREE
      IF(JOPTION(LMIXED))THEN
        DO 1391 KL=1,NLAY
        WRITE(LUN(L),1390)'HITLIST LAYER',KL
        DO 1392 IL=1,NLS(KL)
        WRITE(LUN(L),1393)(HITLIST(J,IL,KL),J=1,1+NLA1(IL,KL))
 1392   CONTINUE
 1391   CONTINUE
 1390   FORMAT(1X,A,I2)
 1393   FORMAT(1X,8I9)
      ENDIF
      IF(JOPTION(LFCOLLI))THEN
        WRITE(LUN(L),213)'LEFTOVER',0.0, AMAWBEAM,
     +   (LEFTOVER(J)/NDELZ,J=1,NBINF-1)
      ENDIF
      IF(JOPTION(LFNCHAN))THEN
        WRITE(LUN(L),213)'NCHANNEL',CRITXDG,CRITYDG,
     +   (NCHANNEL(J)/NDELZ,J=1,NBINF-1)
      ENDIF
 213  FORMAT(//1X,A8,2F8.2/,(1X,10I7))
      IF(JOPTION(LFSCORE))THEN
C   Normalize to 'energy loss per collision, in eV'
        FNORM= 1.0E6/(NTRX*NDELZ)
        WRITE(LUN(L),*) 'SNUCLEAR'
        WRITE(LUN(L),878) (SNUCLEAR(IZ)*FNORM, IZ=0,NBINF-1)
        WRITE(LUN(L),*) 'SCORE'
        WRITE(LUN(L),878) (SCORE(IZ)*FNORM, IZ=0,NBINF-1)
        WRITE(LUN(L),*) 'SVALENCE'
        WRITE(LUN(L),878) (SVALENCE(IZ)*FNORM, IZ=0,NBINF-1)
      ENDIF
  878 FORMAT(1X,7F11.6)
      IF(JOPTION(LFEFINAL))THEN
        WRITE(LUN(L),883) TFMIN, EBIN, TFMAX
        WRITE(LUN(L),884) (EFINAL(IZ), IZ=0,NEBINA)
  883   FORMAT(1X,'Final energy distr: Tmin, step, Tmax',3(1X,F8.5))
  884   FORMAT(1X,10I7)
      ENDIF
      IF(JOPTION(LFVELO))THEN
        WRITE(LUN(L),885) ANGPMDG
      ENDIF
  885 FORMAT(1X,'velocity distr, range:',F7.2)
C OUTPUT TO LUN(3)=7 :
      IF(LUN(L).EQ.7)THEN
        THIK=ZIMAX
        WID=ZIFWHM
        IF(THIK .gt. 999999.0)THIK=999999.0
        IF(WID .gt. 999999.0)WID=999999.0
        CALL THETAPHI(THETA,PHI,K)
        WRITE(FLKOP,3)'FLUX ',VERSION,TIME,THETA,PHI,
     3   T0,Z1,Z2(1,1),THIK,WID,NDELZ,UD(1),
     4   XUN(1),YUN(1),TEMP,NTRX,LATTICE(1)
  3     FORMAT(A5,A27,A20,3F8.3,2F8.1,2F8.0,I6,F8.4,2F8.5,F8.2,I6,I4)
        IF (JOPTION(LFFLUX)) THEN
          NFX1= NFX
          NFY1= NFY
        ELSE
          NFX1=0
          NFY1=0
        ENDIF
        CALL PUTFLX(NBIN,NBINF,NFX1,NFY1,NVA,ML,NLA,NLAY,WEIGHTDE,
     6              FLKOP,RAAK,DAV,DDAV,FLUX,VFINAL)
      ENDIF
      CALL SYNCF(LUN(L))
  30  CONTINUE
      END
C ------------------------------------------------- SUBROUTINE FFIELD
      SUBROUTINE FFIELD
CALCULATES THE VECTOR FIELD FX,FY DUE TO SURROUNDING STRINGS IN
C THE CENTER OF EACH OF THE SUBCELLS IN WHICH THE UNIT CELL IS
C DIVIDED FOR THE FLUX CALCULATION.
C RESULTS ARE STORED IN ARRAYS FFX AND FFY.
      IMPLICIT NONE
#include "FLUX7.COM"
C FFX,FFY: components of the force field of the surrounding atom rows.
      REAL FFX, FFY
      COMMON /COMFIELD/ FFX(NF*NF,NLAYER),FFY(NF*NF,NLAYER)
C SETS ZZA ZZD ATF XUNA YUNA U2A U2ASQ of common block COMPOTE

      INTEGER KL, JXJY, IL, NPOS, J, JX,JY
      INTEGER ILL, IL1, IL2
      REAL YJ, XJ, FX, FY
      do 80 KL=1,NLAY
      DO 10 JXJY=  1,NFX*NFY
      FFX(JXJY,KL)=0.0
      FFY(JXJY,KL)=0.0
   10 CONTINUE
      DO 70 ILL=1,NLS(KL)
      IL1=ILL
C     IF(ILL.GT.NLA0(KL))IL1= ILL-1 ???
      DO 60 IL2=1,NLA1(IL1,KL)
      IL=ILMAP(IL1,KL)-1+IL2
      NPOS=NABUR(ILL,KL)
      ZZA=ZZA2(IL,KL)
      ZZD=ZZD2(IL,KL)
      ATF=RSCREEN(IL,KL)
      XUNA=XUN(KL)/ATF
      YUNA=YUN(KL)/ATF
      U2A=U2(IL,KL)/ATF
      U2ASQ=U2A*U2A
      DO 40 J=NPOS,1,-1
      JXJY=0
      DO 30 JY=1,NFY
      YJ = ( FLOAT(JY) - 0.5 ) / FLOAT(NFY)
      DO 20 JX=1,NFX
      JXJY=JXJY+1
      XJ = ( FLOAT(JX) - 0.5 ) / FLOAT(NFX)
      CALL FORCE(FX,FY,XJ,YJ,XPOS(J,ILL,KL),YPOS(J,ILL,KL))
      FFX(JXJY,KL)=ATOMFRAC(IL,KL)*FFX(JXJY,KL)+FX
      FFY(JXJY,KL)=ATOMFRAC(IL,KL)*FFY(JXJY,KL)+FY
  20  CONTINUE
  30  CONTINUE
  40  CONTINUE
  60  CONTINUE
  70  CONTINUE
  80  CONTINUE
** debugging output
*     write(6,90)NFX,NFY
*     write(6,100)(FFX(JX,1)*1E6,JX=1,JXJY)
*     write(6,100)(FFY(JX,1)*1E6,JX=1,JXJY)
* 90  FORMAT(//'FFX and FFY ... NFX,NFY =',2I5)
*100  FORMAT(/(6F12.7))
**
      END
C ------------------------------------------------- SUBROUTINE TEST
** check if a number is positive and in the range of floating point nrs
      SUBROUTINE TEST(POSNUM,STRING)
      IMPLICIT NONE
      REAL POSNUM
      CHARACTER*(*) STRING

      IF(POSNUM .GT. 0 .and. 1.0/POSNUM .ne. 0)RETURN
      WRITE(6,*)'TEST ',STRING,' error:',POSNUM
      STOP 'TEST'
      END
C ----------------------------------------------- SUBROUTINE EXPORTXYZ
      SUBROUTINE EXPORTXYZ(LOGLU,TRUEX,TRUEY,TRUEZ, PX,PY)
      IMPLICIT NONE
      INTEGER LOGLU
      REAL TRUEX,TRUEY,TRUEZ, PX,PY

#include "FLUX7.COM"
      CHARACTER*5 FORMX, FORMY, FORMZ
      CHARACTER*64 FORMXY

c choose format such that always 3 decimals are printed
      IF(ABS(TRUEX).LT.99.5) THEN
c +99.xxx needs 8 positions
        FORMX='F8.3'
      ELSE IF(ABS(TRUEX).LT.9999.5) THEN
c +9999.xxx needs 10 positions
        FORMX='F10.3'
      ELSE
c enough room for numbers up to +9999999.xxx i.e. 10**7
        FORMX='F13.3'
      ENDIF
      IF(ABS(TRUEY).LT.99.5) THEN
        FORMY='F8.3'
      ELSE IF(ABS(TRUEY).LT.9999.5) THEN
        FORMY='F10.3'
      ELSE
        FORMY='F13.3'
      ENDIF
      IF(ABS(TRUEZ).LT.99.5) THEN
        FORMZ='F8.3'
      ELSE IF(ABS(TRUEZ).LT.9999.5) THEN
        FORMZ='F10.3'
      ELSE
        FORMZ='F13.3'
      ENDIF
      FORMXY= '(''*xyz'','//FORMX//','//FORMY//','//FORMZ//',2F8.4)'

      WRITE(LOGLU,FORMXY) TRUEX,TRUEY,TRUEZ, PX,PY
      END
C --------------------------------------------------- SUBROUTINE WRT5
      SUBROUTINE WRT5(LU, STRING,X1,X2,X3,X4,X5)
      IMPLICIT NONE
      INTEGER LU
      CHARACTER*12 STRING
      REAL X1,X2,X3,X4,X5

      REAL X,Y(5)
      INTEGER N1, N2, K
      CHARACTER*6 FORMX
      CHARACTER*34 FORMXY

      Y(1)=X1
      Y(2)=X2
      Y(3)=X3
      Y(4)=X4
      Y(5)=X5

      FORMXY= '(A'
      N1=2
      do 999 K=1,5
      X= ABS(Y(K))
      IF(X.LT.99.5) THEN
        FORMX=',F10.5'
        N2=6
      ELSE IF(X.LT.9999.5) THEN
        FORMX=',F10.3'
        N2=6
      ELSE
        FORMX=',F12.1'
        N2=6
      ENDIF
      FORMXY= FORMXY(1:N1)//FORMX(1:N2)
      N1= N1+N2
  999 CONTINUE
      FORMXY= FORMXY(1:N1)//')'
      WRITE(LU,FORMXY) STRING,X1,X2,X3,X4,X5
      END
C --------------------------------------------------- SUBROUTINE DPTEST
** check validity of a double precision number
      logical function dptest(name,dpn)
      IMPLICIT NONE
      character*(*) name
      double precision dpn
#include "FLUX7.COM"

      dptest= .FALSE.
      if(dpn.ne.0.0)then
        if(dpn.ne.dpn .or. 1.0/dpn .eq. 0.0) then
          write(6,*)'dptest: ',name, ' is not a valid number:',dpn
          IF(LUOPEN(2).eq.1)THEN
          write(LUN(2),*)'dptest: ',name,' is not a valid number:',dpn
          ENDIF
          dptest=.TRUE.
        endif
      else if(dpn.gt. 1.0e30 .or. dpn.lt. -1.0e30)then
        write(6,*)'dptest: overflow danger variable ',name,dpn
        IF(LUOPEN(2).eq.1)THEN
          write(LUN(2),*)'dptest: overflow danger variable ',name,dpn
        ENDIF
        dptest=.TRUE.
      endif
      end
