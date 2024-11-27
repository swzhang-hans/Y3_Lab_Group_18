      PROGRAM GEOWING
C
C**   BY R.H.J. WILLDEN
C**   JUNE 1998
C
C**   PRODUCES A COMPLETE 3-DIMENSIONAL WING MESH TO BE USED BY PANEL
C**   PROGRAMS, IN PARTICULAR 'combocl.f'
C
      PARAMETER (NDSZ=100,NASZ=81,NRSZ=20,NSZ=3000,NSSZ=51)
C
      REAL PI,RD(NDSZ,2),RA(NASZ,2),LEN(2,NASZ),L(NASZ),THETA,
     &     DTHETA,R,ERROR,MODLEN1,MODLEN2,ROOTC,SSPAN,RC(3),
     &     RN(NSZ,3),MODL,FACTOR,TEMP,H,HMAX,RR(NRSZ,2),A(3),B(3),
     &     MODA,MODB,NACANO,NACAM,NACAP,NACAT,MODUPPER,MODLOWER,SAME,
     &     TAPER,LESWEEP,DIHEDRAL,TWIST,ROOTINC,XLE(NSSZ),YLE(NSSZ),
     &     ZLE(NSSZ),CHORD(NSSZ),ALPHA(NSSZ),BETA,GAMMA,RF(3),S,
     &     A1,A2,TEBI
      INTEGER ND,NTE,NP,PARAM,SP,NNTIP,NNWING,NPTIP,NPWING,PNL(NSZ,4),
     &     WPNL(NSSZ,2),NWP,OVER,NMID,NRAD,K1,K2,DATATYPE,N
      CHARACTER *1 CHOICE
      CHARACTER *20 NAME
      CHARACTER *5 SWEEP
C
      PI = 2*ACOS(0.)
      ERROR = 1.0E-4
      FACTOR = 4./3.
      SAME = 1.0E-12
C
C**   AEROFOIL DATA INPUT
C
      PRINT *,'DO YOU WISH TO USE A 4 DIGIT NACA AEROFOIL SECTION?'
 10   CONTINUE
      PRINT *,'(Y OR N)'
      READ (*,20) CHOICE
 20   FORMAT (A1)
      IF ((CHOICE.EQ.'y').OR.(CHOICE.EQ.'Y')) THEN
         PRINT *,'PLEASE ENTER THE 4 DIGIT NACA DESIGNATION'
         READ *,NACANO
         NACAM = INT(NACANO/1000)
         NACAP = INT((NACANO-NACAM*1000)/100)
         NACAT = INT(NACANO-NACAM*1000-NACAP*100)
         NACAM = NACAM/100
         NACAP = NACAP/10
         NACAT = NACAT/100
         WRITE (*,30) INT(NACAM*100),INT(NACAP*100),INT(NACAT*100)
 30      FORMAT (' YOU HAVE SELECTED ',I2,'% CAMBER AT ',I2,'% CHORD FRO
     &M THE LEADING EDGE',/,' AND A ',I2,'% THICKNESS DISTRIBUTION')
         DATATYPE = 1
      ELSEIF ((CHOICE.EQ.'n').OR.(CHOICE.EQ.'N')) THEN
         OPEN (UNIT=1,FILE='aerofoil.dat',STATUS='OLD')
         REWIND (1)
         READ (1,40) NAME
 40      FORMAT (A20)
         READ (1,*) ND
         DO 50 I = 1,ND
            READ (1,*) (RD(I,J),J=1,2)
 50      CONTINUE
         CLOSE (1)
         WRITE (*,60) NAME
 60      FORMAT (' DISCRETE AEROFOIL DATA FOR ',A20,/,' HAS BEEN READ IN
     & FROM DATA FILE aerofoil.dat')
         CALL STANDARDISE(RD,ND,NTE)
         DATATYPE = 2
      ELSE
         GOTO 10
      ENDIF
C
C**   HOW MANY PANELS ON EACH SIDE OF THE AEROFOIL?
C
      PRINT *,'HOW MANY PANELS ON EACH SIDE OF THE AEROFOIL?'
      READ *,NP
C
C**   USE A COSINE DISTRIBUTION TO DETERMINE THE PERCENTAGE OF THE
C**   SURFACE OCCUPIED BY EACH PANEL
C
      DTHETA = PI/NP
      THETA = 0.
      R = 0.5
      DO 100 I = 1,NP
         THETA = THETA+DTHETA
         L(I) = R-0.5*COS(THETA)
         R = 0.5*COS(THETA)
 100  CONTINUE
C
C**   NOW ADJUST EACH PANEL'S RELATIVE LENGTH SO THAT NO PANEL EXCEEDS 
C**   FACTOR*(ITS NEIGHBOURS' LENGTHS).  THEN NORMALISE THE PANEL
C**   LENGTHS TO BE PERCENTAGES OF THE AEROFOIL SURFACE.
C
      TEMP = NP/2
      NMID = NINT(TEMP)
      DO 110 I = NMID,2,-1
         TEMP = L(I)/L(I-1)
         IF (TEMP.GT.FACTOR) THEN
            L(I-1)=L(I)/FACTOR
         ENDIF
 110  CONTINUE
C
      DO 120 I = 1,NMID
         L(NP+1-I) = L(I)
 120  CONTINUE
C
      MODL = 0.
      DO 130 I = 1,NP
         MODL = MODL+L(I)
 130  CONTINUE
C
      DO 140 I = 1,NP
         L(I) = L(I)/MODL
 140  CONTINUE
C
C**   FIT UPPER SURFACE NODE POINTS BY ITERATION.  IF RA(NP+1,*) IS
C**   WITHIN ERROR OF THE TRAILING EDGE THEN ADJUST RA(NP+1,*) ONTO THE
C**   TRAILING EDGE POSITION AND CONTINUE
C
      MODLOWER = 0.8
      MODUPPER = 1.2
      MODLEN1 = 0.5*(MODLOWER+MODUPPER)
      DO 275 I = 1,NP
         LEN(1,I) = L(I)*MODLEN1
 275  CONTINUE
      PARAM = 1
      RA(1,1) = 0.
      RA(1,2) = 0.
 300  CONTINUE
      IF (DATATYPE.EQ.1) THEN
         CALL NACA4(RA,LEN,NACAM,NACAP,NACAT,PARAM,NP,OVER)
      ELSEIF (DATATYPE.EQ.2) THEN
         CALL AERONODE(RD,RA,LEN,PARAM,NP,NTE,ND,OVER)
      ENDIF
      IF (OVER.EQ.2) THEN
         MODUPPER = MODLEN1
      ELSE
         R = SQRT((RA(NP+1,1)-1.)**2+RA(NP+1,2)**2)
         IF ((R.LE.ERROR).OR.((MODUPPER-MODLOWER).LT.SAME)) GOTO 350
         MODLOWER = MODLEN1
      ENDIF
C
      MODLEN1 = 0.5*(MODUPPER+MODLOWER)
      DO 325 I = 1,NP
         LEN(1,I) = L(I)*MODLEN1
 325  CONTINUE
      GOTO 300
C
 350  CONTINUE
      RA(NP+1,1) = 1.
      RA(NP+1,2) = 0.
C
C**   NOW FIT LOWER SURFACE NODES IN THE SAME MANNER, EXCEPT THAT 
C**   RA(2*NP+1,*) MUST BE WITHIN ERROR OF THE LEADING EDGE
C
      MODLOWER = 0.8
      MODUPPER = 1.2
      MODLEN2 = 0.5*(MODLOWER+MODUPPER)
      DO 375 I = 1,NP
         LEN(2,I) = L(I)*MODLEN2
 375  CONTINUE
      PARAM = 2
 400  CONTINUE
      IF (DATATYPE.EQ.1) THEN
         CALL NACA4(RA,LEN,NACAM,NACAP,NACAT,PARAM,NP,OVER)
      ELSEIF (DATATYPE.EQ.2) THEN
         CALL AERONODE(RD,RA,LEN,PARAM,NP,NTE,ND,OVER)
      ENDIF
      IF (OVER.EQ.2) THEN
         MODUPPER = MODLEN2
      ELSE
         R = SQRT(RA(2*NP+1,1)**2+RA(2*NP+1,2)**2)
         IF ((R.LE.ERROR).OR.((MODUPPER-MODLOWER).LT.SAME)) GOTO 450
         MODLOWER = MODLEN2
      ENDIF
C
      MODLEN2 = 0.5*(MODUPPER+MODLOWER)
      DO 425 I = 1,NP
         LEN(2,I) = L(I)*MODLEN2
 425  CONTINUE
      GOTO 400
C
 450  CONTINUE
      RA(2*NP+1,1) = 0.
      RA(2*NP+1,2) = 0.
C
C**   AEROFOIL NODE POSITIONS NOW KNOWN
C**   NOW CALCULATE THE TRAILING EDGE BISECTOR
C
      TEMP = RA(NP,2)/(1-RA(NP,1))
      A1 = ATAN(TEMP)
      TEMP = RA(NP+2,2)/(1-RA(NP+2,1))
      A2 = ATAN(TEMP)
      TEBI = -0.5*(A1+A2)
      TEBI = TEBI*180/PI
      PRINT *,'TRAILING EDGE BISECTOR ANGLE = ',TEBI
C
C**   NOW DETERMINE POSITION VECTORS OF ALL WING NODES, RN(I,*)
C**   FIRST READ IN WING GEOMETRY PARAMETERS
C
      PRINT *,'ENTER ROOT CHORD'
      READ *,ROOTC
      PRINT *,'ENTER (FLAT) WING SEMI-SPAN (EXCLUDING WING TIPS)'
      READ *,SSPAN
      PRINT *,'ENTER TAPER RATIO'
      READ *,TAPER
      S = SSPAN*ROOTC*(1+TAPER)
      PRINT *,'ENTER nccls FOR NO CENTRE CHORD LINE SWEEP'
      PRINT *,'OR sweep TO APPLY LEADING EDGE SWEEP'
 485  CONTINUE
      READ (*,490) SWEEP
 490  FORMAT (A5)
      IF ((SWEEP.EQ.'sweep').OR.(SWEEP.EQ.'SWEEP')) THEN
         PRINT *,'ENTER LEADING EDGE SWEEP ANGLE (IN DEGREES)'
         READ *,LESWEEP
         LESWEEP = LESWEEP*PI/180
      ELSEIF ((SWEEP.EQ.'nccls').OR.(SWEEP.EQ.'NCCLS')) THEN
         LESWEEP = ATAN(ROOTC*(1-TAPER)/(2*SSPAN))
         PRINT *,'BASED UPON THE TAPER RATIO, THE LEADING EDGE SWEEP'
         PRINT *,'ANGLE HAS BEEN CALCULATED AS',LESWEEP*180/PI,'DEGREES'
      ELSE
         PRINT *,'nccls OR sweep?'
         GOTO 485
      ENDIF
      PRINT *,'ENTER DIHEDRAL ANGLE (IN DEGREES)'
      READ *,DIHEDRAL
      DIHEDRAL = DIHEDRAL*PI/180
      PRINT *,'ENTER TWIST ANGLE (IN DEGREES)'
      READ *,TWIST
      TWIST = TWIST*PI/180
      PRINT *,'ENTER ROOT INCIDENCE (IN DEGREES)'
      READ *,ROOTINC
      ROOTINC = ROOTINC*PI/180
      PRINT *,'HOW MANY PANELS IN THE SPANWISE DIRECTION PER SEMI-SPAN'
      PRINT *,'(EXCLUDING WING TIPS)?'
      READ *,SP
C
C**   USE A COSINE DISTRIBUTION FOR THE SPANWISE DISTRIBUTION OF PANELS
C
      DTHETA = PI/(2*SP)
      THETA = 0.
      R = 0.5
      DO 510 I = 1,2*SP
         THETA = THETA+DTHETA
         L(I) = R-0.5*COS(THETA)
         R = 0.5*COS(THETA)
 510  CONTINUE
C
C**   AS BEFORE,ADJUST EACH PANEL'S RELATIVE WIDTH SO THAT ALL PANELS
C**   DO NOT EXCEED FACTOR*(THEIR NEIGHBOUR'S WIDTH).  THEN ADJUST THE
C**   PANEL WIDTHS SO THAT THEY FIT THE FULL SPAN, USING THE 
C**   PROPORTIONS CALCULATED.
C
      NMID = INT((2*SP+1)/2)
      DO 520 I = NMID,2,-1
         TEMP = L(I)/L(I-1)
         IF (TEMP.GT.FACTOR) THEN
            L(I-1)=L(I)/FACTOR
         ENDIF
 520  CONTINUE
C
      DO 530 I = 1,NMID
         L(2*SP+1-I) = L(I)
 530  CONTINUE
C
      MODL = 0.
      DO 540 I = 1,2*SP
         MODL = MODL+L(I)
 540  CONTINUE
C
      DO 550 I = 1,2*SP
         L(I) = L(I)*2*SSPAN/MODL
 550  CONTINUE
C
C**   NOW DETERMINE HOW MANY PANELS ARE REQUIRED FOR THE SEMI-CIRCULAR
C**   WING TIPS, SUCH THAT THE WIDTH OF THESE PANELS ARE SMALLER THAN
C**   THAT OF THEIR MAIN BODY NEIGHBOURS
C
      HMAX = 0.
      DO 560 I = 1,NP-1
         H = SQRT((RA(I+1,2)-RA(2*NP+1-I,2))**2+(RA(I+1,1)-RA(2*NP+1-I,
     &        1))**2)
         IF (H.GT.HMAX) THEN
            HMAX = H
         ENDIF
 560  CONTINUE
C
      NRAD = INT(((PI*HMAX*ROOTC)/(2*L(1)))+1)
      IF (NRAD.LE.2) THEN
         NRAD = 3
      ENDIF
      DTHETA = PI/NRAD
      THETA = 0.
      DO 570 K = 1,NRAD-1
         THETA = THETA+DTHETA
         RR(K,1) = SIN(THETA)
         RR(K,2) = COS(THETA)
 570  CONTINUE
      NNTIP = (NRAD-1)*(NP-1)
C
C**   CALCULATE ALL LOCAL SECTION PROPERTIES
C
      DO 600 I = 1,2*SP+1
         YLE(I) = -SSPAN
         DO 590 J = 1,I-1
            YLE(I) = YLE(I)+L(J)
 590     CONTINUE
         CHORD(I) = ROOTC*(TAPER+(1-TAPER)*(1-(ABS(YLE(I))/SSPAN)))
         ALPHA(I) = ROOTINC+TWIST*(ABS(YLE(I))/SSPAN)
         XLE(I) = ABS(YLE(I))*TAN(LESWEEP)
         ZLE(I) = ABS(YLE(I))*SIN(DIHEDRAL)
         YLE(I) = YLE(I)*COS(DIHEDRAL)
 600  CONTINUE
C
      SSPAN = SSPAN*COS(DIHEDRAL)
C
C**   WING BODY NODES (I.E. ALL NODES EXCLUDING WING TIPS)
C
      K2 = NNTIP
      DO 630 I = 1,2*SP+1
         DO 620 J = 1,2*NP
            K2 = K2+1
            BETA = ATAN(RA(J,2)/(0.25-RA(J,1)))
            IF (RA(J,1).GT.(0.25)) THEN
               BETA = BETA+PI
            ENDIF
            GAMMA = BETA+ALPHA(I)
            R = CHORD(I)*(SQRT(RA(J,2)**2+(0.25-RA(J,1))**2))
            RN(K2,1) = XLE(I)+(CHORD(I)/4)-R*COS(GAMMA)
            RN(K2,2) = YLE(I)
            RN(K2,3) = ZLE(I)+R*SIN(GAMMA)
 620     CONTINUE
 630  CONTINUE
      NNWING = K2-NNTIP
C
C**   WING TIP NODES (-VE Y)
C
      K2 = 0
      DO 690 I = 1,NP-1
         RC(1) = 0.5*(RN(NNTIP+I+1,1)+RN(NNTIP+2*NP+1-I,1))
         RC(2) = -SSPAN
         RC(3) = 0.5*(RN(NNTIP+I+1,3)+RN(NNTIP+2*NP+1-I,3))
         B(1) = RN(NNTIP+I+1,1)-RN(NNTIP+2*NP+1-I,1)
         B(2) = 0.
         B(3) = RN(NNTIP+I+1,3)-RN(NNTIP+2*NP+1-I,3)
         H = SQRT(B(1)**2+B(3)**2)
         DO 650 K = 1,3
            B(K) = B(K)/H
 650     CONTINUE
         DO 660 K = 1,2
            RF(K) = (RC(K)-RN(NNTIP+1,K))/TAPER
            A(K) = RC(K)-RF(K)
 660     CONTINUE
         MODA = SQRT(A(1)**2+A(2)**2)
         A(1) = A(1)/MODA
         A(2) = A(2)/MODA
         A(3) = 0.
         DO 680 K1 = 1,NRAD-1
            K2 = K2+1
            MODA = RR(K1,1)
            MODB = RR(K1,2)
            DO 670 K = 1,3
               RN(K2,K) = RC(K)+(MODA*A(K)+MODB*B(K))*(H/2)
 670        CONTINUE
 680     CONTINUE
 690  CONTINUE
C
C**   WING TIP NODES (+VE Y)
C
      N = NNTIP+NNWING
      K2 = NNTIP+NNWING
      DO 740 I = 1,NP-1
         RC(1) = 0.5*(RN(N-2*NP+I+1,1)+RN(N+1-I,1))
         RC(2) = SSPAN
         RC(3) = 0.5*(RN(N-2*NP+I+1,3)+RN(N+1-I,3))
         B(1) = RN(N-2*NP+I+1,1)-RN(N+1-I,1)
         B(2) = 0.
         B(3) = RN(N-2*NP+I+1,3)-RN(N+1-I,3)
         H = SQRT(B(1)**2+B(3)**2)
         DO 700 K = 1,3
            B(K) = B(K)/H
 700     CONTINUE
         DO 710 K = 1,2
            RF(K) = (RC(K)-RN(N-2*NP+1,K))/TAPER
            A(K) = RC(K)-RF(K)
 710     CONTINUE
         MODA = SQRT(A(1)**2+A(2)**2)
         A(1) = A(1)/MODA
         A(2) = A(2)/MODA
         A(3) = 0.
         DO 730 K1 = 1,NRAD-1
            K2 = K2+1
            MODA = RR(K1,1)
            MODB = RR(K1,2)
            DO 720 K = 1,3
               RN(K2,K) = RC(K)+(MODA*A(K)+MODB*B(K))*(H/2)
 720        CONTINUE
 730     CONTINUE
 740  CONTINUE
C
C**   NOW DETERMINE EACH PANELS CONNECTIVITY LIST PNL(I,*)
C
C**   WING TIP PANELS (-VE Y)
C
      K1 = 0
      DO 800 J = 1,NP
         DO 775 I = 1,NRAD
            K1 = K1+1
            PNL(K1,1) = (J-1)*(NRAD-1)+I
            PNL(K1,2) = (J-2)*(NRAD-1)+I
            PNL(K1,3) = (J-2)*(NRAD-1)+I-1
            PNL(K1,4) = (J-1)*(NRAD-1)+I-1
            IF (J.EQ.1) THEN
               PNL(K1,2) = NNTIP+1
               PNL(K1,3) = NNTIP+1
            ELSEIF (J.EQ.NP) THEN
               PNL(K1,1) = NNTIP+NP+1
               PNL(K1,4) = NNTIP+NP+1
            ENDIF
            IF ((I.EQ.1).AND.(J.NE.NP)) THEN
               PNL(K1,4) = NNTIP+J+1
            ELSEIF ((I.EQ.NRAD).AND.(J.NE.NP)) THEN
               PNL(K1,1) = NNTIP+2*NP-J+1
            ENDIF
            IF ((I.EQ.1).AND.(J.NE.1)) THEN
               PNL(K1,3) = NNTIP+J
            ELSEIF ((I.EQ.NRAD).AND.(J.NE.1)) THEN
               PNL(K1,2) = NNTIP+2*NP-J+2
            ENDIF
 775     CONTINUE
 800  CONTINUE
      NPTIP = K1
C
C**   ALL PANELS, EXCLUDING TIPS
C
      DO 850 J = 1,2*SP
         DO 825 I = 1,2*NP
            K1 = K1+1
            PNL(K1,1) = NNTIP+(J-1)*2*NP+I
            PNL(K1,2) = NNTIP+J*2*NP+I
            IF (I.NE.(2*NP)) THEN
               PNL(K1,3) = NNTIP+J*2*NP+I+1
               PNL(K1,4) = NNTIP+(J-1)*2*NP+I+1
            ELSE
               PNL(K1,3) = NNTIP+J*2*NP+1
               PNL(K1,4) = NNTIP+(J-1)*2*NP+1
            ENDIF
 825     CONTINUE
 850  CONTINUE
      NPWING = K1-NPTIP
C
C**   WING TIP PANELS (+VE Y)
C
      DO 900 J = 1,NP
         DO 875 I = 1,NRAD
            K1 = K1+1
            PNL(K1,1) = NNTIP+NNWING+(J-2)*(NRAD-1)+I
            PNL(K1,2) = NNTIP+NNWING+(J-1)*(NRAD-1)+I
            PNL(K1,3) = NNTIP+NNWING+(J-1)*(NRAD-1)+I-1
            PNL(K1,4) = NNTIP+NNWING+(J-2)*(NRAD-1)+I-1
            IF (J.EQ.1) THEN
               PNL(K1,1) = NNTIP+NNWING-2*NP+1
               PNL(K1,4) = NNTIP+NNWING-2*NP+1
            ELSEIF (J.EQ.NP) THEN
               PNL(K1,2) = NNTIP+NNWING-NP+1
               PNL(K1,3) = NNTIP+NNWING-NP+1
            ENDIF
            IF ((I.EQ.1).AND.(J.NE.NP)) THEN
               PNL(K1,3) = NNTIP+NNWING-2*NP+J+1
            ELSEIF ((I.EQ.NRAD).AND.(J.NE.NP)) THEN
               PNL(K1,2) = NNTIP+NNWING+1-J
            ENDIF
            IF ((I.EQ.1).AND.(J.NE.1)) THEN
               PNL(K1,4) = NNTIP+NNWING-2*NP+J
            ELSEIF ((I.EQ.NRAD).AND.(J.NE.1)) THEN
               PNL(K1,1) = NNTIP+NNWING+2-J
            ENDIF
 875     CONTINUE
 900  CONTINUE
C
C**   WAKE CONNECTIVITY LIST
C
      NWP = SP*2
      DO 950 I = 1,NWP
         WPNL(I,1) = NNTIP+(I-1)*2*NP+(NP+1)
         WPNL(I,2) = NNTIP+I*2*NP+(NP+1)
 950  CONTINUE
C
C**   DATA OUTPUT STAGE
C
      OPEN (UNIT=11,FILE='mesh.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=12,FILE='view.dat',STATUS='UNKNOWN')
      OPEN (UNIT=13,FILE='auxiliary.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=21,FILE='section.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=22,FILE='relative.dat',STATUS='UNKNOWN')
      REWIND (11)
C      REWIND (12)
      REWIND (13)
C      REWIND (21)
C      REWIND (22)
C
C**   OUTPUT FOR PANEL METHOD TO 'mesh.dat'
C
      WRITE (11,*) NNTIP*2+NNWING
      DO 1010 I = 1,NNTIP*2+NNWING
         WRITE (11,*) (RN(I,J),J=1,3)
 1010 CONTINUE
      WRITE (11,*) NPTIP*2+NPWING
      DO 1020 I = 1,NPTIP*2+NPWING
         WRITE (11,*) (PNL(I,J),J=1,4)
 1020 CONTINUE
      WRITE (11,*) NWP
      DO 1025 I = 1,NWP
         WRITE (11,*) (WPNL(I,J),J=1,2)
 1025 CONTINUE
      CLOSE (11)
C
C**   OUTPUT FOR TECPLOT VIEWING TO 'view.dat'
C
C      WRITE (12,*) 'VARIABLES = "X", "Y", "Z"'
C      WRITE (12,*) 'ZONE N=',NNTIP*2+NNWING,', E=',NPTIP*2+NPWING
C      WRITE (12,*) 'F=FEPOINT, ET=QUADRILATERAL'
C      DO 1030 I = 1,NNTIP*2+NNWING
C         WRITE (12,*) (RN(I,J),J=1,3)
C 1030 CONTINUE
C      DO 1040 I = 1,NPTIP*2+NPWING
C         WRITE (12,*) (PNL(I,J),J=1,4)
C 1040 CONTINUE
C      CLOSE (12)
C
C**   OUTPUT TO FILE 'auxiliary.dat'
C
      WRITE (13,*) S
      WRITE (13,*) TEBI
      WRITE (13,*) 2*NP
      WRITE (13,*) NPTIP
      CLOSE (13)
C
C**   OUTPUT AEROFOIL SECTION TO 'section.dat' FOR VIEWING BY TECPLOT
C
C      WRITE (21,*) 'VARIABLES = "X", "Z"'
C      WRITE (21,*) 'ZONE I=',NP*2+1,'F=POINT'
C      DO 1050 I = 1,2*NP+1
C         WRITE (21,*) (RA(I,J),J=1,2)
C 1050 CONTINUE
C      CLOSE (21)
C
C**   OUTPUT RELATIVE PANEL LENGTHS (CHORDWISE) TO 'relative.dat'
C
C      DO 1060 I = 1,2*NP
C         TEMP = SQRT((RA(I+1,1)-RA(I,1))**2+(RA(I+1,2)-RA(I,2))**2)
C         WRITE (22,*) I,TEMP,TEMP/R
C         R = TEMP
C 1060 CONTINUE
C      CLOSE (22)
C
C**   OUTPUT TO SCREEN
C
      PRINT *,'NUMBER OF NODES  =',NNTIP*2+NNWING
      PRINT *,'NUMBER OF PANELS =',NPTIP*2+NPWING
      PRINT *,'NUMBER OF PANELS PER TIP =',NPTIP
C
C**   DATA OUTPUT COMPLETE
C
      STOP
      END
C
C*****************************SUBROUTINES*******************************
C
      SUBROUTINE STANDARDISE(RD,ND,NTE)
C
C**   THIS SUBROUTINE STANDARDISES THE DISCRETE AEROFOIL DATA BY MOVING
C**   THE LEADING EDGE TO THE ORIGIN, ROTATING THE TRAILING EDGE TO LIE 
C**   ALONG THE X-AXIS AND SCALING THE AEROFOIL SO THAT THE CHORD IS OF
C**   UNIT LENGTH
C
      PARAMETER (NDSZ=100,NASZ=81,NRSZ=10,NSZ=2000,NSSZ=51)
C
      REAL RD(NDSZ,2),XMOVE,ZMOVE,P,ALPHA,BETA,GAMMA,R,SCALE
      INTEGER ND,NTE
C
C**   MOVE LEADING EDGE TO THE ORIGIN, (PRESUMING THE FIRST DATA POINT
C**   IS THE LEADING EDGE)
C
      XMOVE = RD(1,1)
      ZMOVE = RD(1,2)
      DO 2050 I = 1,ND
         RD(I,1) = RD(I,1)-XMOVE
         RD(I,2) = RD(I,2)-ZMOVE
 2050 CONTINUE
C
C**   DUPLICATE LEADING EDGE POINT AT END OF DISCRETE POINT LIST
C
      RD(ND+1,1) = 0.
      RD(ND+1,2) = 0.
C
C**   DETERMINE THE TRAILING EDGE POINT, NTE
C
      P = 0.
      DO 2100 I = 2,ND
         R = SQRT(RD(I,1)**2+RD(I,2)**2)
         IF (R.GT.P) THEN
            P = R
            NTE = I
         ENDIF
 2100 CONTINUE
C
C**   DETERMINE INCLINATION OF CHORD LINE
C
      ALPHA = ATAN(RD(NTE,2)/RD(NTE,1))
C
C**   ROTATE AEROFOIL SO THAT CHORD LINE LIES ALONG X-AXIS
C
      DO 2150 I = 2,ND
         BETA = ATAN(RD(I,2)/RD(I,1))
         GAMMA = BETA - ALPHA
         R = SQRT(RD(I,1)**2+RD(I,2)**2)
         RD(I,1) = R*COS(GAMMA)
         RD(I,2) = R*SIN(GAMMA)
 2150 CONTINUE
C
C**   SCALE AEROFOIL SO THAT CHORD IS OF UNIT LENGTH
C
      SCALE = RD(NTE,1)
      DO 2200 I = 2,ND
         RD(I,1) = RD(I,1)/SCALE
         RD(I,2) = RD(I,2)/SCALE
 2200 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE AERONODE(RD,RA,LEN,PARAM,NP,NTE,ND,OVER)
C
C**   POSITIONS THE AEROFOIL NODES AROUND THE STANDARDISED AEROFOIL,
C**   USING THE LEN ARRAYS AND LINEAR INTERPOLATION BETWEEN DATA POINTS
C
      PARAMETER (NDSZ=100,NASZ=81,NRSZ=10,NSZ=2000,NSSZ=51)
C
      REAL R,RD(NDSZ,2),RA(NASZ,2),LEN(2,NASZ),D(2),G(2),MODD,MODG,A,
     &     B,C,TEMP,ALPHA
      INTEGER PARAM,I1,I2,NP,NTE,ND,J1,J2,K1,K2,OVER,M
C
C**   IF PARAM = 1 THEN POSITION NODES ON UPPER SURFACE, ELSE IF PARAM
C**   = 2 THEN POSITION NODES ON LOWER SURFACE
C
      IF (PARAM.EQ.1) THEN
         I1 = 1
         I2 = NP
         J1 = 2
         J2 = NTE
      ELSEIF (PARAM.EQ.2) THEN
         I1 = NP+1
         I2 = 2*NP
         J1 = NTE+1
         J2 = ND+1
      ENDIF
C
C**   LOOP THROUGH EACH NODE, I, DETERMINING THE POSITION OF THE I+1^TH
C**   NODE USING LINEAR INTERPOLATION BETWEEN DATA POINTS
C
      OVER = 1
      M = 0
      DO 2550 I = I1,I2
         M = M+1
         DO 2400 J = J1,J2
            R = SQRT((RD(J,1)-RA(I,1))**2+(RD(J,2)-RA(I,2))**2)
            IF (R.GE.LEN(PARAM,M)) GOTO 2425
 2400    CONTINUE
         OVER = 2
         GOTO 2600
 2425    CONTINUE
         J1 = J
         K1 = J-1
         K2 = J
         DO 2450 K = 1,2
            D(K) = RD(K2,K)-RD(K1,K)
            G(K) = RD(K1,K)-RA(I,K)
 2450    CONTINUE
         MODD = SQRT(D(1)**2+D(2)**2)
         MODG = SQRT(G(1)**2+G(2)**2)
         DO 2475 K = 1,2
            D(K) = D(K)/MODD
 2475    CONTINUE
         C = MODG**2-(LEN(PARAM,M)**2)
         B = 2*(G(1)*D(1)+G(2)*D(2))
         TEMP = SQRT(B**2-4*C)
         ALPHA = 0.5*(TEMP-B)
         DO 2500 K = 1,2
            RA(I+1,K) = RA(I,K)+G(K)+ALPHA*D(K)
 2500    CONTINUE
 2550 CONTINUE
C
 2600 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE NACA4(RA,LEN,M,P,T,PARAM,NP,OVER)
C
C**   POSITIONS THE AEROFOIL NODES, BY ITERATION, AROUND A 4 DIGIT NACA
C**   AEROFOIL DESIGNATED BY M,P & T.  THE DISTANCES BETWEEN ADJACENT 
C**   NODES ARE GIVEN BY THE LEN ARRAYS
C
      PARAMETER (NDSZ=100,NASZ=81,NRSZ=10,NSZ=2000,NSSZ=51)
C     
      REAL MARGIN,FACTOR,XC,YC,YT,RA(NASZ,2),LEN(2,NASZ),M,P,T,TTHETA,
     &     STHETA,CTHETA,Q,DIST,LIMIT,XLL,XUL,INTERVAL
      INTEGER PARAM,I1,I2,J,K,OVER,NP
C
      MARGIN = 1.0E-6
      FACTOR = 1.0E3
C
C**   IF PARAM = 1 THEN POSITION THE NODES ON THE UPPER SURFACE, ELSE IF
C**   PARAM = 2 THEN POSITION NODES ON THE LOWER SURFACE
C
      IF (PARAM.EQ.1) THEN
         I1 = 1
         I2 = NP
      ELSEIF (PARAM.EQ.2) THEN
         I1 = NP+1
         I2 = 2*NP
      ENDIF
C
C**   LOOP THROUGH EACH NODE, I, DETERMINING THE POSITION OF THE I+1^TH
C**   NODE BY LOWER & UPPER LIMIT ITERATION
C
      OVER = 1
      J = 0
      XLL = 0.
      XUL = 1.
C
      DO 2800 I = I1,I2
         J = J+1
 2750    CONTINUE
         XC = 0.5*(XLL+XUL)
         IF (XC.LE.P) THEN
            YC = (M/(P**2))*(2*P*XC-XC**2)
            TTHETA = (2*M/(P**2))*(P-XC)
         ELSE
            YC = (M/((1-P)**2))*(1-2*P+2*P*XC-XC**2)
            TTHETA = (2*M/((1-P)**2))*(P-XC)
         ENDIF
         CTHETA = (1+TTHETA**2)**(-0.5)
         STHETA = TTHETA*CTHETA
         YT = (T/0.2)*(0.2969*SQRT(XC)-0.126*XC-0.3516*(XC**2)+
     &        0.2843*(XC**3)-0.1015*(XC**4))
         RA(I+1,1) = XC+((-1)**PARAM)*YT*STHETA
         RA(I+1,2) = YC+((-1)**(PARAM+1))*YT*CTHETA
         Q = SQRT((RA(I+1,1)-RA(I,1))**2+(RA(I+1,2)-RA(I,2))**2)
         DIST = LEN(PARAM,J)-Q
C
         IF (ABS(DIST).GT.MARGIN) THEN
            IF (((DIST.LT.0).AND.(PARAM.EQ.1)).OR.((DIST.GT.0).AND.
     &           (PARAM.EQ.2))) THEN
               XUL = XC
            ELSE
               XLL = XC
            ENDIF
            INTERVAL = XUL-XLL
            IF ((DIST-MARGIN).GT.(FACTOR*INTERVAL)) THEN
               OVER = 2
               GOTO 2850
            ELSE
               GOTO 2750
            ENDIF
         ELSE
            IF (PARAM.EQ.1) THEN
               XLL = XC
               XUL = 1.
            ELSEIF (PARAM.EQ.2) THEN
               XLL = 0.
               XUL = XC
            ENDIF
         ENDIF
C
 2800 CONTINUE
C
 2850 CONTINUE
C
      RETURN
      END