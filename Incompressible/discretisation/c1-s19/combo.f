      PROGRAM COMBO
C
C**   BY R.H.J. WILLDEN
C
C**   3-DIMENSIONAL PANEL METHOD USING TWO SINGULARITY DISTRIBUTIONS;
C     BOTH A VORTEX RING DISTRIBUTION AS WELL AS A CONSTANT STRENGTH
C     DENSITY SOURCE DISTRIBUTION, TO SOLVE THE POTENTIAL FLOW PROBLEM
C     ABOUT AN ARBITRARY NON-LIFTING BODY.  IT UTILISES LU DECOMPOSITION
C     FOR THE MATRIX SOLVING PROCEDURE.  THE SURFACE SPEED AND Cp ARE
C     EVALUATED AT EACH PANEL'S COLLOCATION POINT.
C
      INTEGER SIZE,WSIZE
      PARAMETER (SIZE=1250,WSIZE=40)
      REAL PI,RN(SIZE,3),MODQINF,PITCH,YAW,RHO,QU(3),QINF(3),HEAD,R1(3),
     &     R2(3),R3(3),R4(3),NORM(3),RCOL(3),AREA(SIZE),NOR(SIZE,3),
     &     RC(SIZE,3),G,S,TEMP,AVEC(SIZE,SIZE,3),A(SIZE,SIZE),DQ(3),
     &     Q(3),DQ1(3),DQ2(3),DQ3(3),DVEC(SIZE,WSIZE,3),SIGMA(SIZE),
     &     R(5,3),CVEC(SIZE,SIZE,3),C(SIZE,SIZE),B(SIZE),THETA,
     &     F(3),SVEL(SIZE,3),SPEED(SIZE),CP(SIZE),
     &     SU(3),E(SIZE,SIZE),DX2,DY2,CLW,CDW,
     &     RHS(SIZE),GAMMA(SIZE),GVEC(SIZE),WANGLE,TU(3),DX1,DY1,A1,A2,
     &     SAREA,BETA,MODF,LIFT,DRAG,TRANS,CL,CD,WU(WSIZE,3),MODN,MODR4,
     &     GAMMAW(WSIZE),CLSEC(WSIZE),CHORD1,CHORD2,ANGLE,RELAX,
     &     CHORDAV,GB(WSIZE),GTV(WSIZE),DY(WSIZE),Y(WSIZE),LIFTW,DRAGW
C
      INTEGER NNODE,NPANEL,NWPANEL,PNL(SIZE,4),WPNL(WSIZE,2),N1,N2,N3,
     &     N4,MK,ML,NU,NL,PARAM,NN(4),NUM,I,J,K,L,M,N,J1,K1,L1,J2,K2,L2,
     &     MRING,NMAT,ORDER(SIZE),NSIZE,NPTIP,SP
C
C************************DEFINITION OF VARIABLES************************
C
C     NNODE      : NO. OF NODES.
C     NPANEL     : NO. OF PANELS.
C     RN(j,*)    : POSITION VECTOR OF THE j^TH NODE.
C     PNL(j,*)   : CONNECTIVITY LISTING FOR THE j^TH BODY PANEL. NODE
C                  NUMBERS ARE LISTED CLOCKWISE, AS VIEWED FROM OUTSIDE
C                  OF BODY.
C     WPNL(j,*)  : CONNECTIVITY LISTING FOR THE j^TH WAKE PANEL. NODE
C                  NUMBERS ARE LISTED CLOCKWISE, AS VIEWED FROM ABOVE 
C                  THE WAKE.
C     QINF(*)    : FREE STREAM VELOCITY VECTOR.
C     MODQINF    : MODULUS OF QINF(*).
C     PITCH      : PITCH ANGLE OF BODY, RADIANS.
C     YAW        : YAW ANGLE OF BODY, RADIANS.
C     NOR(j,*)   : OUTWARD NORMAL VECTOR OF THE j^TH BODY PANEL OR
C                  UPWARD NORMAL OF THE j^TH WAKE PANEL.
C     RC(j,*)    : POSITION VECTOR OF THE COLLOCATION POINT OF THE j^TH
C                  PANEL, (BODY AND WAKE).
C     AREA(j)    : AREA OF THE j^TH PANEL (BODY AND WAKE).
C                  (FOR NOR(j,*), RC(j,*) & AREA(j,*), THE BODY PANELS
C                  PRECEED THE WAKE PANELS).
C     R1(*)...   : DUMMY POSITION VECTORS.
C       ...R4(*)
C     N1 ... N4  : DUMMY NODE NUMBERS.
C     NORM(*)    : DUMMY OUTWARD NORMAL VECTOR
C     MODNORM    : MODULUS OF NORM(*)
C     RCOL(*)    : DUMMY COLLOCATION POINT POSITION VECTOR
C     G          : DUMMY VORTEX STRENGTH, UNITY
C     Q(*)       : VELOCITY VECTOR
C     DQ(*)      : INCREMENTAL ADDITIONS TO Q(*)
C     A(i,j)     : INFLUENCE COEFFICIENT; DUE TO UNIT STRENGTH VORTEX 
C                  FILAMENTS OF THE j^TH PANEL AT THE COLLOCATION POINT
C                  OF THE i^TH PANEL
C     B(j)       : -1*COMPONENT OF FREE STREAM VELOCITY RESOLVED NORMAL
C                  TO j^TH PANEL
C     ALPHA(*,*) : LOWER FACTORISED MATRIX
C     BETA(*,*)  : UPPER FACTORISED MATRIX
C     DELTA(*)   : INTERMEDIATE SOLUTION ARRAY
C     GAMMA(*)   : VORTEX FILAMENT STRENGTH ARRAY, FINAL SOLUTION
C     VEL(j,*)   : VELOCITY VECTOR AT THE j^TH NODE
C     SVEL(j,*)  : SURFACE VELOCITY VECTOR AT THE j^TH NODE
C     SPEED(j)   : SURFACE SPEED AT THE j^TH NODE
C     CP(j)      : PRESSURE COEFFICIENT, Cp, AT j^TH NODE
C     NOR1(*),   : DUMMY OUTWARD NORMAL VECTORS
C       NOR2(*)
C
C***********************************************************************
C
      PI = 2*ACOS(0.)
      RELAX = 0.5
C
C**   DATA INPUT FROM 'mesh.dat', 'stream.dat' AND 'auxiliary.dat', AND
C**   PREPARE 'progress.dat' FOR RUNNING COMMENTRY
C
      OPEN (UNIT=1,FILE='mesh.dat',STATUS='OLD')
      OPEN (UNIT=2,FILE='stream.dat',STATUS='OLD')
C      OPEN (UNIT=3,FILE='progress.dat',STATUS='UNKNOWN')
      OPEN (UNIT=4,FILE='auxiliary.dat',STATUS='OLD')
      REWIND (1)
      REWIND (2)
C      REWIND (3)
      REWIND (4)
C
      READ (1,*) NNODE
      DO 25 I = 1,NNODE
         READ (1,*) (RN(I,J), J = 1,3)
 25   CONTINUE
C
      READ (1,*) NPANEL
      DO 50 I = 1,NPANEL
         READ (1,*) (PNL(I,J), J = 1,4)
 50   CONTINUE
C
      READ (1,*) NWPANEL
      DO 75 I = 1,NWPANEL
         READ (1,*) (WPNL(I,J), J=1,2)
 75   CONTINUE
      CLOSE (1)
C
      READ (2,*) MODQINF
      READ (2,*) PITCH
      READ (2,*) YAW
      READ (2,*) RHO
      CLOSE (2)
C
      PITCH = PITCH*PI/180
      YAW = YAW*PI/180
      QINF(1) = MODQINF*COS(PITCH)*COS(YAW)
      QINF(2) = MODQINF*COS(PITCH)*SIN(YAW)
      QINF(3) = MODQINF*SIN(PITCH)
      DO 100 I = 1,3
         QU(I) = QINF(I)/MODQINF
 100  CONTINUE
      HEAD = 0.5*RHO*(MODQINF**2)
C
      READ (4,*) SAREA
      READ (4,*) WANGLE
      READ (4,*) MRING
      READ (4,*) NPTIP
      CLOSE (4)
C
      WANGLE = WANGLE*PI/180
C
C      WRITE (3,*) 'DATA INPUT COMPLETE'
C
C**   DETERMINE THE AREA, OUTWARD NORMAL VECTOR AND COLLOCATION POINT 
C**   (CENTROID) OF EACH BODY PANEL
C
      DO 175 I = 1,NPANEL
         N1=PNL(I,1)
         N2=PNL(I,2)
         N3=PNL(I,3)
         N4=PNL(I,4)
         DO 125 J = 1,3
            R1(J) = RN(N1,J)
            R2(J) = RN(N2,J)
            R3(J) = RN(N3,J)
            R4(J) = RN(N4,J)
 125     CONTINUE
         IF ((N1.EQ.N2).OR.(N1.EQ.N4)) THEN
            CALL GEOTRI(R2,R3,R4,NORM,RCOL,AREA(I))
         ELSEIF (N2.EQ.N3) THEN
            CALL GEOTRI(R1,R3,R4,NORM,RCOL,AREA(I))
         ELSEIF (N3.EQ.N4) THEN
            CALL GEOTRI(R1,R2,R4,NORM,RCOL,AREA(I))
         ELSE
            CALL GEOQUAD(R1,R2,R3,R4,NORM,RCOL,AREA(I))
         ENDIF
         DO 150 J = 1,3
            NOR(I,J) = NORM(J)
            RC(I,J) = RCOL(J)
 150     CONTINUE
 175  CONTINUE
C
C      WRITE (3,*) 'BODY GEOMETRY CALCULATIONS COMPLETE'
C
C**   DETERMINE WAKE PANEL GEOMETRIES: UPWARD NORMAL AND COLLOCATION
C**   POINTS
C
      DO 250 I = 1,NWPANEL
         N1 = WPNL(I,1)
         N2 = WPNL(I,2)
         DO 190 K = 1,3
            R1(K) = RN(N1,K)
            R2(K) = RN(N2,K)
            R3(K) = 0.5*(R1(K)+R2(K))
            R4(K) = R2(K)-R1(K)
 190     CONTINUE
         DX1 = RN(N1,1)-RN(N1-MRING/2,1)
         DY1 = RN(N1,2)-RN(N1-MRING/2,2)
         A1 = ATAN(DY1/DX1)
         DX2 = RN(N2,1)-RN(N2-MRING/2,1)
         DY2 = RN(N2,2)-RN(N2-MRING/2,2)
         A2 = ATAN(DY2/DX2)
         ANGLE = RELAX*(WANGLE+0.5*(A1+A2))+(1-RELAX)*PITCH
         CHORDAV = 0.5*(SQRT(DX1**2+DY1**2)+SQRT(DX1**2+DY1**2))
         WU(I,1) = COS(ANGLE)
         WU(I,2) = 0.
         WU(I,3) = SIN(ANGLE)
         DO 195 K = 1,3
            TU(K) = WU(I,K)
 195     CONTINUE
         MODR4 = SQRT(R4(1)**2+R4(2)**2+R4(3)**2)
         CALL CROSS(NORM,TU,R4)
         MODN = SQRT(NORM(1)**2+NORM(2)**2+NORM(3)**2)
         IF ((CHORDAV*0.01).LT.(MODR4*0.1)) THEN
            TEMP = CHORDAV*0.01
         ELSE
            TEMP = MODR4*0.1
         ENDIF
         DO 200 K = 1,3
            NOR(NPANEL+I,K) = NORM(K)/MODN
            RC(NPANEL+I,K) = R3(K)+TU(K)*TEMP
 200     CONTINUE
         AREA(NPANEL+I) = MODR4*TEMP*2
 250  CONTINUE
C
C      WRITE (3,*) 'WAKE GEOMETRY CACULATIONS COMPLETE'
C
C**   DETERMINE PERTURBATIONS AT ALL COLLOCATION POINTS (BODY AS WELL
C**   AS WAKE PANELS) DUE TO ALL UNIT STRENGTH (BODY) VORTEX RINGS.
C
      G = 1.
      DO 400 I = 1,NPANEL+NWPANEL
         DO 305 M = 1,3
            RCOL(M) = RC(I,M)
 305     CONTINUE
         DO 375 J = 1,NPANEL
            DO 310 K = 1,3
               Q(K) = 0.
 310        CONTINUE
            DO 340 K = 1,4
               L = K+1
               IF (L.GT.4) THEN
                  L = L-4
               ENDIF
               N1 = PNL(J,K)
               N2 = PNL(J,L)
               IF (N1.EQ.N2) GOTO 340
               DO 320 M = 1,3
                  R1(M) = RN(N1,M)
                  R2(M) = RN(N2,M)
 320           CONTINUE
               CALL VORTEX(R1,R2,RCOL,DQ,G)
               DO 330 M = 1,3
                  Q(M) = Q(M) + DQ(M)
 330           CONTINUE
 340        CONTINUE
            DO 350 K = 1,3
                  AVEC(I,J,K) = Q(K)
 350        CONTINUE
 375     CONTINUE
C
C         WRITE (3,*) 'VORTEX RING PERTURBATION CALCULATIONS ON PANEL',I,
C     &'  DONE'
C
 400  CONTINUE
C
C      WRITE (3,*) 'ALL VORTEX RING CAUSED PERTURBATIONS FOUND'
C
C**   LOOP THROUGH ALL WAKE PANELS DETERMINING THEIR INFLUENCES.  FIRSTLY
C**   FIND EACH WAKE PANEL'S UPPER AND LOWER SURFACE BORDERING TRAILING
C**   EDGE PANELS.
C
      DO 480 J = 1,NWPANEL
         N1 = WPNL(J,1)
         N2 = WPNL(J,2)
         DO 410 K = 1,3
            R1(K) = RN(N1,K)
            R2(K) = RN(N2,K)
            TU(K) = WU(J,K)
 410     CONTINUE
         M = 0
         DO 430 I = 1,NPANEL
            IF (M.EQ.2) GOTO 440
            DO 420 K = 1,4
               L = K+1
               IF (L.GT.4) THEN
                  L = L-4
               ENDIF
               N3 = PNL(I,K)
               N4 = PNL(I,L)
               IF ((N1.EQ.N4).AND.(N2.EQ.N3)) THEN
                  NU = I
                  M = M+1
                  GOTO 430
               ELSEIF ((N1.EQ.N3).AND.(N2.EQ.N4)) THEN
                  NL = I
                  M = M+1
                  GOTO 430
               ENDIF
 420        CONTINUE
 430     CONTINUE
 440     CONTINUE
C
C**   NOW DETERMINE THE INFLUENCE OF THE J^TH WAKE PANEL ON ALL OTHER
C**   PANELS, I
C
         DO 470 I = 1,NPANEL+NWPANEL
            DO 450 K = 1,3
               RCOL(K) = RC(I,K)
 450        CONTINUE
            G = 1.
            CALL VORTEX(R1,R2,RCOL,DQ1,G)
            CALL TRAILING(R2,RCOL,TU,DQ2,G)
            G = -1.
            CALL TRAILING(R1,RCOL,TU,DQ3,G)
            DO 460 K = 1,3
               DVEC(I,J,K) = DQ1(K)+DQ2(K)+DQ3(K)
               AVEC(I,NU,K) = AVEC(I,NU,K)+DVEC(I,J,K)
               AVEC(I,NL,K) = AVEC(I,NL,K)-DVEC(I,J,K)
 460        CONTINUE
 470     CONTINUE
C     
C         WRITE (3,*) 'INFLUENCE PERTURBATION CALCULATIONS OF THE',J,'
C     &^TH WAKE PANEL DONE'
C
 480  CONTINUE
C
C      WRITE (3,*) 'ALL WAKE INFLUENCE PERTURBATIONS DETERMINED'
C
C**   DETERMINE ALL A(I,J), I.E. MATRIX A
C
      DO 502 I = 1,NPANEL+NWPANEL
         DO 501 J = 1,NPANEL
           A(I,J) = 0.
           DO 500 K = 1,3
              A(I,J) = A(I,J)+AVEC(I,J,K)*NOR(I,K)
 500       CONTINUE
           A(I,J) = A(I,J)*AREA(I)
 501     CONTINUE
 502  CONTINUE
C
C      WRITE (3,*) 'ALL VORTEX INFLUENCE COEFFICIENTS DETERMINED, MATRIX
C     & A KNOWN'
C
C**   DETERMINE SOURCE STRENGTHS BY ALLOWING EACH SOURCE TO SOLELY
C**   CANCEL WITH THE COMPONENT OF THE FREE STREAM NORMAL TO ITS PANEL
C
      DO 510 I = 1,NPANEL
         SIGMA(I) = -2*(QINF(1)*NOR(I,1)+QINF(2)*NOR(I,2)+
     &              QINF(3)*NOR(I,3))
 510  CONTINUE
C
C      WRITE (3,*) 'SOURCE STRENGTHS KNOWN'
C
C**   DETERMINE ALL INFLUENCE COEFFICIENTS FOR CONSTANT (UNIT) 
C**   STRENGTH SOURCES
C
      S = 1.
      DO 580 I = 1,NPANEL+NWPANEL
         DO 570 J = 1,NPANEL
            DO 515 K = 1,3
               CVEC(I,J,K) = 0.
 515        CONTINUE
            IF (I.EQ.J) GOTO 570
            DO 520 K = 1,4
               NN(K) = PNL(J,K)
 520        CONTINUE
            L = 0
            NUM = 4
            DO 540 K = 1,4
               L = L+1
               N = K+1
               IF (N.GT.4) THEN
                  N = N-4
               ENDIF
               IF (NN(K).EQ.NN(N)) THEN
                  NUM = 3
                  L = L-1
                  GOTO 540
               ELSE
                  DO 530 M = 1,3
                     R(L,M) = RN(NN(K),M)
 530              CONTINUE
               ENDIF
 540        CONTINUE
            DO 550 K = 1,3
               R(NUM+1,K) = RC(I,K)
               RCOL(K) = RC(J,K)
               NORM(K) = NOR(J,K)
 550        CONTINUE
            PARAM = 1
            CALL TRANSFORM(R,RCOL,NORM,DQ,PARAM,NUM)
            CALL SOURCE(R,DQ,S,NUM)
            PARAM = 2
            CALL TRANSFORM(R,RCOL,NORM,DQ,PARAM,NUM)
            DO 560 K = 1,3
               CVEC(I,J,K) = DQ(K)
 560        CONTINUE
 570     CONTINUE
C
C         WRITE (3,*) 'SOURCE INFLUENCE PERTURBATION CALCULATIONS ON
C     & PANEL',I,'  DONE'
C
 580  CONTINUE
C
C      WRITE (3,*) 'ALL SOURCE INFLUENCE PERTURBATIONS FOUND'
C
C**   DETERMINE ALL C(I,J), I.E. MATRIX C
C
      DO 595 I = 1,NPANEL+NWPANEL
         IF (I.GT.NPANEL) THEN
            TEMP = QINF(1)*NOR(I,1)+QINF(2)*NOR(I,2)+QINF(3)*NOR(I,3)
         ELSE
            TEMP = 0.
         ENDIF
         DO 590 J = 1,NPANEL
            C(I,J) = 0.
            DO 585 K = 1,3
                  C(I,J) = C(I,J)+CVEC(I,J,K)*NOR(I,K)
 585        CONTINUE
            TEMP = TEMP+C(I,J)*SIGMA(J)
 590     CONTINUE
         B(I) = -TEMP*AREA(I)
 595  CONTINUE
C
C      WRITE (3,*) 'ALL SOURCE INFLUENCE COEFFICIENTS DETERMINED, MATRIX
C     &C KNOWN'
C      WRITE (3,*) 'VECTOR B KNOWN'
C
C**   SET UP THE MATRIX EQUATION; E(*,*) GAMMA(*) = RHS(*).
C
      DO 665 I = 1,NPANEL
         TEMP = 0.
         DO 660 J = 1,NPANEL+NWPANEL
            TEMP = TEMP+A(J,I)*B(J)
 660     CONTINUE
         RHS(I) = TEMP
 665  CONTINUE
C
      DO 675 I = 1,NPANEL
         DO 675 J = 1,NPANEL
            TEMP = 0.
            DO 670 K = 1,NPANEL+NWPANEL
               TEMP = TEMP+A(K,I)*A(K,J)
 670        CONTINUE
            E(I,J) = TEMP
 675  CONTINUE
C
C**   DETERMINE THE GAMMA ARRAY BY USING GAUSSIAN ELIMINATION WITH
C**   PARTIAL PIVOTING.
C
      CALL GEPP(E,RHS,ORDER,NPANEL,SIZE)
      CALL BACKSOLVE(E,RHS,GAMMA,NPANEL,SIZE)
C
C      WRITE (3,*) 'GAMMA ARRAY KNOWN'
C
C**   NOW DETERMINE THE BOUNDARY CONDITION ENFORCEMENT ERRORS GVEC(*).
C
      DO 685 I = 1,NPANEL+NWPANEL
         TEMP = 0.
         DO 680 J = 1,NPANEL
            TEMP = TEMP+A(I,J)*GAMMA(J)
 680     CONTINUE
         GVEC(I) = TEMP-B(I)
 685  CONTINUE
C
C**   COMPUTE THE WAKE GAMMA ARRAY, GAMMAW(*).  FROM THIS WE MAY ALSO
C**   DETERMINE THE SECTIONAL LIFT COEFFICIENT, CLSEC(*)
C
      DO 700 J = 1,NWPANEL
         GAMMAW(J) = 0.
         M = 0
         N1 = WPNL(J,1)
         N2 = WPNL(J,2)
         DO 695 I = 1,NPANEL
            IF (M.EQ.2) GOTO 699
            DO 690 K = 1,4
               L = K+1
               IF (L.GT.4) THEN
                  L = L-4
               ENDIF
               N3 = PNL(I,K)
               N4 = PNL(I,L)
               IF (N3.EQ.N4) GOTO 690
               IF ((N1.EQ.N4).AND.(N2.EQ.N3)) THEN
                  GAMMAW(J) = GAMMAW(J)+GAMMA(I)
                  M = M+1
                  GOTO 695
               ELSEIF ((N1.EQ.N3).AND.(N2.EQ.N4)) THEN
                  GAMMAW(J) = GAMMAW(J)-GAMMA(I)
                  M = M+1
                  GOTO 695
               ENDIF
 690        CONTINUE
 695     CONTINUE
 699     CONTINUE
         CHORD1 = SQRT((RN(N1,1)-RN(N1-(MRING/2),1))**2+
     &        (RN(N1,3)-RN(N1-(MRING/2),3))**2)
         CHORD2 = SQRT((RN(N2,1)-RN(N2-(MRING/2),1))**2+
     &        (RN(N2,3)-RN(N2-(MRING/2),3))**2)
         CHORDAV = 0.5*(CHORD1+CHORD2)
         CLSEC(J) = 2*GAMMAW(J)/(MODQINF*CHORDAV)
 700  CONTINUE
C
C      WRITE (3,*) 'GAMMA WAKE ARRAY KNOWN'
C      WRITE (3,*) 'SECTIONAL CL DISTRIBUTION KNOWN'
C
C**   CALCULATE THE STRENGTH OF THE BOUND AND TRAILING VORTICIES.
C
      GB(1) = GAMMAW(1)
      GTV(1) = GAMMAW(1)
      DY(1) = 0.5*(RN(WPNL(1,2),2)-RN(WPNL(1,1),2))
      Y(1) = RN(WPNL(1,1),2)
C
      DO 705 I = 2,NWPANEL
         GB(I) = 0.5*(GAMMAW(I)+GAMMAW(I-1))
         GTV(I) = GAMMAW(I)-GAMMAW(I-1)
         Y(I) = RN(WPNL(I,1),2)
         DY(I) = 0.5*(RN(WPNL(I,2),2)-RN(WPNL(I-1,1),2))
 705  CONTINUE
C
      GB(NWPANEL+1) = GAMMAW(NWPANEL)
      GTV(NWPANEL+1) = -GAMMAW(NWPANEL)
      DY(NWPANEL+1) = 0.5*(RN(WPNL(NWPANEL,2),2)-RN(WPNL(NWPANEL,1),2))
      Y(NWPANEL+1) = RN(WPNL(NWPANEL,2),2)
C
C**   CALCULATE THE TOTAL LIFT (USING THE WAKE VORTICIES).
C
      TEMP = 0.
      DO 710 I = 1,NWPANEL+1
         TEMP = TEMP+GB(I)*DY(I)
 710  CONTINUE
      LIFTW = RHO*MODQINF*TEMP
      CLW = LIFTW/(HEAD*SAREA)
C
C**   CALCULATE THE TOTAL INDUCED DRAG (USING THE WAKE VORTICIES).
C
      TEMP = 0.
      DO 720 I = 1,NWPANEL+1
         DO 715 J = 1,NWPANEL+1
            IF (J.EQ.I) GOTO 715
            TEMP = TEMP+((GB(I)*GTV(J)*DY(I))/(Y(I)-Y(J)))
 715     CONTINUE
 720  CONTINUE
      DRAGW = (RHO/(4*PI))*TEMP
      CDW = DRAGW/(HEAD*SAREA)
C
C**   LOOP THROUGH ALL PANELS, I, DETERMINING Q(*) AND HENCE SVEL(I,*)
C**   AND CP(I) AT THE COLLOCATION POINT OF EACH, BY CONSIDERING THE 
C**   CONTIBUTIONS FROM THE FREE STREAM, ALL VORTICIES AND ALL SOURCES
C
      DO 900 I = 1,NPANEL
C
C**   CONTRIBUTION OF THE FREE STREAM
C
         DO 730 K = 1,3
            Q(K) = QINF(K)
 730     CONTINUE
C
C**   ALL SOURCE CONTRIBUTIONS EXCEPT THAT OF THE I^TH PANEL, AS IT
C**   WILL MAKE NO CONTRIBUTION TO SVEL(I,*), SINCE IT IS NORMAL TO THE
C**   I^TH PANEL
C
         DO 750 J = 1,NPANEL
            IF (J.EQ.I) GOTO 750
            DO 740 K = 1,3
               Q(K) = Q(K)+CVEC(I,J,K)*SIGMA(J)
 740        CONTINUE
 750     CONTINUE
C
C**   ALL VORTEX CONTRIBUTIONS INCLUDING THAT OF THE I^TH PANEL.  IF
C**   THE I^TH IS A TRAILING EDGE PANEL AVEC(I,I,*) WILL INCLUDE A 
C**   CONTRIBUTION FROM THE ADJACENT WAKE PANEL, WHICH WILL NOT BE
C**   NORMAL TO THE PANEL, UNLIKE THE 'SELF-INDUCTION' COMPONENT OF
C**   AVEC(I,I,*), AND WILL THEREFORE MAKE A CONTRIBUTION TO SVEL(I,*)
C
         DO 770 J = 1,NPANEL
            DO 760 K = 1,3
               Q(K) = Q(K)+AVEC(I,J,K)*GAMMA(J)
 760        CONTINUE
 770     CONTINUE
C
C**   NOW FIND ALL PANELS THAT BORDER THE I^TH PANEL AND FOR EACH SMEAR 
C**   THE VORTEX STRENGTH RESULTING FROM THE BORDERING FILAMENTS OF THE
C**   I^TH AND ITS NEIGHBOUR BETWEEN THEM.  IF AN EDGE OF THE I^TH PANEL
C**   HAPPENS TO BORDER THE WAKE MESH THEN THIS EDGE NEED NOT BE
C**   CONSIDERED AS THE KUTTA CONDITION FIXES THE CIRCULATION TO BE ZERO
C
         DO 772 K = 1,3
            DQ1(K) = 0.
 772     CONTINUE
         NUM = 4
         DO 830 K1 = 1,4
            L1 = K1+1
            IF (L1.GT.4) THEN
               L1 = L1-4
            ENDIF
            N1 = PNL(I,K1)
            N2 = PNL(I,L1)
            IF (N1.EQ.N2) THEN
               NUM = 3
               GOTO 830
            ENDIF
            DO 775 K = 1,NWPANEL
               N3 = WPNL(K,1)
               N4 = WPNL(K,2)
               IF (((N1.EQ.N4).AND.(N2.EQ.N3)).OR.((N1.EQ.N3).AND.
     &              (N2.EQ.N4))) GOTO 830
 775        CONTINUE
            DO 810 J = 1,NPANEL
               IF (J.EQ.I) GOTO 810
               DO 800 K2 = 1,4
                  L2 = K2+1
                  IF (L2.GT.4) THEN
                     L2 = L2-4
                  ENDIF
                  N3 = PNL(J,K2)
                  N4 = PNL(J,L2)
                  IF (N3.EQ.N4) GOTO 800
                  IF ((N1.EQ.N4).AND.(N2.EQ.N3)) THEN
                     DO 780 K = 1,3
                        R1(K) = RN(N1,K)
                        R2(K) = RN(N2,K)
                        NORM(K) = NOR(I,K)
 780                 CONTINUE
                     CALL SMEAR(R1,R2,NORM,AREA(I),AREA(J),GAMMA(I),
     &                    GAMMA(J),DQ)
                     DO 790 K = 1,3
                        DQ1(K) = DQ1(K)+DQ(K)
 790                 CONTINUE
                     GOTO 830
                  ENDIF
 800           CONTINUE
 810        CONTINUE
 830     CONTINUE
         DO 840 K = 1,3
            Q(K) = Q(K)+(2./NUM)*DQ1(K)
 840     CONTINUE
C
C**   ALL CONTRIBUTIONS NOW INCLUDED
C
C**   USE `SVEL = n^(Q^n)' TO DETERMINE TANGENTIAL (SURFACE) VELOCITY
C**   VECTOR FROM VELOCITY VECTOR
C
         DO 870 K = 1,3
            NORM(K) = NOR(I,K)
 870     CONTINUE
C
         CALL CROSS(DQ,Q,NORM)
         CALL CROSS(Q,NORM,DQ)
C
         DO 880 K = 1,3
            SVEL(I,K) = Q(K)
 880     CONTINUE
C
C**   FINALLY CALCULATE TANGENTIAL (SURFACE) SPEED AND HENCE Cp
C
         SPEED(I) = SQRT(SVEL(I,1)**2+SVEL(I,2)**2+SVEL(I,3)**2)
         CP(I) = 1.0 - ((SPEED(I)/MODQINF)**2)
C         WRITE (3,*) 'SURFACE VELOCITY AND Cp CALCULATION ON PANEL', I,
C     &               ' COMPLETE'
C
 900  CONTINUE
C
C      WRITE (3,*) 'SURFACE VELOCITY AND Cp CALCULATIONS COMPLETE'
C
C**   SECTIONAL CALCULATIONS & OUTPUT
C
      CALL SECTIONAL(RN,PNL,RC,SPEED,CP,NPANEL,MRING,NPTIP,SIZE)
C
C**   SPANWISE OUTPUT
C
      CALL SPANWISE(RN,PNL,RC,CLSEC,NPANEL,MRING,NPTIP,SIZE,WSIZE)
C
C**   THE FOLLOWING ROUTINE CALCULATES THE BODY FORCE, F(*), DUE TO
C     THE INTEGRATED PRESSURE DISTRIBUTION
C
      DO 950 J = 1,3
         TEMP = 0.
         DO 925 I = 1,NPANEL
            TEMP = TEMP+CP(I)*AREA(I)*NOR(I,J)
 925     CONTINUE
         F(J) = -HEAD*TEMP
 950  CONTINUE
C
C**   NOW RESOLVE THESE FORCES FOR LIFT, DRAG AND TRANSVERSE FORCE
C**   ALSO FIND THE LIFT AND DRAG COEFFICIENTS
C
      LIFT = F(3)*COS(PITCH)-F(1)*SIN(PITCH)
      DRAG = F(3)*SIN(PITCH)+F(1)*COS(PITCH)
      TRANS = F(2)
      CL = LIFT/(HEAD*SAREA)
      CD = DRAG/(HEAD*SAREA)
C
C      WRITE (3,*) 'BODY FORCE CALCULATIONS DONE'
C
C**   DATA OUTPUT STAGE
C
C**   USER'S DATA OUTPUT REQUIREMENTS
C
C      OPEN (UNIT=11,FILE='panels.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=12,FILE='gamma.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=13,FILE='sigma.dat',STATUS='UNKNOWN')
      OPEN (UNIT=14,FILE='force.dat',STATUS='UNKNOWN')
C      REWIND (11)
C      REWIND (12)
C      REWIND (13)
      REWIND (14)
C
C**   OUTPUT SPEED, U, V, W AND Cp DISTRIBUTIONS TO 'panels.dat' FOR 
C**   VIEWING BY TECPLOT.
C
C      WRITE (11,*) 'VARIABLES = "X","Y","Z","U","V","W","SPEED","Cp",
C     &     "GAMMA"'
C      WRITE (11,*) 'ZONE N=',NNODE,', E=',NPANEL
C      WRITE (11,*) 'F=FEPOINT, ET=QUADRILATERAL'
C      DO 1010 I = 1,NNODE
C         WRITE (11,*) (RN(I,L),L=1,3),0.,0.,0.,0.,0.,0.
C 1010 CONTINUE
C      DO 1020 I = 1,NPANEL
C         WRITE (11,*) (PNL(I,L),L=1,4)
C 1020 CONTINUE
C
C      WRITE (11,*) 'ZONE I=',NPTIP,',J=1,K=1,F=POINT'
C      DO 1030 I = 1,NPTIP
C         WRITE (11,*) (RC(I,J),J=1,3),(SVEL(I,J),J=1,3),SPEED(I),CP(I),
C     &        GAMMA(I)
C 1030 CONTINUE
C
C      WRITE (11,*) 'ZONE I=',MRING,',J=',(NPANEL-NPTIP*2)/MRING,
C     &',K=1,F=POINT'
C      DO 1040 I = NPTIP+1,NPANEL-NPTIP
C         WRITE (11,*) (RC(I,J),J=1,3),(SVEL(I,J),J=1,3),SPEED(I),CP(I),
C     &        GAMMA(I)
C 1040 CONTINUE
C
C      WRITE (11,*) 'ZONE I=',NPTIP,',J=1,K=1,F=POINT'
C      DO 1050 I = NPANEL-NPTIP+1,NPANEL
C         WRITE (11,*) (RC(I,J),J=1,3),(SVEL(I,J),J=1,3),SPEED(I),CP(I),
C     &        GAMMA(I)
C 1050 CONTINUE
C
C      WRITE (11,*) 'ZONE I=',NWPANEL,',J=1,K=1,F=POINT'
C      DO 1060 I = 1,NWPANEL
C         WRITE (11,*) 0.,RC(NPTIP+(I-1)*MRING+1,2),CLSEC(I),0.,0.,0.,
C     &        0.,0.,0.
C 1060 CONTINUE
C
C      L = NPTIP
C      DO 1080 I = 1,(NPANEL-2*NPTIP)/MRING
C         WRITE (11,*) 'ZONE I=',MRING+1,',J=1,K=1,F=POINT'
C         DO 1070 J = 1,MRING+1
C            L = L+1
C            IF (J.GT.MRING) THEN
C               WRITE (11,*) RC(L-MRING,1),RC(L-MRING,2),-CP(L-MRING),
C     &              0.,0.,0.,0.,0.,0.
C            ELSE
C               WRITE (11,*) RC(L,1),RC(L,2),-CP(L),0.,0.,0.,0.,0.,0.
C            ENDIF
C 1070    CONTINUE
C         L = L-1
C 1080 CONTINUE
C
C      CLOSE (11)
C
C**   OUTPUT THE GAMMA ARRAY TO 'gamma.dat'
C     
C      DO 1100 I = 1,NPANEL
C         WRITE (12,*) GAMMA(I)
C 1100 CONTINUE
C      CLOSE (12)
C
C**   OUTPUT THE SIGMA ARRAY TO 'sigma.dat'
C
C      TEMP = 0.
C      DO 1110 I = 1,NPANEL
C         WRITE (13,*) SIGMA(I)
C         TEMP = TEMP+AREA(I)*SIGMA(I)
C 1110 CONTINUE
C      WRITE (13,*) 'VOLUMETRIC FLOW RATE (SOURCES) =',TEMP
C      CLOSE (13)
C
C**   BODY FORCE DATA AND COEFFICIENTS TO 'force.dat'
C
      WRITE (14,1205)
 1205 FORMAT ('**********RESULTS OF PRESSURE INTEGRATION**********')
C
      WRITE (14,1210)
 1210 FORMAT (/,'BODY FORCES ALIGNED WITH THE X,Y & Z AXES RESPECTIVELY'
     &     )
      DO 1220 I = 1,3
         WRITE (14,*) F(I)
 1220 CONTINUE
C
      WRITE (14,1230)
 1230 FORMAT (/,'LIFT, DRAG & TRANSVERSE FORCES RESPECTIVELY')
      WRITE (14,*) LIFT
      WRITE (14,*) DRAG
      WRITE (14,*) TRANS
C
      WRITE (14,1240)
 1240 FORMAT (/,'LIFT & DRAG COEFFICIENTS RESPECTIVELY')
      WRITE (14,*) CL
      WRITE (14,*) CD
C
      WRITE (14,1245)
 1245 FORMAT (//,'**********RESULTS FROM WAKE CIRCULATIONS**********')
C
      WRITE (14,1250)
 1250 FORMAT (/,'LIFT & DRAG FORCES RESPECTIVELY')
      WRITE (14,*) LIFTW
      WRITE (14,*) DRAGW
C
      WRITE (14,1260)
 1260 FORMAT (/,'LIFT & DRAG COEFFICIENTS RESPECTIVELY')
      WRITE (14,*) CLW
      WRITE (14,*) CDW
C
      CLOSE (14)
C
C**   PROGRAMMER'S DATA OUTPUT REQUIREMENTS
C
C      OPEN (UNIT=21,FILE='geometry.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=22,FILE='matrixA.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=23,FILE='matrixC.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=24,FILE='matrixE.dat',STATUS='UNKNOWN')
C      OPEN (UNIT=25,FILE='gvec.dat',STATUS='UNKNOWN')
C      REWIND (21)
C      REWIND (22)
C      REWIND (23)
C      REWIND (24)
C      REWIND (25)
C
C**   OUTPUT BODY PANEL GEOMETRY TO 'geometry.dat'
C     
C      DO 1410 I = 1,NPANEL
C         WRITE (21,*) 'PANEL;', I
C         WRITE (21,*) 'AREA;', AREA(I)
C         WRITE (21,*) 'NORMALS;', (NOR(I,J),J=1,3)
C         WRITE (21,*) 'COLLOCATION;', (RC(I,J),J=1,3)
C 1410 CONTINUE
C
C      CLOSE (21)
C
C**   OUTPUT MATRIX A(I,J) AND VECTOR B(*) TO 'matrixA.dat'
C
C      WRITE (22,*) NPANEL
C      DO 1420 I = 1,NPANEL
C         WRITE (22,*) I
C         WRITE (22,*) (A(I,L),L=1,NPANEL)
C         WRITE (22,*) B(I)
C 1420 CONTINUE
C      CLOSE (22)
C
C**   OUTPUT MATRIX C(I,J) TO 'matrixC.dat'
C
C      WRITE (23,*) NPANEL
C      DO 1430 I = 1,NPANEL
C         WRITE (23,*) (C(I,L),L=1,NPANEL)
C 1430 CONTINUE
C      CLOSE (23)
C
C**   OUTPUT MATRIX E(I,J) AND VECTOR RHS(*) TO 'matrixE.dat'
C
C      WRITE (24,*) NPANEL
C      DO 1440 I = 1,NPANEL
C         WRITE (24,*) I 
C         WRITE (24,*) (E(I,L),L=1,NPANEL)
C         WRITE (24,*) RHS(I)
C 1440 CONTINUE
C      CLOSE (24)
C
C**   OUTPUT BOUNDARY CONDITION ENFORCEMNET ERRORS TO 'gvec.dat'
C
C      DO 1450 I = 1,NPANEL+NWPANEL
C         WRITE (25,*) GVEC(I)
C 1450 CONTINUE
C      CLOSE (25)
C
C**   DATA OUTPUT COMPLETE
C
C      WRITE (3,*) 'DATA OUTPUT COMPLETE'
C      WRITE (3,*) '**********END**********'
C
      STOP
      END
C
C*****************************SUBROUTINES*******************************
C
      SUBROUTINE GEOQUAD(R1,R2,R3,R4,NORM,RCOL,AREA)
C
C**   CALCULATES THE AREA (AREA), OUTWARD NORAML VECTOR (NORM) AND
C**   COLLOCATION POINT (RCOL) (CENTROID), OF A TRAPEZOID,
C**   VERTICES 1,2,3,4
C
      REAL R1(3),R2(3),R3(3),R4(3),R1CR2(3),R1CR3(3),R1CR4(3),R2CR3(3),
     &     R3CR4(3),ACB(3),MODACB,NORM(3),AR1,AR2,AREA,A1(3),A2(3),
     &     RC1(3),RC2(3),RCOL(3)
C
C**   CALCULATE ALL RELEVANT CROSS PRODUCTS
C
      CALL CROSS(R1CR2,R1,R2)
      CALL CROSS(R1CR3,R1,R3)
      CALL CROSS(R1CR4,R1,R4)
      CALL CROSS(R2CR3,R2,R3)
      CALL CROSS(R3CR4,R3,R4)
C
C**   NORMAL VECTOR
C
      DO 2000 I = 1,3
         ACB(I) = -R1CR2(I)+R1CR4(I)-R2CR3(I)-R3CR4(I)
 2000 CONTINUE
      MODACB = SQRT(ACB(1)**2+ACB(2)**2+ACB(3)**2)
      DO 2025 I = 1,3
         NORM(I) = ACB(I)/MODACB
 2025 CONTINUE
C
C**   AREAS
C
      DO 2050 I = 1,3
         A1(I) = R1CR2(I)-R1CR3(I)+R2CR3(I)
         A2(I) = -R1CR3(I)+R1CR4(I)-R3CR4(I)
 2050 CONTINUE
      AR1 = 0.5*SQRT(A1(1)**2+A1(2)**2+A1(3)**2)
      AR2 = 0.5*SQRT(A2(1)**2+A2(2)**2+A2(3)**2)
      AREA = AR1+AR2
C
C**   COLLOCATION POINT, CENTROID
C
      DO 2075 I = 1,3
         RC1(I) = (R1(I)+R2(I)+R3(I))/3
         RC2(I) = (R1(I)+R3(I)+R4(I))/3
         RCOL(I) = (1/AREA)*(AR1*RC1(I)+AR2*RC2(I))
 2075 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE GEOTRI(R1,R2,R3,NORM,RCOL,AREA)
C
C**   CALCULATES THE AREA (AREA), OUTWARD NORAML VECTOR (NORM) AND
C**   COLLOCATION POINT (RCOL) (CENTROID), OF A TRIANGLE, VERTICES 1,2,3
C
      REAL R1(3),R2(3),R3(3),R1CR2(3),R1CR3(3),R2CR3(3),ACB(3),MODACB,
     &     NORM(3),AREA,RCOL(3)
C
C**   CALCULATE ALL RELEVANT CROSS PRODUCTS
C
      CALL CROSS(R1CR2,R1,R2)
      CALL CROSS(R1CR3,R1,R3)
      CALL CROSS(R2CR3,R2,R3)
C
C**   NORMAL VECTOR
C
      DO 2100 I = 1,3
         ACB(I) = -R1CR2(I)+R1CR3(I)-R2CR3(I)
 2100 CONTINUE
      MODACB = SQRT(ACB(1)**2+ACB(2)**2+ACB(3)**2)
      DO 2125 I = 1,3
         NORM(I) = ACB(I)/MODACB
 2125 CONTINUE
C
C**   AREAS
C
      AREA = 0.5*MODACB
C
C**   COLLOCATION POINT, CENTROID
C
      DO 2150 I = 1,3
         RCOL(I) = (R1(I)+R2(I)+R3(I))/3
 2150 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE CROSS(A,B,C)
C
C**   PERFORMS THE VECTOR OPERATION, A = B^C, AND RETURNS THE VECTOR A
C
      REAL A(3),B(3),C(3)
      INTEGER I,J,K
C
      DO 2200 I = 1,3
         J = I+1
         IF (J.GT.3) THEN
            J = J-3
         ENDIF
         K = I+2
         IF (K.GT.3) THEN
            K = K-3
         ENDIF
         A(I) = B(J)*C(K)-B(K)*C(J)
 2200 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE VORTEX(P1,P2,P,DQ,G)
C
C**   CALCULATES INDUCED VELOCITIES DQ, AT POINT P, DUE TO FINITE LENGTH
C**   STRAIGHT VORTEX FILAMENT BETWEEN POINTS 1 AND 2, POSITION VECTORS
C**   P1, P2.
C
      REAL P1(3),P2(3),P(3),DQ(3),G,R1(3),R2(3),R0(3),CORE,R1CR2(3),
     &     MODR1CR2,MODR1,MODR2,TEMP,R0R1,R0R2
C
      PI = 2*ACOS(0.0)
      CORE = 1.0E-16
C
C**   CALCULATE RELATIVE POSITION VECTORS
C
      DO 2300 I = 1,3
         R1(I) = P(I)-P1(I)
         R2(I) = P(I)-P2(I)
         R0(I) = R1(I)-R2(I)
 2300 CONTINUE
C
C**   CALCULATE CROSS PRODUCT
C
      CALL CROSS(R1CR2,R1,R2)
C
C**   CALCULATE MODULI
C
      MODR1CR2 = SQRT(R1CR2(1)**2+R1CR2(2)**2+R1CR2(3)**2)
      MODR1 = SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
      MODR2 = SQRT(R2(1)**2+R2(2)**2+R2(3)**2)
C
C**   DOES VORTEX FILAMENT COINCIDE WITH P
C
      IF ((MODR1.LT.CORE).OR.(MODR2.LT.CORE).OR.((MODR1CR2**2).LT.CORE))
     &   THEN
         DO 2325 I = 1,3
            DQ(I) = 0.
 2325    CONTINUE
         GOTO 2375
      ENDIF
C
C**   DOT PRODUCTS
C
      R0R1 = R0(1)*R1(1)+R0(2)*R1(2)+R0(3)*R1(3)
      R0R2 = R0(1)*R2(1)+R0(2)*R2(2)+R0(3)*R2(3)
C
C**   INCREMENTAL VELOCITIES
C
      TEMP = (G/(4*PI*(MODR1CR2**2)))*((R0R1/MODR1)-(R0R2/MODR2))
      DO 2350 I = 1,3
         DQ(I) = TEMP*R1CR2(I)
 2350 CONTINUE
C
 2375 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE TRANSFORM(RN,RC,N,DQ,PARAM,NUM)
C
C**   TRANSFORMS THE GLOBAL COORDINATE SYSTEM (INFLUENCING) PANEL NODES
C**   AND POINT OF INFLUENCE (COLLOCATION POINT OF THE INFLUENCED PANEL)
C**   COORDINATES, RN(*,*), INTO THE INFLUENCING PANEL'S LOCAL 
C**   COORDINATE SYSTEM (BEING, LYING IN THE XY PLANE, AT Z = 0), VIA 2
C**   ROTATIONS AND A TRANSLATION.  IT ALSO ALLOWS THE VELOCITY VECTOR 
C**   TO BE ROTATED FROM THE LOCAL SYSTEM BACK TO THE GLOBAL, BY 
C**   TREATING IT AS A POSITION VECTOR.  NOTE THAT RC(*) IS THE 
C**   COLLOCATION POINT OF THE INFLUENCING PANEL, NOT THE INFLUENCED
C**   PANEL.
C
      REAL PI,RN(5,3),RC(3),N(3),M(3),L(3),F(3),G(3),H(3),A(3),B(3),
     &     C(3),MODM,MODS1,MODS2,S(2,3),THETA(2),FS,MODC,MODH,DQ(3),
     &     SMALL
      INTEGER PARAM,I1,I2,I3,J2,NUM
C
      PI = 2*ACOS(0.)
      SMALL = 1.0E-8
C
C**   SET LOOP PARAMETERS TO DEFAULT VALUES, SO THAT THEY MAY BE
C**   SUBSEQUENTLY CHANGED IF DESIRED
C
      I1 = 1
      I2 = 2
      I3 = 1
      J2 = NUM+1
C
C**   DETERMINE ROTATION AXES AND ANGLES
C
      MODM = SQRT(N(2)**2+N(3)**2)
      M(1) = 0.
      M(2) = N(2)/MODM
      M(3) = N(3)/MODM
      L(1) = 0.
      L(2) = 0.
      L(3) = 1.
C
      THETA(1) = ATAN(ABS(N(1))/MODM)
      IF ((N(3)).LT.(0.)) THEN
         THETA(2) = PI+ATAN(ABS(N(2))/N(3))
      ELSE
         THETA(2) = ATAN(ABS(N(2))/N(3))
      ENDIF
C
      IF (N(1).EQ.(0.)) THEN
         I1 = 2
      ELSE
         CALL CROSS(A,M,N)
         MODS1 = SQRT(A(1)**2+A(2)**2+A(3)**2)
         DO 2900 I = 1,3
            S(1,I) = A(I)/MODS1
 2900    CONTINUE
      ENDIF
C
      IF (N(2).EQ.(0.)) THEN
         S(2,1) = 1.
         S(2,2) = 0.
         S(2,3) = 0.
      ELSE
         CALL CROSS(B,L,M)
         MODS2 = SQRT(B(1)**2+B(2)**2+B(3)**2)
         DO 2910 I = 1,3
            S(2,I) = B(I)/MODS2
 2910    CONTINUE
      ENDIF
C
C**   IF PARAM = 1 THEN TRANSFORM GLOBAL TO LOCAL COORDINATES, ELSE IF 
C**   PARAM = 2 THEN TRANSFORM LOCAL VELOCITY VECTOR TO GLOBAL VELOCITY
C**   VECTOR
C
      IF (PARAM.EQ.2) THEN
         THETA(1) = -THETA(1)
         THETA(2) = -THETA(2)
         I1 = 2
         I2 = 1
         I3 = -1
         J2 = 1
      ENDIF
C
C**   ROTATIONS ANTICLOCKWISE ABOUT AXIS S(1,*) THEN ABOUT AXIS S(2,*)
C**   BY ANGLES OF THETA(1) AND THETA(2) RESPECTIVELY.  THE ORDER OF
C**   THESE ROTATIONS IS REVERSED IF THE VELOCITY VECTOR IS BEING 
C**   TRANSFORMED
C
      DO 3100 I = I1,I2,I3
         DO 3050 J = 1,J2
            DO 2950 K = 1,3
               IF (PARAM.EQ.1) THEN
                  F(K) = RN(J,K)-RC(K)
               ELSE
                  F(K) = DQ(K)
               ENDIF
 2950       CONTINUE
            FS = F(1)*S(I,1)+F(2)*S(I,2)+F(3)*S(I,3)
            DO 3000 K = 1,3
               G(K) = FS*S(I,K)
               H(K) = F(K)-G(K)
               A(K) = S(I,K)
 3000       CONTINUE
            CALL CROSS(C,A,H)
            MODH = SQRT(H(1)**2+H(2)**2+H(3)**2)
            IF (MODH.LE.SMALL) GOTO 3050
            MODC = SQRT(C(1)**2+C(2)**2+C(3)**2)
            DO 3025 K = 1,3
               B(K) = H(K)/MODH
               C(K) = C(K)/MODC
               IF (PARAM.EQ.1) THEN
                  RN(J,K) = MODH*(B(K)*COS(THETA(I))-C(K)*
     &                 SIN(THETA(I)))+RC(K)+G(K)
               ELSE
                  DQ(K) = MODH*(B(K)*COS(THETA(I))-C(K)*SIN(THETA(I)))+
     &                 G(K)
               ENDIF
 3025       CONTINUE
 3050    CONTINUE
 3100 CONTINUE
C
C**   LATERAL TRANSLATION FOR COORDINATE TRANSFORMATION ONLY
C
      IF (PARAM.EQ.1) THEN
         DO 3200 J = 1,NUM+1
            RN(J,3) = RN(J,3)-RC(3)
 3200    CONTINUE
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE SOURCE(RN,DQ,S,NUM)
C
C**   DETERMINES THE VELOCITY PERTURBATIONS DQ, AT A POINT RN(NUM+1,*),
C**   DUE TO A QUADRILATERAL/TRIANGULAR (NO. OF VERTICES = NUM) CONSTANT
C**   STRENGTH SOURCE ELEMENT, STRENGTH PER UNIT AREA S.  THE
C**   CALCULATIONS ARE FORMED IN THE PANELS LOCAL FRAME OF REFERENCE.
C
      REAL PI,RN(5,3),DQ(3),S,TEMP,CONST(4),R(4),E(4),H(4),D(4),M(4),
     &     TEMP2
      INTEGER NUM,I,J,K,L
C
      PI = 2*ACOS(0.)
C
C**   DETERMINE CONSTANT ARRAYS; R,E,H,D AND M
C
      DO 3300 I = 1,NUM
         J = I+1
         IF (J.GT.NUM) THEN
            J = J-NUM
         ENDIF
         D(I) = SQRT((RN(J,1)-RN(I,1))**2+(RN(J,2)-RN(I,2))**2)
         M(I) = (RN(J,2)-RN(I,2))/(RN(J,1)-RN(I,1))
         R(I) = SQRT((RN(NUM+1,1)-RN(I,1))**2+(RN(NUM+1,2)-RN(I,2))**2+
     &          RN(NUM+1,3)**2)
         E(I) = (RN(NUM+1,1)-RN(I,1))**2+RN(NUM+1,3)**2
         H(I) = (RN(NUM+1,1)-RN(I,1))*(RN(NUM+1,2)-RN(I,2))
 3300 CONTINUE
C
C**   DETERMINE CONSTANT ARRAY; CONST(*)
C
      DO 3350 I = 1,NUM
         J = I+1
         IF (J.GT.NUM) THEN
            J = J-NUM
         ENDIF
         TEMP = (R(I)+R(J)-D(I))/(R(I)+R(J)+D(I))
         CONST(I) = (1/D(I))*LOG(TEMP)
 3350 CONTINUE
C
C**   DETERMINE U AND V COMPONENTS OF PERTURBATION VELOCITY
C
      DO 3450 K = 1,2
         L = K+1
         IF (L.GT.2) THEN
            L = L-2
         ENDIF
         TEMP = 0.
         DO 3400 I = 1,NUM
            J = I+1
            IF (J.GT.NUM) THEN
               J = J-NUM
            ENDIF
            TEMP = TEMP+((-1.)**L)*(RN(J,L)-RN(I,L))*CONST(I)
 3400    CONTINUE
         DQ(K) = (S/(4*PI))*TEMP
 3450 CONTINUE
C
C**   FINALLY, DETERMINE W COMPONENT OF PERTURBATION VELOCITY
C
      TEMP = 0.
      DO 3500 I = 1,NUM
         J = I+1
         IF (J.GT.NUM) THEN
            J = J-NUM
         ENDIF
C
         TEMP2 = ATAN((M(I)*E(I)-H(I))/(RN(NUM+1,3)*R(I)))
         IF ((TEMP2.GE.0).OR.(TEMP2.LT.0)) THEN
            TEMP = TEMP+TEMP2
         ENDIF
C
         TEMP2 = -ATAN((M(I)*E(J)-H(J))/(RN(NUM+1,3)*R(J)))
         IF ((TEMP2.GE.0).OR.(TEMP2.LT.0)) THEN
            TEMP = TEMP+TEMP2
         ENDIF
C
 3500 CONTINUE
      DQ(3) = (S/(4*PI))*TEMP
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE TRAILING(R,RP,SU,DQ,GA)
C
C**   CALCULATES THE VELOCITY PERTURBATIONS,DQ(*) FELT AT RP(*) DUE TO A
C**   SEMI-INFINITE VORTEX FILAMENT, CIRCULATION GA, TRAILING FROM R(*)
C**   IN THE DIRECTION OF SU(*) TO DOWNSTREAM INFINITY
C
      REAL RP(3),R(3),SU(3),DQ(3),GA,F(3),G(3),H(3),FS,SCH(3),MODSCH,
     &     A(3),MODF,MODH,PI,CBETA,TEMP,CORE
C
      PI = 2*ACOS(0.)
      CORE = 1.0E-16
C
C**   DETERMINE VECTORS F(*), G(*) AND H(*)
C
      DO 3600 I = 1,3
         F(I) = RP(I)-R(I)
 3600 CONTINUE
      FS = F(1)*SU(1)+F(2)*SU(2)+F(3)*SU(3)
      DO 3625 I = 1,3
         G(I) = SU(I)*FS
         H(I) = F(I)-G(I)
 3625 CONTINUE
C
C**   IF MODH IS LESS THAN CORE SET VELOCITY PERTURBATIONS TO ZERO AND
C**   RETURN TO MAIN PROGRAM
C
      MODH = SQRT(H(1)**2+H(2)**2+H(3)**2)
      IF (MODH.LE.CORE) THEN
         DO 3640 I = 1,3
            DQ(I) = 0.
 3640    CONTINUE
         GOTO 3680
      ENDIF
C
C**   DETERMINE THE DIRECTION OF THE VELOCITY PERTURBATION, A(*)
C
      CALL CROSS(SCH,SU,H)
      MODSCH = SQRT(SCH(1)**2+SCH(2)**2+SCH(3)**2)
      DO 3650 I = 1,3
         A(I) = SCH(I)/MODSCH
 3650 CONTINUE
C
C**   DETERMINE COS(BETA) = CBETA
C
      MODF = SQRT(F(1)**2+F(2)**2+F(3)**2)
      CBETA = FS/MODF
C
C**   FINALLY CALCULATE THE MAGNITUDE OF THE PERTURBATION AND HENCE THE
C**   ORTHOGANOL VELOCITY PERTURBATIONS
C
      TEMP = (GA/(4*PI*MODH))*(1+CBETA)
      DO 3675 I = 1,3
         DQ(I) = A(I)*TEMP
 3675 CONTINUE
C
 3680 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE SMEAR(R1,R2,N,AR1,AR2,GAMMA1,GAMMA2,DQ)
C
C**   CALCULATES THE SURFACE VELOCITY DUE TO TWO COAXIAL VORTEX FILAMENTS
C**   THAT ARE ON BORDERING PANELS WHICH ARE COINCIDENT WITH A VERTEX 
C
      REAL R1(3),R2(3),DQ(3),AR1,AR2,GAMMA1,GAMMA2,S(3),T(3),N(3),MODS,
     &     MODT,TEMP
C
C**   CALCULATE VECTOR S(*)
C
      DO 3700 I = 1,3
         S(I) = R2(I) - R1(I)
 3700 CONTINUE
C
C**   DETERMINE NON-UNIT VECTOR T
C
      CALL CROSS(T,S,N)
C
C**   MODULI AND UNIT VECTOR T
C
      MODS = SQRT(S(1)**2+S(2)**2+S(3)**2)
      MODT = SQRT(T(1)**2+T(2)**2+T(3)**2)
      DO 3725 I = 1,3
         T(I) = T(I)/MODT
 3725 CONTINUE
C
C**   VELOCITY PERTURBATIONS
C
      TEMP = ((GAMMA1-GAMMA2)*MODS)/(AR1+AR2)
      DO 3750 I = 1,3
         DQ(I) = T(I)*TEMP
 3750 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE GEPP(A,B,ORDER,N,SIZE)
C
C**   THIS SUBROUTINE PERFROMS GAUSSIAN ELIMINATION WITH PARTIAL
C**   PIVOTING (ROW PIVOTING) TO REDUCE THE MATRIX EQUATION;
C**   A(*,*) X(*) = B(*), MATRIX SIZE IS N, TO A FORM SUCH THAT A(*,*)
C**   IS AN UPPER TRIANGULAR MATRIX ONLY.  EACH ELEMENT OF THE LEADING
C**   DIAGONAL OF THE NEW MATRIX WILL HAVE BE 1.
C
      INTEGER SIZE
      REAL A(SIZE,SIZE),B(SIZE),TEMP,VALUE
      INTEGER N,ORDER(SIZE),PIVOT,MOVE
C
C**   INITIALISE THE ORDER ARRAY SO THAT WE MAY KEEP TRACK OF THE ROWS
C
      DO 3800 I = 1,N
         ORDER(I) = I
 3800 CONTINUE
C
C**   LOOP THROUGH EACH ROW OF THE MATRIX USING EACH ROW FOR A COLUMN
C**   ELIMINATION
C
      DO 3950 I = 1,N
C
C**   DETERMINE THE PIVOT
C
         TEMP = 0.
         DO 3825 J = I,N
            VALUE = ABS(A(J,I))
            IF (VALUE.GT.TEMP) THEN
               PIVOT = J
               TEMP = VALUE
            ENDIF
 3825    CONTINUE
C
C**   SWOP THE PIVOT ROW WITH THAT IN THE ELIMINATION ROW POSITION
C
         MOVE = ORDER(I)
         ORDER(I) = ORDER(PIVOT)
         ORDER(PIVOT) = MOVE
         DO 3850 J = I,N
            TEMP = A(I,J)
            A(I,J) = A(PIVOT,J)
            A(PIVOT,J) = TEMP
 3850    CONTINUE
         TEMP = B(I)
         B(I) = B(PIVOT)
         B(PIVOT) = TEMP
C
C**   NORMALISE THE LEADING DIAGONAL TERM OF THE ELIMINATION ROW
C
         TEMP = A(I,I)
         DO 3875 J = I,N
            A(I,J) = A(I,J)/TEMP
 3875    CONTINUE
         B(I) = B(I)/TEMP
C
C**   ELIMINATE ALL ELEMENTS DIRECTLY BELOW THE LEADING DIAGONAL ELEMENT
C**   OF THE ELIMINATION ROW
C
         DO 3925 K = I+1,N
            TEMP = A(K,I)
            DO 3900 J = I,N
               A(K,J) = A(K,J)-TEMP*A(I,J)
 3900       CONTINUE
            B(K) = B(K)-TEMP*B(I)
 3925    CONTINUE
 3950 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE BACKSOLVE(A,B,X,N,SIZE)
C
C**   SOLVE UPPER TRIANGULAR MATRIX FOR X, THE SOLUTION, BY BACKWARD
C**   SUBSTITUTION
C
      INTEGER N,SIZE
      REAL A(SIZE,SIZE),B(SIZE),X(SIZE),TEMP
C
      DO 6000 I = N,1,-1
         TEMP = 0.0
         DO 5900 J = I+1,N
            TEMP = TEMP+A(I,J)*X(J)
 5900    CONTINUE
         X(I) = B(I)-TEMP
 6000 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE SECTIONAL(RN,PNL,RC,SPEED,CP,NPANEL,MRING,NPTIP,SIZE)
C
C**   CALCULATE THE Cp DISTRIBUTIONS AT THE SPECIFIED SECTIONS
C
      CHARACTER BASE*5,EXTEN*3,OUTPUT*8,FORMSTR*8
      INTEGER NPANEL,MRING,NPTIP,SIZE,PNL(SIZE,4)
      PARAMETER (NSECT=10,MRMAX=80)
      REAL RC(SIZE,3),SPEED(SIZE),CP(SIZE),YSECT(NSECT),
     &     XSECT(NSECT,MRMAX),CPSECT(NSECT,MRMAX),CHORD(NSECT),
     &     RN(SIZE,3),XLE(NSECT)
C
C**   READ IN THE SECTION LOCATIONS
C
      OPEN (UNIT=15,FILE='sect_inp',STATUS='UNKNOWN')
      REWIND (15)
C
      READ (15,*) NSEC
C
      DO IS = 1,NSEC
         READ (15,*) YSECT(IS)
      ENDDO
C
      CLOSE (15)
C
C**   NO.'S OF NODES & PANELS
C
      NS = (NPANEL-2*NPTIP)/MRING
      NST = (NS/2)*MRING+NPTIP
C
C**   LOOP OVER THE SECTIONS
C
      DO IS = 1,NSEC
C
C**   DETERMINE THE NEIGHBOURING PANELS
C
         YSEC = YSECT(IS)
         DO I = 1,NS
            N = NST+((I-1)*MRING)+1
            YTEST = RC(N,2)
            IF (YTEST.GE.YSEC) GOTO 100
         ENDDO
C
 100     CONTINUE
C
C**   INTERPLOATION FACTOR
C     
         IA = I-1
         IB = I
C     
         NA = NST+((IA-1)*MRING)+1
         NB = NST+((IB-1)*MRING)+1
C     
         YA = RC(NA,2)
         YB = RC(NB,2)
C
         SAB = YB-YA
         SAP = YSEC-YA
C
         FACTOR = SAP/SAB
C
C**   LOCAL CHORD
C
         MA1 = PNL(NA,1)
         MB1 = PNL(NB,1)
         MC1 = PNL(NB,2)
         MA2 = MA1+(MRING/2)
         MB2 = MB1+(MRING/2)
         MC2 = MC1+(MRING/2)
C
         CHA = SQRT((RN(MA1,1)-RN(MA2,1))**2+
     &        (RN(MA1,2)-RN(MA2,2))**2)
         CHB = SQRT((RN(MB1,1)-RN(MB2,1))**2+
     &        (RN(MB1,2)-RN(MB2,2))**2)
         CHC = SQRT((RN(MC1,1)-RN(MC2,1))**2+
     &        (RN(MC1,2)-RN(MC2,2))**2)
C
         CHAVA = (CHA+CHB)/2.
         CHAVB = (CHB+CHC)/2.
C
         CHORD(IS) = (1-FACTOR)*CHAVA+FACTOR*CHAVB
C
C**   LOCAL LEADING EDGE
C
         XLEA = RN(MA1,1)
         XLEB = RN(MB1,1)
         XLEC = RN(MC1,1)
C
         XLEAVA = (XLEA+XLEB)/2.
         XLEAVB = (XLEB+XLEC)/2.
C
         XLE(IS) = (1-FACTOR)*XLEAVA+FACTOR*XLEAVB

C**   CHORDWISE LOOP OVER THE SECTION
C
         NA = NA-1
         NB = NB-1
C
         DO J = 1,MRING
C
            NA = NA+1
            NB = NB+1
C
            XSECT(IS,J) = RC(NA,1)*(1-FACTOR)+RC(NB,1)*FACTOR
            CPSECT(IS,J) = CP(NA)*(1-FACTOR)+CP(NB)*FACTOR
C
         ENDDO
C
      ENDDO
C
C**   DATA OUTPUT
C
      BASE = 'sect_'
C
      DO IS = 1,NSEC
C
         L = LEN(EXTEN)
         WRITE (FORMSTR,'(''(I'',I2.2,''.'',I2.2,'')'')') L,L
         WRITE (EXTEN,FORMSTR) IS
C
         OUTPUT = BASE//EXTEN
C
         NUM = 50+IS
c
         OPEN (UNIT=NUM,FILE=OUTPUT,STATUS='UNKNOWN')
         REWIND (NUM)
C
         DO J = 1,MRING+1
            IF (J.NE.(MRING+1)) THEN
               XNORM = (XSECT(IS,J)-XLE(IS))/CHORD(IS)
               WRITE (NUM,*) XNORM,CPSECT(IS,J)
            ELSE
               XNORM = (XSECT(IS,1)-XLE(IS))/CHORD(IS)
               WRITE (NUM,*) XNORM,CPSECT(IS,1)
            ENDIF
         ENDDO
C
         CLOSE (NUM)
C
      ENDDO
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE SPANWISE(RN,PNL,RC,CLSEC,NPANEL,MRING,NPTIP,SIZE,WSIZE)
C
      INTEGER NPANEL,MRING,NPTIP,SIZE,WSIZE,PNL(SIZE,4)
      PARAMETER (NSECT=10,MRMAX=80)
      REAL RC(SIZE,3),RN(SIZE,3),CLSEC(WSIZE)
C
C**   CALCULATE THE SEMI-SPAN
C
      NBODY = NPANEL-(NPTIP*2)
      NS = NBODY/MRING
C
      N1 = NPTIP+1
      N2 = NPTIP+NBODY-MRING+1
      N3 = N1+(MRING/2)-1
      N4 = N2+(MRING/2)-1
C
      M1 = PNL(N1,1)
      M2 = PNL(N2,2)
      M3 = PNL(N3,4)
      M4 = PNL(N4,3)
C
      S1 = RN(M2,2)-RN(M1,2)
      S2 = RN(M4,2)-RN(M3,2)
C
      SPAN = (S1+S2)/2.
      SSPAN = SPAN/2.
C
C**   DATA OUTPUT
C
      OPEN (UNIT=31,FILE='span_load',STATUS='UNKNOWN')
      REWIND (31)
C
      DO IS = 1,NS
         WRITE (31,*) RC(NPTIP+(IS-1)*MRING+1,2)/SSPAN,CLSEC(IS)
      ENDDO
C
      CLOSE (31)
C
      RETURN
      END
