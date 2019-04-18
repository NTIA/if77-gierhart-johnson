      SUBROUTINE ASORP(FK,AO,AW)      
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
   19 FORMAT(5X,'FREQUENCY IS TOO HIGH FOR ABSORPTION TABLE USING VALUES        
     XFOR 100 GHZ')     
      DIMENSION ZX(53),ZW(53),FZ(53)                
      COMMON/PRINT/KL,IN,IOT          
          
      DATA FZ/.10,.15,.205,.30,.325,.35,.40,.55,.70,1.0,1.52,2.0,3.0,3.4        
     F,4.0,4.9,8.3,10.2,15.0,17.0,20.0,22.0,23.0,25.0,26.,30.,32.0,33.,3        
     F5.,37.,38.,40.,42.,43.,44.,47.,48.,51.,54.,58.,59.,60.,61.,62.,63.        
     F,64.,68.,70.,72.,76.,84.,90.,100./            
      DATA ZX/.00019,.00042,.00070,.00096,.0013,.0015,.0018,.0024,.003,.        
     X0042,.005,.007,.0088,.0092,.010,.011,.014,.015,.017,.018,.020,.021        
     X,.022,.024,.027,.030,.032,.035,.040,.044,.050,.060,.070,.090,.100,        
     X.15,.23,.50,1.8,4.0,7.0,15.0,8.0,5.0,3.0,1.7,1.2,.90,.50,.35,.20,.        
     X14,.10/           
      DATA ZW/13*0.0,.0001,.00017,.00034,.0021,.009,.025,.045,.10,.22,0.        
     W20,.16,.15,.11,.14,.10,.099,.098,.0963,.0967,.0981,.0987,.099,.100        
     W,.101,.103,.109,.118,.120,.122,.127,.130,.132,.138,.154,.161,.175,        
     W.20,.25,.34,.56/                
          
      TEN=10.           
      F=FK*.001         
      IF(F.LT..1) F=.1                
      IF(F.GT.100.) GO TO 20          
      DO 10 I=1,53      
      IF(F-FZ(I))12,11,10             
   10 CONTINUE          
      GO TO 20          
   12 IF(I.EQ.1) I=2    
   13 L=I-1             
      A=LOG10(F)       
      B=LOG10(FZ(I))                 
      C=LOG10(FZ(L))                 
      R=(A-C)/(B-C)     
      D=LOG10(ZX(I))                 
      E=LOG10(ZX(L))                 
      AR=(R*(D-E))+E    
      AO=TEN**AR        
      IF(I.LE.13) GO TO 21            
      G=LOG10(ZW(I))                 
      H=LOG10(ZW(L))                 
      WR=(R*(G-H))+H    
      AW=TEN**WR        
      RETURN            
   20 WRITE(IOT,19)     
      AO=.10            
      AW=.56            
      RETURN            
   11 AO=ZX(I)          
      AW=ZW(I)          
      RETURN            
   21 AW=0.0000         
      RETURN            
      END               
      SUBROUTINE COMB(YA,YB)          
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
      COMMON /VARY/ IDUM,IWA,IDUM2,IWB,DUM,YC(35)   
      COMMON/PRINT/KL,IN,IOT          
          
      DIMENSION P(35),YLEV(100),XVAL(100),YA(35),YB(35)           
      DIMENSION XX(35),PXVAL(100)     
          
      FNX(Y0,Y1,Y2,X1,X2)=((Y0-Y2)*(X1-X2)/(Y1-Y2))+X2            
          
      DATA  P/.001,.002,.005,.01,.02,.05,.1,.2,.5,1.,2.,5.,10.,   
     C15.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95.,98.,99.,99.5,99.8,            
     C99.9,99.95,99.98,99.99,99.995,99.998,99.999/  
          
2001  FORMAT (/2X,'THERE IS A FLAT SPOT IN ONE OF THE ARRAYS'/)   
2002  FORMAT (/2X,'Y LEVEL ARRAY IS NOT BIG ENOUGH FOR RANGE OF DATA'/)         
          
      NP=35             
      WA=FLOAT(IWA)     
      WB=FLOAT(IWB)     
      TTW=WA+WB         
      NPM1=NP-1         
      DO 21 I=1,NPM1    
      IF (YA(I+1) .LT. YA(I) .AND. YB(I+1) .LT. YB(I)) GO TO 21   
      WRITE(IOT,2001)                 
      RETURN            
21    CONTINUE          
      AMAX=MAX1(YA(1),YB(1))          
      AMIN=MIN1(YA(NP),YB(NP))        
      DIFF=AMAX-AMIN    
      IF (DIFF .GT. 10.) GO TO 23     
      YLEV(1)=AMAX      
      DO 22 I=2,10      
22    YLEV(I)=YLEV(I-1)-.1*DIFF       
      NL=10             
      GO TO 25          
23    CONTINUE          
      F=1.              
      IF (DIFF .GT. 99.) F=2.         
      YLEV(1)=AINT(AMAX)              
      DO 24 I=2,100     
      YLEV(I)=YLEV(I-1)-1.*F          
      IF (YLEV(I) .LT. AMIN) GO TO 25               
      NL=I              
24    CONTINUE          
      WRITE(IOT,2002)                 
      RETURN            
25    CONTINUE          
       NB=1             
        NA=1            
      DO 40 I=1,NL      
      IF (YA(NA) .LT. YLEV(I)) XVAL(I)=0.           
      IF (YLEV(I).LT. YA(NP)) XVAL(I)=WA*100.       
      IF (YA(NA) .LT. YLEV(I) .OR. YLEV(I) .LT. YA(NP)) GO TO 32  
      XVAL(I)=WA*(FNX(YLEV(I),YA(NA),YA(NA+1),P(NA),P(NA+1)))     
31    IF (YLEV(I+1) .GE. YA(NA+1) .OR. NL .EQ. I) GO TO 32        
      NA=NA+1           
      GO TO 31          
32    CONTINUE          
      IF (YB(NB) .LT. YLEV(I)) GO TO 34             
      IF (YLEV(I) .LT. YB(NP)) XVAL(I)=XVAL(I)+WB*100.            
      IF (YLEV(I) .LT. YB(NP)) GO TO 34             
      XVAL(I)=XVAL(I)+WB*FNX(YLEV(I),YB(NB),YB(NB+1),P(NB),P(NB+1))             
33    IF (YLEV(I+1) .GE. YB(NB+1) .OR. NL .EQ. I) GO TO 34        
      NB=NB+1           
      GO TO 33          
34    CONTINUE          
      XVAL(I)=XVAL(I)/TTW             
      PXVAL(I)=QERFI(XVAL(I)*.01)     
40    CONTINUE          
      DO 50 J=1,NP      
      XX(J)=QERFI(P(J)*.01)           
      DO 41 I=2,NL      
      IF (XX(J) .GT. PXVAL(I)) GO TO 42             
41    CONTINUE          
      I=NL              
42    YC(J)=FNX(XX(J),PXVAL(I-1),PXVAL(I),YLEV(I-1),YLEV(I))      
50    CONTINUE          
      WRITE(IOT,2010) (P(J),YA(J),YB(J),YC(J), J=1,35)            
2010  FORMAT(35(2X,'P=',F8.3,'  Y1=',F8.2,'  Y2=',F8.2,'  Y=',F8.2,/)           
     X///)              
      RETURN            
      END               
      SUBROUTINE CNLUT(A,B,C,N,R,RHO,P,D)           
        implicit real*8(A-H,O-Z)
          
C     ROUTINE FOR MODEL APR 77        
      DIMENSION A(50),B(50),C(50),P(50),D(50),X(100),Y(100)       
      DIMENSION Z(50)                 
          
      IF(A(N).LT.A(1)) GO TO 10       
      DO 11 I=1,N       
   11 X(I)=A(I)         
   12 IF(B(N).LT.B(1)) GO TO 13       
      IF(R.LT.0.) GO TO 14            
   15 DO 16 I=1,N       
   16 Y(I)=B(I)         
   17 DO 18 I=1,N       
      P(I)=C(I)         
      IF(C(I).GT..499.AND.C(I).LT..501) M=I         
   18 CONTINUE          
      Z(M)=X(M)+(R*Y(M))              
      DO 19 I=1,N       
      IF(I.EQ.M) GO TO 19             
      YA=X(I)-X(M)      
      YB=Y(I)-Y(M)      
      YU=SQRT ((YA*YA)+(YB*YB)+(2.*R*RHO*YA*YB))    
      IF(I.LT.M) GO TO 20             
      Z(I)=Z(M)+YU      
      GO TO 19          
   20 Z(I)=Z(M)-YU      
   19 CONTINUE          
      DO 23 I=1,N       
      K=N-I+1           
   23 D(I)=Z(K)         
      RETURN            
   10 DO 21 I=1,N       
      K=N-I+1           
   21 X(I)=A(K)         
      GO TO 12          
   13 IF(R.LT.0.) GO TO 15            
   14 DO 22 I=1,N       
      K=N-I+1           
   22 Y(I)=B(K)         
      GO TO 17          
      END               
      SUBROUTINE DFRAC                
        implicit real*8(A-H,O-Z)
          
C     SUBROUTINE TO COMPUTE DIFFRACTION ATTENUATION               
C     ROUTINE FOR MODEL APR 77        
          
    5 FORMAT(5X, 4F7.1,F8.4,2F8.3)    
    6 FORMAT(5X,10F7.1,F8.4,5F8.3,F7.1)             
    7 FORMAT(5X, 6F7.1,F8.4,5F8.3,F7.1)             
   51 FORMAT(8X,'DL7    DL8    TEC1    TEC2     TE4     AC3    D3     AC        
     X4    D4    AV4    GH7     ARK    AKS ')       
   52 FORMAT(5X,2F7.1,3F8.4,8F7.1)    
   57 FORMAT(8X,'AK3    AK4      D    DK4    GH1    GH2      W      AMD         
     X    AED    SWP    AWD    AK5    DK5')         
   60 FORMAT(8X,'AR3    AR4     D3     D4    AK3    AK4     D     DK4           
     X GH1    GH2      W     AMD     AED    SWP    AWD    AK5    DK5')          
   61 FORMAT(8X,'AR3    AR4     D3     D4      W     AMD     AED')              
   70 FORMAT(4(2X,E15.5))             
   71 FORMAT(10X,'W',14X,'D',14X,'DLS',12X,'DL')    
          
      COMMON/DIFPR/HT,HR,DH,AED,AMD,DLS1,DLS2,IPX,KSC,HLT,HRP,AWD,SWP,II        
      COMMON/PARAM/HTE,HRE,D,DL1,DL2,ENS,A,F,ALAK,TE1,TE2,KD,GAO,GAW            
      COMMON /SEA/ ISK,SCK,TP,JM,EPK,SGM            
      COMMON/PRINT/KP,IN,IOT          
      COMMON/CURVE/PI,RAD,DEG,TWPI,PITW             
          
      REAL*8 K1,K2,K3,K4                
          
      FNC(C)=416.4*(F**THIRD)*(1.607-C)             
      FND(C)=.36278/((C*F)**THIRD*(((E-1.)**2+(X*X))**.25))       
      FNE(C)=C*SQRT (E*E+X*X)         
          
      KP=II             
      S=SGM             
      E=EPK             
      IPOL=IPX-1        
      THIRD=1./3.       
      TWTRD=2./3.       
      H1E=1000.*HTE     
      H2E=HRE*1000.     
      HST=HT-HLT        
      HSR=HR-HLT        
      HL1=(HLT-HRP)     
      HL2=HR-HRP        
      HP1=HL1*1000.     
      S=SGM             
      E=EPK             
      DLS=DLS1+DLS2     
      DL=DL1+DL2        
      TE=TE1+TE2        
      CW=0.9            
      CU=.096790        
      TWA=2.*A          
      X=18000.*S/F      
      A1=DL1*DL1/(2.*HTE)             
      A2=DL2*DL2/(2.*HRE)             
      K1=FND(A1)        
      K2=FND(A2)        
      IF(IPOL.EQ.0) GO TO 3           
      K1=FNE(K1)        
      K2=FNE(K2)        
    3 CONTINUE          
          
C     CALCULATION OF GHBAR AND  W     
          
      B5=1.607-K1       
      B6=1.607-K2       
      GH1=GHBAR(F,A1,B5,K1,DL1,H1E)   
      GH2=GHBAR(F,A2,B6,K2,DL2,H2E)   
      AK3=6.-GH1-GH2    
      IF(D.GE.DLS) GO TO 41           
      IF(D.LE.(CW*DLS)) GO TO 50      
      W=0.5*(1.+COS ((PI*(DLS-D))/(DLS*(1.-CW))))   
C     -----------------------WRITE STATEMENTS-------------------- 
      IF(II.GT.0) GO TO 10            
      WRITE(IOT,71)     
      WRITE(IOT,70)W,D,DLS,DL         
      CALL PAGE(2)      
   10 CONTINUE          
C     ----------------------------------------------------------- 
      IF(W.LT..001) GO TO 45          
          
C     CALCULATION OF ROUNDED EARTH DIFFRACTION      
   42 CONTINUE          
      D3=DL+.5*(A*A/F)**THIRD         
      DL7=DL1           
      DL8=DL2           
      ASSIGN 25 TO JD                 
   30 D4=D3+(A*A/F)**THIRD            
      TX=((A*A/F)**THIRD)/A           
      T3=0.5*TX         
      T4=1.5*TX         
      A3=(D3-DL)/T3     
      A4=(D4-DL)/T4     
      K3=FND(A3)        
      K4=FND(A4)        
      IF(IPOL.EQ.0) GO TO 2           
      K3=FNE(K3)        
      K4=FNE(K4)        
2     CONTINUE          
      B1=FNC(K1)        
      B2=FNC(K2)        
      B3=FNC(K3)        
      B4=FNC(K4)        
      X1=B1*DL7/A1**TWTRD             
      X2=B2*DL8/A2**TWTRD             
      X3=X1+X2+(B3*(D3-DL)/(A3**TWTRD))             
      X4=X1+X2+(B4*(D4-DL)/(A4**TWTRD))             
      IF(K1.GE.1.) K1=.99999          
      IF(X1.GT.200.) GO TO 17         
      IF(K1.LE.00001) GO TO 16        
      XL1=450./ABS (LOG10(K1)**3)    
      IF(X1.GE.XL1) GO TO 16          
      FX1=20.*LOG10(K1)+(2.5*1.E-5*X1*X1/K1)-15.   
   20 IF(K2.GE.1.) K2=.99999          
      IF(X2.GT.200.) GO TO 19         
      IF(K2.LE.00001) GO TO 18        
      XL2=450./ABS (LOG10(K2)**3)    
      IF(X2.GE.XL2) GO TO 18          
      FX2=20.*LOG10(K2)+(2.5*1.E-5*X2*X2/K2)-15.   
   21 GX3=.05751*X3-10.*LOG10(X3)    
      GX4=.05751*X4-10.*LOG10(X4)    
      AC3=GX3-FX1-FX2-20.             
      AC4=GX4-FX1-FX2-20.             
      GO TO JD,(25,26)                
   17 FX1=.05751*X1-(10.*LOG10(X1))                
      IF(X1.GT.2000.) GO TO 20        
      W1=.0134*X1*EXP (-.005*X1)      
      FX1=W1*(40.*LOG10(X1)-117.)+(1.-W1)*FX1      
      GO TO 20          
16    T=40.*LOG10(X1)-117.           
      T1=-117.          
      T2=AMIN1((ABS (T)),(ABS (T1)))                
      FX1=T             
      IF(T2.EQ.ABS (T1)) FX1=T1       
      GO TO 20          
   19 FX2=.05751*X2-(10.*LOG10(X2))                
      IF(X2.GT.2000.) GO TO 21        
      W2=.0134*X2*EXP (-.005*X2)      
      FX2=W2*(40.*LOG10(X2)-117.)+(1.-W2)*FX2      
      GO TO 21          
18    T=40.*LOG10(X2)-117.           
      T1=-117.          
      T2=AMIN1((ABS (T)),(ABS (T1)))                
      FX2=T             
      IF(T2.EQ.ABS(T1)) FX2=T1        
      GO TO 21          
   25 AR3=AC3           
      AR4=AC4           
      DR4=D4            
      DR3=D3            
      AMS=(AR4-AR3)/(D4-D3)           
      AES=AR4-AMS*D4    
      IF(W.GT..999) GO TO 43          
          
C     CALCULATION OF SINGLE KNIFE EDGE WITH GHBAR   
          
   45 CONTINUE          
      IF(HL1.LE.0.) GO TO 43          
      TH1=ATAN ((HST/DL1)-(DL1/TWA))                
      TH=2.*ASIN(CU*SQRT(D/(F*DL1*DL2)))            
      TH5=-(-TH+TH1)    
      ATH5=A*TAN (TH5)                
      DLK5=-ATH5+SQRT (ATH5*ATH5+(HSR*TWA))         
      DK5=DLK5+DL1      
      TE5=ATAN ((-HSR/DLK5)-(DLK5/TWA))             
      TH5=TE1+TE5+(DK5/A)             
      TM5=SQRT ((F*DL1*DLK5)/DK5)     
      V5=5.1658*SIN(TH5*.5)*TM5       
      CALL FRNEL(V5,FV5,PH5)          
      AV5=-20.*LOG10(FV5)            
      AMK5=(AV5-AK3)/(DK5-D)          
      AWK=AK3-(AMK5*D)                
      DLST7=SQRT (HL1*TWA)            
      DLSR7=SQRT (HL2*TWA)            
      DL7=DLST7         
      DL8=DLSR7         
      DL=DL7+DL8        
      DLK4=DL           
      ASSIGN 26 TO JD                 
      A1=(DL7*DL7)/(2.*HL1)           
      A2=(DL8*DL8)/(2.*HL2)           
      K1=FND(A1)        
      K2=FND(A2)        
      IF(IPOL.EQ.0) GO TO 29          
      K1=FNE(K1)        
      K2=FNE(K2)        
   29 TEC1=ATAN ((-HL1/DL7)-(DL7/TWA))              
      TEC2=ATAN ((-HL2/DL8)-(DL8/TWA))              
      TE=TEC1+TEC2      
      D3=DL+.5*(A*A/F)**THIRD         
      GO TO 30          
   26 B7=1.607-K1       
      GH7=GHBAR(F,A1,B7,K1,DL7,HP1)   
      AC7=(AC4-AC3)/(D4-D3)           
      ARS=AC4-AC7*D4    
      ARK=ARS+AC7*DLK4                
      TE4=ATAN (((HLT-HR)/DLK4)-(DLK4/TWA))         
      DK4=DLK4+DL1      
      TH=TE1+TE4+(DK4/A)              
      TM2=SQRT ((F*DL1*DLK4)/DK4)     
      V4=5.1658*SIN(TH*.5)*TM2        
      CALL FRNEL(V4,FV,PH)            
      AV4=-20.*LOG10(FV)             
      AKS=AV4-GH1-GH7+ARK             
      AMKD=(AKS-AK3)/(DK4-D)          
      AEK=AK3-(AMKD*D)                
C     -----------------------WRITE STATEMENTS-------------------- 
      IF(II.GT.0) GO TO 11            
      WRITE(IOT,51)     
      WRITE(IOT,52)DL7,DL8,TEC1,TEC2,TE4,AC3,D3,AC4,D4,AV4,GH7,ARK,AKS          
      CALL PAGE(2)      
   11 CONTINUE          
C     ----------------------------------------------------------- 
      AK4=AEK+DK4*AMKD                
      WK=1.-W           
      AK5=AWK+DK5*AMK5                
      IF(W.LT..001) GO TO 36          
          
C     COMBINATION OF ROUNDED EARTH AND KNIFE EDGE DIFFRACTION     
          
      AT3=(WK*AK3)+(W*(AES+(AMS*D)))                
      AT4=(WK*AK4)+(W*(AES+(AMS*DK4)))              
      AT5=(WK*AK5)+(W*(AES+(AMS*DK5)))              
      AMD=(AT4-AT3)/(DK4-D)           
      AED=AT3-(AMD*D)                 
      SWP=(AT5-AT3)/(DK5-D)           
      AWD=AT3-(SWP*D)                 
C     -----------------------WRITE STATEMENTS-------------------- 
      IF(II.GT.0) RETURN              
      WRITE(IOT,60)     
      WRITE(IOT,6)AR3,AR4,DR3,DR4,AK3,AK4,D,DK4,GH1,GH2,W,AMD,AED,SWP,          
     XAWD,AK5,DK5       
      CALL PAGE(2)      
C     ----------------------------------------------------------- 
      RETURN            
          
   36 AED=AEK           
      AMD=AMKD          
      SWP=AMK5          
      AWD=AWK           
C     -----------------------WRITE STATEMENTS-------------------- 
      IF(II.GT.0) RETURN              
      WRITE(IOT,57)     
      WRITE(IOT,7)AK3,AK4,D,DK4,GH1,GH2,W,AMD,AED,SWP,AWD,AK5,DK5 
      CALL PAGE(2)      
C     ----------------------------------------------------------- 
      RETURN            
          
   41 W=1.              
      GO TO 42          
          
   43 AED=AES           
      AMD=AMS           
      AWD=AES           
      SWP=AMS           
C     -----------------------WRITE STATEMENTS-------------------- 
      IF(II.GT.0) RETURN              
      WRITE(IOT,61)     
      WRITE(IOT,5)AR3,AR4,DR3,DR4,W,AMD,AED         
      CALL PAGE(2)      
C     ----------------------------------------------------------- 
      RETURN            
          
   50 W=0.              
      GO TO 45          
      END               
      FUNCTION FDASP(S)               
        implicit real*8(A-H,O-Z)
          
C     ROUTINE FOR MODEL APR 77        
          
C     K IS BASED ON RATIO OF S TO .990              
C     THIS NAKAGAMA-RICE DIST.  HAS TABLES FROM NORTON 55 IRE PAGE 1360         
C     THE TABLES ARE THE NEGATIVE OF THE KK IRE TABLES BUT ARE CHANGED          
C     BEFORE GOING OUT OF THE ROUTINE               
C     K HAS THE OPPOSITE SIGN OF 101 BUT THE SAME AS THE IRE PAPER              
          
      COMMON/VV/VF(36,17)             
      AVEF(YN,XN,YN1,XN1,T)=(YN1*(T - XN) - YN*(T - XN1))/(XN1 - XN)            
          
      R=-S              
      DO 1 I=1,17       
      IF(R-VF(27,I)) 3,2,1            
    1 CONTINUE          
      I=17              
    2 AK=VF(1,I)        
      GO TO 6           
    3 IF(I.EQ.1) GO TO 2              
      AK=AVEF(VF(1,I-1),VF(27,I-1),VF(1,I),VF(27,I),R)            
    6 FDASP=AK          
      RETURN            
      END               
      SUBROUTINE FRNEL(V,FV,PH)       
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
C     SUBROUTINE TO CALCULATE THE FRESNEL INTEGRAL  
          
      DIMENSION A(11),B(11),G(11),D(11)             
      COMMON/CURVE/PI,RAD,DEG,TWPI,PITW             
      COMPLEX*16 PZ,SZ,CZ                
          
      DATA  A/-1.702E-6,-6.808568854,-5.76361E-4,6.920691902,-1.6898657E        
     X-2,-3.05048566,-7.5752419E-2,8.50663781E-1,-2.5639041E-2,-1.502309        
     X60E-1,3.4404779E-2/             
      DATA B/4.255387524,-9.281E-5,-7.7800204,-9.520895E-3,5.075161298,-        
     X1.38341947E-1,-1.363729124,-4.03349276E-1,7.02222016E-1,-2.1619592        
     X9E-1,1.9547031E-2/              
      DATA  G/-2.4933975E-2,3.936E-6,5.770956E-3,6.89892E-4,-9.497136E-3        
     X,1.1948809E-2,-6.748873E-3,2.4642E-4,2.102967E-3,-1.21793E-3,2.339        
     X39E-4/            
      DATA  D/2.3E-8,-9.351341E-3,2.3006E-5,4.851466E-3,1.903218E-3,-1.7        
     X122914E-2,2.9064067E-2,-2.7928955E-2,1.6497308E-2,-5.598515E-3,8.3        
     X8386E-4/          
          
      IF(V.EQ.0.) GO TO 71            
      IF(V.GE.5.) GO TO 74            
      PT=V*V*.25        
      CPSI=TWPI*(PT-INT (PT))         
      X=V*V*PI*.5       
      IF(X.GT.4.) GO TO 10            
      PX=COS(X)*SQRT(X/4.)            
      PY= SIN (-X)*SQRT (X/4.)        
      SUMX=1.59576914                 
      SUMY=-3.3E-8      
      XN= 1.            
      DO 100 I = 1, 11                
      XN=XN*X/4.        
      SUMX=SUMX+A(I)*XN               
  100 SUMY=SUMY+B(I)*XN               
      SZ=CMPLX(SUMX,SUMY)             
      PZ=CMPLX(PX,PY)                 
      CZ=SZ*PZ          
      C=REAL(CZ)        
      S=AIMAG(CZ)       
      GO TO 30          
   10 PX=COS (X)*SQRT (4./X)          
      PY=SIN (-X)*SQRT (4./X)         
      XN=1.             
      SUMX=0.           
      SUMY=.199471140                 
      DO 200 I = 1, 11                
      XN=XN*4./X        
      SUMX=SUMX+G(I)*XN               
  200 SUMY=SUMY+D(I)*XN               
      SZ=CMPLX(SUMX,SUMY)             
      PZ=CMPLX(PX,PY)                 
      CZ=SZ*PZ          
      C=REAL(CZ)        
      S=AIMAG(CZ)       
      C=C+.5            
      S=S-.5            
   30 S=ABS (S)         
      IF(V.LT.0.) GO TO 70            
      FV=.5*SQRT ((1.-(C+S))**2+(C-S)**2)           
      Y=C-S             
      W=1.-(C+S)        
   75 PH=ATAN2(Y,W)     
      PH=PH-CPSI        
      AP=ABS (PH)       
      APX=AP-TWPI*INT (AP/TWPI)       
      IF(PH.LT.0.) GO TO 37           
      PH=APX            
   39 IF(PH.LT.0.) PH=TWPI+PH         
      RETURN            
   37 PH=-APX           
      GO TO 39          
   71 FV=.5             
      PH=0.             
      GO TO 39          
   74 FV=.22508/V       
      PH=.78539816      
      GO TO 39          
   70 FV=.5*SQRT ((1.+(C+S))**2+(C-S)**2)           
      Y=-(C-S)          
      W=1.+(C+S)        
      GO TO 75          
      END               
      SUBROUTINE GANE(X,ITA,HH,TT,GAV,GAH,GDR,NP)   
        implicit real*8(A-H,O-Z)
****  850606 CORRECTION ADDED         
C     ROUTINE COPY A    
C     ROUTINE FOR MODEL MAY 78        
C           HHPBW IS HALF OF THE HALF-POWER BEAMWIDTH IN DEGREES        JUNE 73 
  551 FORMAT (1H+,38X,A10)            
  552 FORMAT (1H+,38X,A4)             
****  850606            
* 553 FORMAT (1H+,38X,A11)            
  553 FORMAT (1H+,38X,A14)            
  554 FORMAT (1H+,38X,A39)            
  555 FORMAT (1H+,38X,A39)            
  556 FORMAT (1H+,38X,A34)            
  557 FORMAT (1H+,38X,A10)            
          
      DIMENSION DA(8),DG(8)           
      DIMENSION RA(24),RB(24)         
      COMMON/PRINT/KL,IN,IOT          
          
      CHARACTER *10 PA1,PA7,PA2*4,PA3*14,PA4*39,PA5*39,PA6*34     
          
      DATA PA1/' ISOTROPIC'/          
      DATA PA4/' 4-LOOP ARRAY (COSINE VERTICAL PATTERN)'/         
      DATA PA5/' 8-LOOP ARRAY (COSINE VERTICAL PATTERN)'/         
      DATA PA6/' I OR II (COSINE VERTICAL PATTERN'/               
      DATA PA7/'NARROWBEAM'/          
          
          
      DATA PA2/' DME'/                
      DATA  DA/-6.,0.,2.5,5.,7.5,7.51,14.99,15.0/   
          
      DATA  DG/-8.,-6.,-3.,0.,-3.0,-20.,-20.01,-30./              
          
          
      DATA PA3/' TACAN (RTA-2)'/      
      DATA RA/-90.,-76.,-60.,-54.,-51.5,-48.,-36.,-33.,-30.,-24.,-18.,-1        
     X2.,-9.,-6.,-2.5,0.,3.,8.,12.,24.,36.,57.,84.,90./           
          
      DATA RB/-29.,-22.,-26.5,-27.4,-21.7,-20.,-5.5,-4.2,-3.5,-4.5,-7.3,        
     X-11.8,-10.,-3.5,-1.,4.,6.5,7.4,7.,-1.4,-1.5,-9.5,-4.,-13.0/ 
          
      FNA(FX,FA,FB,FC,FD)=((FX-FB)*(FC-FD)/(FA-FB))+FD            
          
      A=X               
      IFA=ITA           
      HPBW=HH           
      TLT=TT            
      IP=NP             
      GO TO (10,20,30,40,50,60,70),IFA              
          
          
C     --------- GAIN FOR ISOTROPIC ANTENNA ------------------     
          
   10 GAIN=1.           
      GV=1.             
      GH=1.             
      GO TO 1           
C     --------- GAIN FOR DME ANTENNA ------------------           
   20 D=A*57.29577951                 
      IF(D.LE.(-6.))  GO TO 24        
      IF(D.GE.15.00) GO TO 25         
      DO 21 I=1,8       
      IF(D-DA(I))23,22,21             
   21 CONTINUE          
   25 I=8               
   22 GAIN=10.**(DG(I)*.05)           
      GO TO 9           
   23 IF(I.LE.1) GO TO 24             
      L=I-1             
      GD=FNA(D,DA(I),DA(L),DG(I),DG(L))             
      GAIN=10.**(GD*.05)              
      GO TO 9           
   24 I=1               
      GO TO 22          
C     --------- GAIN FOR RTA-2 ANTENNA ---------------            
   30 CONTINUE          
      D=(A*57.29577951)               
      IF(D.LE.RA(1)) GO TO 35         
      IF(D.GE.RA(24)) GO TO 36        
      DO 31 I=1,24      
      IF(D-RA(I))33,32,31             
   31 CONTINUE          
   36 I=24              
   32 GAIN=10.**((RB(I)-7.4)*.05)     
      GO TO 9           
   33 IF(I.LE.1) GO TO 35             
      L=I-1             
      RD=FNA(D,RA(I),RA(L),RB(I),RB(L))             
      GAIN=10.**((RD-7.4)*.05)        
      GO TO 9           
   35 I=1               
      GO TO 32          
C     --------- GAIN FOR VOR ANTENNA (COSINE PATTERN) ----------- 
   40 GAIN=1.00*COS (A)               
      IF(GAIN.LT..12589) GAIN=.12589                
      GO TO 9           
C     --------------- GAIN FOR ILS LOCALIZER -------------------  
   50 GAIN=1.00*COS (A)               
      IF(GAIN.LT..12589) GAIN=.12589                
      GO TO 9           
C     --------------- GAIN FOR GLIDE SLOPE ---------------------  
   60 GAIN=1.00*COS (A)               
      IF(GAIN.LT..12589) GAIN=.12589                
      GO TO 9           
C     ----------------- JTAC -------------------------------------------        
   70 D=A*57.29577951                 
      TERM=ABS (D-TLT)                
      IF(HPBW.LE.0.)GO TO 10          
      GAIN=(1.+((TERM/HPBW)**2.5))**(-0.5)          
      GO TO 9           
          
    9 IF(IP-2) 6,7,8    
    6 GH=GAIN           
      GV=0.             
      GO TO 1           
    7 GV=GAIN           
      GH=0.             
      GO TO 1           
    8 GH=GAIN           
      GV=GAIN           
      GO TO 1           
    1 GAV=GV            
      GAH=GH            
      GDR=GAIN          
      RETURN            
          
      ENTRY ANTNA       
      A=X               
      IFA=ITA           
      HPBW=HH           
      TLT=TT            
      IP=NP             
      GO TO (501,502,503,504,505,506,507),IFA       
  550 GAIN=0.           
      RETURN            
  501 WRITE(IOT,551)PA1               
      GO TO 550         
  502 WRITE(IOT,552)PA2               
      GO TO 550         
  503 WRITE(IOT,553)PA3               
      GO TO 550         
  504 WRITE(IOT,554)PA4               
      GO TO 550         
  505 WRITE(IOT,555)PA5               
      GO TO 550         
  506 WRITE(IOT,556)PA6               
      GO TO 550         
  507 WRITE(IOT,557)PA7               
      GO TO 550         
        END             
      FUNCTION GHBAR (F,A,B,AK,DHOR,HE)             
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
    6 FORMAT(5X,'GHBAR=',F6.1,'HB=',F10.4,' K=',F7.4,' B=',F7.4)  
      COMMON/PRINT/KP,IN,IOT          
      WG=2.             
      PIG=3.141592654                 
      HB=2.2325*B*B*(F*F/A)**.33333333*(.001*HE)    
      IF(HB.GE.2.5) GO TO 10          
      IF(AK.GT..05) GO TO 11          
      IF(HB.LT..3) GO TO 12           
      GHBAR=-6.5-1.67*HB+6.8*LOG10(HB)             
   13 IF(AK.LE..01) GO TO 2           
      GHB=GHBAR         
   11 IF(HB.LT..25) GO TO 14          
      GHT=-5.9-1.9*HB+6.6*LOG10(HB)                
   15 IF(AK.GT..05) GO TO 16          
      GHBAR=GHT-(GHT-GHB)*((.05-AK)/.04)            
2     CONTINUE          
      FRE=300.*SQRT(.2997925*DHOR/F)                
      IF(HE.LE.FRE) GO TO 3           
      IF(HE.GE.(WG*FRE)) GO TO 4      
      GW=.5*(1.+COS (PIG*(HE-FRE)/FRE))             
      GHBAR=GHBAR*GW    
      GO TO 3           
    4 GHBAR=0.          
    3 IF(KP.GT.0) GO TO 5             
      IF(AK.GT..1) WRITE(IOT,6)GHBAR,HB,AK,B        
      IF(AK.LT..01) WRITE(IOT,6)GHBAR,HB,AK,B       
      IF(HB.LT..01) WRITE(IOT,6)GHBAR,HB,AK,B       
      IF(HB.GT.100.) WRITE(IOT,6)GHBAR,HB,AK,B      
    5 RETURN            
   10 GHBAR=-6.6-.013*HB-2.*LOG10(HB)              
      GO TO 2           
   12 GHBAR=1.2-13.5*HB+15.*LOG10(HB)              
      GO TO 13          
   14 GHT=-13.9+24.1*HB+3.1*LOG10(HB)              
      GO TO 15          
   16 GHB=GHT           
      IF(HB.LT.0.1) GO TO 17          
      GHT=-4.7-2.5*HB+7.6*LOG10(HB)                
   18 GHBAR=GHT-(GHT-GHB)*((.1-AK)/.05)             
      GO TO 2           
   17 GHT=-13.          
      GO TO 18          
      END               
      SUBROUTINE PAGE(N)              
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
    4 FORMAT(1H1)       
    6 FORMAT('  PAGE',I4)             
      SAVE              
          
      COMMON/EGAP/IP,LN               
      COMMON/PRINT/KL,IN,IOT          
          
      IF(N)10,11,12     
   10 IP=0              
   11 IP=IP+1           
      LN=1              
      WRITE(IOT,6)IP    
   13 RETURN            
   12 LN=LN+N           
      IF(LN-53)13,14,14               
   14 WRITE(IOT,4)      
      GO TO 11          
      END               
      FUNCTION QERFI(Q)               
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
C        THE INVERSE OF QERF, GIVES THE STANDARD NORMAL DEVIATE AS A            
C           FUNCTION OF THE COMPLEMENTARY PROBABILITY             
C        TRUNCATED AT 0.00000001 AND 0.99999999     
C        APPROXIMATION DUE TO C. HASTINGS, JR.      
C        MAX ERROR 4.5E-4             
          
      DATA C0,C1,C2/2.515516698,0.802853,0.010328/  
      DATA D1,D2,D3/1.432788,0.189269,0.001308/     
          
      T=Q               
      IF(T .GT. 0.5)     T=1.-T       
      T=AMAX1(T,0.00000001)           
      T=SQRT(-LOG(T)*2.)             
      QERFI=T-((C2*T+C1)*T+C0)/(((D3*T+D2)*T+D1)*T+1.)            
      IF(Q .GT. 0.5)     QERFI=-QERFI               
      RETURN            
      END               
      SUBROUTINE QVARB(KL)             
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
C        COMPUTES THE PARAMETERS WHICH DEFINE THE CURVES OF VARIABILITY,        
C        CALL QVAR INITIALIZES CONSTANTS FOR AN ENTIRE RADIAL     
C        CALL QQVAR(S) COMPUTES THE CONSTANTS OF /PROPV/ WHICH    
C           CORRESPOND TO THE DISTANCE S            
          
      SAVE              
      COMMON/VARY/KLM,MX1,KLM2,MX2,FKE,AD(35)       
          
*      EQUIVALENCE(KLIM,CLIM)          
          
      DIMENSION P(35),ZP(35)          
      DIMENSION XSP1(7),XSP2(7),XSP3(7),XSM1(7),XSM2(7),XSM3(7)   
      DIMENSION XV1(7),XV2(7),XV3(7),BSP1(7),BSP2(7),BSM1(7),BSM2(7)            
      DIMENSION BV1(7),BV2(7),BSD1(7),BZD1(7)       
          
      DATA BV1/-9.67,-0.62,1.26,-9.21,-0.62,-0.39,3.15/           
      DATA BV2/12.7,9.19,15.5,9.05,9.19,2.86,857.9/               
      DATA XV1/144.9,228.9,262.6,84.1,228.9,141.7,2222./          
      DATA XV2/190.3,205.2,185.2,101.1,205.2,315.9,164.8/         
      DATA XV3/133.8,143.6,99.8,98.6,143.6,167.4,116.3/           
      DATA BSM1/2.13,2.66,6.11,1.98,2.68,6.86,8.51/               
      DATA BSM2/159.5,7.67,6.65,13.11,7.16,10.38,169.8/           
      DATA XSM1/762.2,100.4,138.2,139.1,93.7,187.8,609.8/         
      DATA XSM2/123.6,172.5,242.2,132.7,186.8,169.6,119.9/        
      DATA XSM3/94.5,136.4,178.6,193.5,133.5,108.9,106.6/         
      DATA BSP1/2.11,6.87,10.08,3.68,4.75,8.58,8.43/              
      DATA BSP2/102.3,15.53,9.60,159.3,8.12,13.97,8.19/           
      DATA XSP1/636.9,138.7,165.3,464.4,93.2,216.0,136.2/         
      DATA XSP2/134.8,143.7,225.7,93.1,135.9,152.0,188.5/         
      DATA XSP3/95.6,98.6,129.7,94.2,113.4,122.7,122.9/           
      DATA BSD1/1.224,.801,1.420,1.000,1.224,1.473,1.473/         
      DATA BZD1/1.282,2.161,1.282,0.,1.282,1.282,1.282/           
      DATA P/.00001,.00002,.00005,.0001,.0002,.0005,.001,.002,.005,.01,         
     X.02,.05,.10,.15,.20,.30,.40,.50,.60,.70,.80,.85,.90,.95,.98,.99,          
     X.995,.998,.999,.9995,.9998,.9999,.99995,.99998,.99999/      
          
      AFIT(C1,C2,X1,X2,X3)=((DE/X1)**2)/(1.+((DE/X1)**2))*        
     1     (C1+C2/(1.+((DE-X2)/X3)**2))             
      QF(F)=SIN(LOG10(F/200.)*5.)    
          
*      CLIM=S
      KLIM=KL
      IF(KLIM .LE. 0 .OR. KLIM .GT. 7) KLIM=5       
      CV1=BV1(KLIM)     
      CV2=BV2(KLIM)     
      YV1=XV1(KLIM)     
      YV2=XV2(KLIM)     
      YV3=XV3(KLIM)     
      CSM1=BSM1(KLIM)                 
      CSM2=BSM2(KLIM)                 
      YSM1=XSM1(KLIM)                 
      YSM2=XSM2(KLIM)                 
      YSM3=XSM3(KLIM)                 
      CSP1=BSP1(KLIM)                 
      CSP2=BSP2(KLIM)                 
      YSP1=XSP1(KLIM)                 
      YSP2=XSP2(KLIM)                 
      YSP3=XSP3(KLIM)                 
      CSD1=BSD1(KLIM)                 
      ZD=BZD1(KLIM)     
      GP=1.             
      GM=1.             
      IF(FKE.LT.60.) GO TO 10         
      IF(KLIM.EQ.2) GO TO 20          
      IF(KLIM.EQ.4) GO TO 30          
      IF(KLIM.EQ.5) GO TO 20          
      GO TO 10          
   20 IF(FKE.GT.1500.) GO TO 21       
      GP=0.18*QF(FKE) +1.06           
      IF(KLIM.EQ.5) GO TO 22          
      GO TO 10          
   30 IF(FKE.GT.1500.) GO TO 21       
      IF(FKE.LT.200.) GO TO 10        
      GP=0.10*QF(FKE) +1.02           
      GO TO 10          
   21 GP=0.93           
      IF(KLIM.EQ.5) GO TO 23          
      GO TO 10          
   22 GM=0.13*QF(FKE)+1.04            
      GO TO 10          
   23 GM=0.92           
   10 CONTINUE          
      DO 15 I=1,35      
   15 ZP(I)=QERFI(P(I))               
      RETURN            
          
      ENTRY QQVARB(S)                 
          
      DE=S              
      VMD=AFIT(CV1,CV2,YV1,YV2,YV3)   
      SGTM=AFIT(CSM1,CSM2,YSM1,YSM2,YSM3)*GM        
      SGTP=AFIT(CSP1,CSP2,YSP1,YSP2,YSP3)*GP        
      SGTD=SGTP*CSD1    
      YD=SGTP*ZD        
      DO 11 I=1,35      
      YT=VMD            
      Z=ZP(I)           
      IF(Z) 1,11,2      
    1 YT=SGTM*Z+YT      
      GO TO 11          
    2 Z1=Z-ZD           
      IF(Z1) 3,3,4      
    3 YT=SGTP*Z+YT      
      GO TO 11          
    4 YT=SGTD*Z1+YD+YT                
   11 AD(I)=YT          
      RETURN            
      END               
      SUBROUTINE RAIN                 
        implicit real*8(A-H,O-Z)
****  850402 CORRECTIONS INCLUDED     
C     ROUTINE FOR MODEL SEPT 78       
      DIMENSION RR(6,10),RP(24),PP(3,24),FX(12),QRR(11),QRX(12,11)              
      COMMON/DROPS/RNE(35),IZ,STS,FQZ,DT,HA1,HA2,EFA,AN           
      COMMON/CURVE/PI,RAD,DEG,TWPI,PI2              
          
      DATA  FX/2.,3.,4.3,5.,6.,7.5,10.,15.,20.,30.,60.,100./      
      DATA PP/3*1.0,.96,.916,.835,.94,.865,.73,.93,.85,.70,.915,.83,.66,        
     X00.899,.797,.61,.888,.775,.58,.878,.755,.564,.865,.730,.54,.86,.72        
     X,.53,.85,.704,.51,.84,.683,.493,.833,.670,.480,.824,.65,.47,.818,.        
     X64,.452,.813,.627,.44,.805,.61,.422,.798,.592,.408,.79,.575,.392,.        
     X78,.565,.378,.77,.55,.357,.765,.54,.35,.76,.528,.325,.758,.52,.31/        
      DATA  QRR/0.0,0.25,1.25,2.5,5.,12.5,25.,50.,100.,150.,200./ 
      DATA QRX/16*0.0,.0004,.0008,.002,.006,.01,.03,.16,.25,0.,.0004,0.0        
     X008,.001,.002,.004,.01,.04,.09,.21,.76,1.29,.0003,.0007,.002,.003,        
     X.005,.01,.03,.10,.20,.45,1.43,2.19,.0006,.001,.003,.006,.01,.03,.0        
     X7,.23,.43,.93,2.63,3.68,.001,.003,.009,.02,.03,.08,.24,.71,1.18,2.        
     X43,5.46,7.08,.002,.007,.02,.04,.09,.22,.60,1.53,2.49,4.87,9.86,11.        
     X7,.005,.01,.05,.10,.24,.59,1.45,3.28,5.15,9.59,17.,19.6,.01,.03,.1        
     X2,.26,.64,1.55,3.43,6.77,10.4,18.4,29.4,33.7,.02,.05,.21,.47,1.13,        
     X2.71,5.49,10.2,15.7,26.8,40.9,46.8,.04,.08,.34,.73,1.80,4.10,8.10,        
     X14.5,22.0,34.,56.,61./          
      DATA RP/10.,15.,18.,20.,23.,27.,30.,33.,37.,40.,44.,50.,55.,60.,6         
     X5.,70.,80.,90.,100.,115.,140.,155.,185.,200./               
      DATA RR/0.17,.25,.31,.54,.75,1.,.62,.98,1.54,2.07,2.7,3.35, 
     X   1.8,3.1,4.8,6.2,7.8,9.4,3.2,5.4,8.8,11.7,14.,17.,5.1,9.6,14.5,1        
     X9.,23.5,28.5,8.2,17.,25.,33.,40.,48.,11.3,22.8,34.,44.5,54.,67.,14        
     X.6,29.5,43.,57.,68.,84.,18.8,37.8,56.,73.,91.,112.,24.,44.,64.,86.        
     X,110.,160./       
          
      FRA(XX,RA,RB)=(XX-RB)/(RA-RB)   
      FRB(XR,RH,RL)=(XR*(RH-RL))+RL   
****  850402            
*     FNA(FX,FA,FB,FC,FD)=((FX-FB)*(FC-FD)/(FA-FB))+FD            
      FNA(FY,FA,FB,FC,FD)=((FY-FB)*(FC-FD)/(FA-FB))+FD            
          
      RNP=.5            
      F=FQZ*.001        
      DX=DT             
      H1=HA1            
      H2=HA2            
      A=EFA             
      BA=AN             
      IF(STS.GT.20.) STS=20.          
      TE=STS            
      IF(H1.GT.H2) GO TO 80           
      HS=H1             
      HL=H2             
   86 AT=PI2+BA         
      ANUM=HS*SIN (AT)                
      H=TE+A            
      IF(HL.LE.H) GO TO 83            
      IF(H.LT.HS) GO TO 81            
      AS=ASIN (ANUM/H)                
      AE=PI-(AT+AS)     
      IF(BA.GT.1.5620) GO TO 84       
      IF(AE.EQ.0.) GO TO 84           
      DX=(HS*SIN(AE))/SIN(AS)         
      GO TO 82          
   84 DX=H-HS           
   82 IF(DX.GE.TE) DX=TE              
      IF(IZ.GT.6) GO TO 18            
      IF(F.LT.2.) GO TO 23            
      IF(F.GT.100.) GO TO 28          
      DO 24 KI=1,12     
      IF(F.LE.FX(KI)) GO TO 25        
   24 CONTINUE          
      KI=12             
   25 IF(KI.EQ.1) KI=2                
      KL=KI-1           
      QFX=FRA(F,FX(KI),FX(KL))        
   29 IF(TE.GT.5.01) GO TO 20         
      K=1               
   22 DO 10 I=1,10      
      RX=RR(IZ,I)       
      IF(RX.LE.10.) GO TO 11          
      IF(RX.GE.200.) GO TO 12         
      DO 13 N=1,24      
      IF(RX.LE.RP(N)) GO TO 14        
   13 CONTINUE          
      N=24              
   14 IF(N.EQ.1) N=2    
      L=N-1             
      PD=FNA(RX,RP(N),RP(L),PP(K,N),PP(K,L))        
      GO TO 17          
   12 N=24              
      PD=PP(K,N)        
      GO TO 17          
   11 PD=1.             
   17 CONTINUE          
      PD=PD*RX          
      DO 15 N=1,11      
      IF(PD.LT.QRR(N)) GO TO 16       
   15 CONTINUE          
      N=11              
   16 L=N-1             
      RAT=FRA(PD,QRR(N),QRR(L))       
      IF(F.GE.100.) GO TO 41          
      AH=FRB(RAT,QRX(KI,N),QRX(KI,L))               
      AL=FRB(RAT,QRX(KL,N),QRX(KL,L))               
      AX=FRB(QFX,AH,AL)               
      GO TO 42          
   41 AX=FRB(RAT,QRX(KI,N),QRX(KI,L))               
   42 J=25+I            
      RNE(J)=AX*DX      
   10 CONTINUE          
      DO 27 I=1,25      
   27 RNE(I)=0.0        
      RETURN            
   20 IF(TE.GT.10.01) GO TO 21        
      K=2               
      GO TO 22          
   21 K=3               
      GO TO 22          
   23 DO 26 I=1,35      
   26 RNE(I)=0.0        
      RETURN            
   28 QFX=1.            
      KI=12             
      GO TO 29          
   80 HS=H2             
      HL=H1             
      BA=-(AN+(DT/A))                 
      GO TO 86          
   81 IF(AT.GT.PI2) GO TO 85          
      HC=HS*SIN (AT)    
      IF(H.LE.HC) GO TO 85            
      DX=2.*H*SIN(ACOS(HC/H))         
      GO TO 82          
   83 DX=DT             
      GO TO 82          
   85 DX=0.             
      GO TO 23          
   18 RNE(18)=DX*RNP    
      RETURN            
      END               
