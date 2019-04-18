      SUBROUTINE RDEMS(ARG,ADG,IMI,SEC)             
        implicit real*8(A-H,O-Z)
****  850402 CHANGES INCLUDED         
C     ROUTINE FOR MODEL APR 77        
C     SUBROUTINE TO CHANGE RADIANS TO DEGREES, MINUTES AND SECONDS              
          
****  850402 CBUF ADDED               
      CHARACTER CBUF*3                
      COMMON/CURVE/PI,RAD,DEG,TWPI,PITW             
      DE=ABS(ARG)*DEG                 
      IDE=INT(DE)       
      AMINT=60.*(DE-FLOAT (IDE))      
      IMI=INT (AMINT)                 
      SEC=(AMINT-FLOAT (IMI))*60.     
      IF(SEC.GT.59.99995) GO TO 9     
    7 IF(IMI.GT.59) GO TO 8           
    6 IF (IDE .EQ. 0) GO TO 10        
      IF (ARG .LT. 0) IDE=-IDE        
****  850402            
*  12 ENCODE (3,11,ADG) IDE           
   12 WRITE(CBUF(1:3),11) IDE         
      READ (CBUF(1:3),'(A3)') ADG     
   11 FORMAT (I3)       
      RETURN            
   10 IF (ARG .GT. 0.) GO TO 12       
****  850402            
*     ENCODE (3,13,ADG) IDE           
      WRITE(CBUF(1:3),13) IDE         
      READ (CBUF(1:3),'(A3)') ADG     
   13 FORMAT (' -',I1)                
      RETURN            
    9 SEC=0.            
      IMI=IMI+1         
      GO TO 7           
    8 IDE=IDE+1         
      IMI=0             
      GO TO 6           
      END               
      SUBROUTINE RECPX(XI,FK,IR,NP,MS,DK,GH,GV,R,PIC,RLM)         
        implicit real*8(A-H,O-Z)
          
C       ---NOTE---   THIS ANGLE IS LIKE THE FORMULATION IN TN 101 AND IS        
C     PI-C              
C     ROUTINE FOR MODEL APR 77        
C             THIS INCLUDES THE CIRCULAR POLARIZATION             
          
      DIMENSION SSH(10)               
          
      COMMON /SEA/ ISK,SCK,TP,JM,EPK,SGM            
      COMMON/CURVE/PI,RAD,DEG,TWPI,PI2              
          
      DATA SSH/0.,0.02,0.11,0.25,0.46,0.76,1.2,2.,3.,3.3/         
          
      SI=XI             
      DH=DK             
      IC=0              
      TWLD=2.095841232E-2*FK*(-1.)    
      MP=NP             
      IF(SI.LE.0.) GO TO 301          
      IF(SI.GE.PI2) GO TO 300         
      SISI=SIN (SI)     
      COSI=COS (SI)     
      IF(SISI.LE.0.) GO TO 15         
      SQSI=SQRT (SISI)                
   16 IF (IR .EQ. 1 .OR. IR .EQ. 5) GO TO 28        
      IF(MS.GT.0) GO TO 19            
      IF(JM.GT.0) GO TO 41            
      IF(DH.LE.4.) GO TO 17           
      SH=.78*DH*EXP (-.5*(DH**.25))   
   18 EXDH=EXP (TWLD*SH*SISI)         
      DX=(SH*SISI*FK/299.7925)        
      IF(DX.GT.0.3) GO TO 32          
      IF(DX.GE.0.1237) GO TO 33       
      IF(DX.GT.0.0739) GO TO 34       
      IF(DX.GE.0.00325) GO TO 35      
      PD=946.*DX*DX+0.01              
   36 CONTINUE          
   25 IF(MP-2) 10,11,20               
   10 ASSIGN 12 TO N    
      GO TO 6           
   11 CONTINUE          
      ASSIGN 13 TO N    
    6 X=(18000.*SGM)/FK               
      TRM=EPK  -(COSI*COSI)           
      TUPS=SQRT ((TRM*TRM)+(X*X))+TRM               
      P=SQRT (TUPS*.5)                
      GO TO N,(12,13)                 
   12 Q=X/(2.*P)        
      DENOM=(P*P)+(Q*Q)               
      B=1./DENOM        
      AM=(2.*P)/DENOM                 
      RS=(1.+(B*SISI*SISI)-(AM*SISI))/(1.+(B*SISI*SISI)+(AM*SISI))              
      R=SQRT (RS)*GH    
      RLM=R*PD          
      IF(MS.LT.1) R=R*EXDH            
      TOP=-Q            
      BOT=SISI-P        
      TRA=ATAN2(TOP,BOT)              
      TOP =Q            
      BOT=SISI+P        
      GO TO 14          
   13 Q=X/(2.*P)        
      DENOM=(P*P)+(Q*Q)               
      B=((EPK*EPK)+(X*X))/DENOM       
      AM=(2.*((P*EPK  )+(Q*X)))/DENOM               
      RS=(1.+(B*SISI*SISI)-(AM*SISI))/(1.+(B*SISI*SISI)+(AM*SISI))              
      R=SQRT (RS)*GV    
      RLM=R*PD          
      IF(MS.LT.1) R=R*EXDH            
      TOP=(X*SISI)-Q    
      BOT=(EPK  *SISI)-P              
      TRA=ATAN2(TOP,BOT)              
      TOP=(X*SISI)+Q    
      BOT=(EPK  *SISI)+P              
   14 TRB=ATAN2(TOP,BOT)              
      PIC=TRA-TRB       
      IF(IC-1) 51,22,23               
   15 SQSI=0.           
      GO TO 16          
   17 SH=.39*DH         
      GO TO 18          
   19 SH=0.             
      EXDH=1.           
      PD=.01            
      GO TO 25          
   20 IC=1              
      MP=MP-1           
      GO TO 11          
   21 TER= (RLV*RLV)+(RLH*RLH)+(2.*RLV*RLH*COS (PH-PV))           
      RLM=SQRT (TER)/2.               
   51 IF(RLM.LT.0.01) RLM=0.01        
      RETURN            
   22 IC=2              
      RV=R              
      PV=PIC            
      MP=MP-1           
      RLV=RLM           
      GO TO 10          
   23 IC=0              
      RH=R              
      PH=PIC            
      RLH=RLM           
      TER=   ((RV*RV)+(RH*RH)+(2.*RV*RH*COS (PH-PV)))             
      IF(TER.LE.0.) GO TO 30          
      R=SQRT (TER)/2.                 
   31 TOP=(RH*SIN (PH))+(RV*SIN (PV))               
      BOT=(RH*COS (PH))+(RV*COS (PV))               
      IF(BOT.EQ.0.) GO TO 24          
       PC=ATAN2(TOP,BOT)              
      PIC=PC            
      GO TO 21          
   24 PIC=0.0           
      GO TO 21          
   30 R=0.0             
      GO TO 31          
   32 PD=(0.875*EXP(-3.88*DX)+0.01)   
      GO TO 36          
   33 PD=(-1.06*DX)+0.601             
      GO TO 36          
   34 PD=0.45+SQRT (.000843-(DX-.1026)**2)          
      GO TO 36          
   35 PD=6.15*DX        
      GO TO 36          
   28 IF(JM.GT.0) GO TO 41            
      IF(ISK.LT.0) ISK=0              
      IF(ISK.GT.9) ISK=9              
      SH=SSH(ISK+1)     
      GO TO 18          
   41 SH=SCK            
      GO TO  18         
  300 SI=PI2            
      SISI=1.           
      COSI=0.           
      SQSI=1.           
      GO TO 16          
  301 SI=0.             
      SISI=0.           
      COSI=1.           
      SQSI=0.           
      GO TO 16          
      END               
      FUNCTION RYTRA(TT)              
        implicit real*8(A-H,O-Z)
          
C     ROUTINE FOR MODEL APR 77        
      SAVE              
      COMMON/RYTC/ENS,HC,HA,HS,D      
      DIMENSION A(25),RI(25),EN(25),H(25),TEI(25),R(25)           
          
      DATA H/0.00,.01,.02,.05,.1,.2,.305,.5,.7,1.,1.524,2.,3.048,5.,7.,1        
     X0.,20.,30.480,50.,70.,90.,110.,225.,350.,475./              
          
C     ----------SETING UP ARRAY OF REFRACTIVITY------------------ 
          
      DN=-7.32*EXP (0.005577*ENS)     
      CE=LOG(ENS/(ENS+DN))           
      AZ=6370.          
      DUM=0.0           
      AS=AZ+HS          
      DO 10 I=1,25      
      EN(I)=EXP (-CE*H(I) )*ENS*1.E-6               
      RI(I)=1.+EN(I)    
      R(I)=   AZ+H(I)+HS              
   10 CONTINUE          
      DO 20 I=2,25      
      K=I-1             
      DN2N=LOG(RI(I))-LOG(RI(K))    
      DR2R=LOG(R(I))-LOG(R(K))      
      A(I)=DN2N/DR2R    
   20 CONTINUE          
      RYTRA=DUM         
      TT=0.             
      RETURN            
          
C     ----------ENTRANCE FOR TRACING RAY------------------------- 
          
      ENTRY TRCRY(TT)                 
      TE=TT             
      RC=   AZ+HC+HS    
      RA=   AZ+HA+HS    
      ENC=  +1.E-6*ENS*EXP (-CE*HC)   
      RIC=1.+ENC        
      ENA=  +1.E-6*ENS*EXP (-CE*HA)   
      RIA=1.+ENA        
      BALL=0.           
      ATE=TE            
      IF(TE.GE.0.) GO TO 41           
      IF(R(1).EQ.RC) GO TO 73         
      X=R(1)/(2.*RC)    
      Z=(RC-R(1))/R(1)                
      W=(EN(1)-ENC)/RIC               
      TEG=-2.*ASIN (SQRT (X*(Z-W)))   
      GO TO 72          
   73 TEG=0.0           
   72 IF(TE.LT.TEG) TE=TEG            
      ATE=ABS (TE)      
      IF(TE.GE.0.) GO TO 41           
      DO 70 I=2,25      
       Y=2.*(SIN (0.5*ATE))**2        
      Z=(R(I)-RC)/RC    
      W=(ENC-EN(I))*COS (ATE)/RI(I)   
      X=Y+Z-W           
      IF(X.LT.0.0) GO TO 70           
      CT=SQRT (0.5*RC*X/R(I))         
      IF(CT.LE.1.) GO TO 60           
   70 CONTINUE          
   60 CT=2.*ASIN (CT)                 
      BALL=2.*CT*(-A(I)/(A(I)+1.))    
      TEI(I)=CT         
      NK=I+1            
      DO 80 I=NK,25     
      RT=R(I)           
      RIT=RI(I)         
      IF(RT-RC)62,62,61               
   61 RT=RC             
      RIT=RIC           
   62 L=I-1             
      X=RI(L)*R(L)/(RIT*RT)           
      TEI(I)=ACOS (COS (TEI(L))*X)    
      X=2.*(-A(I))/(A(I)+1.)          
      BALL=BALL+(TEI(I)-TEI(L))*X     
      NLA=I             
      IF(RT.EQ.RC) GO TO 40           
   80 CONTINUE          
   40 CONTINUE          
      IF(NLA.LT.2) NLA=2              
      LL=NLA-1          
      TEI(LL)=ATE       
      DO 90 I=NLA,25    
      LC=I-1            
      RT=R(I)           
      RIT=RI(I)         
      ENT=EN(I)         
      IF(RT-RA)47,47,46               
   46 RIT=RIA           
      RT=RA             
      ENT=ENA           
   47 X=RC/(2.*RT)      
      Y=2.*(SIN (0.5*ATE))**2         
      Z=(RT-RC)/RC      
      W=(ENC-ENT)*COS (ATE)/RIT       
      TEI(I)=2.*ASIN (SQRT (X*(Y+Z-W)))             
      X= -A(I)/(A(I)+1.)              
      BALL=BALL+((TEI(I)-TEI(LC))*X)                
      TEA=TEI(I)        
      IF(R(I).GT.RA) GO TO 100        
   90 CONTINUE          
      X=RI(25)*R(25)/RA               
      TEA=ACOS (COS (TEA)*X)          
  100 CA=(TEA-TE+BALL)                
      D=AS*CA           
      DN=D*.5399568034                
      CT=COS (BALL)     
      ST=SIN (BALL)     
      TNT=TAN (TEA)     
      Y=RIA/RIC         
      X=(CT-ST*TNT-Y)/(Y*TAN (TE)-ST-CT*TNT)        
      X=ATAN (X)        
      CX=TE-X           
      RYTRA=TEA         
      RETURN            
   41 DO 85 NL=2,25     
      IF(RC.LE.R(NL)) GO TO 86        
   85 CONTINUE          
      NL=25             
   86 NLA=NL            
      GO TO 40          
      END               
      SUBROUTINE SCAT                 
        implicit real*8(A-H,O-Z)
****  850402 CORRECTIONS INCLUDED     
C     ROUTINE FOR MODEL APR 77        
      COMMON/PARAM/HE(2),D,DL(2),ENS,EFRTH,FREK,ALAM,THE(2),KD,GAO,GAW          
      COMMON/SCTBL/HS(2),ALSC,TWEND,THRFK,HL(2),THETA,HTP,AA,REW  
      COMMON/PRINT/KL,IN,IOT          
          
      DIMENSION RE(2),DZ(2),THA(2),RQ(2)            
      DATA GMA/157.E-6/               
          
C     LOSSES IN THE TROPOSCAT MODE    
C     VALID ONLY ON BEYOND-THE-HORIZON PATHS        
C     ASSUMES EFRTH HAS BEEN OBTAINED STRAIGHT FORWARDLY FROM ENS 
C          AND THAT IT IS LARGER THAN THE GEOGRAPHIC RADIUS.      
          
      GME=1./EFRTH      
          
C    THA ARE RAY SLOPES AT HORIZON    
          
      THA(1)=0.         
      THA(2)=0.         
      THETA=0.          
      HC=AMAX1(HL(1),HL(2))-HTP       
      IF(KD.LE.1) GO TO 3             
      DO 1 J=1,2        
      Q0=DL(J)*GME      
      THA(J)=((HL(J)-HS(J))/(HS(J)+EFRTH)+2.*SIN(Q0*0.5)**2)/SIN(Q0)            
    1 THETA=THA(J)+THETA              
      IF(THETA)2,3,3    
    2 Q0=DL(1)/(DL(1)+DL(2))          
      THA(1)=THA(1)-THETA*Q0          
      THA(2)=THA(2)-THETA*(1.-Q0)     
      THETA=0.          
    3 DZ(2)=D-DL(1)-DL(2)             
      IF(DZ(2).GT.0.) GO TO 10        
          
C     THE SINGLE KNIFE-EDGE IS A SPECIAL CASE       
          
      DZ(1)=0.          
      DZ(2)=0.          
      IF(THETA)900,900,20             
          
C     THE CROSS-OVER POINT.  EXPONENTIAL ATMOSPHERE.              
          
   10 Q0=DZ(2)*GME      
      THETA=Q0+THETA    
      DZ(1)=((Q0*0.5+THA(2))*DZ(2)+HL(2)-HL(1))/THETA             
      DZ(2)=DZ(2)-DZ(1)               
      DN=GMA-GME        
      IF(DN.LE.0.) GO TO 9011         
      HH=ENS*1E-6/DN    
      THETA=0.          
      DO 18 J=1,2       
      X=DZ(J)           
      P=THA(J)          
      Z0=HL(J)-HTP      
      Z1=(GME*X*0.25+P)*X*0.5+Z0      
      Z2=(GME*X*0.5+P)*X+Z0           
      Q0=GMA-DN/EXP(AMIN1(35.,Z0/HH))               
      I=2               
   11 Q1=GMA-DN/EXP(AMIN1(35.,Z1/HH))               
      Q2=GMA-DN/EXP(AMIN1(35.,Z2/HH))               
      Z2=((Q0+2.*Q1)*X/6.+P)*X+Z0     
      I=I-1             
      IF(I)13,13,12     
   12 Z1=((7.*Q0+6.*Q1-Q2)*X/48.+P)*X*0.5+Z0        
      GO TO 11          
          
C     THA BECOME RAY SLOPES AT CROSS-OVER           
C     RE ARE CROSS-OVER HEIGHT ESTIMATES            
          
   13 THA(J)=(Q0+4.*Q1+Q2)*X/6.+P     
      THETA=THA(J)+THETA              
      RE(J)=Z2          
   18 CONTINUE          
      IF(THETA.LE.0.) GO TO 9012      
      X=(RE(2)-RE(1))/THETA           
      DZ(1)=DZ(1)+X     
      DZ(2)=DZ(2)-X     
      IF(DZ(1).LE.0..OR.DZ(2).LE.0.) GO TO 9013     
      HC=(THA(1)*RE(2)+THA(2)*RE(1))/THETA          
   20 WN=FREK/0.0477    
      Q2=2.*WN*THETA    
      DST=DZ(1)         
      DSR=DZ(2)         
      DS=DST+DSR        
          
C     DZ BECOME SLANT RANGES FROM ANTENNA TO CROSS-OVER           
C     DQ IS TOTAL SLANT RANGE         
C     RE BECOME NORMALIZED ANTENNA HEIGHTS          
          
      DQ=0.             
      DO 21 J=1,2       
      X=(HS(J)-HL(J))**2+4.*(EFRTH +HS(J))*(EFRTH+HL(J))          
     X *SIN(DL(J)*0.5/EFRTH)**2       
      DZ(J)=SQRT(X)+DZ(J)             
      DQ=DZ(J)+DQ       
   21 RE(J)=Q2*HE(J)    
      Q2=(200.E-6*ENS-60.E-3)*ENS+6.6               
      Q1=(5.67E-6*ENS-2.32E-3)*ENS+0.031            
      Q1=1.+Q1/EXP(AMIN1(35.,(HC/4.0)**6))          
      Q0=Q1*0.1424      
      ET=Q0*THETA*DQ*0.5              
      SS=(DZ(1)-DZ(2))/DQ             
          
C         ALSC IS ATTENUATION         
          
      ALSC=83.1         
     C     +4.3429*LOG(THETA**3*WN/DQ)             
     C      +8.6859*Q0*HC-17.3718*LOG(Q1)          
     C      -Q2/(1.+(HC/3.6)**2)      
     C      +WSCAT(ET,SS,RE(1),RE(2))               
          
          
C     ----CALCULATION OF OXYGEN AND WATER VAPOR RAYS ------------ 
          
   37 EC1=HS(1)-HTP+EFRTH             
      EC2=HS(2)-HTP+EFRTH             
      HET=HC+EFRTH      
      HER=HC+EFRTH      
      IF(EC1.GT.HET) GO TO 34         
      EX1=EC1           
      HX1=HET           
      AG1=THE(1)        
   35 IF(EC2.GT.HER) GO TO 36         
      EX2=EC2           
      HX2=HER           
      AG2=THE(2)        
   33 IF(DS.GT..001) GO TO 30         
      CALL SORB(EX1,HX1,EFRTH,DL(1),AG1,RQ)         
      REO=RQ(1)         
      REW=RQ(2)         
      CALL SORB(EX2,HX2,EFRTH,DL(2),AG2,RQ)         
      REO=REO+RQ(1)     
      REW=REW+RQ(2)     
   32 AA=GAO*REO+GAW*REW              
      RETURN            
   30 DAT=DL(1)+DST     
      DAR=DL(2)+DSR     
      CALL SORB(EX1,HX1,EFRTH,DAT,AG1,RQ)           
      REO=RQ(1)         
      REW=RQ(2)         
      CALL SORB(EX2,HX2,EFRTH,DAR,AG2,RQ)           
      REO=REO+RQ(1)     
      REW=REW+RQ(2)     
      GO TO 32          
   34 EX1=HET           
      HX1=EC1           
      AG1=-ATAN(THA(1))               
      GO TO 35          
   36 EX2=HER           
      HX2=EC2           
      AG2=-ATAN(THA(2))               
      GO TO 33          
          
C              ERROR RETURN           
          
  900 CONTINUE          
      ALSC=0.           
      THETA=0.          
      DST=0.            
      DSR=0.            
      DS=0.             
      GO TO 37          
9999  FORMAT(5X,A6,8E15.6)            
****  850402            
*011  WRITE(IOT,9999)3HDN=,DN         
9011  WRITE(IOT,9999) 'DN = ',DN      
      GO TO 9010        
****  850402            
*012  WRITE(IOT,9999)3HTH=,THETA      
9012  WRITE(IOT,9999) 'TH = ',THETA   
      GO TO 9010        
****  850402            
*013  WRITE(IOT,9999)3HDZ=,DZ         
9013  WRITE(IOT,9999) 'DZ = ',DZ      
      GO TO 9010        
****  850402            
*010  WRITE(IOT,9999)6HPARAM=,HE,D,DL,THE           
9010  WRITE(IOT,9999) 'PARAM=',HE,D,DL,THE          
*     WRITE(IOT,9999)4HTHA=,THA,THETA               
      WRITE(IOT,9999) 'THA = ',THA,THETA            
*     WRITE(IOT,9999)3HHL=,HL,HTP     
      WRITE(IOT,9999) 'HL = ',HL,HTP                
*     WRITE(IOT,9999)3HDN=,DN,GME,RE,X              
      WRITE(IOT,9999) 'HDN = ',DN,GME,RE,X          
      STOP 'BAD SCAT'                 
          
      END               
      SUBROUTINE SCNTL                
        implicit real*8(A-H,O-Z)
          
C     ROUTINE FOR MODEL APR 77        
          
      COMMON/SOI/A,F,IOS,IPK,V(35)    
      DIMENSION VN(35,6)              
      DIMENSION VL(4)                 
          
      DATA (VN(I,1),I=1,35)/35*0.0/   
      DATA (VN(I,2),I=1,35)/4.4,4.2,3.9,3.6,3.4,3.05,2.8,2.55,2.15,1.85,        
     21.55,1.1,0.8,0.6,0.4,0.2,0.05,0.0,-0.15,-0.35,-0.6,-0.7,-0.9,-1.3,        
     2-1.8,-2.15,-2.55,-3.05,-3.5,-4.0,-4.6,-5.0,-5.5,-6.2,-6.6/  
      DATA (VN(I,3),I=1,35)/5.8,5.5,5.0,4.7,4.4,3.95,3.65,3.35,2.9,2.55,        
     32.2,1.65,1.2,0.95,0.75,0.45,0.2,0,-0.3,-.55,-.9,-1.1,-1.4,-2.0,-2.        
     36,-3.1,-3.65,-4.4,-5.,-5.6,-6.5,-7.2,-8.,-9.2,-10.1/        
      DATA (VN(I,4),I=1,35)/10.8,10.2,9.2,8.2,7.6,6.5,5.8,5.2,4.35,3.8,3        
     4.25,2.5,1.8,1.4,1.1,0.7,0.3,0,-0.35,-0.7,-1.2,-1.5,-1.9,-2.6,-3.5,        
     4-4.2,-4.95,-6.05,-7.0,-7.9,-9.1,-10.0,-11.0,-12.4,-13.4/    
      DATA (VN(I,5),I=1,35)/13.7,13.,12.,11.1,10.3,9.2,8.4,7.5,6.35,5.45        
     5,4.55,3.4,2.55,2.0,1.6,0.98,0.5,0,-0.5,-1.0,-1.8,-2.4,-3.1,-4.4,-6        
     5.0,-7.2,-8.5,-10.,-11.1,-12.2,-14.,-15.4,-16.7,-18.4,-19.5/ 
      DATA (VN(I,6),I=1,35)/14.5,13.8,13.,12.4,11.7,10.7,10.,9.5,8.2,7.5        
     6,6.55,5.2,4.0,3.2,2.6,1.6,0.75,0.0,-0.8,-1.7,-2.95,-3.8,-5.0,-7.05        
     6,-9.8,-12.2,-15.,-19.1,-22.2,-25.8,-30.2,-34.,-37.2,-42.5,-46./           
      DATA VL/0.,4.3,58.,69.7/        
      AVEF(YN,XN,YN1,XN1) = (YN1*(T - XN) - YN*(T - XN1))/(XN1 - XN)            
          
      T=A               
      IF(IOS.LT.0) GO TO 20           
      IS1=IOS+1         
      IF(IS1.LT.1) IS1=1              
      IF(IS1.GT.6) IS1=6              
   19 DO 4 J=1,35       
      V(J)=VN(J,IS1)    
    4 CONTINUE          
   17 IF(IPK.GT.0) GO TO 10           
    8 RETURN            
   10 IF(T.LE.17..OR.T.GE.52.) GO TO 2              
      IF(T.LT.24.) GO TO 11           
      IF(T.GT.45.) GO TO 12           
      PK=2.             
      GO TO 6           
    2 PK=1.             
    6 X=(136./F)**PK    
      DO 7 J=1,35       
    7 V(J)=V(J)*X       
      GO TO 8           
   11 PK=1.+((T-17.)/7.)              
      GO TO 6           
   12 PK=1.+ABS ((52.-T)/7.)          
      GO TO 6           
   20 IF(T.LE.VL(1).OR.T.GE.VL(4)) GO TO 21         
      IF(T.LT.VL(2)) GO TO 22         
      IF(T.LE.VL(3)) GO TO 23         
      L=3               
      M=4               
      K=1               
      N=2               
   18 DO 24 J=1,35      
   24 V(J)=AVEF(VN(J,K),VL(L),VN(J,N),VL(M))        
      GO TO 17          
   21 IS1=2             
      GO TO 19          
   22 L=1               
      M=2               
      K=2               
      N=1               
      GO TO 18          
   23 IS1=1             
      GO TO 19          
      END               
      SUBROUTINE SIG(KSC,FK)          
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
      COMMON /SEA/ ISK,SCK,TP,JM,EPK,SGM            
      DIMENSION SG(7),EP(7)           
      DIMENSION ST(3,2),SM(3,2),TE(3,2)             
      DATA EP/81.,25.,15.,4.,81.,5.,1./             
      DATA SG/5.,.02,.005,.001,.010,.010,10.E+06/   
      DATA ((SM(I,J),I=1,3),J=1,2)/2.70E10,3.7E10,4.7E10,1.E8,1.E8,1.E8/        
      DATA ((ST(I,J),I=1,3),J=1,2)/75.,72.,69.,88.,84.,80./       
      DATA ((TE(I,J),I=1,3),J=1,2)/1.69E-11,1.21E-11,9.2E-12,1.87E-11,1.        
     X36E-11,1.01E-11/                
          
      E0=4.9            
      IR=KSC            
      I=IR              
      IF(IR.EQ.1) GO TO 26            
      IF(IR.EQ.5) GO TO 27            
      EPK=EP(I)         
      SGM=SG(I)         
   18 RETURN            
   26 K=1               
      GO TO 28          
   27 K=2               
   28 IF(TP.LT.9.99) GO TO 38         
      IF(TP.GT.10.01) GO TO 39        
      J=2               
   40 ES=ST(J,K)        
      T=TE(J,K)         
      SMI=SM(J,K)       
      FC=FK*10.E5       
      ALM=299.7925/FK                 
      T1=FC*T*6.28318531              
      T2=ES-E0          
      T3=1.+(T1*T1)     
      T4=2.*SMI/FC      
      T5=60.*ALM        
      EPK=(T2/T3)+E0    
      SGM=((T2*T1/T3)+T4)/T5          
      GO TO 18          
   38 J=1               
      GO TO 40          
   39 J=3               
      GO TO 40          
      END               
      SUBROUTINE SORB(H1,H2,A,R0,CA,RE)             
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
      COMMON/CURVE/PI,RAD,DEG,TWPI,PI2              
      DIMENSION RE(2),TE(2),H(2)      
      TE(1)=3.25        
      TE(2)=1.36        
      BA=CA             
      IF(H1.GT.H2) GO TO 10           
      HS=H1             
      HL=H2             
   11 AT=PI2+BA         
      ANUM=HS*SIN (AT)                
      DO 22 K=1,2       
      H(K)=TE(K)+A      
      IF(HL.LE.H(K)) GO TO 83         
      IF(H(K).LT.HS) GO TO 81         
      AS=ASIN(ANUM/H(K))              
      AE=PI-(AT+AS)     
      IF(BA.GT.1.5620) GO TO 24       
      IF(AE.EQ.0.) GO TO 24           
      RE(K)=(HS*SIN (AE))/SIN (AS)    
      GO TO 22          
   81 IF(AT.GT.PI2) GO TO 85          
      HC=HS*SIN(AT)     
      IF(H(K).LE.HC) GO TO 85         
      RE(K)=2.*H(K)*SIN(ACOS(HC/H(K)))              
      GO TO 22          
   83 RE(K)=R0          
      GO TO 22          
   85 RE(K)=0.          
      GO TO 22          
   24 RE(K)=H(K)-HS     
   22 CONTINUE          
      RETURN            
   10 HS=H2             
      HL=H1             
      BA=-(CA+(R0/A))                 
      GO TO 11          
      END               
      FUNCTION TABLE(XINT)            
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
C     ENTER TNTER WITH DELTA R AND GET SI           
C     ENTER DNTER WITH DELTA R AND GET DISTANCE     
C     ENTER SNTER WITH DISTANCE AND GET SI          
          
      SAVE              
      COMMON/SPLIT/L1,L2,N,X(140),Y(140),D6(140),XS(55),XD(55),XR(55),YS        
     X(55),YD(55),YR(55),L3,ZS(25),ZD(25),ZR(25)    
      COMMON/PRINT/KL,IN,IOT          
      DIMENSION AS(110),AD(110),AR(110)             
C     -----------------SET UP ARRAY-------------------------------------        
      DUM=0.            
      CALL TRMES(XS,XD,XR,L1,YS,YD,YR,L2,AS,AD,AR,L5)             
      CALL TRMES(AS,AD,AR,L5,ZS,ZD,ZR,L3,Y,X,D6,N)  
      M=N               
      TABLE=DUM         
      RETURN            
          
  101 FORMAT(31H OUT OF RANGE FOR INTERPOLATION)    
          
      ENTRY TNTER(XINT)               
      IF(XINT-X(1))7,1,2              
    1 YINT=Y(1)         
      TABLE=YINT        
      RETURN            
    2 K=1               
    3 IF(XINT-X(K+1))6,4,5            
    4 YINT=Y(K+1)       
      TABLE=YINT        
      RETURN            
    5 K=K+1             
      IF(M-K)8,8,3      
    6 YINT=((XINT-X(K))*(Y(K+1)-Y(K))/(X(K+1)-X(K)))+Y(K)         
      TABLE=YINT        
      RETURN            
    7 WRITE(IOT,101)    
      TABLE=Y(1)        
      RETURN            
    8 WRITE(IOT,101)    
      TABLE=Y(M)        
      RETURN            
          
      ENTRY DNTER(XINT)               
      IF(XINT-X(1))17,11,12           
   11 TABLE=D6(1)       
      RETURN            
   12 K=1               
   13 IF(XINT-X(K+1))16,14,15         
   14 TABLE=D6(K+1)     
      RETURN            
   15 K=K+1             
      IF(M-K)18,18,13                 
   16 TABLE=((XINT-X(K))*(D6(K+1)-D6(K))/(X(K+1)-X(K)))+D6(K)     
      RETURN            
   17 WRITE(IOT,101)    
      TABLE=D6(1)       
      RETURN            
   18 WRITE(IOT,101)    
      TABLE =D6(M)      
      RETURN            
          
      ENTRY SNTER(XINT)               
      IF(XINT-D6(1))32,31,37          
   31 TABLE=Y(1)        
      RETURN            
   32 K=1               
   33 IF(XINT-D6(K+1))35,34,36        
   34 TABLE=Y(K+1)      
      RETURN            
   35 K=K+1             
      IF(M-K)38,38,33                 
   36 TABLE=((XINT-D6(K))*(Y(K+1)-Y(K))/(D6(K+1)-D6(K)))+Y(K)     
      RETURN            
   37 WRITE(IOT,101)    
      TABLE=Y(1)        
      RETURN            
   38 WRITE(IOT,101)    
      TABLE =Y(M)       
      RETURN            
      END               
      SUBROUTINE TIMBK(NTB,DE,G1,G9)                
        implicit real*8(A-H,O-Z)
      DIMENSION C1(3,11),C2(3,11),C3(3,11),CN1(3,11),CN2(3,11),CN3(3,11)        
     X,FM(3,11),FIN(3,11),CC(35),Z(3)               
      COMMON/VARY/KLM,MX1,KLM2,MX2,FKE,AD(35)       
          
      DATA C1/1.09E-4,1.72E-4,2.11E-4,1.04E-5,1.05E-5,0.00,       
     C2.02E-4,3.64E-5,3.47E-4,1.22E-2,1.84E-4,1.35E-6,2.58E-4,3.80E-4,          
     C1.05E-6,3.84E-3,1.81E-3,2.04E-4,7.95E-3,3.19E-3,8.00E-4,1.70E-4,          
     C1.64E-6,3.63E-4,4.47E-3,7.42E-4,1.18E-5,2.46E-4,3.45E-6,1.40E-3,          
     C5.25E-4,2.93E-4,1.63E-4/        
          
      DATA C2/1.21E-6,6.39E-8,3.44E-17,4.28E-8,7.00E-13,0.00,1.45E-6,           
     C3.74E-9,3.76E-14,9.81E-6,2.22E-6,5.02E-18,3.41E-6,4.76E-6,5.02E-18        
     C,4.22E-5,5.82E-6,6.61E-18,3.76E-5,2.51E-6,3.91E-16,7.93E-7, 
     C1.43E-7,1.80E-23,1.66E-5,5.55E-5,6.72E-17,1.74E-7,1.25E-8,1.79E-34        
     C,1.57E-6,3.78E-8,1.81E-25/      
          
      DATA C3/8.29E-8,2.93E-10,1.73E-4,3.51E-8,7.64E-9,0.00,4.27E-8,            
     C3.53E-7,5.42E-4,1.09E-8,3.65E-16,1.32E-7,2.01E-11,8.39E-17, 
     C4.14E-9,7.76E-9,6.37E-13,2.82E-9,3.19E-8,5.03E-9,1.20E-5,1.29E-7,         
     C3.14E-7,1.55E-5,2.06E-8,4.37E-8,1.65E-6,1.27E-8,7.50E-7,1.05E-5,          
     C4.70E-7,1.02E-7,8.12E-6/        
          
      DATA CN1/2.28,2.10,1.67,2.71,2.59,0.0,2.15,2.40,1.60,1.36,2.09,           
     C2.80,2.05,1.92,2.70,1.57,1.67,1.87,1.47,1.60,1.68,2.19,3.08,              
     C1.65,1.55,1.84,2.40,2.11,2.87,1.27,1.97,2.00,1.80/          
          
      DATA CN2/2.29,2.79,6.52,2.91,4.80,0.0,2.28,3.28,5.30,2.00,2.29,           
     C6.74,2.25,2.19,6.74,1.76,2.15,6.67,1.76,2.27,5.94,2.37,2.66,8.91,         
     C1.90,1.69,6.32,2.64,3.07,13.23,2.31,2.88,9.59/              
          
      DATA CN3/3.26,4.24,1.82,3.41,3.68,0.0,3.37,2.94,1.58,3.58,6.82,           
     C3.08,4.78,7.10,3.70,3.66,5.38,3.76,3.40,3.69,2.25,3.18,3.03,2.36,         
     C3.48,3.28,2.61,3.62,2.82,2.51,2.90,3.15,2.32/               
          
      DATA FM/9.6,8.2,1.2,9.15,7.05,0.0,9.4,7.8,1.3,10.8,8.0,5.2,8.0,           
     C6.6,2.4,9.6,8.4,5.2,11.2,10.0,7.1,9.5,8.6,1.95,9.98,8.25,5.1,9.37,        
     C7.92,1.05,10.0,8.2,3.0/         
          
      DATA FIN/2.8,2.4,0.5,2.8,2.8,0.0,2.8,2.2,0.6,5.5,4.0,4.0,4.0,3.3,         
     C1.8,5.2,4.1,4.2,5.5,4.4,5.6,3.0,2.6,0.8,5.1,4.0,4.0,2.8,2.45,0.5,         
     C5.4,3.2,1.9/      
          
      DATA CC/3.815,3.69,3.456,3.332,3.168,2.918,2.731,2.528,2.2,1.9507,        
     X1.7166,1.3265,1.,.8087,.6567,.4092,.1976,0.,.1976,.4092,.6567,.808        
     X7,1.,1.2835,1.6025,1.815,2.0098,2.2458,2.4112,2.5675,2.7622,2.9025        
     X,3.0357,3.2052,3.3279/          
          
      IF(NTB.GT.10) GO TO 8           
      J=NTB             
   12 DO 13 I=1,3       
      X=FIN(I,J)+((FM(I,J)-FIN(I,J))*EXP(-C2(I,J)*DE**CN2(I,J)))  
   13 Z(I)=(((C1(I,J)*DE**CN1(I,J))-X)*EXP(-C3(I,J)*DE**CN3(I,J)))+X            
      X10=Z(1)*G1       
      X90=-Z(2)*G9      
      DO 10 I=1,17      
   10 AD(I)=CC(I)*X10+Z(3)            
      AD(18)=Z(3)       
      DO 9 I=19,35      
    9 AD(I)=CC(I)*X90+Z(3)            
      RETURN            
    8 J=NTB-10          
      IF(NTB.GT.18)J=11               
      GO TO 12          
      END               
      SUBROUTINE TRMES(A,B,C,NA,R,S,T,NR,X,Y,Z,N)   
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
      DIMENSION A(*),B(*),C(*),R(*),S(*),T(*),X(*),Y(*),Z(*)      
      I=1               
      J=1               
      N=0               
    4 N=N+1             
      IF(A(I)-R(J))9,7,8              
    9 X(N)=A(I)         
      Y(N)=B(I)         
      Z(N)=C(I)         
      I=I+1             
      IF(I-NA)4,4,5     
    8 X(N)=R(J)         
      Y(N)=S(J)         
      Z(N)=T(J)         
      J=J+1             
      IF(J-NR)4,4,3     
    7 X(N)=A(I)         
      Y(N)=B(I)         
      Z(N)=C(I)         
      I=I+1             
      J=J+1             
      IF(I-NA)11,11,10                
   10 IF(J-NR)5,5,12    
   11 IF(J-NR)4,4,3     
    5 LI=J              
      DO 16 LE=LI,NR    
      N=N+1             
      X(N)=R(LE)        
      Y(N)=S(LE)        
      Z(N)=T(LE)        
   16 CONTINUE          
      GO TO 12          
    3 LI=I              
      DO 18 LE=LI,NA    
      N=N+1             
      X(N)=A(LE)        
      Y(N)=B(LE)        
      Z(N)=C(LE)        
   18 CONTINUE          
   12 RETURN            
      END               
      SUBROUTINE TMESH(A,NA,R,NR,X,N)               
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
      DIMENSION A(*),R(*),X(*)        
      I=1               
      J=1               
      N=0               
    4 N=N+1             
      IF(A(I)-R(J))9,7,8              
    9 X(N)=A(I)         
      I=I+1             
      IF(I-NA)4,4,5     
    8 X(N)=R(J)         
      J=J+1             
      IF(J-NR)4,4,3     
    7 X(N)=A(I)         
      I=I+1             
      J=J+1             
      IF(I-NA)11,11,10                
   10 IF(J-NR)5,5,12    
   11 IF(J-NR)4,4,3     
    5 LI=J              
      DO 16 LE=LI,NR    
      N=N+1             
      X(N)=R(LE)        
   16 CONTINUE          
      GO TO 12          
    3 LI=I              
      DO 18 LE=LI,NA    
      N=N+1             
      X(N)=A(LE)        
   18 CONTINUE          
   12 RETURN            
      END               
      SUBROUTINE VARIB(DE)            
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
    5 FORMAT(' ')       
  766 FORMAT(20X,'WITH CLIMATE ',A30)               
  821 FORMAT(20X,'WITH TIME BLOCK',I2)              
  822 FORMAT(20X,'WITH SUMMER TIME BLOCK')          
  823 FORMAT(20X,'WITH WINTER TIME BLOCK')          
  824 FORMAT(20X,'WITH ALL HOURS TIME BLOCK')       
      DIMENSION AD1(35),AD2(35)       
      CHARACTER*30 VKM(7)             
      SAVE              
      COMMON/VARY/KLM,MX1,KLM2,MX2,FKE,AD(35)       
      COMMON/PRINT/KL,IN,IOT          
          
      DATA VKM/'1 EQUATORIAL','2 CONTINENTAL SUBTROPICAL',        
     X'3 MARITIME SUBTROPICAL','4 DESERT','6 CONTINENTAL TEMPERATE',            
     X'7A MARITIME TEMPERATE OVERLAND','7B MARITIME TEMPERATE OVERSEA'/         
      F=FKE             
      K=0               
      ENTRY VARLB(DE)                 
      IF(KLM.LE.0) GO TO 21           
      IF(KLM.GT.0.AND.KLM.LT.8) GO TO 22            
      IF(KLM.GT.8.AND.KLM.LT.20) GO TO 23           
      KLM=0             
      GO TO 21          
   20 KLM=0             
      GO TO 21          
   21 IF(K.GT.0) GO TO 30             
      IF(F.GT.1600.) GO TO 4          
      QG1=(.21*SIN (5.22*LOG10(F/200.)))+1.28      
      QG9=(.18*SIN (5.22*LOG10(F/200.)))+1.23      
      IF (MX1 .EQ. 0) GO TO 6         
      QG12=QG1          
      QG92=QG9          
    6 CONTINUE          
      ASSIGN 31 TO J    
      K=1               
      GO TO 30          
   22 IF(K.GT.0) GO TO 24             
      ASSIGN 32 TO J    
      K=1               
      IF (MX1 .EQ. 0) CALL QVARB(KLM)               
   24 WRITE(IOT,766)VKM(KLM)          
      GO TO 30          
   23 IF(K.GT.0) GO TO 25             
      IF(F.GT.1600.) GO TO 14         
      QG1=(.21*SIN (5.22*LOG10(F/200.)))+1.28      
      QG9=(.18*SIN (5.22*LOG10(F/200.)))+1.23      
      IF (MX1 .EQ. 0) GO TO 16        
      QG12=QG1          
      QG92=QG9          
   16 CONTINUE          
      NTB=KLM-10        
      ASSIGN 33 TO J    
      K=1               
   25 CONTINUE          
      IF(NTB.GT.0.AND.NTB.LT.9) WRITE(IOT,821)NTB   
      IF(KLM.EQ.19) WRITE(IOT,824)    
      IF(KLM.EQ.10) WRITE(IOT,823)    
      IF(KLM.EQ.9) WRITE(IOT,822)     
      GO TO 30          
   30 WRITE(IOT,5)      
      WRITE(IOT,5)      
      RETURN            
      ENTRY VARDT(DE)                 
      GO TO J,(31,32,33)              
   31 DEE=DE            
      CALL VZD(DEE,QG1,QG9,AD)        
      IF (MX1 .EQ. 0) GO TO 62        
      DO 60 I=1,35      
60    AD1(I)=AD(I)      
      CALL VZD(DEE,QG12,QG92,AD2)     
      CALL COMB(AD1,AD2)              
62    CONTINUE          
      RETURN            
   32 DEE=DE            
      IF (MX1 .NE. 0) GO TO 39        
      CALL QQVARB(DEE)                
      RETURN            
39    CONTINUE          
      CALL QVARB(KLM)                 
      CALL QQVARB(DEE)                
      DO 40 I=1,35      
40    AD1(I)=AD(I)      
      CALL QVARB(KLM2)                
      CALL QQVARB(DEE)                
      DO 41 I=1,35      
41    AD2(I)=AD(I)      
      CALL COMB(AD1,AD2)              
      RETURN            
   33 DEE=DE            
      CALL TIMBK(KLM,DEE,QG1,QG9)     
      IF (MX1 .EQ. 0) GO TO 52        
      DO 50 I=1,35      
50    AD1(I)=AD(I)      
      CALL TIMBK(KLM2,DEE,QG12,QG92)                
      DO 51 I=1,35      
51    AD2(I)=AD(I)      
      CALL COMB(AD1,AD2)              
52    CONTINUE          
      RETURN            
    4 QG1=1.05          
      QG9=1.05          
      IF (MX1 .EQ. 0) GO TO 6         
      QG12=1.05         
      QG92=1.05         
      GO TO 6           
   14 QG1=1.05          
      QG9=1.05          
      IF (MX1 .EQ. 0) GO TO 16        
      QG12=1.05         
      QG92=1.05         
      GO TO 16          
      END               
      SUBROUTINE VZD(DE,G1,G9,A)      
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
      DIMENSION B(35)                 
      DIMENSION C1(3),C2(3),C3(3),CN1(3),CN2(3),CN3(3),FM(3),FIN(3),Z(3)        
     1,Y(35),A(50)      
          
C     MIXED--ALL YEAR TIME BLOCK YS AND CONTINENTAL V(50)         
          
      DATA C1/2.93E-4,5.25E-4,1.59E-5/              
      DATA C2/3.78E-8,1.57E-6,1.56E-11/             
      DATA C3/1.02E-7,4.70E-7,2.77E-8/              
      DATA CN1/2.00,1.97,2.32/        
      DATA CN2/2.88,2.31,4.08/        
      DATA CN3/3.15,2.90,3.25/        
      DATA FIN/3.2,5.4,0.0/           
      DATA FM/8.2,10.0,3.9/           
          
      DO 13 I=1,3       
      X=FIN(I)+((FM(I)-FIN(I))*EXP (-C2(I)*DE**CN2(I)))           
   13 Z(I)=(((C1(I)*DE**CN1(I))-X)*EXP (-C3(I)*DE**CN3(I)))+X     
      Y(13)=-Z(1)*G9    
      Y(23)=Z(2)*G1     
      Y(1)=3.3279*Y(13)               
      Y(2)=3.2052*Y(13)               
      Y(3)=3.0357*Y(13)               
      Y(4)=2.9025*Y(13)               
      Y(5)=2.7622*Y(13)               
      Y(6)=2.5675*Y(13)               
      Y(7)=2.4112*Y(13)               
      Y(8)=2.2458*Y(13)               
      Y(9)=2.0098*Y(13)               
      Y(10)=1.8150*Y(13)              
      Y(11)=1.6025*Y(13)              
      Y(12)=1.2835*Y(13)              
      Y(14)=0.8087*Y(13)              
      Y(15)=0.6567*Y(13)              
      Y(16)=0.4092*Y(13)              
      Y(17)=0.1976*Y(13)              
      Y(18)=0.0000      
      Y(19)=0.1976*Y(23)              
      Y(20)=0.4092*Y(23)              
      Y(21)=0.6567*Y(23)              
      Y(22)=0.8087*Y(23)              
      Y(24)=1.3265*Y(23)              
      Y(25)=1.7166*Y(23)              
      Y(26)=1.9507*Y(23)              
      Y(27)=2.2000*Y(23)              
      Y(28)=2.5280*Y(23)              
      Y(29)=2.7310*Y(23)              
      Y(30)=2.9180*Y(23)              
      Y(31)=3.1680*Y(23)              
      Y(32)=3.3320*Y(23)              
      Y(33)=3.4560*Y(23)              
      Y(34)=3.6900*Y(23)              
      Y(35)=3.8150*Y(23)              
      DO 18 I=1,35      
      KN=36-I           
      B(I)=Y(I)+Z(3)    
      A(KN)= B(I)       
   18 CONTINUE          
      RETURN            
      END               
      FUNCTION WSCAT(H,S,R1,R2)       
        implicit real*8(A-H,O-Z)
      DIMENSION R(2)    
      R(1)=R1           
      R(2)=R2           
      C=12.             
      D=H               
      S1=S              
      S2=S1**2          
      B=1.-S2           
      A=B**2*D+S2*8.+6.               
      DO 10 I=1,2       
      S1=-S1            
      Q1=R(I)            
      IF(Q1) 11,11,12    
11    R(I)=0.           
      GO TO 10          
12    R(I)=1.4142/Q1     
      Q=Q1**2            
      U=(1.-S1)**2*H    
      V=U**2+Q          
      U=U**2/V          
      V=Q/V             
      A=(1.+S1)*8.*U*V+A              
      B=(U*2.+1.)*B     
      C=(R(I)+1.)**2*C                
      D=D/V             
10    CONTINUE 
      A=(B*2.+A)*D      
      Q=R(1)*R(2)       
      IF(Q .NE. 0.) Q=Q*2./(R(1)+R(2))              
      A=C/(Q+1.)+A      
      WSCAT=LOG(A)*4.3429 
      RETURN            
      END               
      SUBROUTINE YIKK(T,PV,V)         
        implicit real*8(A-H,O-Z)
C     ROUTINE FOR MODEL APR 77        
          
C     THIS NAKAGAMA-RICE DIST.  HAS TABLES FROM NORTON 55 IRE PAGE 1360         
C     THE TABLES ARE THE NEGATIVE OF THE KK IRE TABLES BUT ARE CHANGED          
C     BEFORE GOING OUT OF THE ROUTINE               
C     K HAS THE OPPOSITE SIGN OF 101 BUT THE SAME AS THE IRE PAPER              
          
      DIMENSION P(35),PV(50),V(50)    
      COMMON/VV/VF(36,17)             
      DATA P/.00001,.00002,.00005,.0001,.0002,.0005,.001,.002,.005,.01,         
     X.02,.05,.10,.15,.20,.30,.40,.50,.60,.70,.80,.85,.90,.95,.98,.99,          
     X.995,.998,.999,.9995,.9998,.9999,.99995,.99998,.99999/      
      AVEF(YN,XN,YN1,XN1) = (YN1*(T - XN) - YN*(T - XN1))/(XN1 - XN)            
          
      DO 1 I = 1,14     
      IF(T - VF(1,I)) 3,2,1           
    1 CONTINUE          
      I = 14            
    2 DO 4 J = 1,35     
      V(J) = VF(J+1,I)                
    4 PV(J) = P(J)      
      GO TO 6           
    3 IF(I.EQ.1) GO TO 2              
      DO 5  J = 1,35    
      V(J) = AVEF(VF(J+1,I-1),VF(1,I-1),VF(J+1,I),VF(1,I))        
    5 PV(J) =  P(J)     
    6 DO 7 J=1,35       
    7 V(J)=-V(J)        
      RETURN            
      END               
