      SUBROUTINE AALOS                
        implicit real*8(A-H,O-Z)
C     L-O-S ROUTINE FOR ATOA          
C     ROUTINE FOR MODEL APR 77        
          
    5 FORMAT(1H )       
  760 FORMAT(1X,F7.2,13F7.1,F6.1,2F5.1,2F6.1)       
  766 FORMAT(' D',A5,' FREE SPACE 50%    5%    95%    90%    99%   99.9%        
     X 99.99%  .01%   .1%     1%    10%  CA DEG   PL   AA   AY     K            
     XDEE')             
          
      DIMENSION ALM(12),QR(35),SID(24),SPGRD(3),RE(2),BD(35)      
      DIMENSION CFK(3),CMK(3),CFM(3),CKM(3),CKN(3)  
      DIMENSION GLD(8),D1(200),D2(200),D3(200),VD(35)             
      DIMENSION HTX(2),Z(2),TEA(2),DA(2),HPR(2)     
      DIMENSION P(35),QC(35),QA(35),PQA(35),PQK(35),QK(35),PQC(35)              
      DIMENSION XCON(5),NTM(5),XDON(5),NDM(5)       
          
      COMMON/DIFPR/HT ,HR ,DH,AED,SLP,DLST,DLSR,IPL,KSC,HLT,HRP,AWD,SWP         
     X,II               
      COMMON/DROPS/RS(35),IZ,STS,FQX,DX,HA1,HA2,EFA,AN            
      COMMON/EGAP/IP,LN               
      COMMON/PARAM/HTE,HRE,D,DLT,DLR,ENS,EFRTH,FREK,ALAM,TET,TER,KD,GAO,        
     XGAW               
      COMMON/PLTD/NU(8),BX(200,8),BY(200,8)         
      COMMON /SEA/ ISK,SCK,TP,JM,EPK,SGM            
      COMMON/SELCT/RQD,JA,JB,DQM,DBY(12),SELT(12),JG              
      COMMON/SIGHT/DCW,HCW,DMAX,DML,DOZ,IK,EAC,H2,ICC,HFC,PRH,DSL1,PIRP,        
     XQG1,QG9,PFY(200,4),KK,ZH,RDHK,ILB,EAL,IFA,IAA,T1T,HLPBW,T2T,H2PBW,        
     XJT,JS,H1,IJ,JC,NPL,JO,RY2,RY1   
      COMMON/SOI/THO,FS,IOS,IPK,QS(35)              
      COMMON/SPLIT/L1,L2,N,X(140),Y(140),D6(140),XS(55),XD(55),XR(55),YS        
     X(55),YD(55),YR(55),L3,ZS(25),ZD(25),ZR(25)    
      COMMON/VARY/KLM,MX1,KLM2,MX2,FKE,AD(35)       
      COMMON/CURVE/PI,RAD,DEG,TWPI,PI2              
      COMMON/PRINT/KL,IN,IOT          
          
      CHARACTER UND(3)*5              
          
      COMPLEX AT1,AT2                 
          
      DATA ALM/-6.2,-6.15,-6.08,-6.0,-5.95,-5.88,-5.8,-5.65,-5.35,-5.0,-        
     X4.5,-3.7/         
      DATA  CFK/.001,.0003048,.0003048/             
      DATA  CFM/1.,.3048,.3048/       
      DATA  CKM/1000.,3280.839895,3280.839895/      
      DATA  CKN/1.,.6213711922,.5399568034/         
      DATA  CMK/1.,1.609344,1.852/    
      DATA  GLD/0.,.1,.2,.3,.4,.5,.75,1./           
      DATA  NDM/10,19,8,30,5/         
      DATA  NTM/10,40,50,50,0/        
      DATA  P/.00001,.00002,.00005,.0001,.0002,.0005,.001,.002,.005,.01         
     X,.02,.05,.10,.15,.20,.30,.40,.50,.60,.70,.80,.85,.90,.95,.98,.99,         
     X.995,.998,.999,.9995,.9998,.9999,.99995,.99998,.99999/      
      DATA  SID/.2,.5,.7,1.,1.2,1.5,1.7,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,10.        
     X,20.,45.,70.,80.,85.,88.,89./   
      DATA  SPGRD/0.,.06,.1/          
      DATA UND/' KM  ',' S MI',' N MI'/             
      DATA  XCON/1.,5.,10.,10.,0./    
      DATA  XDON/.1,1.,5.,1.,1./      
          
      FNA(FX,FA,FB,FC,FD)=((FX-FB)*(FC-FD)/(FA-FB))+FD            
          
      BSPI=.3183098862                
      ALIM=3.           
      TWHT=2.*HTE       
      AFP=32.45+20.*ALOG10(FREK)      
      F=FREK            
      JG=0              
      NCT=0             
      NOC=0             
      DUM=0.            
      ALA2=ALAM/2.      
      TWPILA=TWPI/ALAM                
      DTR0=ALAM/6.      
      CPI2=1.56         
      ASPA=0.25         
      ASPB=0.25         
      ASPC=ASPA*ASPB*(6.E-8)*F        
      ERTH =6370.       
      A0=ERTH           
      EFN=EFRTH         
      PKL=((3.*PI)/(ALAM))            
      WRITE(IOT,766)UND(IK)           
      CALL PAGE(1)      
      IF(ICC.GT.0) NOC=1              
      IF(NOC.LE.0) GO TO 502          
      RCW=DCW*.5        
      BTC=ATAN (HFC/RCW)              
      EPG=EPK           
      SGG=SGM           
      CALL SIG(ICC,F)                 
      EPC=EPK           
      SGC=SGM           
      EPK=EPG           
      SGM=SGG           
      ABTC=ABS (BTC)    
      R1C=RCW/COS (BTC)               
      SQVT=SQRT (2.*R1C/ALAM)         
      HDI=HTE-HFC       
      TWHC=2.*HFC       
  503 CONTINUE          
      L2=0              
      L1=0              
      N=0               
          
C     ------SETTING UP OF TABLE OF SI, DELTA R AND DISTANCE------ 
          
      TLE1=ACOS(EFN/(HTE+EFN))        
      TLE2=ACOS(EFN/(HRE+EFN))        
      TLE1=-TLE1        
      TLE2=-TLE2        
      YY1=-RY1          
      YY2=-RY2          
      TWT1=(YY1-TLE1)/(EFN-A0)        
      TWT2=(YY2-TLE2)/(EFN-A0)        
      LE=11             
      LK=0              
   61 LK=LK+1           
      IF(LK.GT.LE) GO TO 21           
      IF(LK.LT.4) GO TO 120           
      LB=13-LK          
      GRD=FLOAT (LB)    
      APDR=ALAM/GRD     
  121 IF(APDR.LE.0.) GO TO 122        
      IF(APDR.GT.TWHT) GO TO 21       
      SI=ASIN (APDR/TWHT)             
      ASSIGN 65 TO KR                 
      GO TO 66          
   65 L1=L1+1           
      XS(L1)=SI         
      XD(L1)=DR         
      XR(L1)=D          
      IF(APDR.LE.0.) GO TO 122        
      SI=SQRT (APDR/(2.*DLST))        
      IF(SI.GT.PI2) SI=PI2            
      ASSIGN 123 TO KR                
      GO TO 66          
  123 L2=L2+1           
      YS(L2)=SI         
      YD(L2)=DR         
      YR(L2)=D          
      GO TO 61          
   21 CONTINUE          
      IF(ILB.LE.0) GO TO 162          
      LA=0              
  150 LA=LA+1           
      IF(LA.GT.10) GO TO 162          
      GND=FLOAT (LA)    
      LG=0              
  151 LG=LG+1           
      IF(LG.GT.4) GO TO 150           
      GO TO (155,156,157,158), LG     
****  850606            
*CDC  ADD CDC GO TO ERROR ROUTINE     
*     CALL GOTOER
****  940914
*MS   FOR MICROSOFT USE DOMAIN ERROR
      GRD=SQRT(-1.)
  155 GRD=(4.*GND-1.)/4.              
      GO TO 159         
  156 GRD=GND           
      GO TO 159         
  157 GRD=(4.*GND+1.)/4.              
      GO TO 159         
  158 GRD=(2.*GND+1.)/2.              
      GO TO 159         
  159 APDR=GRD*ALAM     
      IF(APDR.GT.TWHT) GO TO 162      
      SI=ASIN (APDR/TWHT)             
      IF(SI.GT.PI2) SI=PI2            
      ASSIGN 152 TO KR                
      GO TO 66          
  152 L1=L1+1           
      XS(L1)=SI         
      XD(L1)=DR         
      XR(L1)=D          
      SI=SQRT (APDR/(2.*DLST))        
      ASSIGN 153 TO KR                
      GO TO 66          
  153 L2=L2+1           
      YS(L2)=SI         
      YD(L2)=DR         
      YR(L2)=D          
      GO TO 151         
  162 L3=0              
      LK=0              
   67 LK=LK+1           
      IF(LK.GT.24) GO TO 913          
      SI=SID(LK)*RAD    
      ASSIGN 124 TO KR                
      GO TO 66          
  124 L3=L3+1           
      ZS(L3)=SI         
      ZD(L3)=DR         
      ZR(L3)=D          
      GO TO 67          
  913 CONTINUE          
      SI=PI2            
      L3=L3+1           
      ZS(L3)=SI         
      ZD(L3)=TWHT       
      ZR(L3)=0.         
*GH    CALL TABLE(DUM)                 
      DUM=TABLE(DUM)

C     ----USING TABLE TO OBTAIN STRATIGIC DISTANCE POINTS-------- 
          
      DA6=0.            
      DZR=0.            
      LR=0              
      LA=0              
   70 LA=LA+1           
      IF(LA.GT.LE) GO TO 25           
      IF(LA.LT.4) GO TO 88            
      LB=13-LA          
      GRD=FLOAT (LB)    
      DR=ALAM/GRD       
      LD=LD+1           
      IF(DR.GT.TWHT) GO TO 25         
   86 CONTINUE          
      D=DNTER(DR)       
      IF(LB.EQ.6) DA6=D               
      IF(D.GT.DML) GO TO 70           
      LR=LR+1           
      D1(LR)=D          
      GO TO 70          
   25 CONTINUE          
      IF(DLT.GE.DOZ.OR.DOZ.GE.DML) GO TO 81         
      IF(DOZ.LT.DA6.AND.DA6.LT.DML) GO TO 83        
      DZR=DOZ           
   84 IF(ILB.LE.0) GO TO 163          
      LA=0              
  172 LA=LA+1           
      IF(LA.GT.10) GO TO 163          
      GND=FLOAT (LA)    
      LG=0              
  173 LG=LG+1           
      IF(LG.GT.4) GO TO 172           
      GO TO (165,166,167,168), LG     
****  850606            
*CDC  ADD CDC GO TO ERROR ROUTINE     
*     CALL GOTOER
****  940914
*MS   FOR MICROSOFT USE DOMAIN ERROR
      GRD=SQRT(-1.)
  165 GRD=(4.*GND-1.)/4.              
      GO TO 169         
  166 GRD=GND           
      GO TO 169         
  167 GRD=(4.*GND+1.)/4.              
      GO TO 169         
  168 GRD=(2.*GND+1.)/2.              
      GO TO 169         
  169 DR=GRD*ALAM       
      IF(DR.GT.TWHT) GO TO 163        
      D=DNTER(DR)       
      IF(D.GT.DML) GO TO 172          
      LR=LR+1           
      D1(LR)=D          
      GO TO 173         
  163 CONTINUE          
      DXZ=DZR           
      D=DZR             
   13 SI=SNTER(D)       
      ASSIGN 12 TO KR                 
      GO TO 66          
   12 IF(D.GE.DZR) GO TO 14           
      D=DXZ+.001        
      IF (D .GE. DML) GO TO 10        
      DXZ=D             
      GO TO 13          
   10 D=D-.001          
   14 DZR=D             
      WRITE(IOT,750)DZR               
  750 FORMAT(10X,'D ZERO=',F8.3)      
      IF(LR)164,164,154               
  154 D=D1(LR)          
      SILIM=SNTER(D)    
      DO 11 LA=1,LR     
      LV=LR+1-LA        
   11 D3(LA)=D1(LV)     
      D2(1)=DZR         
      CALL TMESH(D2,1,D3,LR,D1,L5)    
  160 LR=0              
          
C     ------------------- LOOPING --------------------------------------        
      IF(JC.GT.0) GO TO 904           
      SPD=1.            
      NSP=0             
  800 NSP=NSP+1         
      IF(NSP.GT.5) GO TO 107          
      MZS=NTM(NSP)      
      IF(MZS.LE.0) GO TO 107          
      MXS=0             
  801 MXS=MXS+1         
      IF(MXS.GT.MZS) GO TO 914        
      D=SPD*CMK(IK)     
      IF(D.GT.DML) GO TO 107          
      LR=LR+1           
      D3(LR)=D          
      SPD=SPD+XCON(NSP)               
      GO TO 801         
  914 CONTINUE          
      SPD=SPD-XCON(NSP)               
      NPP=NSP+1         
      IF(NPP.GT.5) GO TO 107          
      IF(XCON(NPP).EQ.0.) GO TO 107   
      IF(NPP.EQ.0) GO TO 107          
      IXD=INT (SPD/XCON(NPP))         
      SPD=(XCON(NPP)*FLOAT (IXD))+XCON(NPP)         
      GO TO 800         
  107 CONTINUE          
      CALL TMESH(D1,L5,D3,LR,D2,LX)   
      IF(NOC.LE.0) GO TO 75           
          
C     ---------CALCULATION OF COUNTERPOISE STRATIGIC POINTS------ 
      LR=0              
      LK=0              
  600 LK=LK+1           
      IF(LK.GT.13) GO TO 29           
      IF(LK.LT.9) GO TO 601           
      FLK=LK-8          
      LG=0              
  603 LG=LG+1           
      IF(LG.GT.4) GO TO 600           
      FLG=LG            
      GND=((4.*FLK)+FLG)/4.           
  602 APDR=GND*ALAM     
      IF(APDR.GT.TWHC) GO TO 29       
      SI=ASIN (APDR/TWHC)             
      ICPT=1            
      ASSIGN 40 TO KR                 
      GO TO 66          
   40 CONTINUE          
      IF(D.GT.DML) GO TO 604          
      LR=LR+1           
      D3(LR)=D          
  604 IF(LK.LT.9) GO TO 600           
      GO TO 603         
   29 CONTINUE          
      CLIM=D3(LR)       
      CCIM=D3(1)        
      DO 69 I=1,LR      
      LV=LR+1-I         
   69 D1(I)=D3(LV)      
      CALL TMESH(D1,LR,D2,LX,D3,LK)   
  134 LV=0              
  129 LV=LV+1           
      IF(LV.GT.LK) GO TO 111          
      ICPT=0            
      SI=SNTER(D3(LV))                
      ASSIGN 28 TO KR                 
          
C     ----------------RAY OPTICS GEOMETRY------------------------ 
          
   66 CSSI=COS (SI)     
      SNSI=SIN (SI)     
      AK0=EFN/A0        
      ZE=(1./AK0)-1.    
      AKE=1./(1.+(ZE*CSSI))           
      AEFT=A0*AKE       
      DHE=EAC*(AKE-1.)/(AK0-1.)       
      DHE1=EAL*(AKE-1.)/(AK0-1.)      
      HL1=H1-DHE1       
      HTX(1)=HL1-HRP    
      HL=H2-DHE         
      HTX(2)=HL-HRP     
      IF(ICPT.GT.0) GO TO 77          
      A=AEFT            
   78 CONTINUE          
      LC=0              
   62 LC=LC+1           
      IF(LC.GT.2) GO TO 915           
      Z(LC)=A+HTX(LC)                 
      TEA(LC)=ACOS (A*CSSI/Z(LC))-SI                
      DA(LC)=Z(LC)*SIN (TEA(LC))      
      IF(SI.GT.1.56) GO TO 63         
      HPR(LC)=DA(LC)*TAN (SI)         
      GO TO 62          
  915 CONTINUE          
      DTX=ABS (Z(1)-Z(2))             
      TH=TEA(1)+TEA(2)                
      D=AEFT*TH         
      IF(D.LT.0.) D=0.                
      IF(SI.GT.CPI2) GO TO 64         
      AFA=ATAN ((HPR(2)-HPR(1))/(DA(1)+DA(2)))      
      R0=(DA(1)+DA(2))/COS (AFA)      
      R12=(DA(1)+DA(2))/CSSI          
      IF(R0.LT.DTX) R0=DTX            
   68 CA=AFA-TEA(1)     
      ABA=-AFA-TEA(2)                 
      DR=4.*HPR(1)*HPR(2)/(R0+R12)    
      DNM=D*CKN(IK)     
      GO TO KR,(65,28,123,124,40,152,153,12)        
          
C     ----------------------------------------------------------- 
          
   28 IF(D.LT.0.01) GO TO 129         
      IF(D.GT.DML) GO TO 111          
      THFS=(TH*AEFT)/A0               
      ZFS1=6370.+H1-HRP               
      ZFS2=6370.+H2-HRP               
      TRM1=(ZFS2-ZFS1)                
      TRM1=TRM1*TRM1    
      TRM2=SIN(THFS*.5)               
      TRM2=TRM2*TRM2    
      RFS=SQRT(TRM1+(4.*ZFS1*ZFS2*TRM2))            
      DTX=ABS(ZFS2-ZFS1)              
      IF(RFS.LT.DTX) RFS=DTX          
      ALFS=AFP+20.*ALOG10(RFS)        
      BA=CA+TWT1*(AEFT-A0)            
      ACA=ABA+TWT2*(AEFT-A0)          
      PFS=PIRP-ALFS     
      IF(JS.GT.0)GO TO 842            
  843 CALL GANE(ACA,IAA,H2PBW,T2T,GADV,GADH,GAD,IPL)              
      IF(JT.GT.0)GO TO 844            
  845 CALL GANE(BA,IFA,HLPBW,T1T,GFVD,GFDH,GFD,IPL)               
      GOD=GAD*GFD       
      IF (IPL .GE. 3) GOD=0.5*(GFVD*GADV+GFDH*GADH)               
      GPD=20.*ALOG10(GOD)             
      IF(DH.LE.0.) GO TO 42           
      DHD=DH*(1.-(0.8*EXP (-0.02*D)))*1000.         
   44 CALL SORB(Z(1),Z(2),A,R0,CA,RE)               
      AA=GAO*RE(1)+GAW*RE(2)          
      IF(ILB.GT.0.AND.SI.LE.SILIM) GO TO 35         
      IF( DR.GE.ALA2) GO TO 34        
      IF( DR.LE.DTR0) GO TO 26        
      FDR=(1.1-(0.9*COS (PKL*( DR-DTR0))))*.5       
   43 CONTINUE          
      TA=-(TEA(1)+SI)                 
      GA=TA+TWT1*(AEFT-A0)            
      CALL GANE(GA,IFA,HLPBW,T1T,GFV,GFH,GFG,IPL)   
      AHA=-(TEA(2)+SI)                
      AGA=AHA+TWT2*(AEFT-A0)          
      CALL GANE(AGA,IAA,H2PBW,T2T,GAV,GAH,GAG,IPL)  
      GOG=GAG*GFG       
      GV=GAV*GFV        
      GH=GAH*GFH        
      CALL RECPX(SI,F,KSC,IPL,0,DHD,GH,GV,RG,PIC,RDG)             
      REC=0.0           
      IF(TAN (SI).LT.(1./(TWPILA*A))) GO TO 111     
      IF(TAN (SI).GE.0.1) GO TO 54    
      RR12=DA(1)*DA(2)/(CSSI*CSSI)    
      TD1=(2.*RR12/R12)*(1./A)        
      TD4=1.+(SNSI*SNSI)              
      TD5=TD1*TD1       
      TD2=1.+(TD1*TD4/SNSI)+TD5       
      DV=SQRT (1./TD2)                
   55 TD3=R0/R12        
      IF(TD3.GT.1.) TD3=1.            
      REG=RG*DV*TD3     
      RLG=REG           
      IF(NOC.LE.0) GO TO 500          
          
C     ---------CALCULATION OF COUNTERPOISE CONTRIBUTION---------- 
      TEG=ABTC-ABS (SI+TEA(1))        
      TEG=ABS (TEG)     
      VFGD=2.*SIN (TEG*.5)*SQVT       
      IF(ABS (GA).LT.ABTC) VFGD=-VFGD               
      CALL FRNEL(VFGD,FPGD,PHIG)      
      REG=REG*FPGD      
      RDG=RDG*FPGD      
      TRM3=PHIG+(PI2*VFGD*VFGD)       
      IF(D.LT.CLIM.OR.D.GT.CCIM) GO TO 146          
      SIC=CA            
      SIT1=-SIC         
      TEC=ABS (BTC-CA)                
      DARC=2.*HFC*SIN (CA)            
      CALL GANE(SIT1,IFA,HLPBW,T1T,GCV,GCH,GFC,IPL)               
      GOC=GFC*GAG       
      GFCV=GCV*GAV      
      GFCH=GCH*GAH      
      VFCP=2.*SIN (TEC*.5)*SQVT       
      IF(ABS (CA).GT.ABTC) VFCP=-VFCP               
      CALL FRNEL(VFCP,FPCP,PHIC)      
      EPK=EPC           
      SGM=SGC           
      CALL RECPX(SIC,F,ICC,IPL,1,DHD,GFCH,GFCV,RLC,PICC,RDC)      
      EPK=EPG           
      SGM=SGG           
      REC=RLC*FPCP      
      EXPC=(TWPILA*DARC)+PICC+(PHIC+(PI2 *VFCP*VFCP))             
      ATRM=REC*COS (EXPC)             
      BTRM=-REC*SIN (EXPC)            
      AT1=CMPLX(ATRM,BTRM)            
  147 CONTINUE          
          
C     --------------CALCULATION OF LOBING CONTRIBUTION----------- 
      IF(SI.GT.SILIM) GO TO 135       
      EXPG=(TWPILA*DR)+PIC+TRM3       
      ATRM=REG*COS (EXPG)             
      BTRM=-REG*SIN (EXPG)            
          
C     -------------------SUMMATION OF TERMS---------------------- 
  136 AT2=CMPLX(ATRM,BTRM)            
      IF (SI .LE. SILIM .AND. ILB .LE. 0) GO TO 137               
  138 WRL=CABS(GOD+AT1+AT2)           
      WR=WRL*WRL+(.0001*GOD)          
  139 PR=10.*ALOG10(WR)               
      IF(D.LE.DZR) GO TO 148          
      IF(DZR.LT.0.) GO TO 145         
      PL=FNA(D,DML,DZR,PRH,PZ)        
      WL=10.**(.1*PL)                 
  149 CONTINUE          
          
C     ----------------LONG-TERM POWER FADING--------------------- 
          
      PL=PL-GPD         
      IF(D.LE.0.) GO TO 38            
      IF(D-DSL1) 301,301,302          
  301 DEE=(130.*D)/DSL1               
      GO TO 303         
  302 DEE=130.+D-DSL1                 
      GO TO 303         
  303 CALL VARDT(DEE)                 
      IF(CA.LE.0.) GO TO 32           
      IF(CA.GE.1.) GO TO 33           
      FTH=.5-BSPI*(ATAN (20.*ALOG10(32.*CA)))       
      IF(FTH.LE.0.0) GO TO 33         
   52 AL10=PL+(AD(13)*FTH)            
      AY=AL10-ALIM      
      IF(AY.LT.0.) AY=0.              
      IF(AY.GT.10.) AY=10.            
   53 IF(ILB.GT.0.AND.SI.LE.SILIM) GO TO 22         
      DO 31 K=1,35      
      VD(K)=AD(K)*FTH-AY              
      BD(K)=PL+VD(K)    
   31 CONTINUE          
      DO 50 K=1,12      
      ALLM=-ALM(K)      
      IF(BD(K).GT.ALLM) BD(K)=ALLM    
   50 CONTINUE          
   24 CONTINUE          
          
C     -------------VALUES PUT INTO PLOTTING ARRAY---------------- 
      NCT=NCT+1         
      THD=(D/6370.)*DEG               
      THO=THD           
      IF(JC.GT.0) GO TO 101           
      DO 104 NN=1,8     
  104 BX(NCT,NN)=DNM    
  102 CONTINUE          
      IF(KK.GT.1) GO TO 20            
      IF(JO.GT.0) GO TO 92            
      IF(IZ.GT.0) GO TO 93            
   23 PGS=PFS+GPD       
      BY(NCT,1)=PGS     
      BY(NCT,2)=PGS+BD(18)-AA         
      BY(NCT,3)=PGS+BD(12)-AA         
      BY(NCT,4)=PGS+BD(24)-AA         
      BY(NCT,5)=PGS+BD(23)-AA         
      BY(NCT,6)=PGS+BD(26)-AA         
      BY(NCT,7)=PGS+BD(29)-AA         
      BY(NCT,8)=PGS+BD(32)-AA         
      PFY(NCT,1)=PGS+BD(4)-AA         
      PFY(NCT,2)=PGS+BD(7)-AA         
      PFY(NCT,3)=PGS+BD(10)-AA        
      PFY(NCT,4)=PGS+BD(13)-AA        
      WRITE(IOT,760)DNM,(BY(NCT,LZ),LZ=1,8),(PFY(NCT,MW),MW=1,4),THD,           
     XPL,AA,AY,BK,DEE                 
      CALL PAGE(1)      
      IF (JA .GT. 0) GO TO 400        
  404 CONTINUE          
      GO TO 129         
  111 CONTINUE          
      NU(1)=NCT         
      RETURN            
  400 IF (DNM .GE. RQD .AND. JB .GT. 1) GO TO 401   
      DQM=DNM           
      JB=2              
      DO 402 I=1,8      
  402 DBY(I)=BY(NCT,I)                
      DO 403 I=1,4      
      JE=I+8            
  403 DBY(JE)=PFY(NCT,I)              
      GO TO 404         
  401 RAT=(RQD-DQM)/(DNM-DQM)         
      DO 405 I=1,8      
  405 SELT(I)=RAT*(BY(NCT,I)-DBY(I))+DBY(I)         
      DO 406 I=1,4      
      JE=I+8            
  406 SELT(JE)=RAT*(PFY(NCT,I)-DBY(JE))+DBY(JE)     
      WRITE(IOT,760) RQD,SELT         
      JG=1              
      RETURN            
C     -------------------RETURN TO MAIN PROGRAM------------------ 
          
   15 FAY=1.            
      GO TO 17          
   16 FAY=0.1           
      GO TO 17          
          
C     ------------------TROPOSHERIC MULTIPATH-------------------- 
   20 DO 30 I=1,35      
      PQA(I)=P(I)       
      QA(I)=BD(I)-PL    
   30 CONTINUE          
      IF(AY.LE.0.) GO TO 15           
      IF(AY.GE.9.) GO TO 16           
      FAY=(1.1+(0.9*COS((AY/9.)*PI)))/2.            
   17 CONTINUE          
      RSP=REG*FDR*FAY                 
      IF(RE(2).LE.0.) GO TO 45        
      RK=-10.*ALOG10(ASPC*(RE(2)**3))               
      ACK=FDASP(RK)     
      WA=10.**(.1*ACK)                
   46 RST=((RSP*RSP+RDG*RDG)/(GOD*GOD))+WA          
      IF(RST.LE.0.) GO TO 37          
      BK  =+10.*ALOG10(RST)           
      IF(BK.LT.-40.) BK=-40.          
   47 CALL YIKK(BK,PQK,QK)            
      RDHK=BK           
      CALL CNLUT(QA,QK,PQA,35,+1.,0.,PQC,QC)        
      IF(JO.GT.0) GO TO 90            
      IF(IZ.LE.0) GO TO 91            
      DO 98 I=1,35      
   98 QA(I)=QC(I)       
      GO TO 94          
          
   37 BK=-40.           
      GO TO 47          
          
C     -------------------LOBING MODE----------------------------- 
   22 AY=0.             
      TLIM=+20.*ALOG10(GOD+REG+REC)-GPD+(AD(18)*FTH)              
      BLIM=-80.         
      DO 36 K=1,35      
      VD(K)=AD(K)*FTH                 
      BD(K)=PL+VD(K)-AA               
      IF(BD(K).GT.TLIM) BD(K)=TLIM    
      IF(BD(K).LT.BLIM) BD(K)=BLIM    
      BD(K)=BD(K)+AA    
   36 CONTINUE          
      GO TO 24          
          
C     ----------------- MISCELLANEOUS STATEMENTS -----------------------        
   26 FDR=0.1           
      GO TO 43          
   32 FTH=1.0           
      GO TO 52          
   33 FTH=0.0           
      AY=0.0            
      GO TO 53          
   34 FDR=1.            
      GO TO 43          
   35 FDR=0.            
      GO TO 43          
   38 DEE=0.            
      GO TO 303         
   42 DHD=0.0           
      GO TO 44          
   45 WA=.0001          
      GO TO 46          
   54 DV=1.             
      GO TO 55          
   63 HPR(LC)=HTX(LC)                 
      GO TO 62          
   64 AFA=SI            
      R0=AMAX1(DTX,D)                 
      R12=HTX(1)+HTX(2)               
      GO TO 68          
   75 DO 74 LK=1,LX     
   74 D3(LK)=D2(LK)     
      LK=LX             
      LR=LX             
      GO TO 134         
   77 HTX(1)=HFC        
      HTX(2)=HTX(2)-HDI               
      A=AEFT+HDI        
      ICPT=0            
      GO TO 78          
   81 IF(DLT.GT.DA6.OR.DA6.GT.DML) GO TO 82         
   83 DZR=DA6           
      GO TO 84          
   82 DZR=DLT           
      GO TO 84          
   88 GRD=SPGRD(LA)     
      DR=ALAM*GRD       
      LD=LD+1           
      GO TO 86          
          
C     ---------------- COMBINING DISTRIBUTIONS -------------------------        
   92 DO 95 I=1,35      
      QC(I)=BD(I)-PL    
   95 PQA(I)=P(I)       
   90 CALL SCNTL        
      CALL CNLUT(QC,QS,PQA,35,+1.,0.,PQC,QA)        
      IF(IZ.GT.0) GO TO 94            
      DO 99 I=1,35      
   99 BD(I)=QA(I)+PL    
      GO TO 23          
   93 DO 96 I=1,35      
      QA(I)=BD(I)-PL    
   96 PQA(I)=P(I)       
   94 DX=R0             
      HA1=Z(1)          
      HA2=Z(2)          
      AN=CA             
      EFA=A             
      CALL RAIN         
      IF(IZ.GT.6) GO TO 199           
      DO 87 I=1,35      
   87 QR(I)=-RS(I)      
      CALL CNLUT(QA,QR,PQA,35,+1.,0.,PQC,QC)        
   91 DO 97 I=1,35      
   97 BD(I)=QC(I)+PL    
      GO TO 23          
  101 DO 103 NN=1,8     
  103 BX(NCT,NN)=THD    
      GO TO 102         
          
C     ----------------- MISCELLANEOUS STATEMENTS -----------------------        
  120 GRD=SPGRD(LK)     
      APDR=ALAM*GRD     
      GO TO 121         
  122 SI=0.             
      DR=0.             
      D=DLST+DLSR       
      GO TO 123         
  135 ATRM=0.           
      BTRM=0.           
      GO TO 136         
  137 WRL=CABS(GOD+AT2)               
      IF (WRL .LE. GOD)GO TO 138      
      WRL=CABS(GOD+AT1)               
      WR=WRL*WRL+(.0001*GOD)          
      GO TO 139         
  164 D1(1)=DZR         
      L5=1              
      SILIM=0.          
      GO TO 160         
  500 TRM3=0.0          
  146 ATRM=0.           
      BTRM=0.           
      AT1=CMPLX(ATRM,BTRM)            
      RLC=0.0           
      GO TO 147         
  145 DZR=0.            
  148 PL=PR             
      PZ=PR             
      WL=WR             
      GO TO 149         
  199 PL=PL-RS(18)      
      DO 198 I=1,35     
  198 BD(I)=QA(I)+PL    
      GO TO 23          
  502 BTC=0.            
      SQVT=0.           
      HDI=HTE           
      GO TO 503         
  601 GND=GLD(LK)       
      GO TO 602         
  842 T2T=ACA*DEG       
      GO TO 843         
  844 T1T=BA*DEG        
      GO TO 845         
          
C     -------------- LOOPING FOR CENTRAL ANGLE -------------------------        
  904 SDD=.1            
      NSP=0             
  900 NSP=NSP+1         
      IF(NSP.GT.5) GO TO 107          
      MZS=NDM(NSP)      
      IF(MZS.LE.0) GO TO 107          
      MXS=0             
  901 MXS=MXS+1         
      IF(MXS.GT.MZS) GO TO 950        
      D=SDD*6370.*RAD                 
      SPD=D*CKN(IK)     
      IF (D .GT. DML) GO TO 107       
      LR=LR+1           
      D3(LR)=D          
      SDD=SDD+XDON(NSP)               
      GO TO 901         
  950 SDD=SDD-XDON(NSP)               
      NPP=NSP+1         
      IF (NPP .GT. 5) GO TO 107       
      IF(XDON(NPP).EQ.0.) GO TO 107   
      IF(NPP.EQ.0) GO TO 107          
      IXD=INT (SDD/XDON(NPP))         
      SDD=(XDON(NPP)*FLOAT (IXD))+XDON(NPP)         
      GO TO 900         
C     ------------------------------------------------------------------        
          
      END
