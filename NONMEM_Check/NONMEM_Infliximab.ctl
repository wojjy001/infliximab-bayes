$PROBLEM - INFLIXIMAB REFERENCE MODEL
;UNITS ARE MG, L (UG/ML) AND DAYS
;MOULD, 2012

;SET TO READ CURRENTPATIENTDF2 OUTPUT AS A CSV.FILE VIA EXTRACTNMDATA.R
$INPUT  ID EVENTNUM EVENTID TIME AMT RATE SS II DV WEIGHT SEX ADA ALB CLIOFF V1IOFF QIOFF V2IOFF

$DATA ..\John_Standard_Loading_currentpatientdf2.csv IGNORE=C

$SUBROUTINE  ADVAN6 TOL=5

$MODEL  COMP=(CENTRAL,DEFDOS,DEFOBS), COMP=(PERIPH)

$PK
   ;POPULATION CLEARANCE - FROM MOULD PAPER
   CLPOP=THETA(1)
   V1POP=THETA(2)
   QPOP =THETA(3)
   V2POP=THETA(4)

   WTONCL=THETA(5)
   WTONV1=THETA(6)
   WTONQ =THETA(7)
   WTONV2=THETA(8)

   ALBONCL=THETA(9)
   ADAONCL=THETA(10)

   ;POPULATION GROUP CLEARANCE - GIVES PRED
   CLGRP=CLPOP*(WEIGHT/70)**WTONCL*(ALB/4)**(ALBONCL)*(1+ADAONCL*ADA)
   V1GRP=V1POP*(WEIGHT/70)**WTONV1
   QGRP= QPOP *(WEIGHT/70)**WTONQ
   V2GRP=V2POP*(WEIGHT/70)**WTONV2


   ;INDIVIDUAL GROUP CLEARANCE - AS FITTED IN BAYES ESTIMATION STEP GIVEN PARAMETER PRIORS - GIVES IPRED
   CL  = CLGRP*EXP(ETA(1))
   V1  = V1GRP*EXP(ETA(2))
   Q   = QGRP *EXP(ETA(3))
   V2  = V2GRP*EXP(ETA(4))

   S1 = V1                        ; SCALE FACTOR FOR V

$DES
    DADT(1) =  -Q/V1*A(1) +Q/V2*A(2) -CL/V1*A(1)           ;CENTRAL COMPARTMENT
    DADT(2) =   Q/V1*A(1) -Q/V2*A(2)                       ;PERIPHERAL COMPARTMENT

$THETA
       (0,0.294,)    ;POPCL
       (0.1,3.33,)   ;POPV1
       (0,0.0719,)   ;POPQ
       (0,1.14,)     ;POPV2

       0.614   ;WTONCL    #POWER MODEL - REFERENCE WEIGHT 70KG
       0.691   ;WTONV1
       1.10    ;WTONQ
       0.59    ;WTONV2

       -1.17   ;ALBONCL    #POWER MODEL - REFERENCE 4
        0.257  ;ADAONCL    #OFFSET IS 1+0.257 WHEN ADA = 1

;set these to zero so that IPRED is GROUP PREDICTION
;$OMEGA
;0 FIX   ;CLBSV  0.327^2
;0 FIX     ;V1BSV 0.150^2
;0 FIX       ;QBSV  1.10^2
;0 FIX  ;V2BSV  0.799^2

;$SIGMA
;0 FIX  ;RUVPROP 0.419^2

$OMEGA
0.106929   ;CLBSV  0.327^2
0.0225     ;V1BSV 0.150^2
1.21       ;QBSV  1.10^2
0.638401   ;V2BSV  0.799^2

$SIGMA
0.175561  ;RUVPROP 0.419^2


$ERROR
   A1 = A(1)
   A2 = A(2)
   DV2 = A(1)/V1
   IPRED=F
   CL2 = CL
   V12 = V1
   Y = F*(1+ERR(1))


$SIMULATION (1234567)  ONLYSIM SUBPROBLEMS=5000

;$ESTIMATION METHOD=1 INTERACTION MAXEVALS=0 POSTHOC NOABORT NSIG=3 SIGL=9

;$COVARIANCE  UNCONDITIONAL SIGL=12 PRINT=E

$TABLE ID TIME AMT RATE A1 A2
 CLPOP V1POP QPOP V2POP
 CLGRP V1GRP QGRP V2GRP
 CL V1 Q V2
 CL2 V12 DV2 IPRED WEIGHT SEX ADA ALB
NOPRINT ONEHEADER FILE=*.fit
