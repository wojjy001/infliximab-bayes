$PROBLEM - INFLIXIMAB REFERENCE MODEL
;Reference: 

$INPUT ID	TIME PPVCL PPVV1 PPVQ PPVV2 ALB WT ADA ERRPRO AMT CMT EVID RATE DV

$DATA nonmem_simulation_input.csv IGNORE = C

$SUBROUTINE ADVAN6 TOL=5

$MODEL COMP=(CENT,DEFDOS,DEFOBS), COMP=(PERI)

$PK
	;Population parameters
	POPCL = THETA(1)
	POPV1 = THETA(2)
	POPQ = THETA(3)
	POPV2 = THETA(4)

	;Covariate effects
	WT_CL = THETA(5)	;Effect of weight on clearance
	ALB_CL = THETA(9)	;Effect of albumin on clearance
	ADA_CL = THETA(10) ;Effect of anti-drug antibodies on clearance
	WT_V1 = THETA(6)	;Effect of weight on V1
	WT_Q = THETA(7)	;Effect of weight on Q
	WT_V2 = THETA(8)	;Effect of wegiht on V2

	D1 = 0.08333333	;Infusion into the central compartment, specified as duration = 2 hours (converted to days)
	
	;Anti-drug antibodies
	ADACOV = 1	;No anti-drug antibodies
	IF (ADA.EQ.1) ADACOV = 1+ADA_CL	;Anti-drug antibodies present
	
	;Individual parameter values
	CL = POPCL*((WT/70)**WT_CL)*((ALB/4)**ALB_CL)*ADACOV*exp(PPVCL)
	V1 = POPV1*((WT/70)**WT_V1)*exp(PPVV1)
	Q = POPQ*((WT/70)**WT_Q)*exp(PPVQ)
	V2 = POPV2*((WT/70)**WT_V2)*exp(PPVV2)
	
	S1 = V1	;Scale factor for V
	
	;Dummy ETA values (because these have already been simulated and are found in the input dataset)
	DUM1 = ETA(1)
	DUM2 = ETA(2)
	DUM3 = ETA(3)
	DUM4 = ETA(4)

$DES
	DADT(1) = -Q/V1*A(1) +Q/V2*A(2) -CL/V1*A(1)	;Central compartment
	DADT(2) = Q/V1*A(1) -Q/V2*A(2) ;Peripheral compartment

$THETA
	(0,0.294,)	;TVCL
	(0.1,3.33,)	;TVV1
	(0,0.0719,)	;TVQ
	(0,1.14,)	;TVV2

	0.614	;WTONCL	#Power model - referenced to 70 kg
	0.691	;WTONV1
	1.10	;WTONQ
	0.59	;WTONV2

	-1.17	;ALBONCL #Power model - referenced to 4 g/dL
	0.257	;ADAONCL #Dichotomous proportional model - offset is 1+0.257 when ADA = 1

$OMEGA
	0.106929	;CLBSV  0.327^2
	0.0225	;V1BSV 0.150^2
	1.21	;QBSV  1.10^2
	0.638401	;V2BSV  0.799^2

$SIGMA
	0.175561	;RUVPROP 0.419^2

$ERROR
	A1 = A(1)
	A2 = A(2)
	IPRE = F
	Y = F*(1+ERRPRO)
	DUM5 = ERR(1)	;Dummy error value (already simulated as ERRPRO)

$SIMULATION (1234567) ONLYSIM SUBPROBLEMS = 1

;$ESTIMATION METHOD = 1 INTERACTION MAXEVALS = 0 POSTHOC NOABORT NSIG = 3 SIGL = 9

;$COVARIANCE UNCONDITIONAL SIGL = 12 PRINT = E

$TABLE ID TIME AMT RATE A1 A2 CL V1 Q V2 IPRE WT ADA ALB PPVCL PPVV1 PPVQ PPVV2 NOPRINT ONEHEADER FILE = *.fit
