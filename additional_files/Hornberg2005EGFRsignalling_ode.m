%%% Version: 19 Jan 2023
%%%
%%% dX  = <MODELNAME>_ode(t,X,par,model)
%%%
%%% This function defines the system of model ODEs 
%%% 
%%% Input : t           time
%%%         X           state vector
%%%         par         parameter vector
%%%         model       model structure containing the index structure of
%%%                     the model
%%%                   
%%% Output : dX         right hand side of ODEs
%%%
%%% References:
%%%    + Hornberg, J et al.
%%%      'Control of MAPK signaling: from complexity to what really matters'
%%%      Oncogene, Volume 24, 2005
%%%    + Schoeberl, Birgit; Eichler-Johnson, Claudia; Gilles, Ernst Dieter;
%%%      Mueller, Gertrazd
%%%      'Computational modeling of the dynamics of the MAP 
%%%      kinase casade activated by surface and internalized 
%%%      EGF receptors'
%%%      Nature Biotechnology Volume 20, 2002
%%%    + Corresponding Model from BioModels Database:
%%%      http://www.ebi.ac.uk/biomodels-main/BIOMD0000000019
%%%
%%% 
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function dX = Hornberg2005EGFRsignalling_ode(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize rhs vector
dX = 0*X;

%%% -----------------------------------------------------------------------
%%% specify system of ODEs 

%%% Equ 1-13: epidermal growth factor (EGF) and its receptor (EGFR) receptor dynamic with internalisation
%%% 
r1 = par(I.k1) * X(I.EGF) * X(I.EGFR) - par(I.kd1) * X(I.EGF_EGFR);
r2 = par(I.k2) * X(I.EGF_EGFR) * X(I.EGF_EGFR) - par(I.kd2) * X(I.EGF_EGFR_2);
r3 = par(I.k3) * X(I.EGF_EGFR_2) - par(I.kd3) * X(I.EGF_EGFRa_2);
r6 = par(I.k6) * X(I.EGFR) - par(I.kd6) * X(I.EGFRi);
r7 = par(I.k7) * X(I.EGF_EGFRa_2) - par(I.kd7) * X(I.EGF_EGFRia_2);
r8 = par(I.k8) * X(I.EGF_EGFRa_2) * X(I.GAP) - par(I.kd8) * X(I.EGF_EGFRa_2_GAP);
r10 = par(I.k10b) * X(I.EGFRi) * X(I.EGFi) - par(I.kd10) * X(I.EGF_EGFRi);

r11 = par(I.k11) * X(I.EGF_EGFRi) * X(I.EGF_EGFRi) - par(I.kd11) * X(I.EGF_EGFRi_2);
r12 = par(I.k12) * X(I.EGF_EGFRi_2) - par(I.kd12) * X(I.EGF_EGFRia_2);
r13 = par(I.k13) - par(I.kd13) * X(I.EGFR);
r14 = par(I.k14) * X(I.EGF_EGFRia_2) * X(I.GAP) - par(I.kd14) * X(I.EGF_EGFRia_2_GAP);
r60 = par(I.k60) * X(I.EGFRi) - par(I.kd60) * X(I.EGFRideg); 
r61 = par(I.k61) * X(I.EGFi) - par(I.kd61) * X(I.EGFideg);
r62 = par(I.k62) * X(I.EGF_EGFRia_2)- par(I.kd62) * X(I.EGF_EGFRia_2_deg);


dX(I.EGF) = 0;
dX(I.EGFR) = r13 - r1 - r6 ;
dX(I.EGF_EGFR) = r1 + -2*r2;
dX(I.EGF_EGFR_2) = r2 - r3;	
dX(I.EGF_EGFRa_2) = r3 - r7 - r8;

dX(I.EGFRi) = r6 - r10 - r60;
dX(I.EGFRideg) =  r60;
dX(I.EGFi) = -r10 - r61;
dX(I.EGF_EGFRi) = r10 - r11;
dX(I.EGF_EGFRi_2) =  r11 - r12;
dX(I.EGF_EGFRia_2) = r7 + r12 - r14 - r62;
dX(I.EGFideg) = r61;

%%% -----------------------------------------------------------------------
%%% Equ 14: EGF:EGFRa:2:GAP start of both pathways
r16 = par(I.k16) * X(I.Grb2) * X(I.EGF_EGFRa_2_GAP) - par(I.kd16) * X(I.EGF_EGFRa_2_GAP_Grb2);
r22 = par(I.k22) * X(I.Shc) * X(I.EGF_EGFRa_2_GAP) - par(I.kd22) * X(I.EGF_EGFRa_2_GAP_Shc);
r32 = par(I.kd32) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.k32) * X(I.EGF_EGFRa_2_GAP) * X(I.Shca_Grb2_Sos);
r34 = par(I.kd34) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k34) * X(I.EGF_EGFRa_2_GAP) * X(I.Grb2_Sos);
r37 = par(I.kd37) * X(I.EGF_EGFRa_2_GAP_Shca) - par(I.k37) * X(I.EGF_EGFRa_2_GAP) * X(I.Shca);
r39 = par(I.kd39) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) - par(I.k39) * X(I.EGF_EGFRa_2_GAP) * X(I.Shca_Grb2);
r102 = par(I.k102) * X(I.EGF_EGFRa_2_GAP) - par(I.kd102) * X(I.EGF_EGFRia_2_GAP);

dX(I.EGF_EGFRa_2_GAP) = r8 - r16 - r22 + r32 + r34 + r37 + r39 - r102;

%%% -----------------------------------------------------------------------
%%% Equ 15-22: for Pathway without Shc
r4 = par(I.k4) * X(I.EGF_EGFRa_2_GAP_Grb2) * X(I.Prot) - par(I.kd4) * X(I.EGF_EGFRa_2_GAP_Grb2_Prot);
r5 = par(I.kd5) * X(I.EGF_EGFRa_2_GAP_Grb2_Prot) - par(I.k5) * X(I.Proti) * X(I.EGF_EGFRia_2_GAP_Grb2);
r9 = par(I.k9) * X(I.EGF_EGFRa_2_GAP_Grb2) - par(I.kd9) * X(I.EGF_EGFRia_2_GAP_Grb2); 
r17 = par(I.k17) * X(I.Sos) * X(I.EGF_EGFRa_2_GAP_Grb2) - par(I.kd17) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
r18 = par(I.k18) * X(I.Ras_GDP) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.kd18) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP);
r19 = par(I.kd19) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) - par(I.k19) * X(I.Ras_GTP) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
r20 = par(I.k20) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) * X(I.Ras_GTPa) - par(I.kd20) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP);
r21 = par(I.kd21) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) - par(I.k21) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) * X(I.Ras_GDP);
r105 = par(I.k105) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.kd105) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
r106 = par(I.k106) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) * X(I.Prot) - par(I.kd106) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot);
r107 = par(I.kd107) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) - par(I.k107) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) * X(I.Proti);
r108 = par(I.k108) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) - par(I.kd108) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP);
r109 = par(I.k109) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) * X(I.Prot) - par(I.kd109) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot);
r110 = par(I.kd110) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) - par(I.k110) *  X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) * X(I.Proti);
r111 = par(I.k111) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) - par(I.kd111) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP);
r112 = par(I.k112) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) * X(I.Prot) - par(I.kd112) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot);
r113 = par(I.kd113) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot)- par(I.k113) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) * X(I.Proti);
r126 = par(I.k126) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) * X(I.ERK_PP) - par(I.kd126) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP);
r143 = par(I.kd143) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) - par(I.k143) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) * X(I.ERK_PP);

dX(I.EGF_EGFRa_2_GAP_Grb2) = -r4 - r9 + r16 - r17;
dX(I.EGF_EGFRa_2_GAP_Grb2_Prot) = r4 - r5;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos) = r17 - r18 + r19 - r20 + r21 - r34 - r105 - r106 - r126;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) = r126 - r143;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = r106 - r107;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = r18 - r19 - r108 - r109;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = r20 - r21 - r111 - r112;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = r109 - r110;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = r112 - r113;
dX(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) = r143;

%%% -----------------------------------------------------------------------
%%% Equ 23-27: for Pathway without Shc internalised forms
r63 = par(I.k63) * X(I.EGF_EGFRia_2_GAP) * X(I.Grb2) - par(I.kd63) * X(I.EGF_EGFRia_2_GAP_Grb2);
r64 = par(I.k64) * X(I.Sos) * X(I.EGF_EGFRia_2_GAP_Grb2) - par(I.kd64) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
r65 = par(I.k65) * X(I.Ras_GDP) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.kd65) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP);
r66 = par(I.kd66) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) - par(I.k66) * X(I.Rasi_GTP) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
r67 = par(I.k67) * X(I.Rasi_GTPa) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.kd67) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP);
r68 = par(I.kd68) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) - par(I.k68) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) * X(I.Ras_GDP);
r69 = par(I.k69) * X(I.Shc) * X(I.EGF_EGFRia_2_GAP) - par(I.kd69) * X(I.EGF_EGFRia_2_GAP_Shc);
r79 = par(I.kd79) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.k79) * X(I.EGF_EGFRia_2_GAP) * X(I.Shca_Grb2_Sos);
r80 = par(I.kd80) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k80) * X(I.EGF_EGFRia_2_GAP) * X(I.Grb2_Sos);
r81 = par(I.kd81) * X(I.EGF_EGFRia_2_GAP_Shca) - par(I.k81) * X(I.EGF_EGFRia_2_GAP) * X(I.Shca);
r82 = par(I.kd82) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2) - par(I.k82) * X(I.EGF_EGFRia_2_GAP) * X(I.Shca_Grb2);
r127 = par(I.k127) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) * X(I.ERKi_PP) - par(I.kd127) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP);
r146 = par(I.kd146) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) - par(I.k146) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) * X(I.ERKi_PP);
r132 = par(I.k132) * X(I.EGF_EGFRia_2_GAP) - par(I.kd132) * X(I.EGF_EGFRia_2_GAP_deg);
r133 = par(I.k133) * X(I.EGF_EGFRia_2_GAP_Grb2) - par(I.kd133) * X(I.EGF_EGFRia_2_GAP_Grb2_deg);
r134 = par(I.k134) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.kd134) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg);
r135 = par(I.k135) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) - par(I.kd135) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg);
r136 = par(I.k136) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) - par(I.kd136) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg);

dX(I.EGF_EGFRia_2_GAP) =  r14 - r63 - r69 + r79 + r80 + r81 + r82 + r102 - r132;
dX(I.EGF_EGFRia_2_GAP_Grb2) =  r5 + r9 + r63 - r64 - r133;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos) =  r64 - r65 + r66 - r67 + r68 - r80 + r105 + r107 - r127 - r134;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = r65 - r66 + r108 + r110 - r135;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = r67 - r68 + r111 + r113 - r136;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) = r127 - r146;
dX(I.EGF_EGFRia_2_GAP_Grb2_deg) =  + r133;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) = + r134 + r146;
dX(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg) = + r135 + r136;
dX(I.EGF_EGFRia_2_GAP_deg) = + r132 ;
dX(I.EGF_EGFRia_2_deg) = r62;
%%% ----------------------------------------------------------------------
%%% Equ 28-31: Shc dynamic (activation and complex formation with Grb, Sos)
%%%
r33 = par(I.kd33) * X(I.Shca_Grb2_Sos) - par(I.k33) * X(I.Shca) * X(I.Grb2_Sos);
r36 = par(I.k36) * X(I.Shca) -  par(I.kd36) * X(I.Shc);
r38 = par(I.k38) * X(I.Grb2) * X(I.Shca) - par(I.kd38) * X(I.Shca_Grb2);
r40 = par(I.k40) * X(I.Sos) * X(I.Shca_Grb2) - par(I.kd40) * X(I.Shca_Grb2_Sos);

dX(I.Shc) = -r22 + r36 - r69;
dX(I.Shca) = r33 - r36 + r37 - r38 + r81;
dX(I.Shca_Grb2) = r38 + r39 - r40 + r82;
dX(I.Shca_Grb2_Sos) =  r32 - r33 + r40 + r79;

%%% ----------------------------------------------------------------------
%%% Equ 32-41: Pathway including src homology and collagen protein (Shc) 
r23 = par(I.k23) * X(I.EGF_EGFRa_2_GAP_Shc) - par(I.kd23) * X(I.EGF_EGFRa_2_GAP_Shca);
r24 = par(I.k24) * X(I.Grb2) * X(I.EGF_EGFRa_2_GAP_Shca) - par(I.kd24) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2);
r25 = par(I.k25) * X(I.Sos) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) - par(I.kd25) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
r26 = par(I.k26) * X(I.Ras_GDP) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.kd26) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP);
r27 = par(I.kd27) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) - par(I.k27) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) * X(I.Ras_GTP);
r30 = par(I.k30) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) * X(I.Ras_GTPa) - par(I.kd30) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP);
r31 = par(I.kd31) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) - par(I.k31) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) * X(I.Ras_GDP);
r41 = par(I.k41) * X(I.Grb2_Sos) * X(I.EGF_EGFRa_2_GAP_Shca) - par(I.kd41) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
r103 = par(I.k103) * X(I.EGF_EGFRa_2_GAP_Shc) - par(I.kd103) * X(I.EGF_EGFRia_2_GAP_Shc);
r104 = par(I.k104) * X(I.EGF_EGFRa_2_GAP_Shca) - par(I.kd104) * X(I.EGF_EGFRia_2_GAP_Shca);
r114 = par(I.k114) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) - par(I.kd114) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);
r115 = par(I.k115) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) * X(I.Prot) - par(I.kd115) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot);
r116 = par(I.kd116) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) - par(I.k116) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2) * X(I.Proti);
r117 = par(I.k117) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.kd117) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
r118 = par(I.k118) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) * X(I.Prot) - par(I.kd118) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot);
r119 = par(I.kd119) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) - par(I.k119) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) * X(I.Proti);
r120 = par(I.k120) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) - par(I.kd120) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP);
r121 = par(I.k121) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) * X(I.Prot) - par(I.kd121) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot);
r122 = par(I.kd122) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) - par(I.k122) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) * X(I.Proti);
r123 = par(I.k123) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) - par(I.kd123) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP);
r124 = par(I.k124) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) * X(I.Prot) - par(I.kd124) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot);
r125 = par(I.kd125) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot)- par(I.k125) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) * X(I.Proti);
r128 = par(I.k128) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) * X(I.ERK_PP) - par(I.kd128) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP);
r144 = par(I.kd144) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) - par(I.k144) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) * X(I.ERK_PP);

dX(I.EGF_EGFRa_2_GAP_Shc) = r22 - r23 - r103;
dX(I.EGF_EGFRa_2_GAP_Shca) =  r23 - r24 - r37 - r41 - r104;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2) = r24 - r25 - r39 - r114 - r115;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = r115 - r116;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) =  r25 - r26 + r27 - r30 + r31 - r32 + r41 - r117 - r118 - r128;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) = r128 - r144;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = r118 - r119;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = r26 - r27 - r120 - r121;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = r30 - r31 - r123 - r124;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = r121 - r122;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = r124 - r125;
dX(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) = r144;
%%% -----------------------------------------------------------------------
%%% Equ. 42-43: Prot and Proti
r15 = par(I.k15) * X(I.Proti) - par(I.kd15) * X(I.Prot);

dX(I.Prot) = - r4 + r15 - r106 - r109 - r112 - r115 - r118 - r121 - r124;
dX(I.Proti) =  r5 - r15 + r107 + r110 + r113 + r116 + r119 + r122 + r125;

%%% ----------------------------------------------------------------------
%%% Equ 44-49: Pathway including Shc internalised forms
r70 = par(I.k70) * X(I.EGF_EGFRia_2_GAP_Shc) - par(I.kd70) * X(I.EGF_EGFRia_2_GAP_Shca);
r71 = par(I.k71) * X(I.Grb2) * X(I.EGF_EGFRia_2_GAP_Shca) - par(I.kd71) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);
r72 = par(I.k72) * X(I.Sos) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2) - par(I.kd72) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
r73 = par(I.k73) * X(I.Ras_GDP) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.kd73) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP);
r74 = par(I.kd74) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) - par(I.k74) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) * X(I.Rasi_GTP);

r77 = par(I.k77) * X(I.Rasi_GTPa) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.kd77) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP);
r78 = par(I.kd78) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) - par(I.k78) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) * X(I.Ras_GDP);

r83 = par(I.k83) * X(I.Grb2_Sos) * X(I.EGF_EGFRia_2_GAP_Shca) - par(I.kd83) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
r129 = par(I.k129) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) * X(I.ERKi_PP) - par(I.kd129) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP);
r147 = par(I.kd147) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) - par(I.k147) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) * X(I.ERKi_PP);
r137 = par(I.k137) * X(I.EGF_EGFRia_2_GAP_Shc) - par(I.kd137) * X(I.EGF_EGFRia_2_GAP_Shc_deg);
r138 = par(I.k138) * X(I.EGF_EGFRia_2_GAP_Shca) - par(I.kd138) * X(I.EGF_EGFRia_2_GAP_Shc_deg);
r139 = par(I.k139) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2) - par(I.kd139) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg);
r140 = par(I.k140) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.kd140) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg);
r141 = par(I.k141) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) - par(I.kd141) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg);
r142 = par(I.k142) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) - par(I.kd142) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg);


dX(I.EGF_EGFRia_2_GAP_Shc) =  r69 - r70 + r103 - r137;
dX(I.EGF_EGFRia_2_GAP_Shca) = r70 - r71 - r81 - r83 + r104 - r138;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2) =  r71 - r72 - r82 + r114 + r116 - r139;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = r72 - r73 + r74 - r77 + r78 - r79 + r83 + r117 + r119 - r129 - r140;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = r73 - r74 + r120 + r122 - r141;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = r77 - r78 + r123 + r125 - r142;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) = r129 - r147;
dX(I.EGF_EGFRia_2_GAP_Shc_deg) = + r137 + r138;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_deg) = + r139;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) =  + r140 + r147;
dX(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg) = + r141 + r142;
%%% ----------------------------------------------------------------------
%%% Equ 50-55: Common state variables in both pathways: 
%%% GTPase activating protein (GAP),
%%% growth factor receptor-bound protein (Grb2),
%%% son of sevenless (Sos),
%%% Grb2_Sos, 
%%% guanosine triphosphate of Ras (Ras:GTP), 
%%% guanosine diphosphate of Ras (Ras:GDP)
r28 = par(I.k28) * X(I.Ras_GTP) * X(I.Raf) - par(I.kd28) * X(I.Raf_Ras_GTP);
r35 = par(I.kd35) * X(I.Grb2_Sos) - par(I.k35) * X(I.Sos) * X(I.Grb2);
r130 = par(I.k130) * X(I.Sos) * X(I.ERK_PP) - par(I.kd130) * X(I.Sos_ERK_PP);  
r131 = par(I.k131) * X(I.Sos) * X(I.ERKi_PP) - par(I.kd131) * X(I.Sos_ERKi_PP);        
r145 = par(I.kd145) * X(I.Sos_ERK_PP) - par(I.k145) * X(I.Sosi) * X(I.ERK_PP);
r148 = par(I.kd148) * X(I.Sos_ERKi_PP) - par(I.k148) * X(I.Sosi) * X(I.ERKi_PP);

dX(I.GAP) = - r8 - r14;
dX(I.Grb2) = - r16 - r24 + r35 - r38 - r63 - r71;
dX(I.Sos) = - r17 - r25 + r35 - r40 - r64 - r72 - r130 - r131;
dX(I.Sos_ERK_PP) = r130 - r145; 
dX(I.Sos_ERKi_PP) = r131 - r148;
dX(I.Sosi) = r145 + r148;
dX(I.Grb2_Sos) = r33 + r34 - r35 - r41 + r80 - r83;
dX(I.Ras_GTP) = r19 + r27 - r28;
dX(I.Ras_GDP) = - r18 + r21 - r26 + r31 - r65 + r68 - r73 + r78;

%%% -----------------------------------------------------------------------
%%% Equ 56-66: Raf and Ras dynamic with internalised forms
%%% Rapidly accelerated fibrosarcoma (Raf)
%%%
r29 = par(I.kd29) * X(I.Raf_Ras_GTP) - par(I.k29) * X(I.Ras_GTPa) * X(I.Rafa);
r42 = par(I.k42) * X(I.Phosphatase1) * X(I.Rafa) - par(I.kd42) * X(I.Rafa_Pase);
r43 = par(I.kd43) * X(I.Rafa_Pase)- par(I.k43) * X(I.Raf) * X(I.Phosphatase1);
r44 = par(I.k44) * X(I.MEK) * X(I.Rafa) - par(I.kd44) * X(I.MEK_Rafa);
r45 = par(I.kd45) * X(I.MEK_Rafa)- par(I.k45) * X(I.MEK_P) * X(I.Rafa);
r46 = par(I.k46) * X(I.MEK_P) * X(I.Rafa) - par(I.kd46) * X(I.MEK_P_Rafa);
r47 = par(I.kd47) * X(I.MEK_P_Rafa)- par(I.k47) * X(I.MEK_PP) * X(I.Rafa);
r75 = par(I.k75) * X(I.Rasi_GTP) * X(I.Raf) - par(I.kd75) * X(I.Rafi_Rasi_GTP);
r76 = par(I.kd76) * X(I.Rafi_Rasi_GTP) - par(I.k76) * X(I.Rasi_GTPa) * X(I.Rafia);
r84 = par(I.k84) * X(I.Phosphatase1) * X(I.Rafia) - par(I.kd84) * X(I.Rafia_Pase);
r85 = par(I.kd85) * X(I.Rafia_Pase) - par(I.k85) * X(I.Rafia) * X(I.Phosphatase1);
r86 = par(I.k86) * X(I.MEK) * X(I.Rafia) - par(I.kd86) * X(I.MEKi_Rafia);
r87 = par(I.kd87) * X(I.MEKi_Rafia) - par(I.k87) * X(I.MEKi_P) * X(I.Rafia);
r88 = par(I.k88) * X(I.Rafia) * X(I.MEKi_P) - par(I.kd88) * X(I.MEKi_P_Rafia);
r89 = par(I.kd89) * X(I.MEKi_P_Rafia) - par(I.k89) * X(I.MEKi_PP) * X(I.Rafia);

dX(I.Raf) = - r28 +  r43 - r75 +  r85;
dX(I.Ras_GTPa) = - r20 +  r29 - r30;
dX(I.Raf_Ras_GTP) =  r28 - r29;
dX(I.Rafa) = r29 - r42 - r44 + r45 - r46 + r47;
dX(I.Rafa_Pase) = r42 - r43;
dX(I.Phosphatase1) = -r42 + r43 - r84 + r85;

dX(I.Rasi_GTP) =  r66 +  r74 - r75;
dX(I.Rafi_Rasi_GTP) = r75 - r76;
dX(I.Rasi_GTPa) = - r67 +  r76 - r77;
dX(I.Rafia) = r76 - r84 - r86 +  r87 - r88 +  r89;
dX(I.Rafia_Pase) = r84 - r85;

%%% -----------------------------------------------------------------------
%%% Equ 67-80: MEK dynamic with internalised forms
%%% mitogen-activated protein kinase (MEK)
%%%
r48 = par(I.k48) * X(I.MEK_PP) * X(I.Phosphatase2) - par(I.kd48) * X(I.MEK_PP_Pase2);
r49 = par(I.kd49) * X(I.MEK_PP_Pase2)- par(I.k49) * X(I.MEK_P) * X(I.Phosphatase2); 
r50 = par(I.k50) * X(I.Phosphatase2) * X(I.MEK_P) - par(I.kd50) * X(I.MEK_P_Pase2);
r51 = par(I.kd51) * X(I.MEK_P_Pase2)- par(I.k51) * X(I.MEK) * X(I.Phosphatase2);
r52 = par(I.k52) * X(I.ERK) * X(I.MEK_PP) - par(I.kd52) * X(I.ERK_MEK_PP);
r53 = par(I.kd53) * X(I.ERK_MEK_PP)- par(I.k53) * X(I.ERK_P) * X(I.MEK_PP);
r54 = par(I.k54) * X(I.MEK_PP) * X(I.ERK_P) - par(I.kd54) * X(I.ERK_P_MEK_PP);
r55 = par(I.kd55) * X(I.ERK_P_MEK_PP)- par(I.k55) * X(I.ERK_PP) * X(I.MEK_PP);

r90 = par(I.k90) * X(I.MEKi_PP) * X(I.Phosphatase2) - par(I.kd90) * X(I.MEKi_PP_Pase2i);
% ERROR: reaction r91 does NOT involve I.MEK_PP but I.MEK_P;
% Luckily, par(I.k91)=0, so the error has no impact
% r91 = par(I.kd91) * X(I.MEKi_PP_Pase2i)- par(I.k91) * X(I.MEKi_PP) * X(I.Phosphatase2); % original wrong eq. 
r91 = par(I.kd91) * X(I.MEKi_PP_Pase2i)- par(I.k91) * X(I.MEKi_P) * X(I.Phosphatase2); % corrected 13.11.2021 by HUI
r92 = par(I.k92) * X(I.Phosphatase2) * X(I.MEKi_P) - par(I.kd92) * X(I.MEKi_P_Pase2i);
r93 = par(I.kd93) * X(I.MEKi_P_Pase2i) - par(I.k93) * X(I.MEK) * X(I.Phosphatase2);
r94 = par(I.k94) * X(I.ERK) * X(I.MEKi_PP) - par(I.kd94) * X(I.ERKi_MEKi_PP);
r95 = par(I.kd95) * X(I.ERKi_MEKi_PP)- par(I.k95) * X(I.ERKi_P) * X(I.MEKi_PP);
r96 = par(I.k96) * X(I.MEKi_PP) * X(I.ERKi_P) - par(I.kd96) * X(I.ERKi_P_MEKi_PP);
r97 = par(I.kd97) * X(I.ERKi_P_MEKi_PP) - par(I.k97) * X(I.ERKi_PP) * (I.MEKi_PP);


dX(I.MEK) = - r44 + r51 - r86 + r93;
dX(I.MEK_Rafa) = r44 - r45;
dX(I.MEK_P) = r45 - r46 + r49 - r50; 
dX(I.MEK_P_Rafa) = r46 - r47;
dX(I.MEK_PP) =  r47 - r48 - r52 + r53 - r54 + r55;
dX(I.MEK_PP_Pase2) =  r48 - r49;
dX(I.Phosphatase2) = - r48 + r49 - r50 + r51 - r90 + r91 - r92 + r93;
dX(I.MEK_P_Pase2) = r50 - r51;

dX(I.MEKi_Rafia) =  r86 - r87;
dX(I.MEKi_P) =  r87 - r88 + r91 - r92; 
dX(I.MEKi_P_Rafia) = r88 - r89;
dX(I.MEKi_PP) = r89 - r90 - r94 + r95 - r96 + r97;
dX(I.MEKi_PP_Pase2i) = r90 - r91;
dX(I.MEKi_P_Pase2i) = r92 - r93;

%%% -----------------------------------------------------------------------
%%% Equ 81-94: Extracellular signal-regulated kinase (ERK) dynamic with internalised forms
%%%
r56 = par(I.k56) * X(I.ERK_PP) * X(I.Phosphatase3) - par(I.kd56) * X(I.ERK_PP_Pase3);
r57 = par(I.kd57) * X(I.ERK_PP_Pase3) - par(I.k57) * X(I.ERK_P) * X(I.Phosphatase3);
r58 = par(I.k58) * X(I.Phosphatase3) * X(I.ERK_P) - par(I.kd58) * X(I.ERK_P_Pase3);
r59 = par(I.kd59) * X(I.ERK_P_Pase3) - par(I.k59) * X(I.ERK) * X(I.Phosphatase3);
r98 = par(I.k98) * X(I.ERKi_PP) * X(I.Phosphatase3) - par(I.kd98) * X(I.ERKi_PP_Pase3i);
r99 = par(I.kd99) * X(I.ERKi_PP_Pase3i) - par(I.k99) * X(I.ERKi_P) * X(I.Phosphatase3);
r100 = par(I.k100) * X(I.Phosphatase3) * X(I.ERKi_P) - par(I.kd100) * X(I.ERKi_P_Pase3i);
r101 = par(I.kd101) * X(I.ERKi_P_Pase3i) - par(I.k101) * X(I.ERK) * X(I.Phosphatase3);


dX(I.ERK) = - r52 + r59 - r94 + r101;
dX(I.ERK_MEK_PP) = r52 - r53;
dX(I.ERK_P) = r53 - r54 + r57 - r58;
dX(I.ERK_P_MEK_PP) = r54 - r55;
dX(I.ERK_PP) = r55 - r56 - r126 - r128 - r130 + r143 + r144 + r145;
dX(I.Phosphatase3) = - r56 + r57 - r58 + r59 - r98 + r99 - r100 + r101;
dX(I.ERK_PP_Pase3) =  r56 - r57;
dX(I.ERK_P_Pase3) =  r58 - r59;

dX(I.ERKi_P) = r95 - r96 + r99 - r100;
dX(I.ERKi_MEKi_PP) = r94 - r95;
dX(I.ERKi_P_MEKi_PP) =  r96 - r97;
dX(I.ERKi_PP) = r97 - r98 - r127 - r129 - r131 + r146 + r147 + r148;
dX(I.ERKi_PP_Pase3i) =  r98 - r99;
dX(I.ERKi_P_Pase3i) = r100 - r101;

dX(I.AUC_ERK_PP) = X(I.ERK_PP);

%%% -----------------------------------------------------------------------

end

