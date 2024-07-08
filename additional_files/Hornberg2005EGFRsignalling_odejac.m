%%% Version: 19 Jan 2023
%%%
%%% jac  = <MODELNAME>_odejac(t,X,par,model)
%%%
%%% This function defines the jacobian of the ode system 
%%%
%%% Input   t           time
%%%         X           state vector
%%%         par         parameter vector
%%%         model       model structure containing the index structure of
%%%                     the model
%%%                   
%%%
%%% Output  jac         Jacobian of the right hand side with respect to  
%%%                     the state variables
%%%
%%% Reference:
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
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function DF = Hornberg2005EGFRsignalling_odejac(~,X,par,model)

%%% assign model indexing
I  = model.I;

%%% initialize the jacobian
% DF  = spalloc(I.nstates,I.nstates,3*I.nstates);
sz = size(X);
if sz(1) == 1
    DF = repmat(0*X, [length(X) 1]);
elseif sz(2) == 1
    DF = repmat(0*X, [1 length(X)]);
else
    disp(sz)
    error('Wrong size of X for jac fun.')
end

%%% -----------------------------------------------------------------------
%%% define jacobian

%%% -----------------------------------------------------------------------
%%% Equ 1-13: EGF and EGFR receptor dynamic with internalisation
%%%

DF(I.EGF,:) = 0;
DF(I.EGFR,I.EGF) = - par(I.k1)  * X(I.EGFR);
DF(I.EGFR,I.EGFR) = - par(I.k1) * X(I.EGF) - par(I.k6) - par(I.kd13);
DF(I.EGFR,I.EGF_EGFR) = par(I.kd1);
DF(I.EGFR,I.EGFRi) = par(I.kd6);

DF(I.EGF_EGFR,I.EGF) =  par(I.k1)  * X(I.EGFR);
DF(I.EGF_EGFR,I.EGFR) = par(I.k1) * X(I.EGF);
DF(I.EGF_EGFR,I.EGF_EGFR) = - par(I.kd1)- 4* par(I.k2) * X(I.EGF_EGFR);
DF(I.EGF_EGFR,I.EGF_EGFR_2) =  2*par(I.kd2);

DF(I.EGF_EGFR_2,I.EGF_EGFR) = 2 * par(I.k2) * X(I.EGF_EGFR);
DF(I.EGF_EGFR_2,I.EGF_EGFR_2) = - par(I.kd2) - par(I.k3);
DF(I.EGF_EGFR_2,I.EGF_EGFRa_2) = par(I.kd3);

DF(I.EGF_EGFRa_2,I.EGF_EGFR_2) = par(I.k3);
DF(I.EGF_EGFRa_2,I.EGF_EGFRa_2) = - par(I.kd3) - par(I.k7) - par(I.k8) * X(I.GAP);
DF(I.EGF_EGFRa_2,I.EGF_EGFRa_2_GAP) = par(I.kd8) ;
DF(I.EGF_EGFRa_2,I.EGF_EGFRia_2) = par(I.kd7);
DF(I.EGF_EGFRa_2,I.GAP) = - par(I.k8) * X(I.EGF_EGFRa_2);

DF(I.EGFRi,I.EGFR) = par(I.k6) ;
DF(I.EGFRi,I.EGFi) = - par(I.k10b) * X(I.EGFRi);
DF(I.EGFRi,I.EGFRi) = - par(I.kd6) -  par(I.k10b) * X(I.EGFi) - par(I.k60);
DF(I.EGFRi,I.EGFRideg) = par(I.kd60);
DF(I.EGFRi,I.EGF_EGFRi) = par(I.kd10);

DF(I.EGFRideg,I.EGFRi) = par(I.k60);
DF(I.EGFRideg,I.EGFRideg) = - par(I.kd60);

DF(I.EGFi,I.EGFRi) = - par(I.k10b)  * X(I.EGFi);
DF(I.EGFi,I.EGFi) = - par(I.k10b) * X(I.EGFRi) - par(I.k61);
DF(I.EGFi,I.EGFideg) = par(I.kd61);
DF(I.EGFi,I.EGF_EGFRi) = par(I.kd10);

DF(I.EGFideg,I.EGFi) = par(I.k61);
DF(I.EGFideg,I.EGFideg) = - par(I.kd61);

DF(I.EGF_EGFRi,I.EGFRi) = par(I.k10b) * X(I.EGFi);
DF(I.EGF_EGFRi,I.EGFi) = par(I.k10b) * X(I.EGFRi);
DF(I.EGF_EGFRi,I.EGF_EGFRi) = - par(I.kd10) - 2* par(I.k2) * X(I.EGF_EGFRi);
DF(I.EGF_EGFRi,I.EGF_EGFRi_2) = par(I.kd2);

DF(I.EGF_EGFRi_2,I.EGF_EGFRi) = 2* par(I.k2) * X(I.EGF_EGFRi);
DF(I.EGF_EGFRi_2,I.EGF_EGFRi_2) = - par(I.kd2) - par(I.k12);
DF(I.EGF_EGFRi_2,I.EGF_EGFRia_2) = par(I.kd12);

DF(I.EGF_EGFRia_2,I.EGF_EGFRa_2) = par(I.k7);
DF(I.EGF_EGFRia_2,I.EGF_EGFRi_2) = par(I.k12) ;
DF(I.EGF_EGFRia_2,I.EGF_EGFRia_2) = -par(I.kd7) - par(I.kd12) - par(I.k14) * X(I.GAP) - par(I.k62) ;
DF(I.EGF_EGFRia_2,I.EGF_EGFRia_2_deg) = par(I.kd62);
DF(I.EGF_EGFRia_2,I.GAP) = -par(I.k14) * X(I.EGF_EGFRia_2);
DF(I.EGF_EGFRia_2,I.EGF_EGFRia_2_GAP) = par(I.kd14);

DF(I.EGF_EGFRia_2_deg,I.EGF_EGFRia_2) =  par(I.k62);
DF(I.EGF_EGFRia_2_deg,I.EGF_EGFRia_2_deg) = - par(I.kd62) ;

%%% -----------------------------------------------------------------------
%%% Equ 14: EGF:EGFRa:2:GAP start of both pathways

DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2) = par(I.k8) * X(I.GAP);
DF(I.EGF_EGFRa_2_GAP,I.GAP) = par(I.k8) * X(I.EGF_EGFRa_2);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP) = - par(I.kd8) - par(I.k16) * X(I.Grb2)...
    -par(I.k22) * X(I.Shc) - par(I.k32)  * X(I.Shca_Grb2_Sos)...
    - par(I.k34) * X(I.Grb2_Sos) - par(I.k37) * X(I.Shca)...
    - par(I.k39) * X(I.Shca_Grb2) - par(I.k102);
DF(I.EGF_EGFRa_2_GAP,I.Grb2) = - par(I.k16) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP,I.Shc) = - par(I.k22) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Grb2) = par(I.kd16);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Shc) = par(I.kd22);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd32);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.kd34);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRia_2_GAP) = par(I.kd102);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Shca) = par(I.kd37);
DF(I.EGF_EGFRa_2_GAP,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.kd39);
DF(I.EGF_EGFRa_2_GAP,I.Shca_Grb2_Sos) = - par(I.k32) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP,I.Grb2_Sos) = - par(I.k34) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP,I.Shca) = - par(I.k37) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP,I.Shca_Grb2) = - par(I.k39) * X(I.EGF_EGFRa_2_GAP);
%%% -----------------------------------------------------------------------
%%% Equ 15-22: for Pathway without Shc

DF(I.EGF_EGFRa_2_GAP_Grb2,I.EGF_EGFRa_2_GAP_Grb2) = - par(I.k4) * X(I.Prot) - par(I.k9) - par(I.kd16) - par(I.k17) * X(I.Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.EGF_EGFRia_2_GAP_Grb2) = par(I.kd9);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.Prot) = - par(I.k4) * X(I.EGF_EGFRa_2_GAP_Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.EGF_EGFRa_2_GAP_Grb2_Prot) = par(I.kd4);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.Grb2) = par(I.k16) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.EGF_EGFRa_2_GAP) = par(I.k16) * X(I.Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.Sos) = - par(I.k17) * X(I.EGF_EGFRa_2_GAP_Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.kd17);

DF(I.EGF_EGFRa_2_GAP_Grb2_Prot,I.EGF_EGFRa_2_GAP_Grb2) = par(I.k4) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Grb2_Prot,I.Prot) = par(I.k4) * X(I.EGF_EGFRa_2_GAP_Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2_Prot,I.Proti) = par(I.k5) * X(I.EGF_EGFRia_2_GAP_Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2_Prot,I.EGF_EGFRia_2_GAP_Grb2) = par(I.k5) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Grb2_Prot,I.EGF_EGFRa_2_GAP_Grb2_Prot) = - par(I.kd4) - par(I.kd5);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2) = par(I.k17) * X(I.Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Sos) = par(I.k17) * X(I.EGF_EGFRa_2_GAP_Grb2);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.kd17) - par(I.k18) * X(I.Ras_GDP) ...
    - par(I.k19) * X(I.Ras_GTP) - par(I.k20) * X(I.Ras_GTPa) ...
    - par(I.k21) * X(I.Ras_GDP) - par(I.kd34)...
    - par(I.kd105) - par(I.k106) * X(I.Prot)...
    - par(I.k126) * X(I.ERK_PP) ;
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.ERK_PP) = - par(I.k126) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) = par(I.kd126);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Ras_GDP) = - par(I.k18) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k21) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd18) + par(I.kd19);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Ras_GTP) = - par(I.k19) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Ras_GTPa) = - par(I.k20) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd20) + par(I.kd21);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Grb2_Sos) = par(I.k34) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP) = par(I.k34) * X(I.Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.k105);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.Prot) = - par(I.k106) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) ;
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = par(I.kd106);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.k126) * X(I.ERK_PP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP,I.ERK_PP) = par(I.k126) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) + par(I.k143) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) = - par(I.kd126) - par(I.kd143);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) = par(I.k143) * X(I.ERK_PP);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.k106) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot,I.Prot) = par(I.k106) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot,I.Proti) = par(I.k107) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.k107) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = - par(I.kd106) - par(I.kd107);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.Ras_GDP) = par(I.k18)* X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.k18) * X(I.Ras_GDP) + par(I.k19) * X(I.Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = - par(I.kd18) - par(I.kd19) - par(I.k108) - par(I.k109) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.Ras_GTP) = par(I.k19) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) =  par(I.kd108);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.Prot) = - par(I.k109) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = par(I.kd109);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.Ras_GTPa) = par(I.k20) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.k20) * X(I.Ras_GTPa) + par(I.k21) * X(I.Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = - par(I.kd20) - par(I.kd21) - par(I.k111) - par(I.k112) * X(I.Prot) ;
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.Ras_GDP) = par(I.k21) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd111);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.Prot) = - par(I.k112) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = par(I.kd112);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) =  par(I.k109) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot,I.Prot) = par(I.k109) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot,I.Proti) = par(I.k110) *  X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = par(I.k110) *  X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = - par(I.kd110) - par(I.kd109);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = par(I.k112) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot,I.Prot) = par(I.k112) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot,I.Proti) = par(I.k113) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.k113) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = - par(I.kd112) - par(I.kd113);

DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg,I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) = - par(I.k143) * X(I.ERK_PP);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg,I.ERK_PP) = - par(I.k143) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg);
DF(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg,I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) = par(I.kd143);
%%% -----------------------------------------------------------------------
%%% Equ 23-27: for Pathway without Shc internalised forms

DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2) = par(I.k14) * X(I.GAP);
DF(I.EGF_EGFRia_2_GAP,I.GAP) = par(I.k14) * X(I.EGF_EGFRia_2);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP) =  - par(I.kd14) - par(I.k63) * X(I.Grb2)...
    - par(I.k69) * X(I.Shc) - par(I.k79)* X(I.Shca_Grb2_Sos)...
    - par(I.k80) * X(I.Grb2_Sos) - par(I.k81) * X(I.Shca)...
    - par(I.k82) * X(I.Shca_Grb2) - par(I.kd102) - par(I.k132);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_deg) = par(I.kd132);
DF(I.EGF_EGFRia_2_GAP,I.Grb2) = - par(I.k63) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Grb2) = par(I.kd63);
DF(I.EGF_EGFRia_2_GAP,I.Shc) = -par(I.k69) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Shc) = par(I.kd69);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd79);
DF(I.EGF_EGFRia_2_GAP,I.Shca_Grb2_Sos) = - par(I.k79) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.kd80);
DF(I.EGF_EGFRia_2_GAP,I.Grb2_Sos) = - par(I.k80) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Shca) = par(I.kd81);
DF(I.EGF_EGFRia_2_GAP,I.Shca) = - par(I.k81) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.kd82);
DF(I.EGF_EGFRia_2_GAP,I.Shca_Grb2) = - par(I.k82) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP,I.EGF_EGFRa_2_GAP) = par(I.k102) ;

DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRa_2_GAP_Grb2_Prot) = par(I.kd5);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.Proti) = - par(I.k5) * X(I.EGF_EGFRia_2_GAP_Grb2);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRa_2_GAP_Grb2) =  par(I.k9);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRia_2_GAP) = par(I.k63) * X(I.Grb2);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.Grb2) = par(I.k63) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRia_2_GAP_Grb2) = - par(I.kd63) -par(I.k64) * X(I.Sos) - par(I.k5) * X(I.Proti) - par(I.kd9) - par(I.k133);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRia_2_GAP_Grb2_deg) = par(I.kd133);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.Sos) = - par(I.k64) * X(I.EGF_EGFRia_2_GAP_Grb2);
DF(I.EGF_EGFRia_2_GAP_Grb2,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.kd64);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Sos) = par(I.k64)* X(I.EGF_EGFRia_2_GAP_Grb2);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2) = par(I.k64) * X(I.Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.kd64) - par(I.k65) * X(I.Ras_GDP) ...
    - par(I.k66) * X(I.Rasi_GTP) - par(I.k67) * X(I.Rasi_GTPa)...
    - par(I.k68) * X(I.Ras_GDP) - par(I.kd80)...
    - par(I.kd105) - par(I.k127) * X(I.ERKi_PP) - par(I.k134)...
    - par(I.k107) * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.ERKi_PP) = - par(I.k127) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) = par(I.kd127);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) = par(I.kd134);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Proti) = - par(I.k107) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Ras_GDP) = - par(I.k65) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k68) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd65) + par(I.kd66);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Rasi_GTP) = - par(I.k66) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Rasi_GTPa) = - par(I.k67) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd67) + par(I.kd68) ;
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.Grb2_Sos) = par(I.k80) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRia_2_GAP) = par(I.k80) * X(I.Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.k105);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = par(I.kd107);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.Ras_GDP) = par(I.k65) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.k65) * X(I.Ras_GDP) + par(I.k66) * X(I.Rasi_GTP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = - par(I.kd65) - par(I.kd66) ...
    - par(I.kd108) - par(I.k135) - par(I.k110) * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg) = par(I.kd135);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.Proti) = - par(I.k110) *  X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.Rasi_GTP) = par(I.k66) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = par(I.k108);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = par(I.kd110);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.Rasi_GTPa) = par(I.k67) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.k67) * X(I.Rasi_GTPa) + par(I.k68)* X(I.Ras_GDP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = - par(I.kd67) - par(I.kd68)...
    - par(I.kd111) - par(I.k136)- par(I.k113) * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg) = par(I.kd136);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.Proti) = - par(I.k113) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.Ras_GDP) = par(I.k68) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = par(I.k111);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = par(I.kd113);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos) =  par(I.k127) * X(I.ERKi_PP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP,I.ERKi_PP) = par(I.k127) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) + par(I.k146) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) = - par(I.kd127) - par(I.kd146);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) = par(I.k146) * X(I.ERKi_PP);

DF(I.EGF_EGFRia_2_GAP_Grb2_deg,I.EGF_EGFRia_2_GAP_Grb2_deg) = - par(I.kd133);
DF(I.EGF_EGFRia_2_GAP_Grb2_deg,I.EGF_EGFRia_2_GAP_Grb2) = par(I.k133);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg,I.ERKi_PP) = - par(I.k146) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) = - par(I.kd134)- par(I.k146) * X(I.ERKi_PP);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) = par(I.kd146);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.k134);

DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg) =  - par(I.kd135) - par(I.kd136) ;
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = par(I.k135);
DF(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.k136);

DF(I.EGF_EGFRia_2_GAP_deg,I.EGF_EGFRia_2_GAP_deg) = - par(I.kd132);
DF(I.EGF_EGFRia_2_GAP_deg,I.EGF_EGFRia_2_GAP) = par(I.k132);

%%% ----------------------------------------------------------------------
%%% Equ 28-31: Shc dynamic (activation and complex formation with Grb, Sos)
%%%

DF(I.Shc,I.Shc) = - par(I.k22) * X(I.EGF_EGFRa_2_GAP) - par(I.k69) * X(I.EGF_EGFRia_2_GAP) -  par(I.kd36);
DF(I.Shc,I.EGF_EGFRa_2_GAP) = - par(I.k22) * X(I.Shc);
DF(I.Shc,I.EGF_EGFRa_2_GAP_Shc) = par(I.kd22);
DF(I.Shc,I.Shca) = par(I.k36);
DF(I.Shc,I.EGF_EGFRia_2_GAP) = - par(I.k69) * X(I.Shc);
DF(I.Shc,I.EGF_EGFRia_2_GAP_Shc) = par(I.kd69);

DF(I.Shca,I.Shca_Grb2_Sos) = par(I.kd33);
DF(I.Shca,I.Shc) = par(I.kd36);
DF(I.Shca,I.Shca) = - par(I.k33) * X(I.Grb2_Sos) - par(I.k36) ...
    - par(I.k37) * X(I.EGF_EGFRa_2_GAP) - par(I.k38) * X(I.Grb2) ...
    - par(I.k81) * X(I.EGF_EGFRia_2_GAP);
DF(I.Shca,I.Grb2_Sos) = - par(I.k33) * X(I.Shca);
DF(I.Shca,I.Shca_Grb2_Sos) = par(I.kd33);
DF(I.Shca,I.EGF_EGFRa_2_GAP) = - par(I.k37) * X(I.Shca);
DF(I.Shca,I.EGF_EGFRa_2_GAP_Shca) = par(I.kd37);
DF(I.Shca,I.Grb2) = - par(I.k38) * X(I.Shca);
DF(I.Shca,I.Shca_Grb2) = par(I.kd38);
DF(I.Shca,I.EGF_EGFRia_2_GAP_Shca) = par(I.kd81);
DF(I.Shca,I.EGF_EGFRia_2_GAP) = - par(I.k81) * X(I.Shca);

DF(I.Shca_Grb2,I.Shca) = par(I.k38) * X(I.Grb2);
DF(I.Shca_Grb2,I.Grb2) = par(I.k38) * X(I.Shca);
DF(I.Shca_Grb2,I.Shca_Grb2) = - par(I.kd38) - par(I.k39) * X(I.EGF_EGFRa_2_GAP) ...
    - par(I.k40) * X(I.Sos) - par(I.k82) * X(I.EGF_EGFRia_2_GAP);
DF(I.Shca_Grb2,I.EGF_EGFRa_2_GAP) = - par(I.k39) * X(I.Shca_Grb2);
DF(I.Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.kd39);
DF(I.Shca_Grb2,I.Sos) = - par(I.k40) * X(I.Shca_Grb2);
DF(I.Shca_Grb2,I.Shca_Grb2_Sos) = par(I.kd40);
DF(I.Shca_Grb2,I.EGF_EGFRia_2_GAP) = - par(I.k82) * X(I.Shca_Grb2);
DF(I.Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.kd82);

DF(I.Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP) = - par(I.k32) * X(I.Shca_Grb2_Sos);
DF(I.Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd32);
DF(I.Shca_Grb2_Sos,I.Shca_Grb2_Sos) = - par(I.k32) * X(I.EGF_EGFRa_2_GAP) - par(I.kd33) ...
    - par(I.kd40) - par(I.k79) * X(I.EGF_EGFRia_2_GAP);
DF(I.Shca_Grb2_Sos,I.Shca) = par(I.k33) * X(I.Grb2_Sos);
DF(I.Shca_Grb2_Sos,I.Grb2_Sos) = par(I.k33) * X(I.Shca);
DF(I.Shca_Grb2_Sos,I.Sos) = par(I.k40) * X(I.Shca_Grb2);
DF(I.Shca_Grb2_Sos,I.Shca_Grb2) = par(I.k40) * X(I.Sos);
DF(I.Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP) = - par(I.k79) * X(I.Shca_Grb2_Sos);
DF(I.Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd79);

%%% ----------------------------------------------------------------------
%%% Equ 32-41: Pathway including Shc

DF(I.EGF_EGFRa_2_GAP_Shc,I.Shc) = par(I.k22) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Shc,I.EGF_EGFRa_2_GAP) = par(I.k22) * X(I.Shc);
DF(I.EGF_EGFRa_2_GAP_Shc,I.EGF_EGFRa_2_GAP_Shc) = - par(I.kd22) - par(I.k23) - par(I.k103);
DF(I.EGF_EGFRa_2_GAP_Shc,I.EGF_EGFRa_2_GAP_Shca) = par(I.kd23);
DF(I.EGF_EGFRa_2_GAP_Shc,I.EGF_EGFRia_2_GAP_Shc) = par(I.kd103);

DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRa_2_GAP_Shc) =  par(I.k23) ;
DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRa_2_GAP_Shca) = - par(I.kd23) - par(I.k24) * X(I.Grb2) ...
    - par(I.kd37) - par(I.k41) * X(I.Grb2_Sos) ...
    - par(I.k104);
DF(I.EGF_EGFRa_2_GAP_Shca,I.Grb2) = - par(I.k24) * X(I.EGF_EGFRa_2_GAP_Shca);
DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.kd24);
DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRa_2_GAP) = par(I.k37) * X(I.Shca);
DF(I.EGF_EGFRa_2_GAP_Shca,I.Shca) = par(I.k37) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Shca,I.Grb2_Sos) = - par(I.k41) * X(I.EGF_EGFRa_2_GAP_Shca);
DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd41);
DF(I.EGF_EGFRa_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shca) =  par(I.kd104);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.Grb2) = par(I.k24) * X(I.EGF_EGFRa_2_GAP_Shca);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca) = par(I.k24) * X(I.Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2) = - par(I.kd24) - par(I.k25) * X(I.Sos) ...
    - par(I.kd39) - par(I.k114) - par(I.k115) * X(I.Prot) ;
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.Sos) = - par(I.k25) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd25);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP) = par(I.k39) * X(I.Shca_Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.Shca_Grb2) = par(I.k39) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.kd114);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.Prot) = - par(I.k115) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = par(I.kd115);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.k115) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot,I.Prot) = par(I.k115) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = - par(I.kd115) - par(I.kd116);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.k116) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot,I.Proti) = par(I.k116) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Sos) = par(I.k25) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.k25) * X(I.Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.kd25) - par(I.k26) * X(I.Ras_GDP) ...
    - par(I.k27) * X(I.Ras_GTP) - par(I.k30) * X(I.Ras_GTPa) ...
    - par(I.k31) * X(I.Ras_GDP) - par(I.kd32) ...
    - par(I.kd41) - par(I.k117) - par(I.k118)* X(I.Prot) - par(I.k128) * X(I.ERK_PP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.ERK_PP) = - par(I.k128) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) = par(I.kd128);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Ras_GDP) = - par(I.k26) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.k31) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd26) + par(I.kd27);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Ras_GTP) = - par(I.k27) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd30) + par(I.kd31);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Ras_GTPa) = -par(I.k30) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP) = par(I.k32) * X(I.Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Shca_Grb2_Sos) = par(I.k32) * X(I.EGF_EGFRa_2_GAP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Grb2_Sos) =  par(I.k41) * X(I.EGF_EGFRa_2_GAP_Shca);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca) = par(I.k41) * X(I.Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd117);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.Prot) = - par(I.k118) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = par(I.kd118);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) = - par(I.kd128) - par(I.kd144);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP,I.ERK_PP) = par(I.k128) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) + par(I.k144) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.k128) * X(I.ERK_PP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) = par(I.k144) * X(I.ERK_PP);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,I.Prot) = par(I.k118) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.k118) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = - par(I.kd118) - par(I.kd119);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.k119) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,I.Proti) = par(I.k119) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Ras_GDP) = par(I.k26) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.k26) * X(I.Ras_GDP) + par(I.k27) * X(I.Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = - par(I.kd26) - par(I.kd27) ...
    - par(I.k120) - par(I.k121) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Ras_GTP) = par(I.k27) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd120);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Prot) = - par(I.k121) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = par(I.kd121);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Ras_GTPa) = par(I.k30) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.k30) * X(I.Ras_GTPa) + par(I.k31) * X(I.Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = - par(I.kd30) - par(I.kd31) ...
    - par(I.k123) - par(I.k124) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Ras_GDP) = par(I.k31) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) =  par(I.kd123);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Prot) = - par(I.k124) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = par(I.kd124);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot,I.Prot) = par(I.k121) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.k121) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = - par(I.kd121) - par(I.kd122);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) =  par(I.k122) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot,I.Proti) = par(I.k122) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot,I.Prot) = par(I.k124) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.k124) * X(I.Prot);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = - par(I.kd124) - par(I.kd125);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.k125) * X(I.Proti);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot,I.Proti) = par(I.k125) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP);

DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) = par(I.kd144);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) = - par(I.k144) * X(I.ERK_PP);
DF(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg,I.ERK_PP) = - par(I.k144) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg);
%%% -----------------------------------------------------------------------
%%% Equ. 42-43: Prot and Proti

DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2) = - par(I.k4) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Prot) = par(I.kd4);
DF(I.Prot,I.Prot) = - par(I.k4) * X(I.EGF_EGFRa_2_GAP_Grb2) -par(I.k106) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k109) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) ...
    - par(I.k112) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) - par(I.k115) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) ...
    - par(I.k118) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.k121) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) ...
    - par(I.k124) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) - par(I.kd15);
DF(I.Prot,I.Proti) = par(I.k15);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.k106) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = par(I.kd106);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = - par(I.k109) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = par(I.kd109);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = - par(I.k112) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = par(I.kd112);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2) = - par(I.k115) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = par(I.kd115);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.k118) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = par(I.kd118);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = - par(I.k121) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = par(I.kd121);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = - par(I.k124) * X(I.Prot);
DF(I.Prot,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = par(I.kd124);

DF(I.Proti,I.EGF_EGFRa_2_GAP_Grb2_Prot) =  par(I.kd5);
DF(I.Proti,I.Prot) = par(I.kd15);
DF(I.Proti,I.Proti) = - par(I.k15) - par(I.k5) * X(I.EGF_EGFRia_2_GAP_Grb2) ...
    - par(I.k107) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k110) *  X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP)...
    - par(I.k113) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP)  - par(I.k116) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2) ...
    - par(I.k119) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.k122) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP)...
    - par(I.k122) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) - par(I.k125) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Grb2) = - par(I.k5) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot) = par(I.kd107) ;
DF(I.Proti,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot) = par(I.kd110);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot) = par(I.kd113);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = par(I.kd116);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = par(I.kd119);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = par(I.kd122);
DF(I.Proti,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = par(I.kd125);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.k107) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = - par(I.k110) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = - par(I.k113) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Shca_Grb2) = - par(I.k116) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.k119) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = - par(I.k122) * X(I.Proti);
DF(I.Proti,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = - par(I.k125) * X(I.Proti);
%%% ----------------------------------------------------------------------
%%% Equ 44-49: Pathway including Shc internalised forms

DF(I.EGF_EGFRia_2_GAP_Shc,I.Shc) = par(I.k69) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Shc,I.EGF_EGFRia_2_GAP) = par(I.k69) * X(I.Shc);
DF(I.EGF_EGFRia_2_GAP_Shc,I.EGF_EGFRia_2_GAP_Shc) = - par(I.kd69) - par(I.k70) - par(I.kd103) -  par(I.k137);
DF(I.EGF_EGFRia_2_GAP_Shc,I.EGF_EGFRia_2_GAP_Shc_deg) = par(I.kd137);
DF(I.EGF_EGFRia_2_GAP_Shc,I.EGF_EGFRia_2_GAP_Shca) = par(I.kd70);
DF(I.EGF_EGFRia_2_GAP_Shc,I.EGF_EGFRa_2_GAP_Shc) = par(I.k103);

DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shc) = par(I.k70);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shca) = - par(I.kd70) - par(I.k71) * X(I.Grb2) ...
    - par(I.kd81) - par(I.k83) * X(I.Grb2_Sos) - par(I.kd104) - par(I.k138);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shc_deg) = par(I.kd138);
DF(I.EGF_EGFRia_2_GAP_Shca,I.Grb2) = - par(I.k71) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.kd71);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP) = par(I.k81) * X(I.Shca);
DF(I.EGF_EGFRia_2_GAP_Shca,I.Shca) = par(I.k81) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Shca,I.Grb2_Sos) = - par(I.k83) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd83);
DF(I.EGF_EGFRia_2_GAP_Shca,I.EGF_EGFRa_2_GAP_Shca) = par(I.k104);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.Grb2) = par(I.k71) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca) = par(I.k71) * X(I.Grb2);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2) = - par(I.kd71) - par(I.k72) * X(I.Sos) ...
    - par(I.kd82) - par(I.kd114) - par(I.k139)...
    - par(I.k116) * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2_deg) = par(I.kd139);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.Proti) = - par(I.k116) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.Sos) = - par(I.k72) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd72);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRia_2_GAP) = par(I.k82) * X(I.Shca_Grb2);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.Shca_Grb2) = par(I.k82) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.k114);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot) = par(I.kd116);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.k72) * X(I.Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Sos) = par(I.k72) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.kd72) - par(I.k73) * X(I.Ras_GDP) ...
    - par(I.k74) * X(I.Rasi_GTP) - par(I.k77) * X(I.Rasi_GTPa) ...
    - par(I.k78) * X(I.Ras_GDP) - par(I.kd79) ...
    - par(I.kd83) - par(I.kd117) - par(I.k129) * X(I.ERKi_PP) ...
    - par(I.k140) - par(I.k119)* X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) = par(I.kd140);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.ERKi_PP) = - par(I.k129) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) = par(I.kd129);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Ras_GDP) = - par(I.k73) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.k78) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd73) + par(I.kd74);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Rasi_GTP) = - par(I.k74) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Rasi_GTPa) = - par(I.k77) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd77) + par(I.kd78);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP) =  par(I.k79) * X(I.Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Shca_Grb2_Sos) = par(I.k79) * X(I.EGF_EGFRia_2_GAP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Grb2_Sos) = par(I.k83) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca) = par(I.k83) * X(I.Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.k117);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.Proti) = - par(I.k119) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot) = par(I.kd119);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Ras_GDP) = par(I.k73) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.k73) * X(I.Ras_GDP) + par(I.k74) * X(I.Rasi_GTP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = - par(I.kd73) - par(I.kd74) - par(I.kd120) - par(I.k141)...
    - par(I.k122)  * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Rasi_GTP) = par(I.k74) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.k120);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = par(I.kd122);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.Proti) = - par(I.k122) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg) = par(I.kd141);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Rasi_GTPa) = par(I.k77) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.k77) * X(I.Rasi_GTPa) + par(I.k78) * X(I.Ras_GDP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = - par(I.kd77) - par(I.kd78) ...
    - par(I.kd123) - par(I.k142) ...
    - par(I.k125) * X(I.Proti);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Ras_GDP) = par(I.k78) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.Proti) = - par(I.k125) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.k123);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = par(I.kd125);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg) = par(I.kd142);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.k129) * X(I.ERKi_PP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP,I.ERKi_PP) = par(I.k129) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) + par(I.k147) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) = - par(I.kd129) - par(I.kd147);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) = par(I.k147) * X(I.ERKi_PP);

DF(I.EGF_EGFRia_2_GAP_Shc_deg,I.EGF_EGFRia_2_GAP_Shc_deg) = - par(I.kd137) - par(I.kd138);
DF(I.EGF_EGFRia_2_GAP_Shc_deg,I.EGF_EGFRia_2_GAP_Shc) = par(I.k137);
DF(I.EGF_EGFRia_2_GAP_Shc_deg,I.EGF_EGFRia_2_GAP_Shca) = par(I.k138);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) = par(I.kd147);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) = - par(I.kd139) - par(I.kd140) - par(I.k147) * X(I.ERKi_PP);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg,I.ERKi_PP) = - par(I.k147) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.k139);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.k140);

DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg) = - par(I.kd141) - par(I.kd142);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.k141);
DF(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.k142);

%%% ----------------------------------------------------------------------
%%% Equ 50-55: Common state variables in both pathways: GAP,Grb2,Sos,
%%% Grb2_Sos, Ras:GTP, Ras:GDP


DF(I.GAP,I.EGF_EGFRa_2) = - par(I.k8) * X(I.GAP) ;
DF(I.GAP,I.GAP) = - par(I.k8) * X(I.EGF_EGFRa_2) - par(I.k14) * X(I.EGF_EGFRia_2);
DF(I.GAP,I.EGF_EGFRa_2_GAP) = par(I.kd8);
DF(I.GAP,I.EGF_EGFRia_2) = - par(I.k14) * X(I.GAP);
DF(I.GAP,I.EGF_EGFRia_2_GAP) = par(I.kd14);

DF(I.Grb2,I.Grb2) = - par(I.k16) * X(I.EGF_EGFRa_2_GAP) - par(I.k24) * X(I.EGF_EGFRa_2_GAP_Shca) ...
    - par(I.k35) * X(I.Sos) - par(I.k38) * X(I.Shca) ...
    - par(I.k63) * X(I.EGF_EGFRia_2_GAP) - par(I.k71) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.Grb2,I.EGF_EGFRa_2_GAP) = - par(I.k16) * X(I.Grb2);
DF(I.Grb2,I.EGF_EGFRa_2_GAP_Grb2) = par(I.kd16);
DF(I.Grb2,I.EGF_EGFRa_2_GAP_Shca) = - par(I.k24) * X(I.Grb2);
DF(I.Grb2,I.EGF_EGFRa_2_GAP_Shca_Grb2) = par(I.kd24);
DF(I.Grb2,I.Sos) = - par(I.k35) * X(I.Grb2);
DF(I.Grb2,I.Grb2_Sos) = par(I.kd35);
DF(I.Grb2,I.Shca) = - par(I.k38) * X(I.Grb2);
DF(I.Grb2,I.Shca_Grb2) = par(I.kd38);
DF(I.Grb2,I.EGF_EGFRia_2_GAP) = - par(I.k63) * X(I.Grb2);
DF(I.Grb2,I.EGF_EGFRia_2_GAP_Grb2) = par(I.kd63);
DF(I.Grb2,I.EGF_EGFRia_2_GAP_Shca) = - par(I.k71) * X(I.Grb2);
DF(I.Grb2,I.EGF_EGFRia_2_GAP_Shca_Grb2) = par(I.kd71);

DF(I.Sos,I.Sos) = - par(I.k17) * X(I.EGF_EGFRa_2_GAP_Grb2) - par(I.k25) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2) ...
    - par(I.k35) * X(I.Grb2) - par(I.k40) * X(I.Shca_Grb2) ...
    - par(I.k64) * X(I.EGF_EGFRia_2_GAP_Grb2) - par(I.k72) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2)...
    - par(I.k130) * X(I.ERK_PP) - par(I.k131) * X(I.ERKi_PP);
DF(I.Sos,I.EGF_EGFRa_2_GAP_Grb2) = - par(I.k17) * X(I.Sos);
DF(I.Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.kd17);
DF(I.Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2) = - par(I.k25) * X(I.Sos);
DF(I.Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd25);
DF(I.Sos,I.Grb2) = - par(I.k35) * X(I.Sos);
DF(I.Sos,I.Grb2_Sos) = par(I.kd35);
DF(I.Sos,I.ERK_PP) = - par(I.k130) * X(I.Sos);
DF(I.Sos,I.ERKi_PP) = - par(I.k131) * X(I.Sos);
DF(I.Sos,I.Sos_ERK_PP) = par(I.kd130);
DF(I.Sos,I.Sos_ERKi_PP) = par(I.kd131);
DF(I.Sos,I.Shca_Grb2) = - par(I.k40) * X(I.Sos);
DF(I.Sos,I.Shca_Grb2_Sos) = par(I.kd40);
DF(I.Sos,I.EGF_EGFRia_2_GAP_Grb2) = - par(I.k64) * X(I.Sos);
DF(I.Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.kd64);
DF(I.Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2) = - par(I.k72) * X(I.Sos);
DF(I.Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd72);

DF(I.Sos_ERK_PP,I.Sos) = par(I.k130) * X(I.ERK_PP);
DF(I.Sos_ERK_PP,I.ERK_PP) = par(I.k130) * X(I.Sos) + par(I.k145) * X(I.Sosi);
DF(I.Sos_ERK_PP,I.Sosi) = par(I.k145) * X(I.ERK_PP);
DF(I.Sos_ERK_PP,I.Sos_ERK_PP) = - par(I.kd130) - par(I.kd145);

DF(I.Sos_ERKi_PP,I.Sos) =  par(I.k131) * X(I.ERKi_PP);
DF(I.Sos_ERKi_PP,I.ERKi_PP) = par(I.k131) * X(I.Sos) + par(I.k148) * X(I.Sosi);
DF(I.Sos_ERKi_PP,I.Sos_ERKi_PP) = - par(I.kd131) - par(I.kd148);

DF(I.Sosi,I.Sosi) = - par(I.k145) * X(I.ERK_PP) - par(I.k148) * X(I.ERKi_PP);
DF(I.Sosi,I.ERK_PP) = - par(I.k145) * X(I.Sosi);
DF(I.Sosi,I.ERKi_PP) = - par(I.k148) * X(I.Sosi);
DF(I.Sosi,I.Sos_ERK_PP) = par(I.kd145);
DF(I.Sosi,I.Sos_ERKi_PP) = par(I.kd148);

DF(I.Grb2_Sos,I.Grb2_Sos) = - par(I.k33) * X(I.Shca) - par(I.k34) * X(I.EGF_EGFRa_2_GAP) ...
    - par(I.kd35) - par(I.k41) * X(I.EGF_EGFRa_2_GAP_Shca) ...
    - par(I.k80) * X(I.EGF_EGFRia_2_GAP) - par(I.k83) * X(I.EGF_EGFRia_2_GAP_Shca);
DF(I.Grb2_Sos,I.Shca) = - par(I.k33)* X(I.Grb2_Sos);
DF(I.Grb2_Sos,I.Shca_Grb2_Sos) = par(I.kd33);
DF(I.Grb2_Sos,I.EGF_EGFRa_2_GAP) = - par(I.k34) * X(I.Grb2_Sos);
DF(I.Grb2_Sos,I.EGF_EGFRa_2_GAP_Grb2_Sos) = par(I.kd34);
DF(I.Grb2_Sos,I.Sos) = par(I.k35) * X(I.Grb2);
DF(I.Grb2_Sos,I.Grb2) = par(I.k35) * X(I.Sos);
DF(I.Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca) = - par(I.k41) * X(I.Grb2_Sos);
DF(I.Grb2_Sos,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = par(I.kd41);
DF(I.Grb2_Sos,I.EGF_EGFRia_2_GAP) = - par(I.k80) * X(I.Grb2_Sos);
DF(I.Grb2_Sos,I.EGF_EGFRia_2_GAP_Grb2_Sos) = par(I.kd80);
DF(I.Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca) = - par(I.k83) * X(I.Grb2_Sos);
DF(I.Grb2_Sos,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = par(I.kd83);

DF(I.Ras_GTP,I.Ras_GTP) = - par(I.k19) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k27) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) ...
    - par(I.k28) * X(I.Raf);
DF(I.Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.k19) * X(I.Ras_GTP);
DF(I.Ras_GTP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd19);
DF(I.Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.k27) * X(I.Ras_GTP);
DF(I.Ras_GTP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd27);
DF(I.Ras_GTP,I.Raf) = -par(I.k28) * X(I.Ras_GTP);
DF(I.Ras_GTP,I.Raf_Ras_GTP) = par(I.kd28);

DF(I.Ras_GDP,I.Ras_GDP) = - par(I.k18) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k21) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) ...
    - par(I.k26) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) - par(I.k31) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) ...
    - par(I.k65) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k68) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos)...
    - par(I.k73) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) - par(I.k78) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.k18) * X(I.Ras_GDP) - par(I.k21) * X(I.Ras_GDP);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd18);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd21);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.k26) * X(I.Ras_GDP) - par(I.k31) * X(I.Ras_GDP);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd26);
DF(I.Ras_GDP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd31);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.k65) * X(I.Ras_GDP) - par(I.k68) * X(I.Ras_GDP);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd65);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd68);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.k73) * X(I.Ras_GDP) - par(I.k78) * X(I.Ras_GDP);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd73);
DF(I.Ras_GDP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd78);
%%% -----------------------------------------------------------------------
%%% Equ 56-66: Raf and Ras dynamic with internalised forms
%%%

DF(I.Raf,I.Raf) = - par(I.k28) * X(I.Ras_GTP) - par(I.k75) * X(I.Rasi_GTP) - par(I.k43) * X(I.Phosphatase1);
DF(I.Raf,I.Ras_GTP) = - par(I.k28) * X(I.Raf);
DF(I.Raf,I.Raf_Ras_GTP) = par(I.kd28);
DF(I.Raf,I.Rafa_Pase) = par(I.kd43);
DF(I.Raf,I.Rasi_GTP) = - par(I.k75) * X(I.Raf);
DF(I.Raf,I.Rafi_Rasi_GTP) = par(I.kd75);
DF(I.Raf,I.Rafia_Pase) = par(I.kd85);
DF(I.Raf,I.Rafia) = - par(I.k85) * X(I.Phosphatase1);
DF(I.Raf,I.Phosphatase1) = - par(I.k85) * X(I.Rafia) - par(I.k43) * X(I.Raf);

DF(I.Ras_GTPa,I.Ras_GTPa) = - par(I.k20) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k29) * X(I.Rafa) ...
    - par(I.k30) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos);
DF(I.Ras_GTPa,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.k20) * X(I.Ras_GTPa);
DF(I.Ras_GTPa,I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd20);
DF(I.Ras_GTPa,I.Rafa) = - par(I.k29) * X(I.Ras_GTPa);
DF(I.Ras_GTPa,I.Raf_Ras_GTP) = par(I.kd29);
DF(I.Ras_GTPa,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.k30) * X(I.Ras_GTPa);
DF(I.Ras_GTPa,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd30);

DF(I.Raf_Ras_GTP,I.Ras_GTP) = par(I.k28) * X(I.Raf);
DF(I.Raf_Ras_GTP,I.Raf) = par(I.k28) * X(I.Ras_GTP);
DF(I.Raf_Ras_GTP,I.Raf_Ras_GTP) = - par(I.kd28) - par(I.kd29);
DF(I.Raf_Ras_GTP,I.Ras_GTPa) = par(I.k29) * X(I.Rafa);
DF(I.Raf_Ras_GTP,I.Rafa) = par(I.k29) * X(I.Ras_GTPa);

DF(I.Rafa,I.Raf_Ras_GTP) = par(I.kd29);
DF(I.Rafa,I.Ras_GTPa) = - par(I.k29) * X(I.Rafa);
DF(I.Rafa,I.Rafa) = - par(I.k29) * X(I.Ras_GTPa) - par(I.k42) * X(I.Phosphatase1) ...
    - par(I.k44) * X(I.MEK) - par(I.k46) * X(I.MEK_P) - par(I.k45) * X(I.MEK_P)...
    - par(I.k47) * X(I.MEK_PP);
DF(I.Rafa,I.Phosphatase1) = - par(I.k42) * X(I.Rafa);
DF(I.Rafa,I.Rafa_Pase) = par(I.kd42);
DF(I.Rafa,I.MEK) = - par(I.k44) * X(I.Rafa);
DF(I.Rafa,I.MEK_Rafa) = par(I.kd44) + par(I.kd45);
DF(I.Rafa,I.MEK_P) = - par(I.k46) * X(I.Rafa) - par(I.k45) * X(I.Rafa);
DF(I.Rafa,I.MEK_PP) = - par(I.k47) * X(I.Rafa);
DF(I.Rafa,I.MEK_P_Rafa) = par(I.kd46) + par(I.kd47);

DF(I.Rafa_Pase,I.Phosphatase1) = par(I.k42) * X(I.Rafa) + par(I.k43) * X(I.Raf);
DF(I.Rafa_Pase,I.Rafa) = par(I.k42) * X(I.Phosphatase1);
DF(I.Rafa_Pase,I.Raf) = par(I.k43) * X(I.Phosphatase1);
DF(I.Rafa_Pase,I.Rafa_Pase) = - par(I.kd42) - par(I.kd43);

DF(I.Phosphatase1,I.Rafa) = - par(I.k42) * X(I.Phosphatase1);
DF(I.Phosphatase1,I.Phosphatase1) = - par(I.k42) * X(I.Rafa) - par(I.k84) * X(I.Rafia) ...
    - par(I.k85) * X(I.Rafia) - par(I.k43) * X(I.Raf);
DF(I.Phosphatase1,I.Rafa_Pase) = par(I.kd42) + par(I.kd43);
DF(I.Phosphatase1,I.Raf) = - par(I.k43) * X(I.Phosphatase1);
DF(I.Phosphatase1,I.Rafia) = - par(I.k84) * X(I.Phosphatase1) - par(I.k85) * X(I.Phosphatase1);
DF(I.Phosphatase1,I.Rafia_Pase) = par(I.kd84) + par(I.kd85);

DF(I.Rasi_GTP,I.Rasi_GTP) = - par(I.k66) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k74) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) ...
    - par(I.k75) * X(I.Raf);
DF(I.Rasi_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.k66) * X(I.Rasi_GTP);
DF(I.Rasi_GTP,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP) = par(I.kd66);
DF(I.Rasi_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP) = par(I.kd74);
DF(I.Rasi_GTP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.k74) * X(I.Rasi_GTP);
DF(I.Rasi_GTP,I.Raf) = - par(I.k75) * X(I.Rasi_GTP);
DF(I.Rasi_GTP,I.Rafi_Rasi_GTP) = par(I.kd75);

DF(I.Rafi_Rasi_GTP,I.Rasi_GTP) = par(I.k75) * X(I.Raf);
DF(I.Rafi_Rasi_GTP,I.Raf) = par(I.k75) * X(I.Rasi_GTP);
DF(I.Rafi_Rasi_GTP,I.Rafi_Rasi_GTP) = - par(I.kd75) - par(I.kd76);
DF(I.Rafi_Rasi_GTP,I.Rasi_GTPa) = par(I.k76) * X(I.Rafia);
DF(I.Rafi_Rasi_GTP,I.Rafia) = par(I.k76) * X(I.Rasi_GTPa);

DF(I.Rasi_GTPa,I.Rasi_GTPa) = - par(I.k67) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) ...
    - par(I.k76) * X(I.Rafia) - par(I.k77) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos);
DF(I.Rasi_GTPa,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.k67) * X(I.Rasi_GTPa);
DF(I.Rasi_GTPa,I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP) = par(I.kd67);
DF(I.Rasi_GTPa,I.Rafi_Rasi_GTP) = par(I.kd76);
DF(I.Rasi_GTPa,I.Rafia) = - par(I.k76) * X(I.Rasi_GTPa);
DF(I.Rasi_GTPa,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.k77) * X(I.Rasi_GTPa);
DF(I.Rasi_GTPa,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP) = par(I.kd77);

DF(I.Rafia,I.Rafi_Rasi_GTP) = par(I.kd76);
DF(I.Rafia,I.Rasi_GTPa) = - par(I.k76) * X(I.Rafia);
DF(I.Rafia,I.Rafia) = - par(I.k76) * X(I.Rasi_GTPa) - par(I.k84) * X(I.Phosphatase1) ...
    - par(I.k86) * X(I.MEK) - par(I.k88) * X(I.MEKi_P) - par(I.k85) * X(I.Phosphatase1)...
    - par(I.k89) * X(I.MEKi_PP)- par(I.k87) * X(I.MEKi_P);
DF(I.Rafia,I.Phosphatase1) = - par(I.k84) * X(I.Rafia) - par(I.k85) * X(I.Rafia);
DF(I.Rafia,I.Rafia_Pase) = par(I.kd84);
DF(I.Rafia,I.MEK) = - par(I.k86) * X(I.Rafia);
DF(I.Rafia,I.MEKi_Rafia) = par(I.kd86) + par(I.kd87);
DF(I.Rafia,I.MEKi_P) = - par(I.k88) * X(I.Rafia) - par(I.k87) * X(I.Rafia);
DF(I.Rafia,I.MEKi_PP) = - par(I.k89) * X(I.Rafia);
DF(I.Rafia,I.MEKi_P_Rafia) = par(I.kd88) + par(I.kd89);

DF(I.Rafia_Pase,I.Phosphatase1) = par(I.k84) * X(I.Rafia) + par(I.k85) * X(I.Rafia);
DF(I.Rafia_Pase,I.Rafia) = par(I.k84) * X(I.Phosphatase1) + par(I.k85) * X(I.Phosphatase1);
DF(I.Rafia_Pase,I.Rafia_Pase) = - par(I.kd84) - par(I.kd85);


%%% -----------------------------------------------------------------------
%%% Equ 67-80: MEK dynamic with internalised forms
%%%

DF(I.MEK,I.Rafa) = - par(I.k44) * X(I.MEK);
DF(I.MEK,I.MEK) = - par(I.k44) * X(I.Rafa) - par(I.k86) * X(I.Rafia)...
    - par(I.k93) * X(I.Phosphatase2) - par(I.k51) * X(I.Phosphatase2);
DF(I.MEK,I.MEK_Rafa) = par(I.kd44);
DF(I.MEK,I.MEK_P_Pase2) =  par(I.kd51);
DF(I.MEK,I.Rafia) = - par(I.k86) * X(I.MEK);
DF(I.MEK,I.MEKi_Rafia) = par(I.kd86);
DF(I.MEK,I.MEKi_P_Pase2i) = par(I.kd93);
DF(I.MEK,I.Phosphatase2) = - par(I.k93) * X(I.MEK) - par(I.k51) * X(I.MEK);

DF(I.MEK_Rafa,I.MEK) = par(I.k44) * X(I.Rafa);
DF(I.MEK_Rafa,I.Rafa) = par(I.k44) * X(I.MEK) + par(I.k45) * X(I.MEK_P);
DF(I.MEK_Rafa,I.MEK_Rafa) = - par(I.kd44) - par(I.kd45);

DF(I.MEK_P,I.MEK_Rafa) = par(I.kd45);
DF(I.MEK_P,I.MEK_P) = - par(I.k46) * X(I.Rafa) - par(I.k50) * X(I.Phosphatase2) - par(I.k45) * X(I.Rafa) - par(I.k49) * X(I.Phosphatase2);
DF(I.MEK_P,I.Rafa) = - par(I.k46) * X(I.MEK_P) - par(I.k45) * X(I.MEK_P);
DF(I.MEK_P,I.MEK_P_Rafa) = par(I.kd46);
DF(I.MEK_P,I.MEK_PP_Pase2) = par(I.kd49);
DF(I.MEK_P,I.Phosphatase2) = - par(I.k50) * X(I.MEK_P) - par(I.k49) * X(I.MEK_P);
DF(I.MEK_P,I.MEK_P_Pase2) = par(I.kd50);

DF(I.MEK_P_Rafa,I.MEK_P) = par(I.k46) * X(I.Rafa);
DF(I.MEK_P_Rafa,I.Rafa) = par(I.k46) * X(I.MEK_P) +par(I.k47) * X(I.MEK_PP);
DF(I.MEK_P_Rafa,I.MEK_P_Rafa) = - par(I.kd46) - par(I.kd47);
DF(I.MEK_P_Rafa,I.MEK_PP) = par(I.k47) * X(I.Rafa);

DF(I.MEK_PP,I.MEK_P_Rafa) = par(I.kd47);
DF(I.MEK_PP,I.MEK_PP) = - par(I.k48) * X(I.Phosphatase2) - par(I.k52) * X(I.ERK) ...
    - par(I.k54) * X(I.ERK_P) - par(I.k47) * X(I.Rafa)...
    - par(I.k53) * X(I.ERK_P) - par(I.k55) * X(I.ERK_PP);
DF(I.MEK_PP,I.Phosphatase2) = - par(I.k48) * X(I.MEK_PP);
DF(I.MEK_PP,I.MEK_PP_Pase2) = par(I.kd48);
DF(I.MEK_PP,I.Rafa) = - par(I.k47) * X(I.MEK_PP);
DF(I.MEK_PP,I.ERK) = - par(I.k52) * X(I.MEK_PP);
DF(I.MEK_PP,I.ERK_MEK_PP) = par(I.kd52) + par(I.kd53);
DF(I.MEK_PP,I.ERK_P) = - par(I.k54) * X(I.MEK_PP) - par(I.k53) * X(I.MEK_PP);
DF(I.MEK_PP,I.ERK_PP) = - par(I.k55) * X(I.MEK_PP);
DF(I.MEK_PP,I.ERK_P_MEK_PP) = par(I.kd54) + par(I.kd55);

DF(I.MEK_PP_Pase2,I.MEK_PP) = par(I.k48) * X(I.Phosphatase2);
DF(I.MEK_PP_Pase2,I.Phosphatase2) = par(I.k48) * X(I.MEK_PP) + par(I.k49) * X(I.MEK_P);
DF(I.MEK_PP_Pase2,I.MEK_PP_Pase2) = - par(I.kd48) - par(I.kd49);
DF(I.MEK_PP_Pase2,I.MEK_P) = par(I.k49) * X(I.Phosphatase2);

DF(I.Phosphatase2,I.MEK_PP) = - par(I.k48) * X(I.Phosphatase2);
DF(I.Phosphatase2,I.Phosphatase2) = - par(I.k48) * X(I.MEK_PP) - par(I.k50) * X(I.MEK_P) ...
    - par(I.k90) * X(I.MEKi_PP) - par(I.k92) * X(I.MEKi_P)...
    - par(I.k91) * X(I.MEKi_PP) - par(I.k93) * X(I.MEK);
DF(I.Phosphatase2,I.MEK_PP_Pase2) = par(I.kd48) + par(I.k49);
DF(I.Phosphatase2,I.MEK_P_Pase2) = par(I.kd50) + par(I.k51);
DF(I.Phosphatase2,I.MEK_P) = - par(I.k50) * X(I.Phosphatase2);
DF(I.Phosphatase2,I.MEK) = - par(I.k93) * X(I.Phosphatase2);
DF(I.Phosphatase2,I.MEKi_PP) = - par(I.k90) * X(I.Phosphatase2) - par(I.k91) * X(I.Phosphatase2);
DF(I.Phosphatase2,I.MEKi_PP_Pase2i) = par(I.kd90) + par(I.kd91);
DF(I.Phosphatase2,I.MEKi_P) = - par(I.k92) * X(I.Phosphatase2);
DF(I.Phosphatase2,I.MEKi_P_Pase2i) = par(I.kd92) + par(I.kd93);

DF(I.MEK_P_Pase2,I.Phosphatase2) = par(I.k50) * X(I.MEK_P) + par(I.k51) * X(I.MEK);
DF(I.MEK_P_Pase2,I.MEK_P) = par(I.k50) * X(I.Phosphatase2);
DF(I.MEK_P_Pase2,I.MEK) = par(I.k51) * X(I.Phosphatase2);
DF(I.MEK_P_Pase2,I.MEK_P_Pase2) = - par(I.kd50) - par(I.kd51);

DF(I.MEKi_Rafia,I.MEK) = par(I.k86) * X(I.Rafia);
DF(I.MEKi_Rafia,I.MEKi_P) = par(I.k87) * X(I.Rafia);
DF(I.MEKi_Rafia,I.Rafia) = par(I.k86) * X(I.MEK) + par(I.k87) * X(I.MEKi_P);
DF(I.MEKi_Rafia,I.MEKi_Rafia) = - par(I.kd86) - par(I.kd87);

DF(I.MEKi_P,I.MEKi_Rafia) = par(I.kd87);
DF(I.MEKi_P,I.MEKi_P) = - par(I.k88) * X(I.Rafia) - par(I.k92) * X(I.Phosphatase2) - par(I.k87) * X(I.Rafia);
DF(I.MEKi_P,I.Rafia) = - par(I.k88) * X(I.MEKi_P) - par(I.k87) * X(I.MEKi_P);
DF(I.MEKi_P,I.MEKi_P_Rafia) = par(I.kd88);
DF(I.MEKi_P,I.MEKi_PP_Pase2i) = par(I.kd91);
DF(I.MEKi_P,I.MEKi_PP) = - par(I.k91) * X(I.Phosphatase2);
DF(I.MEKi_P,I.Phosphatase2) = - par(I.k92) * X(I.MEKi_P) - par(I.k91) * X(I.MEKi_PP);
DF(I.MEKi_P,I.MEKi_P_Pase2i) = par(I.kd92);

DF(I.MEKi_P_Rafia,I.Rafia) = par(I.k88) * X(I.MEKi_P) + par(I.k89) * X(I.MEKi_PP) ;
DF(I.MEKi_P_Rafia,I.MEKi_P) = par(I.k88) * X(I.Rafia);
DF(I.MEKi_P_Rafia,I.MEKi_PP) = par(I.k89) * X(I.Rafia);
DF(I.MEKi_P_Rafia,I.MEKi_P_Rafia) = - par(I.kd88) - par(I.kd89);

DF(I.MEKi_PP,I.MEKi_P_Rafia) = par(I.kd89);
DF(I.MEKi_PP,I.MEKi_PP) = - par(I.k90) * X(I.Phosphatase2) -  par(I.k94) * X(I.ERK) ...
    - par(I.k96) * X(I.ERKi_P) - par(I.k89) * X(I.Rafia) ...
    - par(I.k95) * X(I.ERKi_P) - par(I.k97) * X(I.ERKi_PP);
DF(I.MEKi_PP,I.Phosphatase2) = - par(I.k90) * X(I.MEKi_PP);
DF(I.MEKi_PP,I.MEKi_PP_Pase2i) =  par(I.kd90);
DF(I.MEKi_PP,I.Rafia) = - par(I.k87) * X(I.MEKi_P) - par(I.k89) * X(I.MEKi_PP);
DF(I.MEKi_PP,I.ERK) = - par(I.k94) * X(I.MEKi_PP);
DF(I.MEKi_PP,I.ERKi_MEKi_PP) = par(I.kd94) + par(I.kd95);
DF(I.MEKi_PP,I.ERKi_P) = - par(I.k96) * X(I.MEKi_PP) - par(I.k95) * X(I.MEKi_PP);
DF(I.MEKi_PP,I.ERKi_PP) = - par(I.k97) * (I.MEKi_PP);
DF(I.MEKi_PP,I.ERKi_P_MEKi_PP) = par(I.kd96) + par(I.kd97);

DF(I.MEKi_PP_Pase2i,I.MEKi_PP) = par(I.k90) * X(I.Phosphatase2) + par(I.k91) * X(I.Phosphatase2);
DF(I.MEKi_PP_Pase2i,I.Phosphatase2) = par(I.k90) * X(I.MEKi_PP) + par(I.k91) * X(I.MEKi_PP);
DF(I.MEKi_PP_Pase2i,I.MEKi_PP_Pase2i) = - par(I.kd90) - par(I.kd91);

DF(I.MEKi_P_Pase2i,I.Phosphatase2) = par(I.k92) * X(I.MEKi_P) + par(I.k93) * X(I.MEK);
DF(I.MEKi_P_Pase2i,I.MEKi_P) = par(I.k92) * X(I.Phosphatase2);
DF(I.MEKi_P_Pase2i,I.MEK) = par(I.k93) * X(I.Phosphatase2);
DF(I.MEKi_P_Pase2i,I.MEKi_P_Pase2i) = - par(I.kd92) - par(I.kd93);

%%% -----------------------------------------------------------------------
%%% Equ 81-94: ERK dynamic with internalised forms
%%%

DF(I.ERK,I.MEK_PP) = - par(I.k52) * X(I.ERK);
DF(I.ERK,I.ERK) = - par(I.k52) * X(I.MEK_PP) - par(I.k94) * X(I.MEKi_PP)...
    - par(I.k59) * X(I.Phosphatase3)- par(I.k101) * X(I.Phosphatase3);
DF(I.ERK,I.ERK_MEK_PP) = par(I.kd52);
DF(I.ERK,I.ERK_P_Pase3) = par(I.kd59);
DF(I.ERK,I.MEKi_PP) = - par(I.k94) * X(I.ERK);
DF(I.ERK,I.ERKi_MEKi_PP) = par(I.kd94);
DF(I.ERK,I.ERKi_P_Pase3i) = par(I.kd101);
DF(I.ERK,I.Phosphatase3) = - par(I.k59) * X(I.ERK) - par(I.k101) * X(I.ERK);

DF(I.ERK_MEK_PP,I.ERK) = par(I.k52) * X(I.MEK_PP);
DF(I.ERK_MEK_PP,I.ERK_P) = par(I.k53) * X(I.MEK_PP);
DF(I.ERK_MEK_PP,I.MEK_PP) = par(I.k52) * X(I.ERK) + par(I.k53) * X(I.ERK_P);
DF(I.ERK_MEK_PP,I.ERK_MEK_PP) = - par(I.kd52) - par(I.kd53);

DF(I.ERK_P,I.ERK_MEK_PP) = par(I.kd53);
DF(I.ERK_P,I.MEK_PP) = - par(I.k54) * X(I.ERK_P) - par(I.k53) * X(I.ERK_P);
DF(I.ERK_P,I.ERK_P) = - par(I.k54) * X(I.MEK_PP) - par(I.k58) * X(I.Phosphatase3) - par(I.k53) * X(I.MEK_PP) - par(I.k57) * X(I.Phosphatase3);
DF(I.ERK_P,I.ERK_P_MEK_PP) = par(I.kd54);
DF(I.ERK_P,I.ERK_PP_Pase3) = par(I.kd57);
DF(I.ERK_P,I.ERK_P_Pase3) = par(I.kd58);
DF(I.ERK_P,I.Phosphatase3) = - par(I.k58) * X(I.ERK_P) - par(I.k57) * X(I.ERK_P);

DF(I.ERK_P_MEK_PP,I.MEK_PP) = par(I.k54) * X(I.ERK_P) + par(I.k55) * X(I.ERK_PP);
DF(I.ERK_P_MEK_PP,I.ERK_P) = par(I.k54) * X(I.MEK_PP);
DF(I.ERK_P_MEK_PP,I.ERK_PP) = par(I.k55) * X(I.MEK_PP);
DF(I.ERK_P_MEK_PP,I.ERK_P_MEK_PP) = - par(I.kd54) - par(I.kd55);

DF(I.ERK_PP,I.ERK_P_MEK_PP) = par(I.kd55);
DF(I.ERK_PP,I.ERK_PP) = - par(I.k56) * X(I.Phosphatase3) - par(I.k55) * X(I.MEK_PP)...
    - par(I.k126) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos) - par(I.k128) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos)...
    - par(I.k130) * X(I.Sos) - par(I.k143) * X(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) ...
    - par(I.k144) * X(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) - par(I.k145) * X(I.Sosi);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos) = - par(I.k126) * X(I.ERK_PP);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP) = par(I.kd126) + par(I.kd143);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos) = - par(I.k128) * X(I.ERK_PP);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP) =  par(I.kd128) + par(I.kd144);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Grb2_Sos_deg) = - par(I.k143) * X(I.ERK_PP);
DF(I.ERK_PP,I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg) = - par(I.k144) * X(I.ERK_PP);
DF(I.ERK_PP,I.Sos) = - par(I.k130) * X(I.ERK_PP);
DF(I.ERK_PP,I.Sosi) = - par(I.k145) * X(I.ERK_PP);
DF(I.ERK_PP,I.Sos_ERK_PP) = par(I.kd130) +  par(I.kd145);
DF(I.ERK_PP,I.MEK_PP) = - par(I.k55) * X(I.ERK_PP);
DF(I.ERK_PP,I.Phosphatase3) = - par(I.k56) * X(I.ERK_PP);
DF(I.ERK_PP,I.ERK_PP_Pase3) = par(I.kd56);

DF(I.Phosphatase3,I.Phosphatase3) = - par(I.k56) * X(I.ERK_PP) - par(I.k58) * X(I.ERK_P) ...
    - par(I.k98) * X(I.ERKi_PP) - par(I.k100) * X(I.ERKi_P)...
    - par(I.k99) * X(I.ERKi_P) - par(I.k101) * X(I.ERK)...
    - par(I.k59) * X(I.ERK) - par(I.k57) * X(I.ERK_P);
DF(I.Phosphatase3,I.ERK_PP) = - par(I.k56) * X(I.Phosphatase3);
DF(I.Phosphatase3,I.ERK_PP_Pase3) = par(I.kd56) + par(I.kd57);
DF(I.Phosphatase3,I.ERK_P) = - par(I.k58) * X(I.Phosphatase3) - par(I.k57) * X(I.Phosphatase3);
DF(I.Phosphatase3,I.ERK) = - par(I.k101) * X(I.Phosphatase3) - par(I.k59) * X(I.Phosphatase3);
DF(I.Phosphatase3,I.ERK_P_Pase3) = par(I.kd58) + par(I.kd59);
DF(I.Phosphatase3,I.ERKi_PP) = - par(I.k98) * X(I.Phosphatase3);
DF(I.Phosphatase3,I.ERKi_PP_Pase3i) = par(I.kd98) + par(I.kd99);
DF(I.Phosphatase3,I.ERKi_P) = - par(I.k100) * X(I.Phosphatase3) - par(I.k99) * X(I.Phosphatase3);
DF(I.Phosphatase3,I.ERKi_P_Pase3i) = par(I.kd100) + par(I.kd101);

DF(I.ERK_PP_Pase3,I.ERK_PP) = par(I.k56) * X(I.Phosphatase3);
DF(I.ERK_PP_Pase3,I.ERK_P) = par(I.k57) * X(I.Phosphatase3);
DF(I.ERK_PP_Pase3,I.Phosphatase3) = par(I.k56) * X(I.ERK_PP) + par(I.k57) * X(I.ERK_P);
DF(I.ERK_PP_Pase3,I.ERK_PP_Pase3) = - par(I.kd56) - par(I.kd57);

DF(I.ERK_P_Pase3,I.Phosphatase3) = par(I.k58) * X(I.ERK_P) + par(I.k59) * X(I.ERK);
DF(I.ERK_P_Pase3,I.ERK_P) = par(I.k58) * X(I.Phosphatase3);
DF(I.ERK_P_Pase3,I.ERK) = par(I.k59) * X(I.Phosphatase3);
DF(I.ERK_P_Pase3,I.ERK_P_Pase3) = - par(I.kd58) - par(I.kd59);

DF(I.ERKi_P,I.ERKi_MEKi_PP) = par(I.kd95);
DF(I.ERKi_P,I.MEKi_PP) = - par(I.k96) * X(I.ERKi_P) - par(I.k95) * X(I.ERKi_P);
DF(I.ERKi_P,I.ERKi_P) = - par(I.k96) * X(I.MEKi_PP) - par(I.k100) * X(I.Phosphatase3) ...
    - par(I.k95) * X(I.MEKi_PP) - par(I.k99) * X(I.Phosphatase3);
DF(I.ERKi_P,I.ERKi_P_MEKi_PP) = par(I.kd96);
DF(I.ERKi_P,I.ERKi_PP_Pase3i) = par(I.kd99);
DF(I.ERKi_P,I.Phosphatase3) = - par(I.k100) * X(I.ERKi_P) - par(I.k99) * X(I.ERKi_P);
DF(I.ERKi_P,I.ERKi_P_Pase3i) = par(I.kd100);

DF(I.ERKi_MEKi_PP,I.ERK) = par(I.k94) * X(I.MEKi_PP);
DF(I.ERKi_MEKi_PP,I.MEKi_PP) = par(I.k94) * X(I.ERK) + par(I.k95) * X(I.ERKi_P);
DF(I.ERKi_MEKi_PP,I.ERKi_P) = par(I.k95) * X(I.MEKi_PP);
DF(I.ERKi_MEKi_PP,I.ERKi_MEKi_PP) = - par(I.kd94) - par(I.kd95);

DF(I.ERKi_P_MEKi_PP,I.MEKi_PP) = par(I.k96) * X(I.ERKi_P) + par(I.k97) * X(I.ERKi_PP);
DF(I.ERKi_P_MEKi_PP,I.ERKi_P) = par(I.k96) * X(I.MEKi_PP);
DF(I.ERKi_P_MEKi_PP,I.ERKi_PP) = par(I.k97) * (I.MEKi_PP);
DF(I.ERKi_P_MEKi_PP,I.ERKi_P_MEKi_PP) = - par(I.kd96) - par(I.kd97);

DF(I.ERKi_PP,I.ERKi_P_MEKi_PP) = par(I.kd97);
DF(I.ERKi_PP,I.ERKi_PP) = - par(I.k98) * X(I.Phosphatase3) - par(I.k97) * (I.MEKi_PP)...
    - par(I.k127) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos) - par(I.k129) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos)...
    - par(I.k131) * X(I.Sos) - par(I.k146) * X(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) ...
    - par(I.k147) * X(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg)  - par(I.k148) * X(I.Sosi)...
    - par(I.k97) * (I.MEKi_PP);
DF(I.ERKi_PP,I.MEKi_PP) = - par(I.k95) * X(I.ERKi_P) - par(I.k97) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos) = - par(I.k127) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP) = par(I.kd127) + par(I.kd146);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos) = - par(I.k129) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP) = par(I.kd129) +  par(I.kd147);
DF(I.ERKi_PP,I.Sos) = - par(I.k131) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.Sos_ERKi_PP) = par(I.kd131) + par(I.kd148);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Grb2_Sos_deg) = - par(I.k146) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg) = - par(I.k147) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.Sosi) = - par(I.k148) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.Phosphatase3) = - par(I.k98) * X(I.ERKi_PP);
DF(I.ERKi_PP,I.ERKi_PP_Pase3i) = par(I.kd98);

DF(I.ERKi_PP_Pase3i,I.ERKi_PP) = par(I.k98) * X(I.Phosphatase3);
DF(I.ERKi_PP_Pase3i,I.ERKi_P) = par(I.k99) * X(I.Phosphatase3);
DF(I.ERKi_PP_Pase3i,I.Phosphatase3) = par(I.k98) * X(I.ERKi_PP) + par(I.k99) * X(I.ERKi_P);
DF(I.ERKi_PP_Pase3i,I.ERKi_PP_Pase3i) = - par(I.kd98) - par(I.kd99);

DF(I.ERKi_P_Pase3i,I.Phosphatase3) = par(I.k100) * X(I.ERKi_P) + par(I.k101) * X(I.ERK);
DF(I.ERKi_P_Pase3i,I.ERKi_P) = par(I.k100) * X(I.Phosphatase3);
DF(I.ERKi_P_Pase3i,I.ERK) = par(I.k101) * X(I.Phosphatase3);
DF(I.ERKi_P_Pase3i,I.ERKi_P_Pase3i) = - par(I.kd100) - par(I.kd101);

DF(I.AUC_ERK_PP,I.ERK_PP) = 1;

%%% -----------------------------------------------------------------------


end
