%%% Version: 19 Jan 2023
%%%
%%% X0 = <MODELNAME>_initialValues(model)
%%%
%%% This function assigns the initial values for all state variables 
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model
%%%
%%% Output : X0         initial values 
%%%
%%%
%%% References:
%%%    - Hornberg, J et al.
%%%     'Control of MAPK signaling: from complexity to what really matters'
%%%     Oncogene, Volume 24, 2005
%%%    - Schoeberl, Birgit; Eichler-Johnson, Claudia; Gilles, Ernst Dieter;
%%%      Mueller, Gertrazd
%%%      'Computational modeling of the dynamics of the MAP 
%%%      kinase casade activated by surface and internalized 
%%%      EGF receptors'
%%%      Nature Biotechnology Volume 20, 2002
%%%    - Corresponding Model from BioModels Database:
%%%      http://www.ebi.ac.uk/biomodels-main/BIOMD0000000019
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function X0 = Hornberg2005EGFRsignalling_initialvalues(model)

%%% assign model indexing
I  = model. I;

%%% initialise state values
X0 = NaN(I.nstates,1);

X0(I.EGF)                                        = 0;
X0(I.EGFR)                                       = 50000.0;
X0(I.EGF_EGFR)                                   = 0.0;
X0(I.EGF_EGFR_2)                                 = 0.0;
X0(I.EGF_EGFRa_2)                                = 0.0;
X0(I.EGFRi)                                      = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Prot)                  = 0.0;
X0(I.EGF_EGFRia_2)                               = 0.0;
X0(I.Proti)                                      = 0.0;
X0(I.EGF_EGFRi)                                  = 0.0;
X0(I.EGF_EGFRi_2)                                = 0.0;
X0(I.Prot)                                       = 81000.0;
X0(I.EGFideg)                                    = 0.0;      
X0(I.GAP)                                        = 12000.0;
X0(I.EGF_EGFRa_2_GAP)                            = 0.0;
X0(I.EGFi)                                       = 0.0;
X0(I.EGF_EGFRia_2_GAP)                           = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2)                      = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos)                  = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP)          = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP)          = 0.0;
X0(I.Grb2)                                       = 11000.0;   
X0(I.EGF_EGFRa_2_GAP_Grb2)                       = 0.0;
X0(I.Sos)                                        = 26300.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos)                   = 0.0;
X0(I.Ras_GDP)                                    = 72000.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP)           = 0.0;
X0(I.Ras_GTP)                                    = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP)           = 0.0;
X0(I.Grb2_Sos)                                   = 40000.0;
X0(I.Shc)                                        = 101000.0;
X0(I.EGF_EGFRa_2_GAP_Shc)                        = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca)                       = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2)                  = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos)              = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP)      = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP)      = 0.0;
X0(I.Shca_Grb2_Sos)                              = 0.0;
X0(I.Shca_Grb2)                                  = 0.0;
X0(I.Shca)                                       = 0.0;
X0(I.Raf)                                        = 40000.0;
X0(I.Raf_Ras_GTP)                                = 0.0;
X0(I.Ras_GTPa)                                   = 0.0;
X0(I.Phosphatase1)                               = 40000.0;
X0(I.Rafa)                                       = 0.0;
X0(I.Rafa_Pase)                                  = 0.0;
X0(I.MEK)                                        = 2.1e7; 
X0(I.MEK_Rafa)                                   = 0.0;
X0(I.MEK_P)                                      = 0.0;
X0(I.MEK_P_Rafa)                                 = 0.0;
X0(I.MEK_PP)                                     = 0.0;
X0(I.MEK_PP_Pase2)                               = 0.0;
X0(I.Phosphatase2)                               = 40000.0; 
X0(I.MEK_P_Pase2)                                = 0.0;
X0(I.ERK)                                        = 2.21e7; 
X0(I.ERK_MEK_PP)                                 = 0.0;
X0(I.ERK_P)                                      = 0.0;
X0(I.ERK_P_MEK_PP)                               = 0.0;
X0(I.ERK_PP)                                     = 0.0;
X0(I.Phosphatase3)                               = 1.0e7;
X0(I.ERK_PP_Pase3)                               = 0.0;
X0(I.ERK_P_Pase3)                                = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shc)                       = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca)                      = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2)                 = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos)             = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP)     = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP)     = 0.0;
X0(I.Rasi_GTP)                                   = 0.0;
X0(I.Rafi_Rasi_GTP)                              = 0.0;
X0(I.Rasi_GTPa)                                  = 0.0;
X0(I.Rafia)                                      = 0.0;
X0(I.Rafia_Pase)                                 = 0.0;
X0(I.MEKi_Rafia)                                 = 0.0;
X0(I.MEKi_P)                                     = 0.0;
X0(I.MEKi_P_Rafia)                               = 0.0;
X0(I.MEKi_PP)                                    = 0.0;
X0(I.MEKi_PP_Pase2i)                             = 0.0;
X0(I.MEKi_P_Pase2i)                              = 0.0;
X0(I.ERKi_MEKi_PP)                               = 0.0;
X0(I.ERKi_P)                                     = 0.0;
X0(I.ERKi_P_MEKi_PP)                             = 0.0;
X0(I.ERKi_PP)                                    = 0.0;
X0(I.ERKi_PP_Pase3i)                             = 0.0;
X0(I.ERKi_P_Pase3i)                              = 0.0;
X0(I.EGFRideg)                                   = 0.0;
X0(I.EGF_EGFRia_2_deg)                           = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_Prot)              = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot)      = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot)      = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Prot)             = 0.0;  
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot)         = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot) = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot) = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP)            = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP)          = 0.0; 
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP)       = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP)     = 0.0;
X0(I.EGF_EGFRa_2_GAP_Grb2_Sos_deg)               = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos_deg)              = 0.0;
X0(I.Sos_ERK_PP)                                 = 0.0;
X0(I.Sos_ERKi_PP)                                = 0.0;
X0(I.Sosi)                                       = 0.0;
X0(I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg)          = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg)         = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_deg)             = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_deg)                  = 0.0;
X0(I.AUC_ERK_PP)                                 = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg)     = 0.0;
X0(I.EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg)          = 0.0;
X0(I.EGF_EGFRia_2_GAP_deg)                       = 0.0;
X0(I.EGF_EGFRia_2_GAP_Shc_deg)                   = 0.0;

Schoeberl_etal_2002 = 0; 
if Schoeberl_etal_2002
    
    beep, pause(1), beep, pause(1), beep,
    X0(I.Grb2) = 5.1e4;
    X0(I.Sos) = 6.63e4;
    X0(I.Ras_GDP)= 1.14e7;
    X0(I.Grb2_Sos) = 0;
    X0(I.Shc) = 1.01e6;
    X0(I.MEK) = 2.2e7;
    X0(I.ERK) = 2.1e7;
   
end

end