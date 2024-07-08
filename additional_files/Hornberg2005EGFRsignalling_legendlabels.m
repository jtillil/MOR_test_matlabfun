%%% Version: 19 Jan 2023
%%%
%%% I  = <MODELNAME>_legendlabels(model)
%%%
%%% This function assigns the better readable labels to states
%%% 
%%% Input:  model       model structure
%%%                   
%%% Output: I           updated strucutre 
%%% 
%%%
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function I = Hornberg2005EGFRsignalling_legendlabels(I)

S = struct('EGF','EGF',...
    'EGFR','EGFR',...
    'EGF_EGFR', 'EGF-EGFR',...
    'EGF_EGFR_2','(EGF-EGFR)_2',...
    'EGF_EGFRa_2','(EGF-EGFR*)_2',...
    'EGFRi','EGFRi',...
    'EGF_EGFRa_2_GAP_Grb2_Prot','(EGF-EGFR*)_2-GAP-Grb2-Prot',...
    'EGF_EGFRia_2','(EGF-EGFRi*)_2',...
    'Proti','Proti',...
    'EGF_EGFRi','EGF-EGFRi',...
    'EGF_EGFRi_2', '(EGF-EGFRi)_2',...
    'Prot','Prot',...
    'EGFideg','EGFideg',...
    'GAP',  'GAP',...
    'EGF_EGFRa_2_GAP','(EGF-EGFR*)_2-GAP',...
    'EGFi','EGFi',...
    'EGF_EGFRia_2_GAP','(EGF-EGFRi*)_2-GAP',...
    'EGF_EGFRia_2_GAP_Grb2', '(EGF-EGFRi*)_2-GAP-Grb2',...
    'EGF_EGFRia_2_GAP_Grb2_Sos', '(EGF-EGFRi*)_2-GAP-Grb2-Sos',...
    'EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GDP','(EGF-EGFRi*)_2-GAP-Grb2-Sos-Ras-GDP',...
    'EGF_EGFRia_2_GAP_Grb2_Sos_Ras_GTP', '(EGF-EGFRi*)_2-GAP-Grb2-Sos-Ras-GTP',...
    'Grb2', 'Grb2',...
    'EGF_EGFRa_2_GAP_Grb2','(EGF-EGFR*)_2-GAP-Grb2',...
    'Sos', 'Sos',...
    'EGF_EGFRa_2_GAP_Grb2_Sos', '(EGF-EGFR*)_2-GAP-Grb2-Sos',...
    'Ras_GDP','Ras-GDP',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP','(EGF-EGFR*)_2-GAP-Grb2-Sos-Ras-GDP',...
    'Ras_GTP','Ras-GTP',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP',  '(EGF-EGFR*)_2-GAP-Grb2-Sos-Ras-GTP',...
    'Grb2_Sos','Grb2-Sos',...
    'Shc', 'Shc',...
    'EGF_EGFRa_2_GAP_Shc', '(EGF-EGFR*)_2-GAP-Shc',...
    'EGF_EGFRa_2_GAP_Shca','(EGF-EGFR*)_2-GAP-Shc*',...
    'EGF_EGFRa_2_GAP_Shca_Grb2','(EGF-EGFR*)_2-GAP-Shc*-Grb2',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-Ras-GDP',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-Ras-GTP',...
    'Shca_Grb2_Sos', 'Shc*-Grb2-Sos',...
    'Shca_Grb2','Shc*-Grb2',...
    'Shca','Shc*',...
    'Raf','Raf',...
    'Raf_Ras_GTP','Raf-Ras-GTP',...
    'Ras_GTPa','Ras-GTP*',...
    'Phosphatase1', 'Phosphatase1',...
    'Rafa','Raf*',...
    'Rafa_Pase','Raf*-Pase',...
    'MEK', 'MEK',...
    'MEK_Rafa','MEK-Raf*',...
    'MEK_P','MEK-P',...
    'MEK_P_Rafa','MEK-P-Raf*',...
    'MEK_PP','MEK-PP',...
    'MEK_PP_Pase2','MEK-PP-Pase2',...
    'Phosphatase2', 'Phosphatase2',...
    'MEK_P_Pase2','MEK-P-Pase2',...
    'ERK','ERK',...
    'ERK_MEK_PP','ERK-MEK-PP',...
    'ERK_P','ERK-P',...
    'ERK_P_MEK_PP','ERK-P-MEK-PP',...
    'ERK_PP', 'ERK-PP',...
    'Phosphatase3', 'Phosphatase3',...
    'ERK_PP_Pase3','ERK-PP-Pase3',...
    'ERK_P_Pase3', 'ERK-P-Pase3',...
    'EGF_EGFRia_2_GAP_Shc', '(EGF-EGFRi*)_2-GAP-Shc',...
    'EGF_EGFRia_2_GAP_Shca', '(EGF-EGFRi*)_2-GAP-Shc*',...
    'EGF_EGFRia_2_GAP_Shca_Grb2',  '(EGF-EGFRi*)_2-GAP-Shc*-Grb2',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos', '(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GDP','(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos-Ras-GDP',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_GTP','(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos-Ras-GTP',...
    'Rasi_GTP','Rasi-GTP',...
    'Rafi_Rasi_GTP','Rafi-Rasi-GTP',...
    'Rasi_GTPa', 'Rasi-GTP*',...
    'Rafia', 'Rafi*',...
    'Rafia_Pase','Rafi*-Pase',...
    'MEKi_Rafia','MEKi-Rafi*',...
    'MEKi_P','MEKi-P',...
    'MEKi_P_Rafia', 'MEKi-P-Rafi*',...
    'MEKi_PP', 'MEKi-PP',...
    'MEKi_PP_Pase2i','MEKi-PP-Pase2i',...
    'MEKi_P_Pase2i','MEKi-P-Pase2i',...
    'ERKi_MEKi_PP','ERKi-MEKi-PP',...
    'ERKi_P','ERKi-P',...
    'ERKi_P_MEKi_PP','ERKi-P-MEKi-PP',...
    'ERKi_PP', 'ERKi-PP',...
    'ERKi_PP_Pase3i', 'ERKi-PP-Pase3i',...
    'ERKi_P_Pase3i', 'ERKi-P-Pase3i',...
    'EGFRideg','EGFRideg',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_Prot','(EGF-EGFR*)_2-GAP-Grb2-Sos-Prot',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP_Prot','(EGF-EGFR*)_2-GAP-Grb2-Sos-Ras-GDP-Prot',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GTP_Prot','(EGF-EGFR*)_2-GAP-Grb2-Sos-Ras-GTP-Prot',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Prot','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Prot',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot',  '(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-Prot',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-Ras-GDP-Prot',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-Ras-GTP-Prot',...
    't_Rafa','Raf*',...
    't_Ras_GTP','Ras-GTP',...
    't_MEK_PP','MEK-PP',...
    't_ERK_PP','ERK-PP',...
    't_SHC_P_t','Total phospho SHC',...
    't_EGF_EGFRa','(EGF-EGFR*)2',...
    'C','C',...
    'EGF_EGFRia_2_deg','(EGF-EGFRi*)2_{deg}',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_ERK_PP','(EGF-EGFR*)_2-GAP-Grb2-Sos-ERK-PP',...
    'EGF_EGFRia_2_GAP_Grb2_Sos_ERKi_PP','(EGF-EGFRi*)_2-GAP-Grb2-Sos-ERKi-PP',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_ERK_PP','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos-ERK-PP',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos_ERKi_PP','(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos-ERKi-PP',...
    'Sos_ERK_PP','Sos-ERK-PP',...
    'Sos_ERKi_PP','Sos-ERKi-PP',...
    'EGF_EGFRa_2_GAP_Grb2_Sos_deg','(EGF-EGFR*)_2-GAP-Grb2-Sos_{deg}',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_deg','(EGF-EGFRi*)_2-GAP-Shc*-Grb2_{deg}',...
    'EGF_EGFRia_2_GAP_Grb2_deg','(EGF-EGFRi*)_2-GAP-Grb2_{deg}',...
    'EGF_EGFRa_2_GAP_Shca_Grb2_Sos_deg','(EGF-EGFR*)_2-GAP-Shc*-Grb2-Sos_{deg}',...
    'EGF_EGFRia_2_GAP_Grb2_Sos_deg','(EGF-EGFRi*)_2-GAP-Grb2-Sos_{deg}',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos_deg','(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos_{deg}',...
    'Sosi','Sosi',...
    'AUC_ERK_PP','AUC ERK-PP',...
    'EGF_EGFRia_2_GAP_Shca_Grb2_Sos_Ras_deg','(EGF-EGFRi*)_2-GAP-Shc*-Grb2-Sos-Ras_{deg}',...
    'EGF_EGFRia_2_GAP_Grb2_Sos_Ras_deg','(EGF-EGFRi*)_2-GAP-Grb2-Sos-Ras_{deg}',...
    'EGF_EGFRia_2_GAP_deg','(EGF-EGFRi*)_2-GAP_{deg}',...
    'EGF_EGFRia_2_GAP_Shc_deg','(EGF-EGFRi*)_2-GAP-Shc_{deg}');

for k = 1:I.nstates
    I.nmstatelegend{k} = S.(I.nmstate{k});
end
