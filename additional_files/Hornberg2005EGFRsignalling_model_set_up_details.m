%%% Version: 19 Jan 2023
%%%
%%% model  = <MODELNAME>_model_set_up_details(model)
%%%
%%% This function specifies all required details that are needed to finally
%%% set of the model
%%% 
%%% Input:  model       model structure (containing the fields model.name 
%%%                     and model.scenario)
%%%                   
%%% Output: model       updated model structure
%%%
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function model = Hornberg2005EGFRsignalling_model_set_up_details(model)

%%% define naming for result files 
model.projectfolder = ''; %%% option to provide an absolute path here
model.savenameroot = [model.projectfolder 'results/' model.name];
if isfield(model, 'scenario')
if ~isempty(model.scenario)
    model.savenameroot = [model.savenameroot '_' model.scenario];
end
end
%%% for testing purposes, reduced number of computations (i.e., timepoints,
%%% states) when computing the indices
model.quicktest = false;

% load indexing to be able to define input, output, conslaws etc
I = feval([model.name '_indexing']);

% define input
model.setup.input      = I.EGF;
model.setup.unit.input = 'M'; 
model.setup.u_ref = 5.0e-8; % with MW(EGF) = 6400g/mol, 50nM = 50nmol/L*6400g/mol = 50*6.4e3 ng/1e3mL = 50*6.4 ng/mL = 320 ng/mL

% define output 
model.setup.output      = I.ERK_PP;
model.setup.unit.output = 'molecules/cell';
model.setup.unit.graphic.transfoutput = @(x) x;
model.setup.unit.graphic.output = 'molecules/cell';  
model.setup.unit.graphic.ylim   = [1e2 5e7]; % in output units

% simulation time span and unit transformation (here: sec to min)
model.setup.tspan     = [0 100*60];  % in [sec]
% model.setup.tspan     = [0 10];  % in [sec]
model.setup.unit.time = 'sec';
model.setup.unit.graphic.transftime = @(x) x/60;
model.setup.unit.graphic.time = 'min';
model.setup.unit.graphic.xlim = [-2 102]; % in output units

% define observed state variables (to be plottet) and plot
model.setup.obsstates = [I.EGFR, I.EGF_EGFR, I.EGF_EGFR_2, I.EGF_EGFRa_2, I.EGF_EGFRa_2_GAP,...
    I.EGF_EGFRa_2_GAP_Grb2_Sos_Ras_GDP, I.EGF_EGFRa_2_GAP_Shca,...
    I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos, I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot,...
    I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GDP_Prot, I.EGF_EGFRia_2_GAP_Shca_Grb2_Sos,...
    I.ERK, I.ERK_P, I.ERK_PP, I.ERK_PP_Pase3, I.ERK_P_MEK_PP, I.ERK_P_Pase3, I.Grb2_Sos,...
    I.MEK, I.MEK_P, I.MEK_PP, I.MEK_P_Rafa, I.MEK_Rafa, I.Phosphatase3, I.Raf, I.Ras_GDP,...
    I.Ras_GTP, I.Ras_GTPa];

% conservation law for this small example was set manually
conlawspec = {'MEK','Ras','Sos','Shc','Prot','Raf','a_Pase','Pase2'};
model.setup.conlaw = Hornberg2005EGFRsignalling_determine_conservation_laws(I,conlawspec); 

% optional: jacobian of ODE and jacobian pattern of extended ODE 
% false = no jacobian provided 
% true = jacobian provided 
model.setup.jacfunprovided = true;

% optional: define your specific legend labels (if false, default values are used)
% use the file 'default_legendlabels.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.legendlabelsprovided = true;

% optional: define your specific line styles (if false, default values are used)
% use the file 'default_state2linestyle.m' in the folder 'generalfiles' as a
% template to define your own legend labels, if desired
model.setup.linestyleprovided = true;
