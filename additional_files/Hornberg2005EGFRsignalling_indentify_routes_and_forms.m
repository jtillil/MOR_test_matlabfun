%%% Version: 19 Jan 2023
%%%
%%% Hornberg2005EGFRsignalling_indentify_routes_and_forms(I)
%%%
%%% define Shc, none-Shc, membrane-bound, internalized etc. states
%%%

function G = Hornberg2005EGFRsignalling_indentify_routes_and_forms(I)

%%% EGFR-Shc states
G.EGFR.all = find(contains(I.nmstate,'EGFR'));
G.GAP_Shc.all = find(contains(I.nmstate,'GAP_Shc'));
G.EGFR_GAP_Shc.all = intersect(G.EGFR.all,G.GAP_Shc.all);

%%% corresponding EGFR-non-Shc states
% delete substrings '_Shca' and '_Shc' from names
dummy = I.nmstate(G.EGFR_GAP_Shc.all);
dummy = replace(dummy,'_Shca','');
dummy = replace(dummy,'_Shc','');
nm_EGFR_GAP_non_Shc_states = dummy;
% determine indices
G.EGFR_GAP_non_Shc.all = NaN(size(nm_EGFR_GAP_non_Shc_states));
for s = 1:length(nm_EGFR_GAP_non_Shc_states)
    statestring = nm_EGFR_GAP_non_Shc_states{s};
    G.EGFR_GAP_non_Shc.all(s) = I.(statestring);
end
%%% some state names appear twice, e.g., when starting from 
%%% ...GAP_Shc and ...GAP_Shca. Remove doublicates
%%% Also, remove all states that end with _GAP and have no Grb2 in it
I_non_Grb2_states = setdiff(1:I.nstates,find(contains(I.nmstate,'Grb2')));
G.EGFR_GAP_non_Shc.all = unique(G.EGFR_GAP_non_Shc.all,'stable');
G.EGFR_GAP_non_Shc.all = setdiff(G.EGFR_GAP_non_Shc.all,I_non_Grb2_states,'stable');

%%% Prot states (note: I.Prot = I_Prot_states(3)
%%% important are 
%%% I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Prot, 
%%% I.EGF_EGFRa_2_GAP_Shca_Grb2_Sos_Ras_GTP_Prot
G.Prot.all = find(contains(I.nmstate,'Prot'));
G.Prot.intern = find(contains(I.nmstate,'Proti'));
G.Prot.membrane = setdiff(G.Prot.all,G.Prot.intern);

%%% internalized forms
G.intern = find(contains(I.nmstate,'i'));

%%% internalized EGFR-forms
G.EGFR.intern = find(contains(I.nmstate,'EGFRi'));
G.EGFR.membrane = setdiff(G.EGFR.all,G.EGFR.intern,'stable');

%%% internalized EGFR-GAP-Shc-forms
G.EGFR_GAP_Shc.intern = intersect(G.EGFR_GAP_Shc.all,G.EGFR.intern,'stable');

%%% membrane EGFR-GAP-Shc-forms
G.EGFR_GAP_Shc.membrane = setdiff(G.EGFR_GAP_Shc.all,G.EGFR_GAP_Shc.intern,'stable');

%%% internalized and membrane EGFR-GAP-non-Shc-forms
G.EGFR_GAP_non_Shc.intern = setdiff(G.EGFR_GAP_non_Shc.all,G.EGFR.intern,'stable');
G.EGFR_GAP_non_Shc.membrane = setdiff(G.EGFR_GAP_non_Shc.all,G.EGFR_GAP_non_Shc.intern,'stable');

%%% internalized Ras-forms
G.Ras.all = find(contains(I.nmstate,'Ras'));
G.Ras.intern = intersect(G.Ras.all,G.intern,'stable');

%%% internalized Raf-forms
G.Raf.all = find(contains(I.nmstate,'Raf'));
G.Raf.intern = intersect(G.Raf.all,G.intern,'stable');

%%% internalized and cytosol MEK-forms
G.MEK.all = find(contains(I.nmstate,'MEK'));
G.MEK.intern = intersect(G.MEK.all,G.intern,'stable');
G.MEK.cytosol = setdiff(G.MEK.all,G.intern,'stable');

%%% internalized and cytosol ERK-forms
G.ERK.all = setdiff(find(contains(I.nmstate,'ERK')),I.AUC_ERK_PP,'stable');
G.ERK.intern = intersect(G.ERK.all,G.intern,'stable');
G.ERK.cytosol = setdiff(G.ERK.all,G.ERK.intern,'stable');

%%% degadration states
G.deg.all = find(contains(I.nmstate,'deg'));
G.deg.intern = intersect(G.deg.all,G.intern,'stable');

%%% some check for consistency
if ~isempty(setxor(G.EGFR.intern,intersect(G.EGFR.all,G.intern)))
    fprintf('\n --> something wrong, since intersect(EGFR,intern) differs from EGFRi\n\n'); beep;
    error('');
end
if ~isempty(setxor(G.Prot.intern,intersect(G.Prot.all,G.intern)))
    fprintf('\n --> something wrong, since intersect(Prot,intern) differs from Proti\n\n'); beep;
    error('');
end

end

