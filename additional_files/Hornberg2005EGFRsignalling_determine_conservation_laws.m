%%% Version: 19 Jan 2023
%%%
%%% call by: Hornberg2005EGFRsignalling_determine_conservation_laws
%%%
%%% Determine the conservation laws of the EGFR signalling model 
%%% by Hornberg et al (2005)
%%% 
%%% Authors: Jane Knoechel & Wilhelm Huisinga
%%%

function conlaw = Hornberg2005EGFRsignalling_determine_conservation_laws(I,conlawspec)

% list (cell array) containing the name(s) of the conserved quantities
L.conlawspec = conlawspec;

% initialise conservation law structure L
L.nconlaws = length(L.conlawspec); % numnber of conlaws
L.statesofconlaw = cell(1,L.nconlaws);
L.statesofallconlaws = [];         % states that belong to some conlaw
L.statesofmultipleconlaws = [];    % states that belong to multiple conlaws
L.valconlaw = NaN(1,L.nconlaws);   % value of conlaw

% initialise some fields that are used later during model order reduction
L.con2conlaw = [];
L.remstatesofconlaw = {};

%%% -----------------------------------------------------------------------
%%% determine conservation laws based on the name(s) of the state
%%% variables that contain specifier (term or phrase) in its name
%%%
for c = 1:length(L.conlawspec)
    
    specifier = L.conlawspec{c};
    
    % determine indices of all states that contain the keyword
    L.statesofconlaw{c} = find(contains(I.nmstate,specifier));
    
    % Phosphatase1 is called a_Pase in a complex, while 
    % Phasphatase2 is called Pase2 in a complex. Manually add the free
    % phosphatses to the corresponding conservation las
    if strcmp(specifier,'a_Pase')
        L.statesofconlaw{c} = [I.Phosphatase1,L.statesofconlaw{c}];
    elseif strcmp(specifier,'Pase2')
        L.statesofconlaw{c} = [I.Phosphatase2,L.statesofconlaw{c}];
    end

    L.statesofmultipleconlaws = [L.statesofmultipleconlaws,intersect(L.statesofallconlaws,L.statesofconlaw{c})];
    L.statesofallconlaws = [L.statesofallconlaws, L.statesofconlaw{c}];
    
end

% remove dublicates in the lists
L.statesofallconlaws = unique(L.statesofallconlaws,'stable');
L.statesofmultipleconlaws = unique(L.statesofmultipleconlaws,'stable');
L.statesofsingleconlaw = setdiff(L.statesofallconlaws,L.statesofmultipleconlaws,'stable');

% order conlaws in increasing number of states
orderconlaws = 1;
if orderconlaws
   
    % number of states involved in the conlaws
    L.nstatesconlaw = NaN(L.nconlaws,1);
    for c = 1:L.nconlaws
        L.nstatesconlaw(c) = length(L.statesofconlaw{c});
    end
    [~,ind] = sort(L.nstatesconlaw);
    
    % sort conlaw quantities 
    L.nstatesconlaw = L.nstatesconlaw(ind);
    L.conlawspec = L.conlawspec(ind);
    L.statesofconlaw = L.statesofconlaw(ind);
    L.valconlaw = L.valconlaw(ind);
    
    
end

%%% assign output 

conlaw.n = L.nconlaws;

for c = 1:length(L.conlawspec)
    
    nmconlaw = L.conlawspec{c};
    conlaw.(nmconlaw).states = L.statesofconlaw{c};

end


%%% list the identifed conservation laws (if verbose = 1)
%%%
verbose = 0;
if verbose
    fprintf('\n\n Identified conservation laws');
    if orderconlaws
        fprintf(' (ordered in increasing number of states)');
    end
    for c = 1:length(L.conlawspec)
        
        specifier = L.conlawspec{c};
        fprintf('\n    %d. %s (%d states): ',c,specifier,length(L.statesofconlaw{c}));
        ndegproducts = sum(contains(I.nmstate(L.statesofconlaw{c}),'deg'));
        if ndegproducts>0
            fprintf('[contains %d deg products]  ',ndegproducts);
        end
        
        fprintf('%s, ',I.nmstate{L.statesofconlaw{c}});
        
    end
    fprintf('\n States only in a single conservation laws (%d states): \n    ',length(L.statesofsingleconlaw));
    fprintf('%s, ',I.nmstate{L.statesofsingleconlaw});
    fprintf('\n States in multiple conservation laws (%d states): \n    ',length(L.statesofmultipleconlaws));
    fprintf('%s, ',I.nmstate{L.statesofmultipleconlaws});
    fprintf('\n\n');
end






