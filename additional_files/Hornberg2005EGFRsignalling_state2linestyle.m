%%% Version: 19 Jan 2023
%%%
%%% call by: Hornberg2005EGFRsignalling_state2linestyle(I)
%%%
%%% This function assigns the line style to the given state variable to 
%%% facilitate plotting and to have always the same style for the specific 
%%% species.
%%%
%%% Citation:
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function state2style = Hornberg2005EGFRsignalling_state2linestyle(I)

%%% define when colormap repeats the colors so that a different line styles 
%%% should be used then
repeat = 7;

styles = {'-','--',':','-.'};
nstyles = length(styles);

% initialize output and assign colors
state2style = cell(I.nstates,1);

for k = 1:I.nstates
    state2style(k) = styles(mod(floor(k/repeat),nstyles)+1);
end

%%% manual changes (if necessary)
state2style{I.Ras_GTP} = ':';
state2style{I.Rafa} = '--';

end
