%% setup
addpath(genpath('./additional_files'))
load("modelEGFR_minimal.mat")

par = model.par;
X = model.X0; % + normrnd(100, 10, [length(model.X0) 1]);
tstart = 0;
I = model.I;

odefun_original = @Hornberg2005EGFRsignalling_ode;
odefun_original_I = @Hornberg2005EGFRsignalling_ode_I;
jacfun_original = @Hornberg2005EGFRsignalling_odejac;

syms t
par_sym = cell2sym(model.I.nmpar(:));
X_sym = cell2sym(model.I.nmstate(:));
odefun_symbolic = odefun_original(t,X_sym,par_sym,model);
odefun_matlabfun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});

jacfun_matlabfun = model.jacfun;

%% time text coded ode function

N = 100000;

tic
for i = 1:N
    res1 = odefun_original_I(tstart, X, par, I);
end
fprintf("file ode function\n")
toc
% 7.3 sec

%% time text coded ode function with supplying I instead of model

tic
for i = 1:N
    res1 = odefun_original_I(tstart, X, par, I);
end
fprintf("\nfile ode function supplying I instead of model (smaller variable to RAM)\n")
toc
% 7.4 sec

%% time symbolically coded matlabfunction ode function

tic
for i = 1:N
    res2 = odefun_matlabfun(X, par);
end
fprintf("\nmatlabfun ode function\n")
toc
% 0.037 sec
% !!! 200 times faster !!!

%% check

% error = res1 - res2;
% 0

%% time jacobian vs fast jacobian

N = 10000;

tic
for i = 1:N
    resr = jacfun_original(tstart, X, par, model);
end
fprintf("\nfile jac function\n")
toc
% N = 10000: 3.5 sec

tic
for i = 1:N
    resr = jacfun_matlabfun(X, par);
end
fprintf("\nmatlabfun jac function\n")
toc
% N = 10000: 0.24 sec
% !!! 10 times faster !!!

%% try pss coding

% syms Sos
% SosID = model.I.Sos;
% SosODE = odefun_symbolic(SosID);
% SosEquation = 0 == SosODE;
% 
% SosSolution = solve(SosEquation, Sos);
% 
% tic
% adapted_odefun_symbolic = odefun_symbolic;
% idx = has(adapted_odefun_symbolic, Sos);
% adapted_odefun_symbolic(idx) = subs(adapted_odefun_symbolic(idx), Sos, SosSolution);
% adapted_odefun_symbolic(SosID) = 0;
% toc
% 
% tic
% adapted_odefun_matlabfun = matlabFunction(odefun_symbolic,'Vars',{X_sym,par_sym});
% toc

%% time adapted symbolically coded matlabfunction ode function

% tic
% for i = 1:N
%     resa = adapted_odefun_matlabfun(X, par);
% end
% toc
% 0.04 sec

%% compute all pss solutions

% tic
% for i = 1:length(X_sym)
%     % disp(X_sym(i))
%     % disp(odefun_symbolic(i))
%     sol = solve(0 == odefun_symbolic(i), X_sym(i));
%     if isempty(sol)
%         sol = 0;
%     end
%     pssSolutions(i) = sol(end);
% end
% toc

% i = 1;
% X_sym(i)
% odefun_symbolic(i)
% solve(0 == odefun_symbolic(i), X_sym(i))

%% try repeated pss coding

% repeated_odefun_symbolic = odefun_symbolic;
% 
% for i = 1:30
%     disp(i)
%     tic
%     % idx = has(repeated_odefun_symbolic, X_sym(i));
%     % repeated_odefun_symbolic(idx) = subs(repeated_odefun_symbolic(idx), X_sym(i), pssSolutions(i));
%     repeated_odefun_symbolic = subs(repeated_odefun_symbolic, X_sym(i), pssSolutions(i));
%     repeated_odefun_symbolic(i) = 0;
%     disp(toc)
% end
% 
% tic
% repeated_odefun_matlabfun = matlabFunction(repeated_odefun_symbolic,'Vars',{X_sym,par_sym});
% toc
% 0.6 - 3 sec

%% time repeated pss symbolically coded matlabfunction ode function

% tic
% for i = 1:N
%     resr = repeated_odefun_matlabfun(X, par);
% end
% toc
% 0.25 - 3 sec (20 pss)
% 6 - 86 sec (30 pss)
