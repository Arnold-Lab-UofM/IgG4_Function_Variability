function [compFormation] = subObjFunc(x, params, fcr)
% SUBCLASS OPTIMIZATION FUNCTION 
%   This serves as the objective function for optimization of subclass distribution as a function of (insert desired question here)
%   INPUT:
%       INSERT HERE LATER
%   OUTPUT
%       INSERT HERE LATER

%% Gathering necessary simulation pieces
% Parameter names
paramNames = ["k_{on} IgG1-RBD","k_{off} IgG1-RBD","k_{on} IgG2-RBD","k_{off} IgG2-RBD",...
     "k_{on} IgG3-RBD","k_{off} IgG3-RBD","k_{on} IgG4-RBD",...
     "k_{off} IgG4-RBD","k_{on} IgG1-FcR","k_{off} IgG1-FcR","k_{on} IgG2-FcR",...
     "k_{off} IgG2-FcR","k_{on} IgG3-FcR","k_{off} IgG3-FcR","k_{on} IgG4-FcR",...
     "k_{off} IgG4-FcR","IgG1 conc","IgG2 conc", "IgG3 conc",...
     "IgG4 conc", "RBD conc", "FcR conc"];

% Complex names
complexes = {'RBD-IgG1', 'RBD-IgG2', 'RBD-IgG3', 'RBD-IgG4','RBD-IgG1-IgG1', 'RBD-IgG1-IgG2','RBD-IgG1-IgG3', 'RBD-IgG1-IgG4', 'RBD-IgG2-IgG2', 'RBD-IgG2-IgG3'...
    'RBD-IgG2-IgG4','RBD-IgG3-IgG3','RBD-IgG3-IgG4','RBD-IgG4-IgG4','FcR-RBD-IgG1-IgG1', 'FcR-RBD-IgG1-IgG2','FcR-RBD-IgG1-IgG3', 'FcR-RBD-IgG1-IgG4', 'FcR-RBD-IgG2-IgG2', 'FcR-RBD-IgG2-IgG3'...
    'FcR-RBD-IgG2-IgG4','FcR-RBD-IgG3-IgG3','FcR-RBD-IgG3-IgG4','FcR-RBD-IgG4-IgG4'};

%% Running Simulation
% Inserting the correct params
updatedParams = params;
updatedParams(17:20) = x;

% Execute with params of interest
[yend, steadystate, complexes]  = Simulate(updatedParams, paramNames, complexes, fcr);
compFormation = yend(33);
end