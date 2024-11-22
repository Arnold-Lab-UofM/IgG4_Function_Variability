function [yend, steadystate, complexes]  = Simulate(params, paramnames, complexes, fcr_ttl)
% SIMULATE Runs the model one time for a given set of parameters and inputs
%   INPUT:
%       params: Vector (len = 22)
%           Array containing the associated IgG subclass concentrations,
%           IgG affinities for antigen and FcR, and experimental
%           concentrations
%       paramnames: Vector strings (len = 22)
%           Names for each of the parameters
%       complexes: Vector strings (len = 24)
%           Names for each of the complexes tracked during the model
%           simulation
%       fcr_ttl: String
%           Name of the Fc receptor being simulated (optional parameter)
%   OUTPUT:
%       yend: Vector (len = 33)
%           The final, steady state values of each of the complexes
%           following the ODE simulation. yend(33) represents the sum of
%           all bound Fc receptor dimers.
%       steadystate: boolean
%           A statement to check that a given simulation has reached steady
%           state
%       complexes: vector strings (len = 33)
%           Names for each of the complexes including the free species and
%           some of the cumulative values

%% Run Normal Simulation
% Initialize Global Variables
global tStep tArray varArray;
tStep = 1; 
tArray = zeros(1e6,1); 
varArray = zeros(1e6,6);  

% Set initial values, time span, and solver options
y0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Initial complex concentrations
tspan = [0 100000];
options = odeset('AbsTol',1e-50,'RelTol',1e-10);

% Run the ODE solver
[t,y] = ode15s(@ODEs,tspan,y0,options,params);

% Post-processing 
tArray = tArray(1:min(tStep,length(y(:,1))));   % truncate unused variables
varArray = varArray(1:min(tStep,length(y(:,1))),:);
ycut = y(1:min(tStep,length(y(:,1))),:);

% free species concentrations over time
g1y=varArray(:,1);
g2y=varArray(:,2);
g3y=varArray(:,3);
g4y=varArray(:,4);
ey=varArray(:,5);
fy=varArray(:,6);

% Adding free species concentrations and complex summations to the
% results matrix
ycut(:,25) = g1y;
ycut(:,26) = g2y;
ycut(:,27) = g3y;
ycut(:,28) = g4y;
ycut(:,29) = ey;
ycut(:,30) = fy;
ycut(:,31) = ycut(:,15)+ycut(:,17)+ycut(:,22);
ycut(:,32) = ycut(:,15)+ycut(:,16)+ycut(:,17)+ycut(:,18)+ycut(:,20)+ycut(:,22)+ycut(:,23);
ycut(:,33) = sum(ycut(:,15:24),2); 

complexes = [complexes cellstr(["IgG1", "IgG2", "IgG3", "IgG4", "env",...
    "FcR", "FcR complexes with IgG1 and IgG3 only", ...
    "FcR complexes including IgG1 and IgG3","All FcR complexes"])];

% Concentrations at the last time point
yend = ycut(end,:);

% Steady state definition and check
steadystate = sum(abs(yend-ycut(end-floor(0.05*length(ycut(:,1))),:)) < 1e-5*yend);
    
end
