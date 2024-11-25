function [params, paramNames, complexNames] = getMonoclonalParams(iggConcs, fcrName)
% GET MONOCLONAL PARAMS Import parameters for a the simulation
%   Input:
%       iggConcs: Array Double [igg1, igg2, igg3, igg4]
%          IgG1, IgG2, IgG3, and IgG4 concentrations
%       fcrName: string {'FcgRIIA-131H', 'FcgRIIA-131R'}
%           Name of the FcR of interest for a given simulation
%   Output:
%       params: Array (len = 22)
%           Array containing all of the parameters required to simulate the
%           model
%       paramNames: Array (len = 22)
%           Self explanatory
%       complexNames: Array (len = 24)
%           Array containing the names of all of the complexes at the end
%           of the model simulation

%% IgG-Ag Affinities
% We assume equivalent affinity across all subclasses for antigen (Kd ~ .22 nM)
agOn = 5.07e-3; % (nm^-1*s^-1)
agOff = 1.127e-3; % (s^-1)

% Copy the parameters for the solver input
f1 = agOn;%RBD-IgG1 
r1 = agOff;%         
f2 = agOn;%RBD-IgG2
r2 = agOff;
f3 = agOn;%RBD-IgG3
r3 = agOff;
f4 = agOn;%RBD-IgG4
r4 = agOff;

%% IgG-FcgR Affinities
% Pick out FcR of interest (note that this defaults to FcgRIIA-131H)
if (fcrName == "FcgRIIIA-158V")
    fcrInd = 1;
elseif (fcrName == "FcgRIIIA-158F")
    fcrInd = 2;
elseif (fcrName == "FcgRIIA-131H")
    fcrInd = 3;
elseif (fcrName == "FcgRIIA-131R")
    fcrInd = 4;
else
    disp("Invalid FcgR name - please select a valid FcR name");
    return;
end

%% Bruhns et al, 2009 IgG-Fc/FcgR Affinity Measurements
% Organization of the Bruhns measurements
FcR_kon = [20   11.7  52   35; % IgG1-Fc kon
           0.7  0.3   4.5  1.0; % IgG2-Fc kon
           98   77    8.9  9.1; % IgG3-Fc kon
           2.5   2    1.7  2.1;]*10^-6; % IgG4-Fc kon
FcR_koff = 1e-2;

% Store the correct FcgR parameters
kOns = FcR_kon(:, fcrInd);

% Copy the parameters for the solver input
f5 = kOns(1);%IgG1-FcR (nM-1)
r5 = FcR_koff;
f6 = kOns(2);%IgG2-FcR
r6 = FcR_koff;
f7 = kOns(3);%IgG3-FcR
r7 = FcR_koff;
f8 = kOns(4);%IgG4-FcR
r8 = FcR_koff;

%% Concentrations
% Antibody measurements (set from the samples)
g1tot    = iggConcs(1);%nM 
g2tot    = iggConcs(2);%nM 
g3tot    = iggConcs(3);%nM 
g4tot    = iggConcs(4);%nM

%% Experimental values
% Total antigen
etot = 2.5; % nM

% Total dimeric Fc receptor
ftot    = 20;% nM 

%% Organizing the final output
params = [f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,f7,r7,f8,r8,...
    g1tot,g2tot, g3tot,g4tot, etot, ftot];

paramNames = ["k_{on} IgG1-RBD","k_{off} IgG1-RBD","k_{on} IgG2-RBD","k_{off} IgG2-RBD",...
     "k_{on} IgG3-RBD","k_{off} IgG3-RBD","k_{on} IgG4-RBD",...
     "k_{off} IgG4-RBD","k_{on} IgG1-FcR","k_{off} IgG1-FcR","k_{on} IgG2-FcR",...
     "k_{off} IgG2-FcR","k_{on} IgG3-FcR","k_{off} IgG3-FcR","k_{on} IgG4-FcR",...
     "k_{off} IgG4-FcR","IgG1 conc","IgG2 conc", "IgG3 conc",...
     "IgG4 conc", "RBD conc", "FcR conc"];

complexNames = {'RBD-IgG1', 'RBD-IgG2', 'RBD-IgG3', 'RBD-IgG4','RBD-IgG1-IgG1', 'RBD-IgG1-IgG2','RBD-IgG1-IgG3', 'RBD-IgG1-IgG4', 'RBD-IgG2-IgG2', 'RBD-IgG2-IgG3'...
    'RBD-IgG2-IgG4','RBD-IgG3-IgG3','RBD-IgG3-IgG4','RBD-IgG4-IgG4','FcR-RBD-IgG1-IgG1', 'FcR-RBD-IgG1-IgG2','FcR-RBD-IgG1-IgG3', 'FcR-RBD-IgG1-IgG4', 'FcR-RBD-IgG2-IgG2', 'FcR-RBD-IgG2-IgG3'...
    'FcR-RBD-IgG2-IgG4','FcR-RBD-IgG3-IgG3','FcR-RBD-IgG3-IgG4','FcR-RBD-IgG4-IgG4'};

end
