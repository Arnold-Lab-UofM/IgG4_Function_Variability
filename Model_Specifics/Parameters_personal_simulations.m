% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

function [params, paramnames, complexes] = Parameters_personal_simulations(g1_F, ...
    g2_F, g3_F, g4_F, g1_F_off, g2_F_off, g3_F_off, g4_F_off, g1tot_ind, g2tot_ind, ...
    g3tot_ind, g4tot_ind)

%% define parameters
colors = [[66, 134, 244]/255; [219, 48, 48]/255; [137, 196, 74]/255; [178, 91, 175]/255; [99, 210, 216]/255; [255, 248, 58]/255; [0 0 0]; [1 1 1]];

complexes = {'RBD-IgG1', 'RBD-IgG2', 'RBD-IgG3', 'RBD-IgG4','RBD-IgG1-IgG1', 'RBD-IgG1-IgG2','RBD-IgG1-IgG3', 'RBD-IgG1-IgG4', 'RBD-IgG2-IgG2', 'RBD-IgG2-IgG3'...
    'RBD-IgG2-IgG4','RBD-IgG3-IgG3','RBD-IgG3-IgG4','RBD-IgG4-IgG4','FcR-RBD-IgG1-IgG1', 'FcR-RBD-IgG1-IgG2','FcR-RBD-IgG1-IgG3', 'FcR-RBD-IgG1-IgG4', 'FcR-RBD-IgG2-IgG2', 'FcR-RBD-IgG2-IgG3'...
    'FcR-RBD-IgG2-IgG4','FcR-RBD-IgG3-IgG3','FcR-RBD-IgG3-IgG4','FcR-RBD-IgG4-IgG4'};

f1 = 5.07e-4;%RBD-IgG1 from IgG vs IgA model
r1 = 3.38e-3;%reduced OM by 2
f2 = 5.07e-4;%RBD-IgG2
r2 = 3.38e-3;
f3 = 5.07e-4;%RBD-IgG3
r3 = 3.38e-3;
f4 = 5.07e-4;%RBD-IgG4
r4 = 3.38e-3;
f5 = g1_F;%IgG1-FcR (nM-1)
r5 = g1_F_off;
f6 = g2_F;%IgG2-FcR
r6 = g2_F_off;
f7 = g3_F;%IgG3-FcR
r7 = g3_F_off;
f8 = g4_F;%IgG4-FcR
r8 = g4_F_off;

g1tot    = g1tot_ind;%nM 
g2tot    = g2tot_ind;%nM 
g3tot    = g3tot_ind;%nM 
g4tot    = g4tot_ind;%nM 
etot    = 25;%nM 
ftot    = 20;%nM 

params = [f1,r1,f2,r2,f3,r3,f4,r4,f5,r5,f6,r6,f7,r7,f8,r8,...
    g1tot,g2tot, g3tot,g4tot, etot, ftot];
 
paramnames = ["k_{on} IgG1-RBD","k_{off} IgG1-RBD","k_{on} IgG2-RBD","k_{off} IgG2-RBD",...
     "k_{on} IgG3-RBD","k_{off} IgG3-RBD","k_{on} IgG4-RBD",...
     "k_{off} IgG4-RBD","k_{on} IgG1-FcR","k_{off} IgG1-FcR","k_{on} IgG2-FcR",...
     "k_{off} IgG2-FcR","k_{on} IgG3-FcR","k_{off} IgG3-FcR","k_{on} IgG4-FcR",...
     "k_{off} IgG4-FcR","IgG1 conc","IgG2 conc", "IgG3 conc",...
     "IgG4 conc", "RBD conc", "FcR conc"];
end