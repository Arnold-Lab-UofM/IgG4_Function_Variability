%%% MONOCLONAL SPIKING SURFACES
%%% This script generates the landscape of responses for 
%%% Author: Robert Theisen, Arnold Lab, University of Michigan, Ann Arbor

%% Housekeeping
clear all;
close all;
clc;

%% Set overarching parameters
% Determining the antigen specific subclass titers for import
ag = "Trimer";
var = "WT";

% Setting the concentration multiplier to proportonally adjust the input
% concentrations and the associated spikes
concMultiplier = 25;
sIgG4Spike = .05*concMultiplier;
lIgG4Spike = .2*concMultiplier;


%% Import pre-converted IgG titer data
dfConverted = importConvertedSpiking("WT", "Trimer", false);

% Calculate the median concentrations as a baseline for the surfaces
baselineMedians = median(dfConverted{:, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]})*concMultiplier;

%% Establish parameter ranges
% Get baseline parameters for both FcgRIIA polymorphisms
[baseParams2AH, paramNames, complexNames] = getMonoclonalParams(baselineMedians, "FcgRIIA-131H");
[baseParams2AR, paramNames, complexNames] = getMonoclonalParams(baselineMedians, "FcgRIIA-131R");

% Define the ranges for IgG1 and IgG4
igg1_range = logspace(log10(.1*baselineMedians(1)), log10(10*baselineMedians(1)), 50);
igg4_range = logspace(log10(.9*baselineMedians(4)), log10(150*baselineMedians(4)), 100);


%% Iteratively call model
% Define data carriers
IIAH_Surface = zeros(50, 100);
IIAR_Surface = zeros(50, 100);

% Iteration
for i = 1:length(igg1_range)
    for j = 1:length(igg4_range)
        % FcgRIIA-H
        fcr = "FcgRIIA-131H";
        tempParams = baseParams2AH;
        tempParams(17) = igg1_range(i);
        tempParams(20) = igg4_range(j);
        [yBase, steadystate, complexes] = Simulate(tempParams, paramNames, complexNames, fcr);
        IIAH_Surface(i, j) = yBase(33);

        % FcgRIIA-R
        fcr = "FcgRIIA-131R";
        tempParams = baseParams2AR;
        tempParams(17) = igg1_range(i);
        tempParams(20) = igg4_range(j);
        [ySpike, steadystate, complexes] = Simulate(tempParams, paramNames, complexNames, fcr);
        IIAR_Surface(i, j) = ySpike(33);
    end
    % Progress reporter
    fprintf("Completed Cycle %i/50\n", i);
end

%% Draw lines equivalent to the baseline and each of the IGG4 spikes
% Base line
igg4_base_line_2ah = zeros(1, length(igg1_range));
igg4_base_line_2ar = zeros(1, length(igg1_range));
igg4_base_value = ones(1, length(igg1_range)) * baseParams2AH(20);

% Low spike line
igg4_spike_line_2ah = zeros(1, length(igg1_range));
igg4_spike_line_2ar = zeros(1, length(igg1_range));
igg4_spike_value = ones(1, length(igg1_range)) * (baseParams2AH(20) + sIgG4Spike);

% High Spike line
igg4_high_spike_line_2ah = zeros(1, length(igg1_range));
igg4_high_spike_line_2ar = zeros(1, length(igg1_range));
igg4_high_spike_value = ones(1, length(igg1_range)) * (baseParams2AH(20) + lIgG4Spike);

for i = 1:length(igg1_range)
    % Initialize
    temp2AHParams = baseParams2AH;
    temp2ARParams = baseParams2AR;
    temp2AHParams(17) = igg1_range(i);
    temp2ARParams(17) = igg1_range(i);

    % Baseline calculation
    [yBase, steadystate, complexes] = Simulate(temp2AHParams, paramNames, complexNames, "FcgRIIA-131H");
    igg4_base_line_2ah(i) = yBase(33);
    [yBase, steadystate, complexes] = Simulate(temp2ARParams, paramNames, complexNames, "FcgRIIA-131R");
    igg4_base_line_2ar(i) = yBase(33);

    % Small Spike
    temp2AHParams(20) = temp2AHParams(20) + sIgG4Spike;
    [yBase, steadystate, complexes] = Simulate(temp2AHParams, paramNames, complexNames, "FcgRIIA-131H");
    igg4_spike_line_2ah(i) = yBase(33);
    temp2ARParams(20) = temp2ARParams(20) + sIgG4Spike;
    [yBase, steadystate, complexes] = Simulate(temp2ARParams, paramNames, complexNames, "FcgRIIA-131R");
    igg4_spike_line_2ar(i) = yBase(33);

    % Large Spike
    temp2AHParams(20) = temp2AHParams(20) + lIgG4Spike - sIgG4Spike;
    [yBase, steadystate, complexes] = Simulate(temp2AHParams, paramNames, complexNames, "FcgRIIA-131H");
    igg4_high_spike_line_2ah(i) = yBase(33);
    temp2ARParams(20) = temp2ARParams(20) + lIgG4Spike - sIgG4Spike;
    [yBase, steadystate, complexes] = Simulate(temp2ARParams, paramNames, complexNames, "FcgRIIA-131R");
    igg4_high_spike_line_2ar(i) = yBase(33);
end

%% Calculating the color data
% By local gradient (x is IgG4, y is IgG1)
[dfdx, dfdy] = gradient(IIAH_Surface);
[dfdx1, dfdy1] = gradient(IIAR_Surface);

%% Plotting
% Colormap loading
red_white_blue = load("rwb_cmap.mat");

% 2AH figure (Figure 3U)
figure;
surf(igg4_range, igg1_range, IIAH_Surface, CData=dfdx);
colormap(red_white_blue.USA_2AH);
clim([-.0025, .0009])
% 2AH lines
hold on;
plot3(igg4_base_value, igg1_range, igg4_base_line_2ah, Color="black", LineWidth=3);
plot3(igg4_spike_value, igg1_range, igg4_spike_line_2ah, Color="black", LineWidth=3);
plot3(igg4_high_spike_value, igg1_range, igg4_high_spike_line_2ah, Color="black", LineWidth=3)
hold off;
% Labels
xlabel("IgG4 Concentration (nM)");
ylabel("IgG1 Concentration (nM)");
zlabel("Complex Formation (nM)");
%title("FcγRIIaH")
set(gca, 'XScale', 'log', 'YScale', 'log');
zlim([0, .25]);


% FcgR2AR figure (Figure 3V)
figure;
surf(igg4_range, igg1_range, IIAR_Surface, CData=dfdx1);
colormap(red_white_blue.USA_2AH);
clim([-.0025, .0009])
% 2AR lines
hold on;
plot3(igg4_base_value, igg1_range, igg4_base_line_2ar, Color="black", LineWidth=3);
plot3(igg4_spike_value, igg1_range, igg4_spike_line_2ar, Color="black", LineWidth=3);
plot3(igg4_high_spike_value, igg1_range, igg4_high_spike_line_2ar, Color="black", LineWidth=3)
hold off;
% Labels
xlabel("IgG4 Concentration (nM)");
ylabel("IgG1 Concentration (nM)");
zlabel("Complex Formation (nM)");
%title("FcγRIIaR")
set(gca, 'XScale', 'log', 'YScale', 'log');
zlim([0, .25]);
