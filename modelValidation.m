%%% MODEL VALIDATION (Figure 3T)
%%% This script is designed to compare the model predicted values to the
%%% experimentally measured values in the multiplex assay.
%%% Author: Robert Theisen, Arnold Lab, University of Michigan, Ann Arbor
%"FcγRIIaH WT Trimer", "FcγRIIaR WT Trimer"

%% Housekeeping
clear all;
close all;
clc;

%% Import IgG subclass titers
% Set Overarching Parameters
concMultiplier = 25;
df_converted = importConvertedSpiking("WT", "Trimer", false);

% Adjust concentrations by multiplier
iggConcs = df_converted{:, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]} .* concMultiplier;


%% Simulate the model for all samples
% FcgRIIA-131H
[pred2ah, actual2ah, spRho2ah] = validateModel(iggConcs, df_converted{:, "FcγRIIaH WT Trimer"}, "FcgRIIA-131H");

% FcgRIIA-131R
[pred2ar, actual2ar, spRho2ar] = validateModel(iggConcs, df_converted{:, "FcγRIIaR WT Trimer"}, "FcgRIIA-131R");

%% Plot the predicted vs actual
% FcgRIIA-131H
figure;
validationPlot(pred2ah, actual2ah, spRho2ah, "FcgRIIA-131H", true);

% FcgRIIA-131R
figure;
validationPlot(pred2ar, actual2ar, spRho2ar, "FcgRIIA-131R", true);
