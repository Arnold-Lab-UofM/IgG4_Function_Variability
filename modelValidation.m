%%% MODEL VALIDATION (Figure 3T)
%%% This script is designed to compare the model predicted values to the
%%% experimentally measured values in the multiplex assay
%%% Author: Robert Theisen, Arnold Lab, University of Michigan, Ann Arbor

%% Housekeeping
clear all;
close all;
clc;

%% Establish validation parameters


% Import data
concMultiplier = 25;
fcrData = "FcγRIIaH WT Trimer"; % ["FcγRIIaH WT Trimer", "FcγRIIaR WT Trimer"]
df_converted = importConvertedSpiking("WT", "Trimer", false);

% Run validation
[pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]} .* concMultiplier, ...
                                        df_converted{:, fcrData}, fcr);
% Plotting
validationPlot(pred, actual, spRho, fcr, true);