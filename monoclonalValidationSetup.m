%%% Monoclonal Validation Setup
%%% Script to run simulations with ABSOLUTE changes in individual IgG
%%% subclass levels in accordance with a booster dose
%%% Author: Robert Theisen      Last updated: November 6th, 2023

%% Housekeeping
clear all;
close all;
clc;

%% Importing parameters
% Pick FcR for Simulations
fcr = 'FcgRIIA-131H';

% Set simulated boosting multipliers
rIgG1 = 1.075;
rIgG2 = 2.43;
rIgG3 = .06;
rIgG4 = 199;

%% Levels for the addition of IgG4
saIgG4 = .05;
laIgG4 = .2;

%% Levels for reduction of IgG3
% scIgG3 = .25;
% mcIgG3 = .50; % roughly consistent with sci immuno (fig 3)
% lcIgG3 = .75;
% vlcIgG3 = .95; % roughly consistent with sci immuno (fig 1)

%% Importing all relevant vaccinee data
% Select options
vacc = "Pfizer";
tp = "Dose 2";
agInterest = "Trimer";
var = "WT";
refConc = 5;

% Run the full conversion
%[~, allMfis, allConverted, allMedians] = importAll(refConc);
[~, ~, d2Iso, d2Converted] = importSpiking(agInterest, var);

% Isolate Dose 2 Pfizer
%[d2Converted, d2Medians] = sortImportAll(allConverted, vacc, tp, agInterest);

% Import Parameters
% Get dose 2 params
[d2Params, paramNames, complexNames] = getBaselineParams([0, 0, 0, 0], fcr);


%% Full Dose 2
nd2Vaccinees = size(d2Converted, 1);
ag = agInterest;

% Dose 2
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Set params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]};
    persd2Params(17:20) = persd2Params(17:20) * 200;

    % Run simulation
    [ybased2, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybased2(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];
end

% Setup tables
d2Pers = table();
d2Pers.IDs = idCarrier;
d2Pers.ComplexFormation = dataCarrier;



%% Small IgG4 Spike
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]};
    %persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};
    persd2Params(20) = saIgG4 + persd2Params(20);
    persd2Params(17:20) = persd2Params(17:20) * 200;

    % Run simulation
    [ybaseb, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybaseb(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
sIgG4Pers = table();
sIgG4Pers.IDs = idCarrier;
sIgG4Pers.ComplexFormation = dataCarrier;

%% Large IgG4 Spike
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]};
    %persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};
    persd2Params(20) = laIgG4 + persd2Params(20);
    persd2Params(17:20) = persd2Params(17:20) * 200;

    % Run simulation
    [ybaseb, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybaseb(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
lIgG4Pers = table();
lIgG4Pers.IDs = idCarrier;
lIgG4Pers.ComplexFormation = dataCarrier;


%% Consolidation
finalTable = table();
finalTable.ID = d2Pers.IDs;
finalTable.Baseline = d2Pers.ComplexFormation;
finalTable.("Small IgG4 Spike") = sIgG4Pers.ComplexFormation;
finalTable.("Large IgG4 Spike") = lIgG4Pers.ComplexFormation;

%% Sort into our smaller vaccinee population
% Get the IDs
[~, ~, dfIds, ~] = importSpiking("Trimer", "WT");

% Isolate in the final output
idMatchup = ismember(finalTable.ID, dfIds.("Sample ID"));
finalTable = finalTable(idMatchup, :);

%% Align Experimental and Predicted Approaches
% Sort the final table
sortedFinalTable = sortrows(finalTable, "ID");

% Split and sort the measured values
measuredBaseline = dfIds(dfIds{:, "IgG4 Cocktail Concentration"} == 0, :);
sortedMeasuredBaseline = sortrows(measuredBaseline, "Sample ID");

measuredLowSpike = dfIds(dfIds{:, "IgG4 Cocktail Concentration"} == .05, :);
sortedMeasuredLowSpike = sortrows(measuredLowSpike, "Sample ID");

measuredHighSpike = dfIds(dfIds{:, "IgG4 Cocktail Concentration"} == .2, :);
sortedMeasuredHighSpike = sortrows(measuredHighSpike, "Sample ID");

% Calculate Differences in Prediction
predictedDiff = sortedFinalTable{:, "Large IgG4 Spike"} - sortedFinalTable{:, "Baseline"};
measuredDiff = sortedMeasuredHighSpike{:, "FcγRIIaH WT Trimer"} - sortedMeasuredBaseline{:, "FcγRIIaH WT Trimer"};

%% Plotting/Analyzing Results
% Plot
scatter(predictedDiff, measuredDiff, 'blue', 'filled');
xlabel("Predicted Difference (nM)");
ylabel("Measured Difference (MFI)");
title("Large Spike - Baseline (3A, Full Conversion)")
xline(0);
yline(0);

% Stats
rho = corr(predictedDiff, measuredDiff, type="Pearson");

% Text
NW = [min(xlim) max(ylim)] + [diff(xlim) -diff(ylim)]*.05;
rhoString = sprintf("%1.3f", rho);
text(NW(1), NW(2), "\rho = " + rhoString, FontSize=18, VerticalAlignment="top", HorizontalAlignment="left");

%{
%% Knock in model (selective addition)
% 1) Small IgG3 cut
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};

    % Modify each IgG subclass concentration by known ratio values (comment
    % out whatever we want to remove)
    %persd2Params(17) = rIgG1*persd2Params(17);
    %persd2Params(18) = rIgG2*persd2Params(18);
    persd2Params(19) = (1 - scIgG3)*persd2Params(19);
    %persd2Params(20) = sIgG4 + persd2Params(20);

    % Run simulation
    [ybased2, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybased2(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
modbPers1 = table();
modbPers1.IDs = idCarrier;
modbPers1.ComplexFormation = dataCarrier;


% 2) Medium IgG3 cut
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};

    % Modify each IgG subclass concentration by known ratio values (comment
    % out whatever we want to remove)
    %persd2Params(17) = rIgG1*persd2Params(17);
    %persd2Params(18) = rIgG2*persd2Params(18);
    persd2Params(19) = (1 - mcIgG3)*persd2Params(19);
    %persd2Params(20) = mIgG4 + persd2Params(20);

    % Run simulation
    [ybased2, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybased2(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
modbPers2 = table();
modbPers2.IDs = idCarrier;
modbPers2.ComplexFormation = dataCarrier;


% 3) Large IgG3 cut
dataCarrier = zeros(nd2Vaccinees, 1);

% Looking into some complex specific information
igg11_ag = zeros(nd2Vaccinees, 1);
igg14_ag = zeros(nd2Vaccinees, 1);
igg44_ag = zeros(nd2Vaccinees, 1);

idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};

    % Modify each IgG subclass concentration by known ratio values (comment
    % out whatever we want to remove)
    %persd2Params(17) = rIgG1*persd2Params(17);
    %persd2Params(18) = rIgG2*persd2Params(18);
    persd2Params(19) = (1 - lcIgG3)*persd2Params(19);
    %persd2Params(20) = lIgG4 + persd2Params(20);

    % Run simulation
    [ybased2, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybased2(33);
    igg11_ag(i) = ybased2(5);
    igg14_ag(i) = ybased2(8);
    igg44_ag(i) = ybased2(14);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
modbPers3 = table();
modbPers3.IDs = idCarrier;
modbPers3.ComplexFormation = dataCarrier;


% 4) Very large IgG3 cut
dataCarrier = zeros(nd2Vaccinees, 1);
idCarrier = [];
for i = 1:nd2Vaccinees
    % Establish params
    persd2Params = d2Params;
    persd2Params(17:20) = d2Converted{i, ["IgG1", "IgG2", "IgG3", "IgG4"]};

    % Modify each IgG subclass concentration by known ratio values (comment
    % out whatever we want to remove)
    %persd2Params(17) = rIgG1*persd2Params(17);
    %persd2Params(18) = rIgG2*persd2Params(18);
    persd2Params(19) = (1 - vlcIgG3)*persd2Params(19);
    %persd2Params(20) = vlIgG4 + persd2Params(20);

    % Run simulation
    [ybased2, steadystate, complexes] = Simulate(persd2Params, paramNames, complexNames, fcr);

    % Record
    dataCarrier(i) = ybased2(33);
    idCarrier = [idCarrier; d2Converted{i, "Sample ID"}];

end

% Setup tables
modbPers4 = table();
modbPers4.IDs = idCarrier;
modbPers4.ComplexFormation = dataCarrier;


%% Consolidate Final Table
final_table = table(d2Converted.("Sample ID"), d2Pers.ComplexFormation, modbPers1.ComplexFormation, modbPers2.ComplexFormation,...
                    modbPers3.ComplexFormation, modbPers4.ComplexFormation);


%% Check on Statistical Differences
% Setup Groups
g1 = d2Pers.ComplexFormation;
g2 = modbPers4.ComplexFormation;
diff = g2 - g1;

% Make table and sort
testTable = table(g1, g2, diff);
testTable = sortrows(testTable,"diff","ascend");
final_table.diff = diff;
final_table = sortrows(final_table,"diff","ascend");

% Loop through to calculate group differences
rankSum = zeros(length(g1), 1);
for i = 2:length(g1)
g1Sub = testTable.g1(1:i);
g2Sub = testTable.g2(1:i);
rankSum(i) = ranksum(g1Sub, g2Sub);
end

%% Plotting Complex Specific Information
%}