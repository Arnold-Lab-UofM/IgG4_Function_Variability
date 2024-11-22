function [ax] = validationPlot(pred, actual, spRho, fcr, logScale, funcAssay)
% VALIDATION PLOT Show the output of a model validation run
%   INPUT:
%       Pred: Vector (len = nVaccinees)
%           Model predicted results
%       Actual: Vector (len = nVaccinees)
%           Measured results (either functional assay or multiplex MFIs)
%       spRho: double
%           Spearman correlation coefficient between pred and actual
%       fcr: string
%           Name of the Fc receptor for graph labeling
%       log: boolean
%           If true, log scales the axes
%       funcAssay: boolean = false
%           Name of the functional assay if its being used for testing.
%           Default argument is false.
%   OUTPUT:
%       ax: MATLAB axis object
%           A handle for the plot object (not required to see the graph)

%% Input checking
isFuncAssay = true;
if (nargin <= 5)
    isFuncAssay = false;
end

%% Create output plot
% Plotting
figure;
scatter(pred, actual, [], 'blue', 'filled');
xlabel("Predicted Complex Formation (nM)");
if isFuncAssay
    ylabel(funcAssay)
else
    ylabel("Measured Complex Formation (MFI)");
end
title(fcr + " Validation");

% Modifications
ax = gca;
if (logScale)
    ax.XScale = 'log';
    ax.YScale = 'log';
end

% Placing the Spearman Rho on the plot
if logScale
    NW = [min(xlim) max(ylim)] + [min(xlim) -diff(ylim)]*.05;
    rhoString = sprintf("%1.3f", spRho);
    text(NW(1), NW(2), "\rho = " + rhoString, FontSize=18, VerticalAlignment="top", HorizontalAlignment="left");
else
    NW = [min(xlim) max(ylim)] + [diff(xlim) -diff(ylim)]*.05;
    rhoString = sprintf("%1.3f", spRho);
    text(NW(1), NW(2), "\rho = " + rhoString, FontSize=18, VerticalAlignment="top", HorizontalAlignment="left");
end


end