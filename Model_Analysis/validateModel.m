function [pred, actual, spRho] = validateModel(subclassConcs, fcrMFI, fcr)
%VALIDATE MODEL Run the model and compare output to the FcR assay
%   INPUT:
%       subclasConcs: array (nVaccinees, 4)
%           Array containing the converted subclass concentrations for each
%           of the vaccinees
%       fcrMFI: vector (nVaccinees, 1)
%           Vector containing the MFIs for the FcR of choise for the
%           validation
%       fcr: string ("FcgRIIA-131H", "FcgRIIIA-158V")
%           Dictates which FcgR we are using for parameter collection
%   OUTPUT:
%       pred: vector (nVaccinees)
%           Model predicted FcR complex formation (in nM)
%       actual: vector (nVaccinees)
%           Measured FcR complex formation values
%       spRho
%           Spearman correlation coefficient of models predictive accuracy


%% Collect Background data
nVaccinees = size(subclassConcs, 1);

%% Output data carriers
pred = zeros(nVaccinees, 1);

%% Iterate
for i = 1:nVaccinees
    % Get personalized parameter set
    [baseline_params, paramnames, complexes] = getMonoclonalParams(subclassConcs(i, :), fcr);

    % Run Simulation
    [ybase, steadystate, complexes] = Simulate(baseline_params, paramnames, complexes, fcr);

    % Check steady state
    if(steadystate == false)
        disp("Steady State Failed");
    end

    % Collect output
    pred(i) = ybase(33);

end

%% Postprocessing
actual = fcrMFI;
spRho1 = corr(pred, actual, type="Spearman");
logSpRho = corr(log10(pred), log10(actual), type="Spearman");
perRho = corr(pred, actual, type="Pearson");
logPerRho = corr(log10(pred), log10(actual), type="Pearson");
%spRho = [spRho1; logSpRho; perRho; logPerRho];
spRho = spRho1;
%fprintf(fcr + " Spearman Rho: %1.3f %1.3f %1.3f %1.3f\n", spRho');

end