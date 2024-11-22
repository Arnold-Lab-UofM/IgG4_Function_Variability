function [pred, actual, spRho] = validateModel(subclassConcs, fcrMFI, fcr)
%VALIDATE MODEL Run the model and compare output to the FcgR assay
%   INPUT:
%       subclassConcs: array (nVaccinees, 4)
%           Array containing the converted subclass concentrations for each
%           of the vaccinees
%       fcrMFI: vector (nVaccinees, 1)
%           Vector containing the MFIs for the FcR of choice for the
%           validation
%       fcr: string ["FcgRIIA-131H", "FcgRIIA-131R"]
%           Dictates which FcgR we are using for parameter collection
%   OUTPUT:
%       pred: vector (nVaccinees)
%           Model predicted FcR complex formation (in nM)
%       actual: vector (nVaccinees)
%           Measured FcR complex formation values (in MFI)
%       spRho: double
%           Spearman rho correlation coefficient of the predicted and
%           actual


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
spRho = corr(pred, actual, type="Spearman");
end