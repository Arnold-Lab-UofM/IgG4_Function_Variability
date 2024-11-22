function [dfIso] = importConvertedSpiking(var, ag, igg4Spike)
% IMPORT CONVERTED SPIKING Imports IgG4 spiking experimental data
% INPUT:
%   ag: string ["Trimer"]
%       String to determine which antigen to pull from the data
%   var: string ["WT"]
%       Variant of choice
%   igg4Spike: boolean
%       If true, includes the data at the spiking levels (added IgG4)
% OUTPUT:
%   df_converted: table
%       Contains the converted IgG subclass concentrations and the
%       associated multiplex assay measurements

%% Import from file
df = readtable("./Data/Samples/IgG4_spiking_nm.csv", "VariableNamingRule", "preserve");

%% Isolate the section of interest
% Determine if we want the spike or not
if igg4Spike
    n = 72;
else
    n = 24;
end

% Isolate the Antibodies
varNames = string(df.Properties.VariableNames);
isoAgIndex = contains(varNames, ag);
isoVarIndex = contains(varNames, var);
dfPartialIso = df(:, isoAgIndex & isoVarIndex);

% Get ID information
background = df(:, ["Var1", "IgG4 added (nM)"]);

%% Pull Final Data
dfIso = [background(1:n, :), dfPartialIso(1:n, :)];
dfIso = renamevars(dfIso, "Var1", "ID");
end



