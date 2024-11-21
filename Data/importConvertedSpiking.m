function [dfIso] = importConvertedSpiking(var, ag, igg4Spike)
% IMPORT CONVERTED SPIKING Imports IgG4 spiking experimental data -
% converted by Carissa
% INPUT:
%   ag: string ("RBD", "Trimer")
%       String to determine which antigen to pull from the data
%   var: string ("WT", "Delta", "BA.2", "BA.5")
%       Variant of choice
%   fcr: string ("FcgR$" where $: I, IIA-131H, IIA-131R, IIIA-158V,
%        IIIA-158F)
%        Determines which Fc receptor MFIs to pull from the data
%   igg4Spike: boolean
%        If true, includes the data at the spiking levels (added IgG4)
% OUTPUT:
%   df_converted: table
%       Same as df_iso but with nM values for IgG subclass concentrations

%% Import from file
df = readtable("./Data/Samples/IgG4_spiking_nm.csv", "VariableNamingRule", "preserve");

% %% Align Fcr input with data
% % Slightly different naming convention in the code as in the data file
% % (notice the gamma character)
% if fcr == "FcgRI"
%     fcrData = "FcγRI";
% elseif fcr == "FcgRIIA-131H"
%     fcrData = "FcγRIIaH";
% elseif fcr == "FcgRIIA-131R"
%     fcrData = "FcγRIIaR";
% elseif fcr == "FcgRIIb"
%     fcrData = "FcγRIIb";
% elseif fcr == "FcgRIIIA-158V"
%     fcrData = "FcγRIIIaV";
% elseif fcr == "FcgRIIIA-158F"
%     fcrData = "FcγRIIIaF";
% else
%     disp("Improper FcR Input - try again");
%     return;
% end

%% Isolate the section of interest
% Determine if we want the spike or not
if igg4Spike
    n = 72;
else
    n = 24;
end

% Isolate the Antibodies (only have WT Trimer for now)
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



