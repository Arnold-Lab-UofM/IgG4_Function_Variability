%%% UNIFIED VALIDATION SCRIPT FOR ALL OF OUR AVAILABLE DATA
%%% This script is designed to run the validation personalized to each data
%%% set and generate the desired output plots. This script contains options
%%% for each data set available in our analysis.
%%% Author: Robert Theisen, Arnold Lab, University of Michigan, Ann Arbor

%% Housekeeping
clear all;
close all;
clc;

%% Choosing Validation
% NOTE: Each selection will also require specification of certain features
% (Vaccine, timepoint, antigen, etc. scroll down to the specific case to
% check

% NOTE ON FULL VALIDATION:
%   1) DO NOT run with baseline involved unless it is specifially necessary
%   2) Stick to the same dilution (Pfizer Dose 2 and Dose 3, AZ dose 3)
%           The current import function skips AZ dose 2 when it is included

valChoice = "Baseline Spike";
% Options: "Preliminary Data", "Functional Assays", "Full Data",
%           "Dual Timepoint", "Dual Platform", "Full Validation"
%           "D2/D3 Monoclonal", "Baseline Spike"

%% Other Master Control options
fcr = "FcgRIIA-131H"; % Selection for parameters to pull and graph titles
% Options: "FcgRIIA-131H", "FcgRIIIA-158V"




%% SWITCH STATEMENT
switch(valChoice)
%%%%%%%%%%%%%%%%%%% Validation with Preliminary Data %%%%%%%%%%%%%%%%%%%%%%
    case "Preliminary Data"
        % Import Data
        vacc = "COR Pfizer"; % "COR Pfizer, "CP Pfizer", "COR AZ", "CP AZ"
        tp = "Booster"; % "Baseline", "Dose 1", "Dose 2", *"6mth Dose 2", *"Booster" -- * = COR Pfizer Only

        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcgR3aV"; % "FcgR2aH", "FcgR3aV"
        [df_full, df_vacc, df_converted, ~] = importMidDilution(vacc, tp);

        % Run the Validation
        [pred, actual, spRho] = validateModel(df_converted{:, ["SARS2 RBD IgG1", "SARS2 RBD IgG2", "SARS2 RBD IgG3", "SARS2 RBD IgG4"]},...
                                                        df_converted{:, "SARS2 RBD " + fcr_Data}, fcr);

        % Create and Display Output
        validationPlot(pred, actual, spRho, fcr);

        % Common customizations of output



%%%%%%%%%%%%%%%%%%% Validation with Functional Assays %%%%%%%%%%%%%%%%%%%%%
    case "Functional Assays"
        % Import Data
        ag = "TRIMER S"; % "RBD", "TRIMER S"
        [df_full, df_interest, df_converted, IgG_idx, meas_idx] = importFuncAssays(ag);

        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcR3aV"; % "FcgR2aH", "FcgR3aV"

        % Selected functional assay -- CHECK BEFORE RUNNING
        functAssay = true; % Selection on whether to look at functional assay
        funcAssayChoice = "ADCC - FcγRIIIa activation (AUC)"; % FcR3: "ADCC (Area under curve)", "ADCC - FcγRIIIa activation (AUC)"
                                                              % FcR2: "ADP cells (%ADP)", "ADP beads (phagocytic score)"
        
        % Run the Validation
        if functAssay
            [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 SARS2 " + ag, "IgG2 SARS2 " + ag, "IgG3 SARS2 " + ag, "IgG4 SARS2 " + ag]},...
                                                  df_full{:, funcAssayChoice}, fcr);
        else
            [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 SARS2 " + ag, "IgG2 SARS2 " + ag, "IgG3 SARS2 " + ag, "IgG4 SARS2 " + ag]},...
                                                  df_converted{:, fcr_Data + " SARS2 " + ag}, fcr);
        end

        % Create and Display Output
        if functAssay
            validationPlot(pred, actual, spRho, fcr, funcAssayChoice);
        else
            validationPlot(pred, actual, spRho, fcr);
        end
        
        % Common customizations of output


%%%%%%%%%%%%%%%%%%%%%% Validation with Ruth's Data %%%%%%%%%%%%%%%%%%%%%%%%
    case "Full Data"
        % Import Data
        vacc = "Pfizer"; % "Pfizer, "AstraZeneca"
        tp = "Dose 2"; % "Baseline", "Dose 2", "Dose 3"
        ag = "RBD WT Sino"; % "NP", "S2", "S1", "RBD WT Sino", "RBD alpha", "RBD beta", "RBD delta", 
                          % "Trimer", "RBD WT WEHI", "RBD BA.2", "RBD BA.5"
        refConc = 1.14; % ug/mL
        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcgRIIIaV"; % "FcgRIIaH", "FcgRIIIaV"
        [df_full, df_vacc, df_converted, ~] = importRuth(vacc, tp, ag, refConc);

        % Run the Validation
        [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 SARS2 " + ag, "IgG2 SARS2 " + ag, "IgG3 SARS2 " + ag, "IgG4 SARS2 " + ag]},...
                                                        df_converted{:, fcr_Data + " SARS2 " + ag}, fcr);

        % Create and Display Output
        validationPlot(pred, actual, spRho, fcr);

        % Common customizations of output

%%%%%%%%%%%%%%%%%%%%%% Dual Timepoint Ruth's Data %%%%%%%%%%%%%%%%%%%%%%%%%
    case "Dual Timepoint"
        % Import Data
        vacc = "Pfizer"; % "Pfizer, "AstraZeneca"
        tps = ["Dose 2", "Dose 3"]; % "Baseline", "Dose 2", "Dose 3"
        ag = "RBD alpha"; % "NP", "S2", "S1", "RBD WT Sino", "RBD alpha", "RBD beta", "RBD delta", 
                          % "Trimer", "RBD WT WEHI", "RBD BA.2", "RBD BA.5"
        refConc = 9; % ug/mL
        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcgRIIIaV"; % "FcgRIIaH", "FcgRIIIaV"
        [df_full, df_vacc, df_converted, ~] = importDualTimepoint(vacc, tps, ag, refConc);

        % Run the Validation
        [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 SARS2 " + ag, "IgG2 SARS2 " + ag, "IgG3 SARS2 " + ag, "IgG4 SARS2 " + ag]},...
                                                        df_converted{:, fcr_Data + " SARS2 " + ag}, fcr);

        % Create and Display Output
        ax = validationPlot(pred, actual, spRho, fcr);

%%%%%%%%%%%%%%%%%%%%%% Dual Platform Ruth's Data %%%%%%%%%%%%%%%%%%%%%%%%%%
    case "Dual Platform"
        % Import Data
        vaccs = ["Pfizer", "AstraZeneca"]; % "Pfizer, "AstraZeneca"
        tp = "Dose 2"; % "Baseline", "Dose 2", "Dose 3"
        ag = "RBD alpha"; % "NP", "S2", "S1", "RBD WT Sino", "RBD alpha", "RBD beta", "RBD delta", 
                          % "Trimer", "RBD WT WEHI", "RBD BA.2", "RBD BA.5"
        refConc = 9; % ug/mL
        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcgRIIIaV"; % "FcgRIIaH", "FcgRIIIaV"
        [df_full, df_vacc, df_converted, ~] = importDualPlatform(vaccs, tp, ag, refConc);

        % Run the Validation
        [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 SARS2 " + ag, "IgG2 SARS2 " + ag, "IgG3 SARS2 " + ag, "IgG4 SARS2 " + ag]},...
                                                        df_converted{:, fcr_Data + " SARS2 " + ag}, fcr);

        % Create and Display Output
        ax = validationPlot(pred, actual, spRho, fcr);

%%%%%%%%%%%%%%%%%%%%%% All of Ruth's Data at Once %%%%%%%%%%%%%%%%%%%%%%%%%
    case "Full Validation"
        % Import Data
        vaccs = ["Pfizer", "AstraZeneca"]; % "Pfizer, "AstraZeneca"
        tps = ["Dose 2", "Dose 3"]; % "Baseline", "Dose 2", "Dose 3"
        ags = ["NP", "S2", "S1", "RBD WT Sino", "RBD alpha", "RBD delta",... 
               "Trimer", "RBD WT WEHI", "RBD BA.2", "RBD BA.5"];
        refConc = 9; % ug/mL
        % Selected FcR for specific data set (ensure this aligns with the
        % master FcR set above -- the notation in the files is different)
        fcr_Data = "FcgRIIIaV"; % "FcgRIIaH", "FcgRIIIaV"

        % Temporary iterator (for diagnosing high vs low MFI groups)
        dataCarrier = zeros(20, 4);
        refConcs = 1:20;
        for i = 1:20
            % Run import function
            [df_full, df_vacc, df_converted, medians] = importAll(refConcs(i), vaccs, tps, ags);

            % Sort for specific group (OPTIONAL)
            df_converted = sortImportAll(df_converted, "Pfizer", "Dose 2", "RBD WT Sino");

            % Run Everything
            [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1", "IgG2", "IgG3", "IgG4"]}, df_converted{:, fcr_Data}, fcr);
            dataCarrier(i, :) = spRho;
            fprintf("Finished Reference Concetration %i/%i\n", i, length(refConcs));
            validationPlot(pred, actual, spRho(1), fcr);
        end
        % Create and Display output
        % ax = validationPlot(pred, actual, spRho, fcr);

    case "D2/D3 Monoclonal"
        % Import data
        concMultiplier = 25;
        fcrData = "FcγRIIIaF WT Trimer"; % ["FcγRIIaH WT Trimer", "FcγRIIaR WT Trimer"]
        dose = 3;
        df_converted1 = importConvertedBooster("WT", "Trimer", 2);
        df_converted2 = importConvertedBooster("WT", "Trimer", 3);
        df_converted = [df_converted1; df_converted2];

        % Run validation
        [pred, actual, spRho] = validateModel(df_converted2{:, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]} .* concMultiplier, ...
                                                df_converted2{:, fcrData}, fcr);
        % Plotting
        validationPlot(pred, actual, spRho, fcr, true);

    case "Baseline Spike"
         % Import data
        concMultiplier = 25;
        fcrData = "FcγRIIaH WT Trimer"; % ["FcγRIIaH WT Trimer", "FcγRIIaR WT Trimer"]
        df_converted = importConvertedSpiking("WT", "Trimer", false);

        % Run validation
        [pred, actual, spRho] = validateModel(df_converted{:, ["IgG1 WT Trimer", "IgG2 WT Trimer", "IgG3 WT Trimer", "IgG4 WT Trimer"]} .* concMultiplier, ...
                                                df_converted{:, fcrData}, fcr);
        % Plotting
        validationPlot(pred, actual, spRho, fcr, true);
    
    otherwise
        disp("Please reselect input - no viable option detected");
end