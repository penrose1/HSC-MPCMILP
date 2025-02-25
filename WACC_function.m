function [WACC] = WACC_function()
    % General
    D = 0.7; % Debt ratio
    E = 0.3; % Equity ratio
    k_D = 0.025; % Interest rate debt financing
    k_E = 0.15; % Required return on equity
    % Infl = 0.015; % Inflation
    CT = 0.25; % Corporate tax
    % PV 
    DPV = 0.9; % Debt ratio of PV 
    EPV = 0.1; % Equity ratio on PV 
    k_DPV = 0.01; % Intarest rate debt financing on PV
    k_EPV = 0.09; % Required return on equity on PV
    % WT
    DWT = 0.8; % Debt ratio of WT
    EWT = 0.2; % Equity ratio on WT 
    k_DWT = 0.01; % Intarest rate debt financing on WT
    k_EWT = 0.11; % Required return on equity on WT

    WACC_G = (D*k_D*(1-CT)+E*k_E)/(D+E);
    WACC_PV = (DPV*k_DPV*(1-CT)+EPV*k_EPV)/(DPV+EPV);
    WACC_WT = (DWT*k_DWT*(1-CT)+EWT*k_EWT)/(DWT+EWT);

    % return a row vector of WACC: [WACC_PV WACC_WT ...]
    WACC = [WACC_PV WACC_WT WACC_G WACC_G WACC_G];
end
