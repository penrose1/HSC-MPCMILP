%%%%   
        % 3:48 AM 16th Thur 2025 
        % Calculates T hours of operation 
        % 
function [opex] = OPEX_Helper(opex_Grid, n)
    WACC = WACC_function();
    WACC = [WACC WACC(5)]; % Add grid WACC. Same as the last 3 elements in WACC vector
    Discounted_WACC = zeros(1,length(WACC));
    opex_PV = 0.025*297000;
    opex_WT = 0.03*90000;
    opex_Batt = 0.015*250;
    opex_Ele = 0.03*1000;
    opex_HHS = 0.01*490; 
    for j = 1:length(WACC)
        Discounted_WACC(j) = 1/power((1+WACC(j)),n); % To be used for discounting opex 
    end
    Cost = [opex_PV opex_WT opex_Batt opex_Ele opex_HHS opex_Grid];
    opex= Discounted_WACC*Cost';
