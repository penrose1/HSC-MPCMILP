%%%%
        % 3:54 AM 26th Jan 2025 Plotalaan
        % Compute CAPEX 
        % C ---> CAPEX(component_i)
        % X ---> size from PSO
function [capex] = CAPEXHelper(X)
    % CAPEX in vector C

    % C(1) capex of PV [euro/ha]
    % C(2) capex of WT [euro/unit]
    % C(3) capex of Battery [euro/kWh]
    % C(4) capex of Electrolyzer [euro/kW]
    % C(5) capex of HHD storage [euro/kg]

    capex_Stack =320; % CAPEX of electrolyzer stack [euro/kW]
    capex_Ele = 1000; % CAPEX of electrolyzer [euro/kW]
    WACC = WACC_function();
    capex_Ele_total = capex_Ele + capex_Stack*(1/(1+WACC(4))^6 + 1/(1+WACC(4))^11); 
    capex = [297000 90000 250 capex_Ele_total 490];
    
    capex = capex*X';