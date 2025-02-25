%%%% 
        % 9:33 AM 18th Jan 2025 Plotolaan 
        % Simulates microgrids real time operations 
        % Adopted from MainProgam_v1.m 

function [] = MainProgramStandaloneMILP_PSO(Grid, dt_f, T_f, Ndays_f, Nh_f, index_f, Hs, CO2, Np, f, nvars, simType, BT)
    dt = dt_f; % Sampling time of 1h
    Ndays = Ndays_f;
    T = T_f; %  Hours in a day
    Nh = Nh_f; % Total hours in Ndays --> Simulation time
    index = index_f;
    %% Electricity prices
    Cg_Ex = Grid.Cg_Ex; % Grid power cost (per kWh)
    Cg_Im = Grid.Cg_Im; % Selling price (per kWh)
       
    %% Battery energy storage system (BSS)
    Qnom_bat = X_in(3)*1; % kWh 
    eta_c = sqrt(0.9);
    eta_d = sqrt(0.9);
    Pb_max = 1000; % kW
    Pb_min = 0; % kW
    A = 1;
    B = (100/Qnom_bat)*[eta_c -1/eta_d]*dt;
    C = 1; 
    D = zeros(1,2);
    SOC_max = 90;
    SOC_min = 10; % Minimum battery SOC (%)
    
    Parameters = [A B C D soc_0_in eta_c eta_d Qnom_bat Pb_max Pb_min SOC_max SOC_min];
    Batt.signals.values = Parameters;
    Batt.time = length(Parameters)';
    clear A B C D SOC_init eta_c eta_d Qnom_bat Parameters;

    %% Generation with stochasticity [per unit]
    % The aim here is to account for stochasticity in generation whereby
    % the prediction generation data (from Tromp's CASE) are used to generate
    % actual generation (real time). The assumption is that, the actual generation is
    % varying normally with the mean of Xpv and Xwind and with the standard deviation
    % of SDpv and SDwind for PV and wind respectively.
    % Generate actual PV profile ( PV_predicted + ~N(Xpv, SDpv)
    Ppv_temp = load('Data/Ppv_data.mat');
    P_pv_prediction= Ppv_temp.PPV(index);
    Ppv_actual = 0.0*randn(Nh+Np-1,1) + P_pv_prediction;
    
    % Generate actual wind profile ( wind_predicted + ~N(Xpv, SDpv)
    Pwind_temp = load('Data/PWind_data.mat');
    % Calculate wind power based on speed 
    % Some ifs 
    P_wind_prediction = Pwind_temp.PWT(index);
    Pwind_actual = 0.0*randn(Nh+Np-1,1) + P_wind_prediction;
    clear Ppv_temp Pwind_temp;
    
    %% Demand profiles ----> No uncertainty is considered for this profile. i.e prediction profile == actual profile
    % 1. H2 for heating (H2Heating_demand.mat contains demand for the whole year)
    % This demand must be served by the plant **ALL TIME 
    H2_heating = load('Data/H2Heating_data.mat'); % In kg
    H2_heating = H2_heating.H2HeDemand(index);
    % H2kWhkg = load('Data/H2kWhkg.mat');
    % P_electrolyzer = H2kWhkg.H2kWhkg*H2_heating/Electrolyzer.eta_e; % From kg to kWh
    HeatingD = H2_heating; % kg to kWh
    % Other demands will be added here later on 
    
    GenerationDemandProfilesPlotting(HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Nh)
    
    init_soc.Battery = SOC_init;
    init_soc.Hs = Hs.SOC_init;

    PSO_main(Grid, dt_f, T_f, Ndays_f, Nh_f, index_f, Hs, CO2, Np, f, nvars, simType, BT);
    
    end 