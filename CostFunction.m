%%%
            % 01:25 AM Thur 16th Jan 2025 
            % Cost function based on the net present cost (NPC)

%%%
function NPC = CostFunction(X)
    disp('X: ')
    disp(X)
    %% Initialization of the optimization loop 
    N = 1;
    Ndays = 10;
    T = 24; % # of hours in a day 
    Nh = T*Ndays; % Total hours in Ndays
    Np = 24;
    Hstart = 0; % Start hour of simulation
    index = Hstart+1:Hstart+Nh+Np-1;
    dt = 1;

    %% Monte Carlo ---> 1st step
    % For N lifetime years generates scenarios using Monte Carlo with SD = 20 % 
    SD = 0.0;
    % PV
    Ppv_one_year_S = load('Data/Ppv_data.mat');
    Ppv_one_year = Ppv_one_year_S.PPV;
    Ppv_N_years = repmat(Ppv_one_year(index), 1, N);
    Ppv_N_years = X(1)*Ppv_N_years.*(1 + SD*randn(Nh+Np-1, N));
    clear Ppv_one_year Ppv_one_year_S; % Save memory 

    % WT
    Pwind_one_year_S = load('Data/PWind_data.mat');
    Pwind_one_year = Pwind_one_year_S.PWT;
    Pwind_N_years = repmat(Pwind_one_year(index), 1, N);
    Pwind_N_years = X(2)*Pwind_N_years.*(1 + SD*randn(Nh+Np-1, N));
    clear Pwind_one_year Pwind_one_year_S; % Save memory 

    % H2_Heating demand 
    H2_Heating_one_year_S = load('Data/H2Heating_data.mat');
    H2_Heating_one_year = H2_Heating_one_year_S.H2HeDemand;
    H2HEating_N_years = repmat(H2_Heating_one_year(index), 1, N);
    H2HEating_N_years = H2HEating_N_years.*(1 + SD*randn(Nh+Np-1, N)); 
    clear H2_Heating_one_year H2_Heating_one_year_S; % Save memory

    
    % Battery energy storage system (BSS)
    Bs = X(3);
    Qnom_bat = Bs*1; % kWh
    eta_c = sqrt(0.9);
    eta_d = sqrt(0.9);
    Pb_max = 100; % kW
    Pb_min = 0; % kW
    A = 1;
    B = (100/Qnom_bat)*[eta_c -1/eta_d]*dt;
    C = 1; 
    D = zeros(1,2);
    SOC_init = 10;
    SOC_max = 90;
    SOC_min = 10; % Minimum battery SOC (%)
    
    Parameters = [A B C D SOC_init eta_c eta_d Qnom_bat Pb_max Pb_min SOC_max SOC_min];
    Batt.signals.values = Parameters;
    Batt.time = length(Parameters)';
    clear A B C D eta_c eta_d Qnom_bat Parameters;

    % Heating storage
    Hs.SOC_max = 100;   % Hydrogen storage SOC max
    Hs.SOC_min = 1e-10; % Hydrogen storage SOC min
    Hs.SOC_init = 10;    % Hydrogen storage SOC initial
    
    Grid.Pgrid_max_Im = 1000;
    Grid.Pgrid_max_Ex = 1000; % kW
    Grid.Pgrid_min_Ex = 0;

    Grid.Cg_Ex = 0.065*ones(Nh+Np-1, 1)+  0.00*randn(Nh+Np-1, 1); % Grid power cost (per kWh)
    Grid.Cg_Im = 0.065*ones(Nh+Np-1, 1) + 0.00*randn(Nh+Np-1, 1); % Grid revenue (per kWh)
    Cost.Cg_ExF = 0.065*ones(Nh+Np-1, 1)+  0.00*randn(Nh+Np-1, 1); % Grid power cost (per kWh)
    Cost.Cg_ImF = 0.065*ones(Nh+Np-1, 1)+  0.00*randn(Nh+Np-1, 1); % Grid power cost (per kWh)


    % CO2 cost
    CO2.PVE = 93;  % gCO2eq per kWh
    CO2.WTE = 14;  % gCO2eq per kWh
    CO2.GRE = 557; % gCO2eq per kWh
    CO2.NGE = 193; % gCO2eq per kWh
    CO2.ER = 0.7;  % gCO2eq reduction factor
    CO2.required = 1; % 1--> CO2 emission is considered, 0--> Not considered
    CO2.ExportedREConsidered = 1; % 1 --> CO2 is sold to grid, 0 otherwise | Carbon credits

    %% Main loop 
    OPEX = zeros(N, 1);
    % Simulation for the N years 
    % CAPEX calculation per year 
    CAPEX = CAPEXHelper(X);
    Capacities.Hs = X(5);
    Capacities.Pnom_El = X(4);
    for y = 1:N

        [gridCost, boolVal] = MPCMILPSimulation_PSO(Nh, Np, T, Ppv_N_years, Pwind_N_years, H2HEating_N_years,Grid, Cost, Batt, dt, Capacities, Hs, CO2, Ndays);
        % Check for infeasibility in the MILP optimizer 
        if boolVal == -1
            OPEX(y)= 1e25;
        else
            OPEX(y) = OPEX_Helper(gridCost, y); % One year is over
        end
    end    

    % Evaluates NPC among N
    NPC = (sum(OPEX) + CAPEX);

end % N years are over







