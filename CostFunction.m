%%%
            % 01:25 AM Thur 16th Jan 2025 
            % Cost function based on the net present cost (NPC)

%%%
function NPC = CostFunction(X)
    disp('X: ')
    disp(X)
    %% Initialization of the optimization loop 
    N = 1;
    Ndays = 364;
    T = 24; % # of hours in a day 
    Nh = T*Ndays; % Total hours in Ndays
    Np = 24;
    Hstart = 0; % Start hour of simulation
    index = Hstart+1:Hstart+Nh+Np-1;
    dt = 1;
    PSO_Penalt = Inf;
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
    Grid.Pgrid_max_Ex = 100; % kW
    Grid.Pgrid_min_Ex = 0;

    Grid.Cg_Ex = 0.065*ones(Nh+Np-1, 1)+  0.001*randn(Nh+Np-1, 1); % Grid power cost (per kWh)
    Grid.Cg_Im = 0.065*ones(Nh+Np-1, 1) + 0.001*randn(Nh+Np-1, 1); % Grid revenue (per kWh)

    opex_G = [Grid.Cg_Im Grid.Cg_Ex];


    % CO2 cost
    CO2.PVE = 93;  % gCO2eq per kWh
    CO2.WTE = 14;  % gCO2eq per kWh
    CO2.GRE = 557; % gCO2eq per kWh
    CO2.NGE = 193; % gCO2eq per kWh
    CO2.ER = 0.7;  % gCO2eq reduction factor
    CO2.required = 1; % 1--> CO2 emission is considered, 0--> Not considered
    CO2.ExportedREConsidered = 1; % 1 --> CO2 is sold to grid, 0 otherwise | Carbon credits
    
    Cg = [Grid.Cg_Im', -0.8*Grid.Cg_Ex'];
    Cb = [-Grid.Cg_Ex', -1*Grid.Cg_Im']; % Favour battery over grid
    Csoc = (1/14)*ones(1,Nh+Np-1);
    Cre = [Grid.Cg_Im' Grid.Cg_Im'];
    Cpel = Grid.Cg_Im';
    
    Cre_b = (1/3)*Grid.Cg_Im';
    Cre_g = (1/3)*Grid.Cg_Ex';
    Cre_Ele = (1/3)*Grid.Cg_Im';
    Cg_b = 1*Grid.Cg_Im';
    Cg_Ele = 1*Grid.Cg_Im';
    Cb_Ele = (1/3)*Grid.Cg_Im';
    Cb_g = (1)*Grid.Cg_Im';
    
    % As the prediction window recedes, more data,(Np-1) are required. --> This
    % will be different with shrinking horizon approach
    n_binary = 2;
    nvars = 21;

    % Cost function
    f = [0*Cg, 0*Cb, -Csoc, zeros(1, n_binary*(Nh+Np-1)), -1*Cre, -Csoc, zeros(1,3*(Nh+Np-1)), -0*Cpel, -Cre_b, -Cre_g, -Cre_Ele, Cg_b, Cg_Ele, -Cb_Ele, Cb_g];
    %% Main loop 
    OPEX = zeros(N, 1);
    % Simulation for the N years 
    % CAPEX calculation per year 
    CAPEX = CAPEXHelper(X);
    Capacities.Hs = X(5);
    Capacities.Pnom_El = X(4);
    for y = 1:N            
           [Pg, boolVar] = MPCMILPSimulation_PSO(Nh, Np, Ppv_N_years, Pwind_N_years, H2HEating_N_years,Grid, Batt, dt, Capacities, Hs, CO2, f, nvars);
            % Check for infeasibility in the MILP optimizer 
            if boolVar == -1
                %  Not a successful MILP operation. Discard X from PSO by heavily penalizes it
                OPEX(y) = PSO_Penalt; % To penlized failured PSO selected size
                
            else
                %  Successful MILP operation 
                % Compute grid cost
                GOp = [opex_G(1+Nh*(y-1):Nh*y, 1); opex_G(1+Nh*(y-1):Nh*y, 2)];
                ExOp = [Pg(1:Nh); Pg(Nh + 1:2*Nh)];
                % OPEX_Helper(sum(GOp'*ExOp), y)
                OPEX(y) = OPEX_Helper(sum(GOp'*ExOp), y); % One year is over
            end    
    end 
     % Evaluates NPC among N
    NPC = (sum(OPEX) + CAPEX);
end % N years are over







