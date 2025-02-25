%%%
            % 01:25 AM Thur 16th Jan 2025 
            % Cost function based on the net present cost (NPC)

%%%
function NPC = CostFunction(X, Config, f)
    disp('X: ')
    disp(X)
    %% Initialization of the optimization loop 
    N = Config.Years;
    Ndays = Config.Ndays; % Number of days simulated
    T = 24; % Resolution 
    Nh = T*Ndays;% 8760; % Number of hours
    index = Config.index;
    soc_0 = Config.SOC_Batt;
    Grid_OPEX = zeros(Ndays, 1);
    OPEX = zeros(N, 1);
    %% Monte Carlo ---> 1st step
    % For N lifetime years generates scenarios using Monte Carlo with SD = 20 % 
    SD = 0.0;
    % PV
    Ppv_one_year_S = load('Data/Ppv_data.mat');
    Ppv_one_year = Ppv_one_year_S.PPV;
    Ppv_N_years = repmat(Ppv_one_year(index), 1, N);
    Ppv_N_years = Ppv_N_years.*(1 + SD*randn(Nh, N));
    clear Ppv_one_year Ppv_one_year_S; % Save memory 

    % WT
    Pwind_one_year_S = load('Data/PWind_data.mat');
    Pwind_one_year = Pwind_one_year_S.PWT;
    Pwind_N_years = repmat(Pwind_one_year(index), 1, N);
    Pwind_N_years = Pwind_N_years.*(1 + SD*randn(Nh, N));
    clear Pwind_one_year Pwind_one_year_S; % Save memory 

    % H2_Heating demand 
    H2_Heating_one_year_S = load('Data/H2Heating_data.mat');
    H2_Heating_one_year = H2_Heating_one_year_S.H2HeDemand;
    H2_N_years = repmat(H2_Heating_one_year(index), 1, N);
    H2_N_years = H2_N_years.*(1 + SD*randn(Nh, N)); 
    clear H2_Heating_one_year H2_Heating_one_year_S; % Save memory

    %% Electricity prices
    % Cg_Im = 0.065*(1 + 0.1*randn(Nh, 1)); % Buying price (per kWh)
    % Cg_Ex = 0.065*(1 + 0.1*randn(Nh, 1)); % Selling price (per kWh)
    opex_G = [0.065*(1 + 0*randn(Nh, 1)) 0.065*(1 + 0.0*randn(Nh, 1))];

    %% Main loop 
    % Simulation for the N years 
    % CAPEX calculation per year 
    CAPEX = CAPEXHelper(X);
    for y = 1:N    
        
        for d = 1:Ndays 
            % Solve for optimal problem for T hours
            % Send profiles and components sizes to optimizer&scheduler (MILP). 
            % The function returns a dispached solution  
            [Xdispached, boolVar, soc_0] = MainProgramMILPCalledByCostFunc(X, soc_0, Ppv_N_years(1+T*(d-1):T*d, y), Pwind_N_years(1+T*(d-1):T*d, y), H2_N_years(1+T*(d-1):T*d, y), dt, T, ...
                opex_G(1+T*(d-1):T*d, 1), opex_G(1+T*(d-1):T*d, 2));
            % OPEX calculation per day
            % Pass only grid power returned from MILP optimizer to compute OPEX per year 

            % ************START****************DECISION VARIABLE ORDER************START****************

                                              %** Continuous **
                                                % P_grid_in(t) 1
                                                % P_grid_out(t) 2
                                                % Pb_charge(t) 3
                                                % Pb_discharge(t) 4
                                                % SOC(t) 5
                        
                                              %** Binary **
                                                % delta_g 6
                                                % delta_b 7
                        
                                              %** Continuous ** newly added
                                                % lambda_pv(t) ---> [0-1] 8
                                                % lambda_w(t)  ---> [0-1] 9
                                                % SOC_Hs(k) 10
                                                % H_El2Hs(k) 11
                                                % H_El2Hd(k) 12
                                                % H_Hs2Hd(k) 13
                                                % Pel_in(k) 14
                                                % P_RE_B(k) 15
                                                % P_RE_g(k) 16
                                                % P_RE_Ele(k) 17
                                                % P_g_B(k) 18
                                                % P_g_Ele(k) 19
                                                % P_B_Ele(k) 20
                                                % P_B_g(k) 21
                                                
            % ************END****************DECISION VARIABLE ORDER************END****************
            % Check for infeasibility in the MILP optimizer 
            if boolVar == -1
                %  Not a successful MILP operation. Discard X from PSO by heavily penalizes it
                Grid_OPEX(d) = Config.PSO_Penalt; % To penlized failured PSO selected size
                
            else
                %  Successful MILP operation 
                % Compute grid cost
                GOp = [opex_G(1+T*(d-1):T*d, 1); opex_G(1+T*(d-1):T*d, 2)];
                ExOp = [Xdispached(1:T); Xdispached(T + 1:2*T)];
                Grid_OPEX(d) = sum(GOp'*ExOp); % OPEX per day % sum(OPEXPer_T_Helper(GOp, ExOp)); 

            end
        end % One year is over
        OPEX(y) = OPEX_Helper(Grid_OPEX);
    end % N years are over

    % Evaluates NPC among N
    NPC = (sum(OPEX) + CAPEX);
end % End function







