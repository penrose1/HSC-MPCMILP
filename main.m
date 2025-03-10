
clear,clc, close all;
sim = 'pso';
simulation.simType = "mpcmilp_and_shrinkingH";
simulation.stochastic = 0;
newPrices = 0;
switch sim
    case 'pso'
        [x,fval,exitflag,output,points] = PSO_main();

    case 'nH_sim'
        
        Np = 24;
        dt = 1;
        Ndays = 1;
        T = 24; % # of hours in a day 
        Nh = T*Ndays; % Total hours in Ndays
        Hstart = T*0; % Start hour of simulation
        index = Hstart+1:Hstart+Nh+Np-1;
        Capacities.Apv = 0;       % Number of PV hectares 
        Capacities.WT = 4;       % Number of wind turbines
        Capacities.Bs = 250;      % Battery kWh 
        Capacities.Pnom_El = 261.888; % Electrolyzer capacity in kW
        Capacities.Hs = 17582.69628119343;      % H2 storage in kg
        
        B.SOC_max = 90;
        B.SOC_min = 10;
        B.SOC_init = 10;
        B.Pb_max = 100;
        B.Pb_min = 0;
        Hs.SOC_max = 100;   % Hydrogen storage SOC max
        Hs.SOC_min = 1e-10; % Hydrogen storage SOC min
        Hs.SOC_init = 10;   % Hydrogen storage SOC initial
        
        Grid.Pgrid_max_Im = 1000;
        Grid.Pgrid_max_Ex = 1000; % kW
        Grid.Pgrid_min_Ex = 0;
        
        
        if newPrices == 1
            simulation.folder = "Plots-varyingDP";
            Grid.Cg_Ex = 0.065*ones(Nh+Np-1, 1)+  0.001*randn(Nh+Np-1, 1); % Grid power cost (per kWh)
            Cg_ExF = 0.9*Grid.Cg_Ex;
            Grid.Cg_Im = 0.065*ones(Nh+Np-1, 1) + 0.001*randn(Nh+Np-1, 1); % Grid revenue (per kWh)
            Cg_ImF = 0.9*Grid.Cg_Im;
            Cg_Ex = Grid.Cg_Ex;
            Cg_Im = Grid.Cg_Im;
            save('Data/Cg_Ex.mat','Cg_Ex')
            save('Data/Cg_Im.mat','Cg_Im')
            sd = 0.15;
            rand_pv_wt = 0+sd*randn((Nh)*Np,1);
            save('Data/rand_pv_wt.mat',"rand_pv_wt");
        else
            simulation.folder = "mpcmilp_and_shrinkingHx-stochastic";
            Cg_Im = load('Data/Cg_Im.mat');
            Grid.Cg_Im = 1*Cg_Im.Cg_Im;
            Cg_ImF = Cg_Im.Cg_Im;
            Cg_Ex = load('Data/Cg_Ex.mat');
            Grid.Cg_Ex = 1*Cg_Ex.Cg_Ex;
            Cg_ExF = Cg_Ex.Cg_Ex;
            Cost.Cg_ExF = Cg_ExF;
            Cost.Cg_ImF = Cg_ImF; % allow to flactuate grid prices separetely
            rand_pv_wt = load('Data/rand_pv_wt.mat','rand_pv_wt');
            rand_pv_wt = rand_pv_wt.rand_pv_wt;
        end
        
        % CO2 cost
        CO2.PVE = 93;  % gCO2eq per kWh
        CO2.WTE = 14;  % gCO2eq per kWh
        CO2.GRE = 557; % gCO2eq per kWh
        CO2.NGE = 193; % gCO2eq per kWh
        CO2.ER = 0.7;  % gCO2eq reduction factor
        CO2.required = 1; % 1--> CO2 emission is considered, 0--> Not considered
        CO2.ExportedREConsidered = 1; % 1 --> CO2 is sold to grid, 0 otherwise | Carbon credits
        
        Cg = [Cg_ImF', -1*Cg_ExF'];
        Cb = [-Cg_ExF', -1*Cg_ImF']; % Favour battery over grid
        Csoc = (1/14)*ones(1,Nh+Np-1);
        Cre = [Cg_ImF' Cg_ImF'];
        Cpel = Cg_ImF';
        
        Cre_b = (1/3)*Cg_ImF';
        Cre_g = (1/3)*Cg_ExF';
        Cre_Ele = (1)*Cg_ImF';
        Cg_b = 1*Cg_ImF';
        Cg_Ele = 1*Cg_ImF';
        Cb_Ele = (1/1)*Cg_ImF';
        Cb_g = (0)*Cg_ImF';
        
        % As the prediction window recedes, more data,(Np-1) are required. --> This
        % will be different with shrinking horizon approach
        n_binary = 2;
        nvars = 21;
        nvars_vector = nvars*Nh;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% VECTOR OF DECISION VERIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% VECTOR OF DECISION VERIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        f.MPCMILP = [0*Cg, 0*Cb, -0*Csoc, zeros(1, n_binary*(Nh+Np-1)), -1*Cre, -0*Csoc, zeros(1,3*(Nh+Np-1)), -1*Cpel, -Cre_b, Cre_g, -Cre_Ele, Cg_b, Cg_Ele, Cb_Ele, Cb_g];
        f.MILP = [0*Cg(1:2*Nh), 0*Cb(1:2*Nh), -Csoc(1:Nh), zeros(1, n_binary*Nh), -1*Cre(1:2*Nh), -Csoc(1:Nh), zeros(1,3*Nh), 0*Cpel(1:Nh), -Cre_b(1:Nh), -Cre_g(1:Nh), -Cre_Ele(1:Nh), Cg_b(1:Nh), Cg_Ele(1:Nh), -Cb_Ele(1:Nh), Cb_g(1:Nh)];
        MainProgramStandaloneMILP(Capacities, Grid, Cost, dt, T, Ndays, Nh, index, Hs, CO2, Np, f, nvars, simulation, B, rand_pv_wt);
        
    otherwise 
        clc, clear, close all;
        disp('Bad selection. No program was loaded.')
        return
end