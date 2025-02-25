function [] = MainProgramStandaloneMILP(Capacities, Grid, dt_f, T_f, Ndays_f, Nh_f, index_f, Hs, CO2, Np, f, nvars, sim, BT)
folder = sim.folder;
dt = dt_f; % Sampling time of 1h
Ndays = Ndays_f;
T = T_f; %  Hours in a day
Nh = Nh_f; % Total hours in Ndays --> Simulation time
index = index_f;
%% Electricity prices
Cg_Ex = Grid.Cg_Ex; % Grid power cost (per kWh)
Cg_Im = Grid.Cg_Im; % Selling price (per kWh)

%% PV power plant 
Apv = Capacities.Apv; % Number of PV hectares 
%% Wind power plant 
WT = Capacities.WT; % Number of wind turbines

%% Battery energy storage system (BESS)
Bs = Capacities.Bs;
Qnom_bat = Bs*1; % kWh
eta_c = sqrt(0.9);
eta_d = sqrt(0.9);
Pb_max = BT.Pb_max; % kW
Pb_min = BT.Pb_min; % kW
A = 1;
B = (100/Qnom_bat)*[eta_c -1/eta_d]*dt;
C = 1; 
D = zeros(1,2);
SOC_init = BT.SOC_init;
SOC_max = BT.SOC_max;
SOC_min = BT.SOC_min; % Minimum battery SOC (%)

Parameters = [A B C D SOC_init eta_c eta_d Qnom_bat Pb_max Pb_min SOC_max SOC_min];
Batt.signals.values = Parameters;
Batt.time = length(Parameters)';
clear A B C D eta_c eta_d Qnom_bat Parameters;

%% Generation with stochasticity 
% The aim here is to account for stochasticity in generation whereby
% the prediction generation data (from Tromp's CASE) are used to generate
% actual generation (real time). The assumption is that, the actual generation is
% varying normally with the mean of Xpv and Xwind and with the standard deviation
% of SDpv and SDwind for PV and wind respectively.
% Generate actual PV profile ( PV_predicted + ~N(Xpv, SDpv)
Ppv_temp = load('Data/Ppv_data.mat');
P_pv_prediction= Apv*Ppv_temp.PPV(index);
Ppv_actual = 0.0*randn(Nh+Np-1,1) + P_pv_prediction;

% Generate actual wind profile ( wind_predicted + ~N(Xpv, SDpv)
Pwind_temp = load('Data/PWind_data.mat');
% Calculate wind power based on speed 
% Some ifs 
P_wind_prediction = WT*Pwind_temp.PWT(index);
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
%% Simulation 
switch sim.simType
    case 'milp'
        % MILP open loop 
        x = MILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im(1:Nh), Cg_Ex(1:Nh), Grid, Batt, dt, Capacities, Hs, CO2, f.MILP);
        %% Ploting/visualization of results 
        % Invoke helper function to plot the result of control setpoints from the
        % controller (MILP)
        PlotingMILP_v5(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);
    
    case 'mpcmilp'
        tic;
        x = MPCMILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f.MPCMILP, nvars);    
        
        elapsedTime = toc;
        % Display elapsed time in milliseconds
        elapsedTime_ms = elapsedTime * 1000;  % Convert seconds to milliseconds
        disp(['Elapsed time: ', num2str(elapsedTime_ms, '%.3f'), ' milliseconds']);
        disp(['Elapsed time: ', num2str(elapsedTime, '%.3f'), ' seconds']);
        disp(['Elapsed time: ', num2str(elapsedTime/60, '%.3f'), ' Minutes']);

        % MCPMILP results 
        PlotingMILP_v5(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);

    case 'both'
        % MILP
        [x1, milp] = MILPSimulation(Nh, Np, T, Ppv_actual(1:Nh), Pwind_actual(1:Nh), HeatingD(1:Nh), Cg_Im(1:Nh), Cg_Ex(1:Nh), Grid, Batt, dt, Capacities, Hs, CO2, f.MILP); 
        % Invoke helper function to plot the result of control setpoints from the controller (MILP)
        PlotingMILP_v5(x1, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);

        % MPC-MILP
        [x, mpcmilp] = MPCMILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f.MPCMILP, nvars);    
        PlotingMILP_v5(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);
        Results_Table(milp, mpcmilp) % Plotting for both MILP and MPC
    case 'multi-Np'
        % Run the same simulation with different prediction horizon 
        % Storage vectors
        NNp = 4;
        % CO2
        PV_co2 = zeros(1,NNp);
        WT_co2 = zeros(1,NNp);
        GR_co2 = zeros(1,NNp);

        % Grid Revenues/Cost
        Im = zeros(1,NNp);
        Ex = zeros(1,NNp);

        % SOC end of day
        SOC_Batt = zeros(1,NNp);
        SOC_HDS = zeros(1,NNp);
        
        % RE usage 
        PVUsage = zeros(1,NNp);
        WTUsage = zeros(1,NNp);
        indx = [24];
        k = 1;
        for PredictionL = indx
            [x, mpcmilp, elapsedTime] = MPCMILPSimulation(Nh, PredictionL, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f.MPCMILP, nvars);
            PlotingMILP_v5_saving(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc,folder, PredictionL);
            PV_co2(k) = mpcmilp.co2.PVE;
            WT_co2(k) = mpcmilp.co2.WTE;
            GR_co2(k) = mpcmilp.co2.GRE;

            Im(k) = mpcmilp.Im;
            Ex(k) = mpcmilp.Ex;

            SOC_Batt(k) = mpcmilp.SOC_Batt;
            SOC_HDS(k) = mpcmilp.SOC_HDS;

            PVUsage(k) =  mpcmilp.PVUsage;
            WTUsage(k) =  mpcmilp.WTUsage;
            k= k+ 1;
        end
        % Display elapsed time in milliseconds
        elapsedTime_ms = elapsedTime * 1000;  % Convert seconds to milliseconds
        disp(['Elapsed time: ', num2str(elapsedTime_ms, '%.3f'), ' milliseconds']);
        disp(['Elapsed time: ', num2str(elapsedTime, '%.3f'), ' seconds']);
        disp(['Elapsed time: ', num2str(elapsedTime/60, '%.3f'), ' Minutes']);

        % Histograms showing the effect of varying prediction horizon 
        annotation.title = 'CO2 Emission';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'PVE , WTE , GRE [gco2perkWh]';
        annotation.label={"3h", "6h", "12h", "24h"};
        Var(1,:) = [PV_co2(1) PV_co2(2) PV_co2(3) PV_co2(4)];
        Var(2,:) = [WT_co2(1) WT_co2(2) WT_co2(3) WT_co2(4)];
        Var(3,:) = [GR_co2(1) GR_co2(2) GR_co2(3) GR_co2(4)];
        Multi_PredictionHorizonPlotting(1:3, annotation, Var, {"PV ", " WT ", " GR"},folder+"/co2.svg");
        clear Var;
        annotation.title = 'Grid usage';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'Pg_{Im} , Pg_{Ex} [kW]';
        annotation.label={"3h", "6h", "12h", "24h"};
        Var(1,:) = [Im(1) Im(2) Im(3) Im(4)];
        Var(2,:) = [Ex(1) Ex(2) Ex(3) Ex(4)];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"P_{Im} ", " P_{Ex}"},folder+"/Griduse.svg");

        clear Var;
        annotation.title = 'Storage SOC';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'SOC_{Batt} , SOC_{HS} [%]';
        annotation.label={"3h", "6h", "12h", "24h"};
        Var(1,:) = [SOC_Batt(1) SOC_Batt(2) SOC_Batt(3) SOC_Batt(4)];
        Var(2,:) = [SOC_HDS(1) SOC_HDS(2) SOC_HDS(3) SOC_HDS(4)];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"SOC_{Batt} ", " SOC_{HDS}"},folder+"/SOC.svg");

        clear Var;
        annotation.title = 'RE usage';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'RE available / RE used [%]';
        annotation.label={"3h", "6h", "12h", "24h"};
        Var(1,:) = [PVUsage(1) PVUsage(2) PVUsage(3) PVUsage(4)];
        Var(2,:) = [WTUsage(1) WTUsage(2) WTUsage(3) WTUsage(4)];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"PV ", " WT"},folder+"/usage.svg");


    case 'multi-Np-shrinking'
        % Run the same simulation with different prediction horizon 
        % Storage vectors
        NNp = 1;
        % CO2
        PV_co2 = zeros(1,NNp);
        WT_co2 = zeros(1,NNp);
        GR_co2 = zeros(1,NNp);

        % Grid Revenues/Cost
        Im = zeros(1,NNp);
        Ex = zeros(1,NNp);

        % SOC end of day
        SOC_Batt = zeros(1,NNp);
        SOC_HDS = zeros(1,NNp);
        
        % RE usage 
        PVUsage = zeros(1,NNp);
        WTUsage = zeros(1,NNp);
        indx = 24;
        k = 1;
        for PredictionL = indx
            [x, mpcmilp, elapsedTime] = MPCMILPSimulationShrinkingHorizon(Nh, PredictionL, T, Ppv_actual(1:Nh), Pwind_actual(1:Nh), HeatingD(1:Nh), Cg_Im(1:Nh), Cg_Ex(1:Nh), Grid, Batt, dt, Capacities, Hs, CO2, f.MILP, nvars, Ndays);

            PV_co2(k) = mpcmilp.co2.PVE;
            WT_co2(k) = mpcmilp.co2.WTE;
            GR_co2(k) = mpcmilp.co2.GRE;

            Im(k) = mpcmilp.Im;
            Ex(k) = mpcmilp.Ex;

            SOC_Batt(k) = mpcmilp.SOC_Batt;
            SOC_HDS(k) = mpcmilp.SOC_HDS;

            PVUsage(k) =  mpcmilp.PVUsage;
            WTUsage(k) =  mpcmilp.WTUsage;
            k= k+ 1;
        end

        % Display elapsed time in milliseconds
        elapsedTime_ms = elapsedTime * 1000;  % Convert seconds to milliseconds
        disp(['Elapsed time: ', num2str(elapsedTime_ms, '%.3f'), ' milliseconds']);
        disp(['Elapsed time: ', num2str(elapsedTime, '%.3f'), ' seconds']);
        disp(['Elapsed time: ', num2str(elapsedTime/60, '%.3f'), ' Minutes']);

        PlotingMILP_v5(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);
        %
        % annotation.title = 'CO2 Emission contribution';
        % annotation.xlabel = 'Time(Day)';
        % annotation.ylabel = 'PVE WTE GRE (gco2perkWh)';
        % annotation.label={"3h", "6h", "12h", "24h"};
        % Var(1,:) = [PV_co2(1) PV_co2(2) PV_co2(3) PV_co2(4)];
        % Var(2,:) = [WT_co2(1) WT_co2(2) WT_co2(3) WT_co2(4)];
        % Var(3,:) = [GR_co2(1) GR_co2(2) GR_co2(3) GR_co2(4)];
        % Multi_PredictionHorizonPlotting(1:3, annotation, Var, {"PV", "WT", "GR"});
        % clear Var;
        % annotation.title = 'CO2 Emission contribution';
        % annotation.xlabel = 'Time(Day)';
        % annotation.ylabel = 'PVE WTE GRE (gco2perkWh)';
        % annotation.label={"3h", "6h", "12h", "24h"};
        % Var(1,:) = [Im(1) Im(2) Im(3) Im(4)];
        % Var(2,:) = [Ex(1) Ex(2) Ex(3) Ex(4)];
        % Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"Im", "Ex"});
        % 
        % clear Var;
        % annotation.title = 'CO2 Emission contribution';
        % annotation.xlabel = 'Time(Day)';
        % annotation.ylabel = 'PVE WTE GRE (gco2perkWh)';
        % annotation.label={"3h", "6h", "12h", "24h"};
        % Var(1,:) = [SOC_Batt(1) SOC_Batt(2) SOC_Batt(3) SOC_Batt(4)];
        % Var(2,:) = [SOC_HDS(1) SOC_HDS(2) SOC_HDS(3) SOC_HDS(4)];
        % Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"SOC_{Batt}", "SOC_{HDS}"});
        % 
        % clear Var;
        % annotation.title = 'CO2 Emission contribution';
        % annotation.xlabel = 'Time(Day)';
        % annotation.ylabel = 'PVE WTE GRE (gco2perkWh)';
        % annotation.label={"3h", "6h", "12h", "24h"};
        % Var(1,:) = [PVUsage(1) PVUsage(2) PVUsage(3) PVUsage(4)];
        % Var(2,:) = [WTUsage(1) WTUsage(2) WTUsage(3) WTUsage(4)];
        % Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"PVUsage", "WTUsage"});

    case 'mpcmilp_and_shrinkingH'
        % Run the same simulation with different prediction horizon 
        % Storage vectors
        NNp = 2;
        % CO2
        PV_co2 = zeros(1,NNp);
        WT_co2 = zeros(1,NNp);
        GR_co2 = zeros(1,NNp);

        % Grid Revenues/Cost
        Im = zeros(1,NNp);
        Ex = zeros(1,NNp);

        % SOC end of day
        SOC_Batt = zeros(1,NNp);
        SOC_HDS = zeros(1,NNp);
        
        % RE usage 
        PVUsage = zeros(1,NNp);
        WTUsage = zeros(1,NNp);

        [~, mpcmilp,~] = MPCMILPSimulationShrinkingHorizon(Nh, Np, T, Ppv_actual(1:Nh), Pwind_actual(1:Nh), HeatingD(1:Nh), Cg_Im(1:Nh), Cg_Ex(1:Nh), Grid, Batt, dt, Capacities, Hs, CO2, f.MILP, nvars, Ndays);

        PV_co2(1) = mpcmilp.co2.PVE;
        WT_co2(1) = mpcmilp.co2.WTE;
        GR_co2(1) = mpcmilp.co2.GRE;

        Im(1) = mpcmilp.Im;
        Ex(1) = mpcmilp.Ex;

        SOC_Batt(1) = mpcmilp.SOC_Batt;
        SOC_HDS(1) = mpcmilp.SOC_HDS;

        PVUsage(1) =  mpcmilp.PVUsage;
        WTUsage(1) =  mpcmilp.WTUsage;

        [~, mpcmilp, ~] = MPCMILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f.MPCMILP, nvars);
        % PlotingMILP_v5_saving(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc,folder, PredictionL);
        PV_co2(2) = mpcmilp.co2.PVE;
        WT_co2(2) = mpcmilp.co2.WTE;
        GR_co2(2) = mpcmilp.co2.GRE;

        Im(2) = mpcmilp.Im;
        Ex(2) = mpcmilp.Ex;

        SOC_Batt(2) = mpcmilp.SOC_Batt;
        SOC_HDS(2) = mpcmilp.SOC_HDS;

        PVUsage(2) =  mpcmilp.PVUsage;
        WTUsage(2) =  mpcmilp.WTUsage;
        % PlotingMILP_v5(x, Nh, T, Ndays, HeatingD(1:Nh), Ppv_actual(1:Nh), Pwind_actual(1:Nh), Cg_Ex(1:Nh), Cg_Im(1:Nh), CO2, init_soc);

        annotation.title = 'CO2 Emission';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'PVE , WTE , GRE [gco2perkWh]';
        annotation.label={"Shrinking", "Fixed"};
        Var(1,:) = [PV_co2(1) PV_co2(2)];
        Var(2,:) = [WT_co2(1) WT_co2(2)];
        Var(3,:) = [GR_co2(1) GR_co2(2)];
        Multi_PredictionHorizonPlotting(1:3, annotation, Var, {"PV ", " WT ", " GR"},folder+"/co2.svg");
        clear Var;
        annotation.title = 'Grid usage';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'Pg_{Im} , Pg_{Ex} [kW]';
        annotation.label={"Shrinking", "Fixed"};
        Var(1,:) = [Im(1) Im(2)];
        Var(2,:) = [Ex(1) Ex(2)];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"P_{Im} ", " P_{Ex}"},folder+"/Griduse.svg");

        clear Var;
        annotation.title = 'Storage SOC';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'SOC_{Batt} , SOC_{HS} [%]';
        annotation.label={"Shrinking", "Fixed"};
        Var(1,:) = [SOC_Batt(1) SOC_Batt(2)];
        Var(2,:) = [SOC_HDS(1) SOC_HDS(2) ];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"SOC_{Batt} ", " SOC_{HDS}"},folder+"/SOC.svg");

        clear Var;
        annotation.title = 'RE usage';
        annotation.xlabel = 'Time [Day]';
        annotation.ylabel = 'RE available / RE used [%]';
        annotation.label={"Shrinking", "Fixed"};
        Var(1,:) = [PVUsage(1) PVUsage(2)];
        Var(2,:) = [WTUsage(1) WTUsage(2)];
        Multi_PredictionHorizonPlotting(1:2, annotation, Var, {"PV ", " WT"},folder+"/usage.svg");


    otherwise
        disp('Invalid Simulation case.')
        return
end
end