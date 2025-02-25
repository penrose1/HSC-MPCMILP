function [Pg, boolVar] = MPCMILPSimulation_PSO(Nh, Np, Ppv_actual, Pwind_actual, HeatingD, Grid, Batt, dt, Capacities, Hs, CO2, f, nvars)


% Call MILP_v4_Callee.m to get optimized values for each h
% Storage location
Pg_Im = zeros(Nh,1);
Pg_Ex = zeros(Nh,1);

% Implementing receding horizon principle
% Demand, generation and controls are computed by using prediction horion data 
% First control action is implemented while others are rejected thereafter,
% the optimization problem is recomputed again

for k = 1:Nh
    fprintf("Start sim hour: %d\n", k);
    ff = [];
    for kf = 1:nvars
        ff = [ff, f((kf-1)*(Nh+Np-1)+k:Np+k-1+(kf-1)*(Nh+Np-1))];
    end 
    [x, boolVar, ~] = MILP_v4_Callee(Ppv_actual(k:Np+k-1), Pwind_actual(k:Np+k-1), HeatingD(k:Np+k-1), ...
        Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np, Capacities, Hs, CO2,ff);
    if boolVar == 1
        Pg_Im(k) = x(1);
        Pg_Ex(k) = x(Np+1);
        Batt.signals.values(7) = x(4*Np+1); % Update the initial SOC
        Hs.SOC_init = x(9*Np+1);
        fprintf("End sim hour: %d\n", k);
    else 
        fprintf("Sim hour: %d failed\n", k);
        Pg = 0;
        return; % On failure, return 0. PSO selects new sizes of the components 
    end
end
Pg = [Pg_Im; Pg_Ex];
end