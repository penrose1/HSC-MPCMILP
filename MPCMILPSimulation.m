function [x,mpcmilp, elapsedTime] = MPCMILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f, nvars)

% Call MILP_v4_Callee.m to get optimized values for each h
% Storage location
Pg_Im = zeros(Nh,1);
Pg_Ex = zeros(Nh,1);
Pb_c = zeros(Nh,1);
Pb_d = zeros(Nh,1);
SOC = zeros(Nh,1);
delta_g = zeros(Nh,1);
delta_b = zeros(Nh,1);
lambda_pv = zeros(Nh,1);
lambda_w = zeros(Nh,1);
SOC_Hs = zeros(Nh,1);
H_El2Hs = zeros(Nh,1);
H_El2Hd = zeros(Nh,1);
H_Hs2Hd = zeros(Nh,1);
Pel = zeros(Nh,1);
P_RE_B = zeros(Nh,1);
P_RE_g =zeros(Nh,1);
P_RE_Ele = zeros(Nh,1);
P_g_B = zeros(Nh,1);
P_g_Ele = zeros(Nh,1);
P_B_Ele = zeros(Nh,1);
P_B_g = zeros(Nh,1);

% Implementing receding horizon principle
% Demand, generation and controls are computed by using prediction horion data 
% First control action is implemented while others are rejected thereafter,
% the optimization problem is recomputed again
elapsedTime = 0.0;
for k = 1:Nh
    ff = [];
    for kf = 1:nvars
        ff = [ff, f((kf-1)*(Nh+Np-1)+k:Np+k-1+(kf-1)*(Nh+Np-1))];
    end 
    tic;
    [x, ~, ~] = MILP_v4_Callee(Ppv_actual(k:Np+k-1), Pwind_actual(k:Np+k-1), HeatingD(k:Np+k-1), ...
        Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np, Capacities, Hs, CO2,ff);
    elapsedTime = elapsedTime+toc;
    Pg_Im(k) = x(1);
    Pg_Ex(k) = x(Np+1);
    Pb_c(k) = x(2*Np+1);
    Pb_d(k) = x(3*Np+1);
    SOC(k) = x(4*Np+1);
    delta_g(k) = x(5*Np+1);
    delta_b(k) = x(6*Np+1);
    lambda_pv(k) = x(7*Np+1);
    lambda_w(k) = x(8*Np+1);
    SOC_Hs(k) = x(9*Np+1);
    H_El2Hs(k) = x(10*Np+1);
    H_El2Hd(k) = x(11*Np+1);
    H_Hs2Hd(k) = x(12*Np+1);
    Pel(k) = x(13*Np+1);
    P_RE_B(k) = x(14*Np+1);
    P_RE_g(k) = x(15*Np+1);
    P_RE_Ele(k) = x(16*Np+1);
    P_g_B(k) = x(17*Np+1);
    P_g_Ele(k) = x(18*Np+1);
    P_B_Ele(k) = x(19*Np+1);
    P_B_g(k) = x(20*Np+1);
    Batt.signals.values(7) = x(4*Np+1); % Update the initial SOC
    Hs.SOC_init = x(9*Np+1);
end
disp('Grid cost')
Im = sum(Cg_Im(1:Nh).*Pg_Im);
disp(Im);
disp('Grid revenue')
Ex = sum(Cg_Ex(1:Nh).*Pg_Ex);
disp(Ex);
disp('MPC-MILP grid expenses')
disp(sum(Cg_Im(1:Nh).*Pg_Im)-sum(Cg_Ex(1:Nh).*Pg_Ex))
x = [Pg_Im; Pg_Ex; Pb_c; Pb_d; SOC; lambda_pv; lambda_w; SOC_Hs; H_El2Hs; H_El2Hd; H_Hs2Hd; Pel; P_RE_B ; P_RE_g; P_RE_Ele; P_g_B; P_g_Ele; P_B_Ele; P_B_g ];
mpcmilp.co2.PVE = sum(Ppv_actual(1:Nh).*lambda_pv(1:Nh))*CO2.PVE;
mpcmilp.co2.WTE = sum(Pwind_actual(1:Nh).*lambda_w(1:Nh))*CO2.WTE;
mpcmilp.co2.GRE = sum((Pg_Im(1:Nh) - Pg_Ex(1:Nh)))*CO2.GRE;
mpcmilp.Im = Im;
mpcmilp.Ex = Ex;
mpcmilp.SOC_Batt = Batt.signals.values(7);
mpcmilp.SOC_HDS = Hs.SOC_init;
mpcmilp.PVUsage = 100*sum(lambda_pv)/Nh;
mpcmilp.WTUsage = 100*sum(lambda_w)/Nh;