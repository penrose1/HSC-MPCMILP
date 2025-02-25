function [x, milp] = MILPSimulation(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Batt, dt, Capacities, Hs, CO2, f)

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
for k = 1:(Nh/Np)
    ff = [];
    for kf = 1:Nh:length(f)
        ff = [ff, f(Np*(k-1)+kf:Np*k -1 +kf)];
    end 
    [x, ~, ~] = MILP_v4_Callee(Ppv_actual(1 + (k-1)*T:k*T), Pwind_actual(1 + (k-1)*T:k*T), ...
        HeatingD(1 + (k-1)*T:k*T), ...
        Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, T, Capacities, Hs, CO2,ff);
    Pg_Im(1+(k-1)*T:T*k) = x(1:T);
    Pg_Ex(1+(k-1)*T:T*k) = x(T+1:2*T);
    Pb_c(1+(k-1)*T:T*k) = x(2*T+1:3*T);
    Pb_d(1+(k-1)*T:T*k) = x(3*T+1:4*T);
    SOC(1+(k-1)*T:T*k) = x(4*T+1:5*T);
    delta_g(1+(k-1)*T:T*k) = x(5*T+1:6*T);
    delta_b(1+(k-1)*T:T*k) = x(6*T+1:7*T);
    lambda_pv(1+(k-1)*T:T*k) = x(7*T+1:8*T);
    lambda_w(1+(k-1)*T:T*k) = x(8*T+1:9*T);
    SOC_Hs(1+(k-1)*T:T*k) = x(9*T+1:10*T);
    H_El2Hs(1+(k-1)*T:T*k) = x(10*T+1:11*T);
    H_El2Hd(1+(k-1)*T:T*k) = x(11*T+1:12*T);
    H_Hs2Hd(1+(k-1)*T:T*k) = x(12*T+1:13*T);
    Pel(1+(k-1)*T:T*k) = x(13*T+1:14*T);
    P_RE_B(1+(k-1)*T:T*k) = x(14*T+1:15*T);
    P_RE_g(1+(k-1)*T:T*k) = x(15*T+1:16*T);
    P_RE_Ele(1+(k-1)*T:T*k) = x(16*T+1:17*T);
    P_g_B(1+(k-1)*T:T*k) = x(17*T+1:18*T);
    P_g_Ele(1+(k-1)*T:T*k) = x(18*T+1:19*T);
    P_B_Ele(1+(k-1)*T:T*k) = x(19*T+1:20*T);
    P_B_g(1+(k-1)*T:T*k) = x(20*T+1:21*T);
    Batt.signals.values(7) = x(5*T); % Update the initial SOC
    Hs.SOC_init = x(10*T);
end
disp('Grid cost')
disp(sum(Cg_Im.*Pg_Im));
disp('Grid revenue')
disp(sum(Cg_Ex.*Pg_Ex));
disp('MILP grid expenses')
disp(sum(Cg_Im.*Pg_Im)-sum(Cg_Ex.*Pg_Ex))
x = [Pg_Im; Pg_Ex; Pb_c; Pb_d; SOC; lambda_pv; lambda_w; SOC_Hs; H_El2Hs; H_El2Hd; H_Hs2Hd; Pel; P_RE_B ; P_RE_g; P_RE_Ele; P_g_B; P_g_Ele; P_B_Ele; P_B_g ];
milp.co2.PVE = sum(Ppv_actual(1:Nh).*lambda_pv)*CO2.PVE;
milp.co2.WTE = sum(Pwind_actual(1:Nh).*lambda_w)*CO2.WTE;
milp.co2.GRE = sum((Pg_Im(1:Nh) - Pg_Ex(1:Nh)))*CO2.GRE;
end