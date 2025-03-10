function [x,mpcmilp, elapsed_time, pv_ret, wt_ret] = MPCMILPSimulationShrinkingHorizon(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Cg_Im, Cg_Ex, Grid, Cost, Batt, dt, Capacities, Hs, CO2, f, nvars, Ndays, rand_pv_wt, stochastic)


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

% the length of horizon is taken as 24H and implemented only for 1 day
elapsed_time = 0.0;
pv_ret = zeros(Nh,1);
wt_ret = zeros(Nh,1);
for n =1:Ndays
    for d = 1:T/Np
        for k = 1:Np
            % Slicing demands and generation data to MILP engine 
            if stochastic == 1
                pv = Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np).*(1+rand_pv_wt(Np*(k-1)+(d-1)*Np*Np + (n-1)*Np*T +1:Np*k-k+1+(d-1)*Np*Np+ (n-1)*Np*T));
                pv_ret(k+(d-1)*Np+(n-1)*T,1) =  pv(1);
                wt = Pwind_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np).*(1+rand_pv_wt(Np*(k-1)+(d-1)*Np*Np + (n-1)*Np*T +1:Np*k-k+1+(d-1)*Np*Np+ (n-1)*Np*T));
                wt_ret(k+(d-1)*Np+(n-1)*T,1) = wt(1);
                hd =  HeatingD(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np );
            else
                % disp(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np)
                pv = Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np);
                pv_ret(k+(d-1)*Np+(n-1)*T,1) =  pv(1);
                wt = Pwind_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np );
                wt_ret(k+(d-1)*Np+(n-1)*T,1) = wt(1);
                hd =  HeatingD(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np );
            end
    
            tic;
            [x, ~, ~] = MILP_v4_Callee(pv, wt, hd, ...
            Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np - k+1, Capacities, Hs, CO2, Grid, Cost, "shrinking", Np - k+1);
            elapsed_time = toc+elapsed_time;

            Pg_Im(k+(d-1)*Np+(n-1)*T) = x(1);
            Pg_Ex(k+(d-1)*Np+(n-1)*T) = x(1*(Np-k+1)+1);
            Pb_c(k+(d-1)*Np+(n-1)*T) = x(2*(Np-k+1)+1);
            Pb_d(k+(d-1)*Np+(n-1)*T) = x(3*(Np-k+1)+1);
            SOC(k+(d-1)*Np+(n-1)*T) = x(4*(Np-k+1)+1);
            delta_g(k+(d-1)*Np+(n-1)*T) = x(5*(Np-k+1)+1);
            delta_b(k+(d-1)*Np+(n-1)*T) = x(6*(Np-k+1)+1);
            lambda_pv(k+(d-1)*Np+(n-1)*T) = x(7*(Np-k+1)+1);
            lambda_w(k+(d-1)*Np+(n-1)*T) = x(8*(Np-k+1)+1);
            SOC_Hs(k+(d-1)*Np+(n-1)*T) = x(9*(Np-k+1)+1);
            H_El2Hs(k+(d-1)*Np+(n-1)*T) = x(10*(Np-k+1)+1);
            H_El2Hd(k+(d-1)*Np+(n-1)*T) = x(11*(Np-k+1)+1);
            H_Hs2Hd(k+(d-1)*Np+(n-1)*T) = x(12*(Np-k+1)+1);
            Pel(k+(d-1)*Np+(n-1)*T) = x(13*(Np-k+1)+1);
            P_RE_B(k+(d-1)*Np+(n-1)*T) = x(14*(Np-k+1)+1);
            P_RE_g(k+(d-1)*Np+(n-1)*T) = x(15*(Np-k+1)+1);
            P_RE_Ele(k+(d-1)*Np+(n-1)*T) = x(16*(Np-k+1)+1);
            P_g_B(k+(d-1)*Np+(n-1)*T) = x(17*(Np-k+1)+1);
            P_g_Ele(k+(d-1)*Np+(n-1)*T) = x(18*(Np-k+1)+1);
            P_B_Ele(k+(d-1)*Np+(n-1)*T) = x(19*(Np-k+1)+1);
            P_B_g(k+(d-1)*Np+(n-1)*T) = x(20*(Np-k+1)+1);
            Batt.signals.values(7) = x(4*(Np-k+1)+1); % Update the initial SOC
            Hs.SOC_init = x(9*(Np-k+1)+1);
        end
    end
    
end

% disp('Grid cost with shrinking horizon')
Im = sum(Cg_Im(1:Nh).*Pg_Im);
% disp(Im);
% disp('Grid revenue with shrinking horizon')
Ex = sum(Cg_Ex(1:Nh).*Pg_Ex);
% disp(Ex);
% disp('MPC-MILP grid expenses with shrinking horizon')
% disp(sum(Cg_Im(1:Nh).*Pg_Im)-sum(Cg_Ex(1:Nh).*Pg_Ex))
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