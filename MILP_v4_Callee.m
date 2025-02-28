function [x, boolVar, fval] = MILP_v4_Callee(P_pv, P_wind, HeatingD, P_grid_max_Im, P_grid_max_Ex, P_grid_min_Ex, Batt, dt, Np, Capacities, Hs, CO2, Grid, Cost, scheme_type, Np_shrink)

% To use this solver for both shrinking and fixed length horizon, create an
% illusion on Np. Set Np manual to be Np_shrink for shrinking horizon.

if scheme_type == "shrinking"
    Np = Np_shrink;
end

% figure;
% plot(P_wind+P_pv,'LineWidth',1.2);
% hold on;
% plot(HeatingD,'LineWidth',1.2);
% legend('$P_{RE}$','$P_{Ele}$','Interpreter','latex');

%%%%%    
    % MPC-MILP - Modified to be a callee from MainProgram.m (1.04 AM Tue 14th Jan 2025)
    % 1106hrs Mon 23rd Dec 2024
    % MPC conconepts are added into the MILP program
%%%%%
% Parameters

%% MAX and MIN impossed by the system
% P_grid_max = 100; % Max grid power import/export (in kW)
Pb_max = Batt.signals.values(11);
Pb_min = Batt.signals.values(12);
SOC_max = Batt.signals.values(13); % Maximum battery SOC
SOC_min = Batt.signals.values(14); % Minimum battery SOC
Mb = Pb_max;
% Mg_max = P_grid_max;
Mg = P_grid_max_Im;


%% Electrolyzer parameters
eta_El = 0.6;
Pnom_El = Capacities.Pnom_El;
%% Hydrogen storage 
Qnom_Hs = Capacities.Hs;
SOC_Hs_min = Hs.SOC_min;
SOC_Hs_max = Hs.SOC_max;
SOC_Hs_init = Hs.SOC_init;

%% Battery additional parameters

eta_c = Batt.signals.values(8);
eta_d = Batt.signals.values(9);
SOC_init = Batt.signals.values(7);
Qnom = Batt.signals.values(10);

%% Decision Variables:

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
% Number of variables
n_continuous = 19 * Np;  % 11 continuous variables per time step (grid import, grid export, charge, discharge, SOC, lambda_pv, lambda_w)
n_binary = 2 * Np; % 2 binary variables for grid import/export and charge/discharge per time step
start_index_binary = 5*Np + 1;
stop_index_binary = 7*Np;
n_vars = n_continuous + n_binary;

%% Inequality Constraints
% Constraints matrix A and vector b
% Inequality Constraints equations for: Pb_charge, Pb_discharge, Pg_Im and
% Pg_Ex
% Decision variable vector, X: [Pg_Im; Pg_Ex; Pb_c; Pb_d; SOC; delta_g; delta_b; lambda_pv; lambda_w]

% 1. Pb_c(t) <= delta_b(t)*Mb
A1 = zeros(Np, n_vars);
A1(:,2*Np + 1:3*Np) = eye(Np); % Set Pb_c diagonal elements to one
A1(:, 6*Np + 1:7*Np) = -Mb*eye(Np); % Set delta_b 
b1 = zeros(Np,1);

% 2. Pd_d(t) <= (1 - delta_b(t))*Mb
A2 = zeros(Np, n_vars);
A2(:, 3*Np + 1:4*Np) = eye(Np); % Set Pb_d diagonal elements to one
A2(:, 6*Np + 1:7*Np) = Mb*eye(Np); % Set delta_b 
b2 = Mb*ones(Np,1);

% 3. Pg_Im(t) <= delta_g(t)*Mg
A3 = zeros(Np, n_vars);
A3(:, 1:Np) = eye(Np); % Set Pg_Im diagonal elements to one
A3(:, 5*Np + 1:6*Np) = -Mg*eye(Np); % Set delta_g 
b3 = zeros(Np,1);

% 4. Pg_Ex(t) <= (1 - delta_g(t))*Mg
A4 = zeros(Np, n_vars);
A4(:, Np + 1:2*Np) = eye(Np); % Set Pg_Ex diagonal elements to one
A4(:, 5*Np + 1:6*Np) = Mg*eye(Np); % Set delta_g 
b4 = Mg*ones(Np,1);

% 5. CO2 emission constraint 
% Ppv(k)*lambda_pv(k)*PVE + Pw(k)*lambda_w(k)*WTE + (Pg_Im(k) - Pg_Ex(k))*GRE <= (1-ER)*(HHD*NGE)

% Create matrix Aug_P_RE to implement P_pv(t)*lambda_pv(t) +
% P_wind(t)*lambda_w(t) --> also to be used later in otgher equations 
% For P_pv(t)
Aug_pv = zeros(Np);
Aug_pw = zeros(Np);
for i = 1:Np
    Aug_pv(i,i) =  P_pv(i);
end

% For P_wind(t)
for i = 1:Np
    Aug_pw(i,i) =  P_wind(i);
end

Aug_P_RE = [Aug_pv Aug_pw]; % To be used later on 
Aug_P_REE = [CO2.PVE*Aug_pv CO2.WTE*Aug_pw]; % Set Ppv(k)*lambda_pv(k)*PVE and Pw(k)*lambda_w(k)*WTE

A5 = zeros(Np, n_vars); % Create matrix Aeq1
A5(:, 1:Np) = CO2.GRE*eye(Np,Np); % Set Pg_Im(t)
A5(:, Np + 1:2*Np) = -CO2.ExportedREConsidered*CO2.GRE*eye(Np,Np); % Set -Pg_Ex(t) --> What happen if ignored?
A5(:, 7*Np +1:9*Np) = Aug_P_REE;
b5 = (1-CO2.ER)*(CO2.NGE*HeatingD );%+ 0*MObility);

% Combined As and bs to form A and b respectively to be used in the intlinprog() routine
A = [A1; A2; A3; A4; CO2.required*A5];
b = [b1; b2; b3; b4; CO2.required*b5];

%% Equality Constraints
% -> Energy balance at each time step
% 1. Pg_Im(t) - Pg_Ex(t) - Pb_c(t) + Pb_d(t) + Ppv(t)*lambda_pv(t) + Pwt(t)*lambda_w(t) - Pel(t) = 0 
Aeq1 = zeros(Np, n_vars); % Create matrix Aeq1
Aeq1(:, 1:Np) = eye(Np,Np); % Set Pg_Im(t)
Aeq1(:, Np + 1:2*Np) = -eye(Np,Np); % Set Pg_Ex(t)
Aeq1(:, 2*Np +1:3*Np) = -eye(Np,Np); % Set Pb_c(t)
Aeq1(:, 3*Np +1:4*Np) = eye(Np,Np); % Set Pb_d(t)
Aeq1(:, 13*Np +1:14*Np) = -eye(Np,Np); % Set Pel(k)
Aeq1(:, 7*Np +1:9*Np) = Aug_P_RE; % Set {lambda_pv(t) & lambda_w(t)}
beq1 = zeros(Np,1); % RHS 

% 2. 0*Pg_Im(t) + 0*Pg_Ex(t) + 100*Pb_c(t)*dt*eta_c/Qnom - 100*Pb_d(t)*dt/(eta_d*Qnom) + SOC(t-1) - SOC(t) + 0*delta_g + 0*delta_b = 0
Aeq2 = zeros(Np, n_vars);
Aeq2(:, 2*Np + 1:3*Np) = (100/Qnom)*dt*eta_c*eye(Np); % Set Pb_c(t)
Aeq2(:, 3*Np + 1:4*Np) = -(100/Qnom)*(dt/eta_d)*eye(Np); % Set Pb_d(t)
Aeq2(:, 4*Np + 1:5*Np) = -eye(Np); % Set SOC(k) to -1

for k = 2:Np
    Aeq2(k, 4*Np + k -1) = 1; % Set elements of SOC(k) to +1. Diagonal elements of a (T-1)x(T-1) matrix
end

beq2 = zeros(Np, 1); % Set all elements of RHS to zero
beq2(1,1) = -SOC_init; % Set the initia SOC to compliment. --> *SOC*

% 3. SOC_Hs evolution (SOC in % )
% SOC_Hs(k-1) + SOC_Hs(k) + 100*(H_El2Hs(k) - H_Hs2Hd(k))/Qnom_Hs = 0
Aeq3 = zeros(Np, n_vars);
Aeq3(:, 10*Np + 1:11*Np) = (100/Qnom_Hs)*dt*eye(Np); % Set H_El2Hs(k)
Aeq3(:, 12*Np + 1:13*Np) = -(100/Qnom_Hs)*dt*eye(Np); % Set H_Hs2Hd(k)
Aeq3(:, 9*Np +1:10*Np) = -eye(Np);

for k =2:Np
    Aeq3(k, 9*Np + k -1) = 1; % Set elements of SOC_Hs(k) to +1. Diagonal elements of a (T-1)x(T-1) matrix
end

beq3 = zeros(Np, 1); % Set all elements of RHS to zero
beq3(1,1) = -SOC_Hs_init; % Set the initia SOC to compliment. --> *SOC*

% 4. Electrolyzer power balance 
% Eel(k)*eta_Ele/kWhkg = H_El2Hs(k) + H_El2Hd(k)
Aeq4 = zeros(Np, n_vars);
Aeq4(:, 10*Np +1:11*Np) = eye(Np); % Set H_El2Hs(k)
Aeq4(:, 11*Np + 1:12*Np) = eye(Np); % Set H_El2Hd(k)
Aeq4(:, 13*Np + 1:14*Np) = -(eta_El)*eye(Np);
beq4 = zeros(Np, 1);

% 5. H2 Heating demand balance 
% H_El2Hd(k) + H_Hs2Hd(k) = HD(k)
Aeq5 = zeros(Np, n_vars);
Aeq5(:, 11*Np + 1:12*Np) = eye(Np); % Set H_El2Hd(k)
Aeq5(:, 12*Np + 1:13*Np) = eye(Np); % Set H_Hs2Hd(k)
beq5 = HeatingD;

% 6 P_B_Ele(k) + P_Re_Ele(k) + P_g_Ele(k) = P_Ele_in(k)
Aeq6 = zeros(Np, n_vars);
Aeq6(:,19*Np + 1:20*Np) = eye(Np); 
Aeq6(:,16*Np + 1:17*Np) = eye(Np); 
Aeq6(:,18*Np + 1:19*Np) = eye(Np); 
Aeq6(:,13*Np + 1:14*Np) = -eye(Np); 
beq6 = zeros(Np,1);

% 7. P_RE(k) = P_RE_B(k) + P_RE_g(k) + P_RE_Ele(k) 
Aeq7 = zeros(Np, n_vars);
Aeq7(:,14*Np + 1:15*Np) = eye(Np); 
Aeq7(:,15*Np + 1:16*Np) = eye(Np); 
Aeq7(:,16*Np + 1:17*Np) = eye(Np);  
Aeq7(:, 7*Np +1:9*Np) = -Aug_P_RE; % Set {lambda_pv(t) & lambda_w(t)}
beq7 = zeros(Np,1);

% 8. Pb_c(k) = P_RE_B(k) + P_g_B(k)
Aeq8 = zeros(Np, n_vars);
Aeq8(:,14*Np + 1:15*Np) = eye(Np); 
Aeq8(:,17*Np + 1:18*Np) = eye(Np);
Aeq8(:,2*Np + 1:3*Np) = -eye(Np); 
beq8 = zeros(Np,1);

% 9. Pb_d(k) = P_B_g(k) + P_B_Ele(k)
Aeq9 = zeros(Np, n_vars);
Aeq9(:,20*Np + 1:21*Np) = eye(Np); 
Aeq9(:,19*Np + 1:20*Np) = eye(Np);
Aeq9(:,3*Np + 1:4*Np) = -eye(Np); 
beq9 = zeros(Np,1);

% 10. Pg_Im(k) = P_g_Ele(k) + P_g_B(k) 
Aeq10 = zeros(Np, n_vars);
Aeq10(:,18*Np + 1:19*Np) = eye(Np); 
Aeq10(:,17*Np + 1:18*Np) = eye(Np);
Aeq10(:,1:Np) = -eye(Np); 
beq10 = zeros(Np,1);

% 11. Pg_Ex(k) = P_B_g(k) P_RE_g(k) 
Aeq11 = zeros(Np, n_vars);
Aeq11(:,20*Np + 1:21*Np) = eye(Np); 
Aeq11(:,15*Np + 1:16*Np) = eye(Np);
Aeq11(:,Np + 1:2*Np) = -eye(Np); 
beq11 = zeros(Np,1);

% Single matrix and vector for Aeq and beq to be used in intlinprog()
Aeq = [Aeq1; Aeq2; Aeq3; Aeq4; Aeq5; Aeq6; Aeq7; Aeq8; Aeq9; Aeq10; Aeq11;];
beq = [beq1; beq2; beq3; beq4; beq5; beq6; beq7; beq8; beq9; beq10; beq11];

%% Upper and Lower bounds
% Lower and upper bounds for continuous variables (all >= 0)

lb = zeros(1, n_vars); % Set all decision variables at once 
lb(2*Np + 1:4*Np) = repmat(Pb_min, 1, 2*Np); % Set Pb_min (same for dis/charge)--> 3rd and 4th block
lb(4*Np + 1:5*Np) = repmat(SOC_min, 1, Np); % Set SOC_min --> 5th block entry 
lb(9*Np + 1:10*Np) = repmat(SOC_Hs_min, 1, Np); % Set SOC_Hs_min --> 10th block entry 
lb(Np + 1:2*Np) = repmat(P_grid_min_Ex, 1, Np); % Set Pg_Ex_min --> 2rd bOClock

ub = [repmat(P_grid_max_Im, 1, Np), repmat(P_grid_max_Ex, 1, Np), repmat(Pb_max, 1, Np), ...
    repmat(Pb_max, 1, Np), repmat(SOC_max, 1, Np), ones(1, n_binary), ones(1,2*Np), ...
    repmat(SOC_Hs_max,1,Np), Inf*ones(1,3*Np), repmat(Pnom_El, 1, Np), Inf*ones(1,7*Np)];

%% Cost function
% % Cost function coefficients
% Cg = [Cg_Im', -Cg_Ex'];
% Cb = [Cg_Ex', Cg_Im']; % Favour battery over grid
% Csoc = zeros(1,Np);
% Csoc(1) = gamma;
% Csoc(Np) = -gamma;
% Cre = [-1*Cg_Im' -1*Cg_Im'];
% Cpel = Cg_Im';
% f = [0*Cg, 0*Cb, -Csoc, zeros(1, n_binary), Cre, zeros(1,4*Np), Cpel, zeros(1,7*Np)];

Cg_Im = Grid.Cg_Im(1:Np)';
Cg_Ex = Grid.Cg_Ex(1:Np)';

Cg_ImF = Cost.Cg_ImF(1:Np)';
Cg_ExF = Cost.Cg_ExF(1:Np)';

% As the prediction window recedes, more data,(Np-1) are required. --> This
% will be different with shrinking horizon approach
Cg = [Cg_Im, -1*Cg_Ex];
Cb = [-Cg_ExF, -1*Cg_ImF]; % Favour battery over grid
Csoc = (1/14)*ones(1,Np);
Cre = [Cg_ImF Cg_ImF];
Cpel = Cg_ImF;

Cre_b = (1/3)*Cg_ImF;
Cre_g = (1/3)*Cg_ExF;
Cre_Ele = (1)*Cg_ImF;
Cg_b = 1*Cg_ImF;
Cg_Ele = 1*Cg_ImF;
Cb_Ele = (1/1)*Cg_ImF;
Cb_g = (0)*Cg_ImF;

f = [0*Cg, 0*Cb, -Csoc, zeros(1, n_binary), -1*Cre, -Csoc, zeros(1, 3*Np), 0*Cpel, -Cre_b, -Cre_g, -Cre_Ele, Cg_b, Cg_Ele, -Cb_Ele, Cb_g];
% [~ ,f_size] = size(f);
% fprintf('Size of f: %d\n', f_size);

%% Configuration and invocation of  intlinprog routine
% Define the optimization problem using 'intlinprog' for MILP

intcon = start_index_binary:stop_index_binary; % Binary decision variables: grid import/export, battery charge/discharge
options = optimoptions('intlinprog', 'Display', 'off', 'AbsoluteGapTolerance', 0.01); % Options
% Solving the MILP with try and catch 
try 
    [x, fval, exitflag, ] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);
    if exitflag == 0 || exitflag == -2 || exitflag == -3 || exitflag == -9
        boolVar = -1;
        disp('failed');
        x = 0; % For infeasible region. i.e. sizes specified by PSO should be discard (Inf cost)
    else 
        boolVar = 1;
    end
    % disp('Founded')
catch
    x = 0; % For infeasible region.{When sizing, the sizes specified by PSO should be discard (inf cost)}
    boolVar = -1;
    disp('catched exception');
end
