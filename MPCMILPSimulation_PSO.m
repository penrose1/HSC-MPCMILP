
% Implementing receding horizon principle
% Demand, generation and controls are computed by using prediction horion data 
% First control action is implemented while others are rejected thereafter,
% the optimization problem is recomputed again


function [gridUse, boolVar] = MPCMILPSimulation_PSO(Nh, Np, T, Ppv_actual, Pwind_actual, HeatingD, Grid, Cost, Batt, dt, Capacities, Hs, CO2, Ndays)
Cg_Im = Grid.Cg_Im;
Cg_Ex = Grid.Cg_Ex;
% Call MILP_v4_Callee.m to get optimized values for each h
% Storage location
Pg_Im = zeros(Nh,1);
Pg_Ex = zeros(Nh,1);

% Implementing receding horizon principle
% Demand, generation and controls are computed by using prediction horion data 
% First control action is implemented while others are rejected thereafter,
% the optimization problem is recomputed again

for n = 1:Ndays % Day d
    fprintf("Day : %d START\n", n);
    for d = 1:T/Np % Which part of the day n are you in?
        for k = 1:Np % Hour k in day d
            % if stochastic == 1
            %     pv = Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1).*(1+rand_pv_wt(Np*(k-1)+1:Np*k));
            % 
            %     wt = Pwind_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1).*(1+rand_pv_wt(Np*(k-1)+1:Np*k));
            % 
            %     hd =  HeatingD(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
            % else
                % pv = Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
                % 
                % wt = Pwind_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
                % 
                % hd =  HeatingD(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
            % end

            pv = Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
                
            wt = Pwind_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);
            
            hd =  HeatingD(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np + k - 1);

            [x, boolVar, ~] = MILP_v4_Callee(pv, wt, hd, ...
                Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np, Capacities, Hs, CO2, Grid, Cost, "fixed", Np);
            
            if boolVar == 1
                Pg_Im(k+(d-1)*Np+(n-1)*T) = x(1);
                Pg_Ex(k+(d-1)*Np+(n-1)*T) = x(Np+1);
                Batt.signals.values(7) = x(4*(Np)+1); % Update the initial SOC
                Hs.SOC_init = x(9*(Np)+1);
            else
                Pg_Im(k+(d-1)*Np+(n-1)*T) = 1e25;
                Pg_Ex(k+(d-1)*Np+(n-1)*T) = 0;
                break;
            end
            if boolVar == -1
                break;
            end
        end
        if boolVar == -1
            break;
        end
    end
    fprintf("Day : %d END   \n", n);
    if boolVar == -1
        break;
    end
end
if boolVar == -1
    gridUse = 1e25;
else
    Im = sum(Cg_Im(1:Nh).*Pg_Im);
    Ex = sum(Cg_Ex(1:Nh).*Pg_Ex);
    gridUse = Im-Ex;
end
end

