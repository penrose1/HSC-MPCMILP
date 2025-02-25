clc, clear;
np =6;
t=24;
ndays = 2;
nh = ndays*t;
c1 =1:24;
c2=25:48;
c3=49:72;
c4=73:96;
c5=97:120;
c6=121:144;

% c1 =1:10;
% c2=11:20;
% c3=21:30;
% c4=31:40;
% c5=41:50;
% c6=51:60;

f=[c1 c2 c3 c4 c5 c6];
nvars=3;

for n =1:ndays
    for d=1:t/np
    for k = 1:np
        ff=[];
        for kf=1:nvars
           ff = [ff f(k+(n-1)*t+(d-1)*np+(kf-1)*nh: np +(n-1)*t+(d-1)*np+(kf-1)*nh)];
            if kf ==3
                % disp('Day: ')
                % disp(n)
                % disp("H: ")
                % disp(k)
                % disp('Var: ')
                disp(k)
                % disp('Starting index:')
                % disp((kf-1)*nh+np*(n-1)+k)
                disp(ff)
            end
        

        end
    end
    end
    
end

% [x, ~, ~] = MILP_v4_Callee(Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np), Pwind_actual(k+(n-1)*T:Np*n), HeatingD(k+(n-1)*T:Np*n), ...
%             Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np-k+1, Capacities, Hs, CO2,ff);
%             elapsed_time = toc+elapsed_time;