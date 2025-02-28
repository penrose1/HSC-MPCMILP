% clc, clear;
np =10;
t=10;
ndays = 2;
nh = ndays*t;
% c1 =1:24;
% c2=25:48;
% c3=49:72;
% c4=73:96;
% c5=97:120;
% c6=121:144;

% c1 =1:10+9;
% c2=11:20+9;
% c3=21:30+9;
% c4=31:40+9;
% c5=41:50+9;
% c6=51:60+9;
c1=ones(1,19);
c2=2*ones(1,19);
c3=3*ones(1,19);
c4=4*ones(1,19);
c5=5*ones(1,19);
c6=6*ones(1,19);
pv = 1:nh+np;
rand_pv_wt = linspace(1,np*nh,np*nh);
f=[c1 c2 c3 c4 c5 c6];
nvars=3;
recede_type = 'shrinking'
switch recede_type
    case 'fixed'
        for n =1:ndays
            for d=1:t/np
                for k = 1:np
                    ff=[];
                    for kf=1:nvars
                       ff = [ff f(k+(n-1)*t+(d-1)*np+(kf-1)*nh: np +(n-1)*t+(d-1)*np+(kf-1)*nh + k - 1)];
                        % if kf ==3
                        %     fprintf('Step #: %d\n', k);
                        %     % disp(ff)
                        % end
                    end
                    disp(pv(k+(n-1)*t+(d-1)*np:np +(n-1)*t+(d-1)*np + k - 1));
                    disp(rand_pv_wt(np*(k-1)+(d-1)*np*np + (n-1)*np*t +1:np*k+(d-1)*np*np+ (n-1)*np*t));
                end
            end
            
        end

    case 'shrinking'
        for n =1:ndays
            for d=1:t/np
                for k = 1:np
                    ff=[];
                    for kf=1:nvars
                       % ff = [ff f(k+(n-1)*t+(d-1)*np+(kf-1)*nh: np +(n-1)*t+(d-1)*np+(kf-1)*nh)];
                        if kf ==3
                            % fprintf('Step #: %d\n', k);
                            % disp(ff)
                        end
                    end
                    disp(pv(k+(n-1)*t+(d-1)*np:np +(n-1)*t+(d-1)*np))
                    disp(rand_pv_wt(np*(k-1)+(d-1)*np*np + (n-1)*np*t +1:np*k-k+1+(d-1)*np*np+ (n-1)*np*t));
                end
            end
            
        end
end




% [x, ~, ~] = MILP_v4_Callee(Ppv_actual(k+(n-1)*T+(d-1)*Np:Np +(n-1)*T+(d-1)*Np), Pwind_actual(k+(n-1)*T:Np*n), HeatingD(k+(n-1)*T:Np*n), ...
%             Grid.Pgrid_max_Im, Grid.Pgrid_max_Ex, Grid.Pgrid_min_Ex, Batt, dt, Np-k+1, Capacities, Hs, CO2,ff);
%             elapsed_time = toc+elapsed_time;