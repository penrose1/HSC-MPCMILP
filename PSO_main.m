function [x,fval,exitflag,output,points] = PSO_main()
%%%%
        % 4:47 AM Thur 16 Jan 2025 
        % Perform PSO 

%%%%
func = @(X)CostFunction(X);
nvars = 5;

% Start timing
tic;

% X(1) Area of PV needed [ha] 
% X(2) Number of WT needed [#]
% X(3) Battery capacity [kWh] 
% X(4) Electrolyser capacity [kW] 
% X(5) Hydrogen storage capacity for heating 

lb = [0 3 200 100 500];
ub = [2 10 1000 3000 20000];
% M = 200;
% min = [1, 1, 50]; 
% max = [4000, 400, 10000];
initialSoln = [0 4 250 262 18000];
options = optimoptions('particleswarm','Display','iter','FunValCheck','on', 'InitialPoints',initialSoln);%, ...
    % 'OutputFcn',@PlotSwarmR3Euclidean);
[x,fval,exitflag,output,points] = particleswarm(func, nvars, lb, ub, options);

elapsedTime = toc;

% Display elapsed time in milliseconds
elapsedTime_ms = elapsedTime * 1000;  % Convert seconds to milliseconds
disp(['Elapsed time: ', num2str(elapsedTime_ms, '%.3f'), ' milliseconds']);
disp(['Elapsed time: ', num2str(elapsedTime, '%.3f'), ' seconds']);
disp(['Elapsed time: ', num2str(elapsedTime/60, '%.3f'), ' Minutes']);
Sound(1);
