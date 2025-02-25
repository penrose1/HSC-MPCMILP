clc, clear;
function stop = myOutputFcn(optimValues, state)
    % persistent solutions fvals  % Declare persistent variables to store solutions and fvals
    
    stop = false;  % Continue optimization by default
    file = fopen('data.txt','a');
    if file == -1
        error('Faile to open');
    end
    X = [optimValues.iteration optimValues.bestfval optimValues.bestx]
    formatSpec = "%.0f  %.4f  %.4f  %.4f\r\n";

fprintf(file,formatSpec, X);
fclose(file);
    % disp(['Iteration:',num2str(optimValues.iteration), ', Best values:', ...
    %     num2str(optimValues.bestfval), ', Best soln:', num2str(optimValues.bestx)]);
end

% Define optimization problem
objectiveFunction = @(x) sum(x.^2);  % Example objective function
numVariables = 2;  % Number of decision variables
lb = [-5, -5];  % Lower bounds
ub = [5, 5];    % Upper bounds

% Create options with custom output function
options = optimoptions('particleswarm', 'OutputFcn', @myOutputFcn);

% Run particleswarm with the custom output function
[bestSolution, bestValue] = particleswarm(objectiveFunction, numVariables, lb, ub, options);
