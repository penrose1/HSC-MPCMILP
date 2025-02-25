% Number of sets per scenario (3 sets per scenario)
num_sets = 3;

% Define the parameters for each set and scenario
Parameter1 = [1.1, 1.2, 1.3; 2.1, 2.2, 2.3];  % 2 scenarios, 3 parameter sets each
Parameter2 = [1.4, 1.5, 1.6; 2.4, 2.5, 2.6];
Parameter3 = [1.7, 1.8, 1.9; 2.7, 2.8, 2.9];

% Define the results for each set and scenario
Result1 = [123.45, 124.56, 125.67; 98.76, 99.87, 100.98];
Result2 = [67.89, 68.90, 69.01; 45.32, 46.43, 47.54];
Result3 = [12.34, 13.45, 14.56; 56.78, 57.89, 58.90];

% Scenario labels
Scenario = {'Scenario 1'; ''; ''; 'Scenario 2'; ''; ''};  % Only label the first row of each scenario

% Timestamps (repeating for each row)
Timestamp = repmat({datestr(now)}, 6, 1);  % 6 rows in total

% Notes (repeating for each scenario)
Notes = {'Notes for Scenario 1'; ''; ''; 'Notes for Scenario 2'; ''; ''};

% Create the table
simulationResults = table(Scenario, Parameter1(:), Parameter2(:), Parameter3(:), ...
    Result1(:), Result2(:), Result3(:), Timestamp, Notes, ...
    'VariableNames', {'Scenario', 'Parameter1', 'Parameter2', 'Parameter3', ...
                      'Result1', 'Result2', 'Result3', 'Timestamp', 'Notes'});

% Display the table
disp(simulationResults);

writetable(simulationResults,'Sim.csv')

file_path = which('Sim.csv');
% Start Excel
Excel = actxserver('Excel.Application');

% Set Excel to be visible (optional)
Excel.Visible = true;

% Open a specific Excel file
Workbook = Excel.Workbooks.Open(file_path);


%%
iterations = 1:5;  % Example iteration numbers
time_taken = [0.1, 0.2, 0.15, 0.12, 0.18];  % Example time taken for each iteration
accuracy = [0.98, 0.96, 0.97, 0.95, 0.99];  % Example accuracy for each iteration

% Create the table
results_table = table(iterations', time_taken', accuracy', 'VariableNames', {'Iteration', 'TimeTaken', 'Accuracy'});

% Display the table
disp(results_table);
% New data for iteration 6
new_iteration = 6;
new_time_taken = 0.14;
new_accuracy = 0.96;

% Create a new row as a table
new_row = table(new_iteration, new_time_taken, new_accuracy, 'VariableNames', {'Iteration', 'TimeTaken', 'Accuracy'});

% Append the new row to the original table
results_table = [results_table; new_row];

disp(results_table);

 % % Specify the URL you want to open in Edge
    % url = 'https://www.example.com';
    % 
    % % Open the URL in Microsoft Edge
    % status = system(['start microsoft-edge:', file_path]);
    % 
    % % Check if the command was successful
