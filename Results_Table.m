function [] = Results_Table(milp, mpcmilp)

    % Define the parameters for each set and scenario
    Parameter1 = [milp.co2.PVE, mpcmilp.co2.PVE];  % 2 scenarios, 3 parameter sets each
    Parameter2 = [milp.co2.WTE, mpcmilp.co2.WTE];
    Parameter3 = [milp.co2.GRE, mpcmilp.co2.GRE];
    
    % Scenario labels
    simType = {'MILP';  'MPCMILP'; };  % Only label the first row of each scenario
    
    % % Timestamps (repeating for each row)
    % Timestamp = repmat({datestr(now)}, 6, 1);  % 6 rows in total
    % 
    % % Notes (repeating for each scenario)
    % Notes = {'Notes for Scenario 1'; ''; ''; 'Notes for Scenario 2'; ''; ''};
    
    % Create the table
    results_table = table(simType, Parameter1(:), Parameter2(:), Parameter3(:),...
        'VariableNames', {'Simulation type', 'PV', 'WT', 'GR'});
    
    % Display the table
    disp(results_table);
    writetable(results_table,'Results.xlsx')
    file_path = which('Results.xlsx');
    % system(['start excel' "file_path"'])
    
    % Start Excel
    Excel = actxserver('Excel.Application');
    
    % Set Excel to be visible (optional)
    Excel.Visible = true;
    
    % Open a specific Excel file
    Workbook = Excel.Workbooks.Open(file_path);
end
    
   
