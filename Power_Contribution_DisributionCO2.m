function Power_Contribution_DisributionCO2(span, annotation, Pagg) 

% Create a bar chart with var number of inputs
figure; % Create a new figure window
bar(span, Pagg, 'grouped'); % Grouped bar chart

% Add labels and title
xlabel(annotation.xlabel);
ylabel(annotation.ylabel);
% title(annotation.title);

% Add a legend
legend(annotation.label, 'Location', 'NorthEast')

% Set x-axis labels for each hour
% xticks(span);
% xticklabels(arrayfun(@(x) sprintf('%02d:00', x), span, 'UniformOutput', false));

% Optionally, add gridlines for better readability
grid on;
