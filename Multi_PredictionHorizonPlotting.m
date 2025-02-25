function Multi_PredictionHorizonPlotting(span, annotation, Var, groupName, namefile) 

% Create a bar chart with var number of inputs
figure; % Create a new figure window
bar(span, Var, 'grouped'); % Grouped bar chart

% Add labels and title
xlabel(annotation.xlabel,'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(annotation.ylabel,'FontName', 'Times New Roman', 'FontSize', 12);
% title(annotation.title);

% Add a legend
legend(annotation.label, 'Location', 'best','FontName', 'Times New Roman', 'FontSize', 12);

% Set x-axis labels for each hour
% xticks(span);
xticklabels( groupName)
% Optionally, add gridlines for better readability
set(gca,'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, namefile);%,"ContentType","vector");
% grid on;
