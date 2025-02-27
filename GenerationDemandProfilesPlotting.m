
function GenerationDemandProfilesPlotting(H2HeatingDemand, PpvGeneration, WindGeneration, T, namefile)
hours = 1:T;
% Create a bar chart with 3 sets of data
figure; % Create a new figure window
bar(hours, [H2HeatingDemand'; PpvGeneration'; WindGeneration'], 'grouped'); % Grouped bar chart

% Add labels and title
xlabel('Hour of the Day');
ylabel('Heating-Hydrogen demand[kg] Power [kW]');
% title('24-Hour HHD, PV and Wind Data');

% Add a legend
legend({'HHD', 'PV', 'Wind'}, 'Location', 'NorthEast');

% Set x-axis labels for each hour
xticks(hours);
xticklabels(arrayfun(@(x) sprintf('%02d:00', x), hours, 'UniformOutput', false));

% Optionally, add gridlines for better readability
% grid on;
figure;
% line plot of generation and deamand profiles 
yyaxis left;
plot(hours, PpvGeneration', 'b-', hours, WindGeneration', 'g-', 'LineWidth',1.2)
ylabel('Available Solar PV and Wind power [kW]')
yyaxis right;
plot(hours, H2HeatingDemand', 'r--','LineWidth',1.2);
ylabel('Hydrogen-Heating Demand [kg]')
xlabel('Time [h]');
legend('PV','WT','HHD')
set(gca,'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, namefile);%,"ContentType","vector");


% % Example data
% x = 1:T;  % x-axis positions
% y1 = [10, 20, 30];  % Data for the first left y-axis
% y2 = [15, 25, 35];  % Data for the second left y-axis
% y3 = [100, 200, 300];  % Data for the right y-axis
% 
% % Create a figure
% figure;
% 
% % Hold the plot to overlay the bars
% hold on;
% 
% % Plot the first two datasets on the left y-axis
% bar(x - 0.1, PpvGeneration', 0.2, 'FaceColor', 'b');  % Dataset 1 (blue)
% bar(x + 0.1, WindGeneration', 0.2, 'FaceColor', 'g');  % Dataset 2 (green)
% 
% % Create the right y-axis using yyaxis
% yyaxis right
% bar(x, H2HeatingDemand, 0.4, 'FaceColor', 'r');  % Dataset 3 (red)
% 
% % Set the left and right y-axis labels
% yyaxis left
% ylabel('RE-generation');  % Label for the left y-axis
% % ylim([0 40]);  % Set limits for the left y-axis
% 
% yyaxis right
% ylabel('Heating-Hydrogen demand');  % Label for the right y-axis
% % ylim([0 350]);  % Set limits for the right y-axis
% 
% % Set axis labels and title
% xlabel('Hyderogen gas[kg]');
% % title('Grouped Bar Graph with Different Scales on Left and Right Y-Axes');
% 
% % Adjust legend
% legend('PV', 'WT', 'HHD', 'Location', 'best');
% 
% figure;
% 
% % Bar width (same for all datasets)
% barWidth = 0.25;
% 
% % Plot the first two datasets on the left y-axis
% hold on;  % Hold the plot to overlay multiple bars
% 
% % Group the first two datasets (y1 and y2) and plot them with the same width
% bar(x - barWidth, PpvGeneration', barWidth, 'FaceColor', 'b');  % Dataset 1 (blue)
% bar(x, WindGeneration', 'FaceColor', 'g');  % Dataset 2 (green)
% 
% % Create the right y-axis using yyaxis
% yyaxis right
% bar(x + barWidth, H2HeatingDemand', barWidth, 'FaceColor', 'r');  % Dataset 3 (red)
% 
% % Set the left and right y-axis labels
% yyaxis left
% ylabel('Left Y-Axis Scale');  % Label for the left y-axis
% % ylim([0 40]);  % Set limits for the left y-axis
% 
% yyaxis right
% ylabel('Right Y-Axis Scale');  % Label for the right y-axis
% % ylim([0 350]);  % Set limits for the right y-axis
% 
% % Set axis labels and title
% xlabel('X-axis Label');
% title('Grouped Bar Graph with Same Size Bars and Different Scales on Left and Right Y-Axes');
% 
% % Adjust legend
% legend('Dataset 1 (Left)', 'Dataset 2 (Left)', 'Dataset 3 (Right)', 'Location', 'best');
