% figure;
% plot(1*P_electrolyzer,'--r','LineWidth',2)
% hold on 
% wt = 8;
% ap = 4;
% plot(P_wind_prediction + P_pv_actual)
% legend('Pele','Pwind+Ppv')

%% 
% Your program code here

% Sound after program finishes
Fs = 44100;  % Sampling frequency (Hz)
t = 0:1/Fs:10;  % Time vector (1 second long)
y = sin(2 * pi * 440 * t);  % 440 Hz tone (A4 note)
% plot(y)
% Play sound
% sound(y, Fs);

%%
% % Example data
% x = 0:0.1:10;               % X-axis values
% y1 = sin(x);                % First Y-axis data
% y2 = cos(x) * 100;          % Second Y-axis data with a different scale
% 
% % Create a plot with two y-axes
% figure;
% 
% % Activate the left y-axis (default)
% yyaxis left
% plot(x, y1, '-b', 'LineWidth', 2);
% ylabel('sin(x)');           % Label for the left y-axis
% axis tight;
% 
% % Activate the right y-axis
% yyaxis right
% plot(x, y2, '-r', 'LineWidth', 2);
% ylabel('cos(x) * 100');     % Label for the right y-axis
% 
% % Add labels and title
% xlabel('X-axis');
% title('Plot with Different Y-Axis Scales');
% 
% % Display the grid
% grid on;

%% Hitogram
% Sample data
data = randn(1000, 1); 

% Create a customized histogram
figure; 
Pele_min = min(P_electrolyzer(1:T));
Pele_max = max(P_electrolyzer(1:T));
span = Pele_min:1:Pele_max;
histogram(P_electrolyzer(1:T), 'BinEdges', span, 'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.5);
title('Customized Histogram');
xlabel('Value');
ylabel('Frequency');

%% Bar chart 
% Example data
products = {'Product A', 'Product B', 'Product C', 'Product D'};
demand = [200, 150, 300, 250];

% Create a bar chart
figure; % Create a new figure window
bar(demand); % Create bar chart
set(gca, 'XTickLabel', products); % Label the x-axis with product names

% Add labels and title
xlabel('Products');
ylabel('Demand');
title('Demand for Different Products');

% Optionally, add gridlines for better readability
grid on;

%% 
% Example data (24-hour data for Load and PV)
hours = 0:23; % 24 hours (0 to 23)
load_data = [150, 160, 170, 180, 190, 200, 220, 240, 250, 260, 270, 290, ...
             300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410]; % Load data (kW)
pv_data = [0, 0, 0, 0, 0, 0, 50, 100, 150, 200, 250, 300, ...
           320, 330, 320, 310, 300, 290, 250, 200, 150, 100, 50, 0]; % PV data (kW)

% Create a bar chart with 2 sets of data
figure; % Create a new figure window
bar(hours, [load_data; pv_data]', 'grouped'); % Grouped bar chart

% Add labels and title
xlabel('Hour of the Day');
ylabel('Power (kW)');
title('24-Hour Load and PV Data');

% Add a legend
legend({'Load', 'PV'}, 'Location', 'NorthEast');

% Set x-axis labels for each hour
xticks(hours);
xticklabels(arrayfun(@(x) sprintf('%02d:00', x), hours, 'UniformOutput', false));

% Optionally, add gridlines for better readability
grid on;

%% 
% Create the bar plot
figure;
y1 = [1 0 3 0 5];
y2 = [-0 -2 -0 -4 -0];

x = [1 2 3 4 5];
h1 = bar(x, y1);
hold on;
h2 = bar(x, y2);
hold off;
% Set color for positive and negative bars
h1.FaceColor = 'flat';
h2.FaceColor = 'flat';
h1.CData = repmat([0, 1, 0], 5, 1);  % Green for positive bars
h2.CData = repmat([1, 0, 0], 5, 1);  % Red for negative bars

% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
title('Bar Plot with Positive and Negative Values');
