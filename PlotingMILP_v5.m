function [] = PlotingMILP_v5(x, Nh, T, Ndays, HDemand, Ppv_available, Pwt_available, Cg_Im, Cg_Ex, CO2, init_soc)
% To be used for plotting results from MILP optimizer and from the
% simulation when MILP's controls are applied (Uncertainty or determistic case) 
% Adapted from PlotingMILP_v4(x, T, P_load, P_pv, C_grid, C_feed_in)
% 2:12 PM 14th Jan 2025 
% ENTRANCE 

%% Extract the results
% Decision variable vector, X: [Pg_Im; Pg_E; Pb_c; Pb_d; SOC; delta_g; delta_b]
Pg_Im = [x(1:Nh)];
Pg_Ex = [x(Nh+1:2*Nh)];
Pb_c = [x(2*Nh+1:3*Nh)];
Pb_d = [x(3*Nh+1:4*Nh)];
SOC = [x(4*Nh+1:5*Nh)];
SOC = [init_soc.Battery; SOC];
lambda_pv = [x(5*Nh+1:6*Nh)];
lambda_w = [x(6*Nh+1:7*Nh)];
SOC_Hs = [init_soc.Hs; x(7*Nh+1:8*Nh)];
H_El2Hs = [x(8*Nh+1:9*Nh)];
H_El2Hd = [x(9*Nh+1:10*Nh)];
H_Hs2Hd = [x(10*Nh+1:11*Nh)];
Pel = [x(11*Nh+1:12*Nh)];
P_RE_B = [x(12*Nh+1:13*Nh)];
P_RE_g = [x(13*Nh+1:14*Nh)];
P_RE_Ele = [x(14*Nh+1:15*Nh)];
P_g_B = [x(15*Nh+1:16*Nh)];
P_g_Ele = [x(16*Nh+1:17*Nh)];
P_B_Ele = [x(17*Nh+1:18*Nh)];
P_B_g = [x(18*Nh+1:19*Nh)];
% Time span 
span = 1:Nh;
span_soc = 0:Nh;

%##################
% Control setpoints
figure;
yyaxis left % For power (kW)

plot(span, Pb_c - Pb_d, '--c*', span, Pel - lambda_pv.*Ppv_available - lambda_w.*Pwt_available, '-r', span, Pel - lambda_pv.*Ppv_available - lambda_w.*Pwt_available - Pg_Im + Pg_Ex - Pb_d + Pb_c, 'g--', ...
    span, Pg_Im -Pg_Ex,'-', 'LineWidth', 2);
ylabel('P_{B}[kW], P_{Ele}- P_{RE}[kW], \Sigma P[kW], P_{G}[kW]');

yyaxis right % For SOC (%)
plot(span_soc, SOC, '--*', 'LineWidth', 2);
ylabel('Battery SOC(%)');

% title('Control setpoints and Battery SOC evolution');
xlabel('Time (h)');


legend('$P_{B}$', '$P_{Ele} - P_{RE}$','$\Sigma P$', ...
    '$P_{G}$', '$Battery \:SOC$', 'Location', 'northwest', 'Interpreter','latex')
grid on

%##################
% Plot curtailed power
figure;
subplot(2,1,1)
plot(span, Pwt_available,'r',span, (1-lambda_w).*Pwt_available,'b--*','LineWidth',1.2);
legend('$Available \:Wind$','$Curtailed \:wind$', 'Interpreter', 'latex')
subplot(2,1,2)
plot(span, Ppv_available,'r',span, (1-lambda_pv).*Ppv_available,'b--*','LineWidth',1.2);
legend('$Available \:PV$','$Curtailed \:PV$', 'Interpreter', 'latex')

%##################
% RE used and Pele
figure;
subplot(2,1,1)
plot(span, Pel,'r',span, (lambda_pv).*Ppv_available + (lambda_w).*Pwt_available,'b--*','LineWidth',1.2);
legend('$P_{Ele}$','$P_{RE}$', 'Interpreter', 'latex')
subplot(2,1,2)
plot(span, Pel,'r',span, Ppv_available + Pwt_available,'b--*','LineWidth',1.2);
legend('$P_{Ele}$','$P_{RE} available$', 'Interpreter', 'latex')

%##################
% H2 storage SOC, and internal setpoints
figure;
subplot(2,1,1)
plot(span,H_El2Hs, 'g-',span,H_El2Hd, 'r-',span,H_Hs2Hd, 'b-', span,HDemand, 'm--*', 'LineWidth', 1.2)
ylabel("H_{El2Hs}, H_{EL2Hd}, H_{Hs2Hd}, HH_{Demand} [kg]")
xlabel('Time [h]')
legend('$H_{El2Hs}$','$H_{EL2Hd}$','$H_{Hs2Hd}$', '$HH_{Demand}$','Interpreter','latex');
subplot(2,1,2)
yyaxis left
plot(span, Pel, span, Ppv_available+Pwt_available, 'LineWidth',1.2)
ylabel('P_{Ele}');
xlabel('Time [h]');
yyaxis right
plot(span_soc, SOC_Hs, 'LineWidth',1.2)
ylabel('SOC_{Hs} (%)')
legend('$P_{Ele}$', '$P_{RE}$', '$SOC_{Hs}$', 'Interpreter','latex')

%##################
% Grid Expenxes/Revenues
figure;
y1 = -Cg_Im.*Pg_Im; % Grid charges 
y2 = Cg_Ex.*Pg_Ex;  % Grid revenues
h1 = bar(span, y1);
hold on;
h2 = bar(span, y2);
hold off;
% Set color for positive and negative bars
h1.FaceColor = 'flat';
h2.FaceColor = 'flat';
h1.CData = repmat([0, 1, 0], Nh, 1);  % Green for positive bars
h2.CData = repmat([1, 0, 0], Nh, 1);  % Red for negative bars

% Add labels and title
xlabel('Time(h)');
y_label = sprintf('Cost/Revenue(Euros)');
ylabel(y_label,'Interpreter','latex');
% title('Grid cost/revenue in Euros','Interpreter','latex');

%##################
% RE distribution
annotation.title = 'RE distribution';
annotation.xlabel = 'Time[h]';
annotation.ylabel = 'P_{RE2B}, P_{RE2G}, P_{RE2Ele} [kW]';
annotation.label={"P_{RE2B}", "P_{RE2G}", "P_{RE2Ele}"};
Pagg = [P_RE_B'; P_RE_g'; P_RE_Ele'];
Power_Contribution_Disribution(span, annotation, Pagg);

%##################
% Electrolyzer contribution
annotation.title = 'Electrolyzer contribution';
annotation.xlabel = 'Time(h)';
annotation.ylabel = 'RE,B,G(kWh)';
annotation.label={"P_{B2ELe}", "P_{G2Ele}", "P_{RE2Ele}"};
Pagg = [P_B_Ele'; P_g_Ele'; P_RE_Ele'];
Power_Contribution_Disribution(span, annotation, Pagg);

%##################
% Battery contribution
annotation.title = 'Battery contribution';
annotation.xlabel = 'Time(h)';
annotation.ylabel = 'RE,G(kWh)';
annotation.label={"P_{RE2B}", "P_{G2B}"};
Pagg = [P_RE_B'; P_g_B'];
Power_Contribution_Disribution(span, annotation, Pagg);

%##################
% Grid contribution
annotation.title = 'Grid contribution';
annotation.xlabel = 'Time(h)';
annotation.ylabel = 'RE,B (kWh)';
annotation.label={"P_{RE2G}", "P_{B2G}"};
Pagg = [P_RE_g'; P_B_g'];
Power_Contribution_Disribution(span, annotation, Pagg);

%##################
% CO2 Emission contribution hourly basis 
annotation.title = 'Hourly CO2 Emission contribution';
annotation.xlabel = 'Time(h)';
annotation.ylabel = 'PVE,WTE,GRE(gco2perkWh)';
annotation.label={"PVE", "WTE", "GRE"};
Pagg = [Ppv_available'.*lambda_pv'*CO2.PVE; Pwt_available'.*lambda_w'*CO2.WTE; (Pg_Im' - CO2.ExportedREConsidered*Pg_Ex')*CO2.GRE];
Power_Contribution_Disribution(span, annotation, Pagg);

%##################
% CO2 Emission contribution daily basis
annotation.title = 'Daily CO2 Emission contribution';
annotation.xlabel = 'Time(Day)';
annotation.ylabel = 'PVE WTE GRE (gco2perkWh)';
annotation.label={"PVE", "WTE", "GRE"};
Pagg1 = zeros(Ndays,3);
for k = 1:Ndays
    Pagg1(k,:) = [sum(Ppv_available((k-1)*T+1:T*k)'.*lambda_pv((k-1)*T+1:T*k)')*CO2.PVE, sum(Pwt_available((k-1)*T+1:T*k)'.*lambda_w((k-1)*T+1:T*k)')*CO2.WTE, sum((Pg_Im((k-1)*T+1:T*k)' - CO2.ExportedREConsidered*Pg_Ex((k-1)*T+1:T*k)'))*CO2.GRE];
end
Power_Contribution_DisributionCO2(1:Ndays, annotation, Pagg1);


