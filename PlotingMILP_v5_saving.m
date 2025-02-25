function [] = PlotingMILP_v5_saving(x, Nh, T, Ndays, HDemand, Ppv_available, Pwt_available, Cg_Im, Cg_Ex, CO2, init_soc, folder, index)
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
% Time span 
span = 1:Nh;
span_soc = 0:Nh;

%%
% close all;
% figure;
% plot(1:T, P_load, 1:T, P_pv, 1:T, P_wind, 'r--*', 'LineWidth', 2);
% ylabel('Power (kW)', 'Interpreter', 'latex');
% xlabel('Time (h)');
% figure;
% yyaxis left % For power (kW)
% plot(1:T, Pb_c, 1:T, -Pb_d, 1:T, P_load - P_pv - P_wind, 'r--*', 1:T, P_load - P_pv - P_wind- Pg_Im + Pg_Ex - Pb_d + Pb_c, 'g--o', ...
%     1:T, -Pg_Ex, 1:T, Pg_Im, 'LineWidth', 2);
% ylabel('Power (kW)', 'Interpreter', 'latex');
% yyaxis right % For SOC (%)
% plot(0:T, SOC, 'LineWidth', 2);
% ylabel('$SOC(\%)$', 'Interpreter', 'latex');
% legend('$Pb_{c}$', '$Pb_{d}$', '$P_{load} - P_{pv} - P_{wind}$','$\Sigma P$', ...
%     '$P^{grid}_{out}$', '$P^{grid}_{in}$', '$SOC$', 'Location', 'northwest', 'Interpreter','latex')
% grid on;
figure;
yyaxis left % For power (kW)

plot(span, Pb_c - Pb_d, '--c*', span, Pel - lambda_pv.*Ppv_available - lambda_w.*Pwt_available, '-r', span, Pel - lambda_pv.*Ppv_available - lambda_w.*Pwt_available - Pg_Im + Pg_Ex - Pb_d + Pb_c, 'g--', ...
    span, Pg_Im -Pg_Ex,'-', 'LineWidth', 2);
ylabel('P_{B}[kW], P_{Ele}- P_{RE}[kW], \Sigma P[kW], P_{G}[kW]','FontName', 'Times New Roman', 'FontSize', 12);

yyaxis right % For SOC (%)
plot(span_soc, SOC, '--*', 'LineWidth', 2);
ylabel('Battery SOC(%)','FontName', 'Times New Roman', 'FontSize', 12);

% title('Control setpoints and Battery SOC evolution');
xlabel('Time [h]','FontName', 'Times New Roman', 'FontSize', 12);


legend('$P_{B}$', '$P_{Ele} - P_{RE}$','$\Sigma P$', ...
    '$P_{G}$', '$Battery \:SOC$','NumColumns',2, 'Location', 'best', 'Interpreter','latex','FontName', 'Times New Roman', 'FontSize', 12);
% grid on
set(gca,'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, folder+"/control_setpoints-"+num2str(index)+".svg")
%% H2 storage SOC, and internal setpoints
% figure;
% yyaxis left;
% plot(span,H_El2Hs, 'g-',span,H_El2Hd, 'r-',span,H_Hs2Hd, 'b-', span,HDemand, 'm--*', 'LineWidth', 1.2)
% ylabel("H_{El2Hs}, H_{EL2Hd}, H_{Hs2Hd}, HH_{Demand} [kg]")
% yyaxis right;
% plot(span_soc,SOC_Hs, "LineWidth",1.2)
% ylabel('$SOC_{Hs} (\%)$')
% legend('$H_{El2Hs}$','$H_{EL2Hd}$','$H_{Hs2Hd}$', '$HH_{Demand}$', '$SOC_{Hs}$','Interpreter','latex');
% xlabel('Time(h)');
% % title('Electrolyzer power contribution and heating storage SOC evolution')
% % Electrolyzer output hydrogen flow 
% figure;
% yyaxis left;
% plot(span,H_El2Hs, 'g-',span,H_El2Hd, 'r-',span,H_Hs2Hd, 'b-', span,HDemand, 'm--*', 'LineWidth', 1.2)
% ylabel("H_{El2Hs}, H_{EL2Hd}, H_{Hs2Hd}, HH_{Demand} [kg]")
% yyaxis right;
% plot(span, Pel, '-', 'LineWidth', 1.2);
% ylabel('P_{Ele}');
% legend('$H_{El2Hs}$','$H_{EL2Hd}$','$H_{Hs2Hd}$', '$HH_{Demand}$', '$P_{Ele}$' ,'Interpreter','latex');
% % title('Hydrogen flow')

figure;
subplot(2,1,1)
plot(span,H_El2Hs, 'g-',span,H_El2Hd, 'r-',span,H_Hs2Hd, 'b-', span,HDemand, 'm--*', 'LineWidth', 1.2)
ylabel("H_{El2Hs}, H_{EL2Hd}, H_{Hs2Hd}, HHD [kg]",'FontName', 'Times New Roman', 'FontSize', 12);
% xlabel('Time [h]','FontName', 'Times New Roman', 'FontSize', 12);
legend('$H_{El2Hs}$','$H_{EL2Hd}$','$H_{Hs2Hd}$', '$HH_{Demand}$','NumColumns',2, 'Interpreter','latex','FontName', 'Times New Roman', 'FontSize', 12, 'Location','best');
subplot(2,1,2)
yyaxis left
plot(span, Pel, span, Ppv_available+Pwt_available, 'LineWidth',1.2)
ylabel('P_{Ele} [kW]','FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time [h]','FontName', 'Times New Roman', 'FontSize', 12);
yyaxis right
plot(span_soc, SOC_Hs, 'LineWidth',1.2)
ylabel('SOC_{Hs} [%]','FontName', 'Times New Roman', 'FontSize', 12);
legend('$P_{Ele}$', '$P_{RE}$', '$SOC_{Hs}$','NumColumns',2,  'Interpreter','latex','FontName', 'Times New Roman', 'FontSize', 12, 'Location','best');
set(gca,'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, folder+"/electrolyzer_setpoints-"+num2str(index)+".svg");