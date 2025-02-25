
%% 
clear,clc, close all;
dt = 1; % Sampling time of 1h
Ndays = 1;
T = 24; % Simulation time 
Nh = T*Ndays; % Total hours in Ndays

Hstart = 400; % Start hour of simulation
index = Hstart+1:Hstart+Nh;
%% Electricity prices
Cg_Ex = 0.065*ones(Nh, 1);% + 0.* randn(NT, 1); % Grid power cost (per kWh)
Cg_im = 0.065*ones(Nh, 1);% + 0. * randn(NT, 1); % Feed-in tariff (per kWh)

%% Wind power plant 
WT = 5; % Number of wind turbines
%% PV power plant 
% [3.046955584757793,1.335548309315877,97.648898886325510]
% [3.128094492574218,1.711850739314564,80.108524379089760]
Apv = 100; % Number of PV hectares 


%% Electrolyzer unit
Electrolyzer.eta_e = 0.6; % Electrolyzer efficiency 

%% Grid 
Grid.Pgrid_max_Im = 10000;
Grid.Pgrid_max_Ex =  000; % kW
Grid.Pgrid_min_Ex = 0;
