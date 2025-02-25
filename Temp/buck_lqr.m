% Define the system parameters
L = 10e-6;  % Inductance (Henries)
C = 100e-6; % Capacitance (Farads)
R = 10;     % Load resistance (Ohms)
Vin = 24;   % Input voltage (Volts)

% State-space matrices
A = [0, -1/L;
     1/C, -1/(R*C)];
B = [Vin/L; 
     0];

% Initial conditions
x0 = [0; 0];  % Initial state [i_L(0), v_out(0)]
xref = [1.2; 12];
Q = [10 0; 0 10];
R = 1;

[K, S, E] = lqr(A, B, Q, R);

fnc = @(t,x) ((A - B*K)*x + B*K*xref);
% Time vector for simulation
t = 0:0.001:0.1;  % Simulate for 0.1 seconds with a step of 1ms

% Duty cycle input as a function of time (e.g., 50% duty cycle)
D = 0.5 * ones(size(t));  % Constant duty cycle for simplicity

% Simulate the system response using ode45
% Define the system of differential equations
dxdt = @(t, x) A * x + B * Vin * D(round(t*1000)+1);  % D(t) applied at each time step

% Solve the system using ode45
[t, x] = ode45(fnc, t, x0);

% Extract the state variables
i_L = x(:, 1);  % Inductor current
v_out = x(:, 2); % Output voltage

% Plot the results
figure;

subplot(3,1,1);
plot(t, i_L, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Inductor Current (A)');
title('Inductor Current vs Time');
grid on;

subplot(3,1,2);
plot(t, v_out, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Output Voltage (V)');
title('Output Voltage vs Time');
grid on;

for k=1:length(x)
    u(k) = -K*x(k,:)';
end
subplot(3,1,3);
plot(t, u, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Control action');
title('Control input vs Time');
grid on;
