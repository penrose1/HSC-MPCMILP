% Define system matrices A and B for a simple system
A = [0 1; 0 0]; % Example: simple double integrator
B = [0; 1];

% Define the reference trajectory (for example, a constant reference)
x_ref = [1; 0]; % constant reference, could be time-varying

% Define weighting matrices Q and R
Q = [10 0; 0 1];  % penalize deviation from the reference (state error)
R = .1;       % penalize the control effort

% Solve the Riccati equation to find the optimal gain matrix K
[K, P, E] = lqr(A, B, Q, R);

% Display the gain matrix
disp('Optimal Gain Matrix K:');
disp(K);
% Define the closed-loop system dynamics
A_cl = A - B * K; % Closed-loop A matrix
B_cl = B;          % Closed-loop B matrix

% Define the initial state (let's start at the origin)
x0 = [0; 0];

% Time vector for simulation
t = 0:0.01:10;  % Simulate for 10 seconds

% Define the system using MATLAB's ODE solver
% The control law u = -K(x - x_ref) will drive the state towards x_ref
ode = @(t, x) (A_cl * x + B_cl * (-K * (- x_ref)));

% Simulate the system using ode45
[t, x] = ode45(ode, t, x0);

% Plot the results
figure;
plot(t, x(:, 1), 'r', 'LineWidth', 2); % State 1 vs time
hold on;
plot(t, x(:, 2), 'b', 'LineWidth', 2); % State 2 vs time
xlabel('Time (s)');
ylabel('States');
legend('x_1', 'x_2');
title('LQR Tracking Problem');
grid on;
