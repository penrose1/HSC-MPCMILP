function stop = PlotSwarmR3Euclidean(options, state, flag)
    stop = false;  % Indicate that the optimization should continue
    
    % Extract particle positions from state
    particles = options.swarm;  % This is where particles are stored
    
    % Plot the particle positions in 3D
    clf;  % Clear the current figure
    plot3(particles(:,1), particles(:,2), particles(:,3), 'ro');  % Plot particles as red circles
    hold on;
    
    % Optional: Plot the objective function's 3D surface (if you want)
    % [X, Y] = meshgrid(-5:0.1:5, -5:0.1:5);
    % Z = X.^2 + Y.^2;  % Example 3D surface for visualization
    % surf(X, Y, Z);  % Plot the surface
    % alpha 0.5;  % Make the surface semi-transparent
    
    % Formatting the plot
    title('Particle Swarm Optimization (3D)');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis([0 105 0 105 0 900]);  % Set axis limits for better visualization
    drawnow;  % Update the plot
end
