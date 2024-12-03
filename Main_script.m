%% Tumor Growth Model with Visualization and Animation
clear; clc;
global alpha beta gamma

% Parameters
alpha = 0.8; 
beta = 0.5; 
gamma = 10;

dx = 1; X = 210; dt = 0.004; T = 16;
c0 = 1;

% Define functions
f = @(c) 0.5 * (1 - tanh(4 * c - 2));   % Function f(c)
h = @(c) 0.5 * f(c);                   % Function h(c)
g = @(c) beta * exp(beta * c);         % Gompertz growth rate, g(c)

% Spatial and temporal grid
x = dx:dx:X; 
Nx = round(X/dx); 
Nt = round(T/dt);

% Initialize cell densities and nutrient concentration
p = exp(-0.1 .* x); % Initial density of proliferating cells
q = zeros(1, Nx);   % Quiescent cells start at 0
n = zeros(1, Nx);   % Necrotic cells start at 0
c = zeros(1, Nx);   % Nutrient concentration

% Arrays to store results
P = zeros(Nt, Nx); 
Q = zeros(Nt, Nx); 
N = zeros(Nt, Nx);

% Store initial conditions
P(1, :) = p; 
Q(1, :) = q; 
N(1, :) = n;

% Time-stepping loop
for k = 1:Nt
    r = p + q; % Total live cell density
    c = (c0 .* gamma ./ (gamma + p)) .* (1 - alpha .* (p + q + n)); % Nutrient concentration
    
    % Set flux terms to zero at boundaries
    u = zeros(1, Nx);
    v = zeros(1, Nx);
    
    % Calculate flux terms for internal points
    for i = 2:Nx-1
        u(i) = ((p(i+1) - p(i-1)) * r(i) * (r(i+1) - r(i-1)) + ...
                4 * p(i) * r(i) * (r(i+1) - 2 * r(i) + r(i-1)) - ...
                p(i) * (r(i+1) - r(i-1))^2) / (2 * (dx * r(i))^2);
        v(i) = ((q(i+1) - q(i-1)) * r(i) * (r(i+1) - r(i-1)) + ...
                4 * q(i) * r(i) * (r(i+1) - 2 * r(i) + r(i-1)) - ...
                q(i) * (r(i+1) - r(i-1))^2) / (2 * (dx * r(i))^2);
    end
    
    % Update cell densities using the finite difference scheme
    nextp = p + dt .* (u + g(c) .* p .* (1 - (p + q + n)) - f(c) .* p);
    nextq = q + dt .* (v + f(c) .* p - h(c) .* q);
    nextn = n + dt .* (h(c) .* q);
    
    % Update current values
    p = nextp;
    q = nextq;
    n = nextn;

    % Store results
    P(k, :) = p; 
    Q(k, :) = q; 
    N(k, :) = n;
end

% Scale the model results to match experimental data
scale = 443.7249; % Scaling factor
P_scaled = P * scale; 
Q_scaled = Q * scale; 
N_scaled = N * scale;

% Compute total live and dead cells
total_live = P_scaled + Q_scaled; % Proliferating + Quiescent
total_dead = N_scaled;           % Necrotic

% Time points for plotting
time_points = 0:dt:T-dt; 

% Experimental data from Nirmala et al.
exp_time = [0, 1, 2, 3, 4, 5, 6]; % Time points
exp_live = [3000, 7015, 10000, 13000, 15500, 18000, 19500]; % Live cell counts
exp_dead = [0, 0, 1500, 3000, 4500, 6000, 7000];            % Dead cell counts

% Plot comparison of live and dead cell counts
figure(4); % New figure for comparison
hold on;

% Plot model predictions
plot(time_points, sum(total_live, 2), 'b-', 'LineWidth', 2); % Live cells
plot(time_points, sum(total_dead, 2), 'r-', 'LineWidth', 2); % Dead cells

% Plot experimental data
plot(exp_time, exp_live, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Experimental live cells
plot(exp_time, exp_dead, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Experimental dead cells

% Customize the plot
xlabel('Time (t)');
ylabel('Cell Count');
title('Comparison of Model Predictions with Experimental Data');
legend('Model: Live Cells', 'Model: Dead Cells', ...
       'Experimental: Live Cells', 'Experimental: Dead Cells', ...
       'Location', 'Northwest');
grid on;

% Plot Proliferating Cells
figure(1);
hold on;
for t = 1:500:Nt
    plot(x, P(t, :), 'LineWidth', 1.2); % Plot at intervals
end
xlabel('Space, x');
ylabel('Proliferating Cells, p(x, t)');
title('Proliferating Cells Over Time');
legend('t=0', 't=2', 't=4', 't=6', 't=8', 't=10', 't=12', 't=14');
axis([0 210 0 1]);
grid on;

% Plot Quiescent Cells
figure(2);
hold on;
for t = 1:500:Nt
    plot(x, Q(t, :), 'LineWidth', 1.2); % Plot at intervals
end
xlabel('Space, x');
ylabel('Quiescent Cells, q(x, t)');
title('Quiescent Cells Over Time');
legend('t=0', 't=2', 't=4', 't=6', 't=8', 't=10', 't=12', 't=14');
axis([0 210 0 0.6]);
grid on;

% Plot Necrotic Cells
figure(3);
hold on;
for t = 1:500:Nt
    plot(x, N(t, :), 'LineWidth', 1.2); % Plot at intervals
end
xlabel('Space, x');
ylabel('Necrotic Cells, n(x, t)');
title('Necrotic Cells Over Time');
legend('t=0', 't=2', 't=4', 't=6', 't=8', 't=10', 't=12', 't=14');
axis([0 210 0 1]);
grid on;

% Animation of tumor growth
figure(6); % New figure for animation
hold on;
axis equal;
axis([-220 220 -220 220]);
aviFile = VideoWriter(fullfile(pwd, 'tumor_growth.avi'));
aviFile.FrameRate = 5; % Adjust frame rate
open(aviFile);

for t = 1:500:Nt
    clf; % Clear figure for new frame
    hold on;
    for i = 1:Nx
        % Proliferating cells
        theta = 2 * pi * rand(1, round(P(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'b.', 'MarkerSize', 5);
        
        % Quiescent cells
        theta = 2 * pi * rand(1, round(Q(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'r.', 'MarkerSize', 5);
        
        % Necrotic cells
        theta = 2 * pi * rand(1, round(N(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'k.', 'MarkerSize', 5);
    end
     % Plot legend (only needs to be set once outside the loop)
    hBlue = plot(NaN, NaN, 'b.', 'MarkerSize', 10); % Blue dot
    hRed = plot(NaN, NaN, 'r.', 'MarkerSize', 10);  % Red dot
    hBlack = plot(NaN, NaN, 'k.', 'MarkerSize', 10); % Black dot
    xlabel('Radial Distance (x)'); % X-axis label
    ylabel('Radial Distance (y)'); % Y-axis label
    title('Radial Distribution of Tumor Cell Types Over Time');
    legend([hBlue, hRed, hBlack], ...
        {'Blue: Proliferating Cells', 'Red: Quiescent Cells', 'Black: Necrotic Cells'}, ...
        'Location', 'northeast');
    drawnow;
    writeVideo(aviFile, getframe(gcf));
end

close(aviFile);
disp(['Animation saved as ', fullfile(pwd, 'tumor_growth.avi')]);
