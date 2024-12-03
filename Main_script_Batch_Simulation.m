%% Tumor Growth Model: Batch Simulations with Visualization and Animation
clear; clc;

global alpha beta gamma

% Parameter ranges for simulations
alpha_values = [0.7, 0.8, 0.9];   % Nutrient availability
beta_values = [0.4, 0.5, 0.6];    % Growth rate
gamma_values = [8, 10, 12];       % Nutrient saturation

% Simulation settings
dx = 1; X = 210; dt = 0.004; T = 16; % Space and time steps
c0 = 1;                              % Initial nutrient concentration

% Initialize results container
sim_results = struct(); % Store all results for later analysis

% Simulation loop: Vary parameters
sim_num = 0; % Simulation counter
for alpha = alpha_values
    for beta = beta_values
        for gamma = gamma_values
            sim_num = sim_num + 1; % Increment simulation counter

            % Define functions (using current beta)
            f = @(c) 0.5 * (1 - tanh(4 * c - 2));   % Function f(c)
            h = @(c) 0.5 * f(c);                   % Function h(c)
            g = @(c) beta * exp(beta * c);         % Gompertz growth rate, g(c)

            % Create spatial and temporal grids
            x = dx:dx:X; 
            Nx = round(X/dx); 
            Nt = round(T/dt);

            % Initialize cell densities and nutrient concentration
            p = exp(-0.1 .* x); % Initial density of proliferating cells
            q = zeros(1, Nx);   % Quiescent cells start at 0
            n = zeros(1, Nx);   % Necrotic cells start at 0
            c = zeros(1, Nx);   % Nutrient concentration

            % Arrays to store results for this simulation
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

            % Store simulation results
            sim_results(sim_num).alpha = alpha;
            sim_results(sim_num).beta = beta;
            sim_results(sim_num).gamma = gamma;
            sim_results(sim_num).P = P;
            sim_results(sim_num).Q = Q;
            sim_results(sim_num).N = N;
            sim_results(sim_num).x = x;
            sim_results(sim_num).time = (0:dt:T-dt);
        end
    end
end

% Save results to file
save('tumor_growth_simulations.mat', 'sim_results');
disp('Simulations complete. Results saved to "tumor_growth_simulations.mat".');

%% Visualize Results for a Single Simulation (First Simulation)
sim = sim_results(1); % Access the first simulation results

% Plot Proliferating Cells
figure(1);
hold on;
for t = 1:500:Nt
    plot(sim.x, sim.P(t, :), 'LineWidth', 1.2); % Plot at intervals
end
xlabel('Space, x');
ylabel('Proliferating Cells, p(x, t)');
title(sprintf('Proliferating Cells (alpha=%.2f, beta=%.2f, gamma=%.2f)', ...
              sim.alpha, sim.beta, sim.gamma));
legend('t=0', 't=2', 't=4', 't=6', 't=8', 't=10', 't=12', 't=14');
grid on;

% Animation of tumor growth
figure(2); % New figure for animation
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
        theta = 2 * pi * rand(1, round(sim.P(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'b.', 'MarkerSize', 5);
        
        % Quiescent cells
        theta = 2 * pi * rand(1, round(sim.Q(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'r.', 'MarkerSize', 5);
        
        % Necrotic cells
        theta = 2 * pi * rand(1, round(sim.N(t, i) * 100)); 
        r = i + randn(1, length(theta)); 
        plot(r .* cos(theta), r .* sin(theta), 'k.', 'MarkerSize', 5);
    end
    drawnow;
    writeVideo(aviFile, getframe(gcf));
end

close(aviFile);
disp(['Animation saved as ', fullfile(pwd, 'tumor_growth.avi')]);

% Load simulation results
load('tumor_growth_simulations.mat');

% Example: Plot proliferating cells for the first simulation
sim = sim_results(1); % Access first simulation
figure;
imagesc(sim.x, sim.time, sim.P); % Plot P(x, t) as a heatmap
xlabel('Space, x');
ylabel('Time, t');
title(sprintf('Proliferating Cells (alpha=%.2f, beta=%.2f, gamma=%.2f)', ...
              sim.alpha, sim.beta, sim.gamma));
colorbar;

% Compare total live cells across simulations
figure;
hold on;
for i = 1:length(sim_results)
    sim = sim_results(i);
    total_live = sum(sim.P, 2) + sum(sim.Q, 2); % Proliferating + Quiescent
    plot(sim.time, total_live, 'DisplayName', ...
         sprintf('alpha=%.2f, beta=%.2f, gamma=%.2f', sim.alpha, sim.beta, sim.gamma));
end
xlabel('Time (t)');
ylabel('Total Live Cells');
title('Comparison of Live Cell Counts Across Simulations');
legend show;
grid on;
