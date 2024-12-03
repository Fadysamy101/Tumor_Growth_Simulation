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
