clc, clear, close all;

node_1_0 = [5; 1];
node_2_0 = [1; 1];

n_max = 500;
dt = 0.2;

%%

load('state_spaces.mat');

%% 

tt = 1 : n_max;

switch_signal = min(1, ones(1, n_max) + square(tt * 2 * pi / 80) .* double(tt >= 100));
switch_signal = timeseries(switch_signal', tt);

u1 = double(tt >= 5);
u1 = timeseries(u1, tt);

u1 = 0 * u1;
u2 = 0 * u1;

%%

% K = -place(ss_nom.A, ss_nom.B, 2 * [0.45, 0.4, 0.49, 0.4]);
% ss_ctr_nom = ss(zeros(4), zeros(4, 4), zeros(2, 4), K);
% 
% K = -place(ss_fault.A, ss_fault.B,  2  * [0.45, 0.4, 0.49, 0.4]);
% ss_ctr_fault = ss(zeros(4), zeros(4, 4), zeros(2, 4), K);

sim_out = sim("net_model_2_nodes_ext_switch.slx");

%% Plots

figure
subplot(2, 1, 1)
hold on
plot(sim_out.tout * dt, sim_out.out_1.Data(:, 1), 'LineWidth', 1.5)
plot(sim_out.tout * dt, sim_out.test_out_1.Data(:, 1), 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 1', 'Test SS - Delta 1')
grid

subplot(2, 1, 2)
hold on
plot(sim_out.tout * dt, sim_out.out_1.Data(:, 2), 'LineWidth', 1.5)
plot(sim_out.tout * dt, sim_out.test_out_1.Data(:, 2), 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 1', 'Test SS - Omega 1')
grid

figure
subplot(2, 1, 1)
hold on
plot(sim_out.tout * dt, sim_out.out_2.Data(:, 1), 'LineWidth', 1.5)
plot(sim_out.tout * dt, sim_out.test_out_2.Data(:, 1), 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 2', 'Test SS - Delta 2')
grid

subplot(2, 1, 2)
hold on
plot(sim_out.tout * dt, sim_out.out_2.Data(:, 2), 'LineWidth', 1.5)
plot(sim_out.tout * dt, sim_out.test_out_2.Data(:, 2), 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 2', 'Test SS - Omega 2')
grid

figure
plot(tt * dt, switch_signal.data)
grid