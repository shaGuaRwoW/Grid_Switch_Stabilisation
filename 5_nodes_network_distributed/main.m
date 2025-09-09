clc, clear, close all;

node_1_0 = pi / 180 * [5, 0]';
node_2_0 = pi / 180 * [-8, 0]';
node_3_0 = pi / 180 * [2, 0]';
node_4_0 = pi / 180 * [-8, 0]';
node_5_0 = pi / 180 * [-3, 0]';

n_max = 200;
dt = 0.2;

%%

load('ss_matrices_A_conns.mat');

% load('controllers.mat');
% load('null_controller.mat')
% load('distributed_sw_stable_ctrl_sparse_slow.mat');
load('distributed_sw_stable_ctrl_no_sparse_fast.mat');

n_T_local = 5;
n_states = 5 * n_T_local;
n_of_systems = 9;

T_A = zeros(n_states, n_states, n_of_systems);
T_B = zeros(n_states, 10, n_of_systems);
T_C = zeros(5, n_states, n_of_systems);
T_D = zeros(5, 10, n_of_systems);

ss_nom = ss(A_matrices(:, :, 9), B, C, 0, 1);

%% Input signals definition

tt = 1 : n_max;

u1 = double(tt >= 120);
u1 = timeseries(u1, tt);

u5 = 0.1 * u1;
u1 = 0 * u1;
u2 = 0 * u1;
u3 = 0 * u1;
u4 = 0 * u1;

w1 = 0 * u1;
w2 = 0 * u2;
w3 = 0 * u3;
w4 = 0 * u4;
w5 = 0 * u5;

% w3 = 0.5 * randn(n_max, 1);
% w3 = sin(2 * pi / 1.8 * tt);
% w3 = timeseries(w3, tt);

%% Network configuration signals definition

% mode_sel_dat = 9 * double(tt >= 0 & tt <= 30) + 5 * double(tt >= 31 & tt <= 150) ...
%     + 9 * double(tt >= 151 & tt <= 350) + 1 * double(tt >= 351 & tt <= n_max);
% mode_sel_dat = randi(9, 1, n_max);
% mode_sel_dat = mod(tt, 9) + 1;
% mode_sel_dat = 9 * ones(1, n_max);
load('random_switch.mat')
mode_sel = timeseries(mode_sel_dat, tt);

% Considered failure modes:
% 1 : line_12 and line_21 down
% 2 : line_24 and line_42 down
% 3 : line_34 and line_43 down
% 4 : line_13 and line_31 down
% 5 : line_15 and line_51 down
% 6 : line_25 and line_52 down
% 7 : line_45 and line_54 down
% 8 : line_35 and line_53 down
% 9 : ALL GOOD

%% Simulate

% sim_out = sim("net_model_5_nodes.slx");
sim_out = sim("net_model_5_nodes_distributed_controllers.slx");

%% Unpack output

delta_1 = sim_out.net_out.Data(:, 1);
delta_2 = sim_out.net_out.Data(:, 3);
delta_3 = sim_out.net_out.Data(:, 5);
delta_4 = sim_out.net_out.Data(:, 7);
delta_5 = sim_out.net_out.Data(:, 9);

omega_1 = sim_out.net_out.Data(:, 2);
omega_2 = sim_out.net_out.Data(:, 4);
omega_3 = sim_out.net_out.Data(:, 6);
omega_4 = sim_out.net_out.Data(:, 8);
omega_5 = sim_out.net_out.Data(:, 10);

ss_delta_1 = sim_out.ss_out.Data(:, 1);
ss_delta_2 = sim_out.ss_out.Data(:, 3);
ss_delta_3 = sim_out.ss_out.Data(:, 5);
ss_delta_4 = sim_out.ss_out.Data(:, 7);
ss_delta_5 = sim_out.ss_out.Data(:, 9);

ss_omega_1 = sim_out.ss_out.Data(:, 2);
ss_omega_2 = sim_out.ss_out.Data(:, 4);
ss_omega_3 = sim_out.ss_out.Data(:, 6);
ss_omega_4 = sim_out.ss_out.Data(:, 8);
ss_omega_5 = sim_out.ss_out.Data(:, 10);

%% Plots

% Deltas

figure
subplot(2, 3, 1)
hold on
plot(sim_out.tout * dt, delta_1, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_delta_1, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 1', 'Test SS - Delta 1')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 2)
hold on
plot(sim_out.tout * dt, delta_2, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_delta_2, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 2', 'Test SS - Delta 2')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 3)
hold on
plot(sim_out.tout * dt, delta_3, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_delta_3, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 3', 'Test SS - Delta 3')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 4)
hold on
plot(sim_out.tout * dt, delta_4, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_delta_4, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 4', 'Test SS - Delta 4')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 5)
hold on
plot(sim_out.tout * dt, delta_5, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_delta_5, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Delta 5', 'Test SS - Delta 5')
set(gca, 'FontSize', 15)
grid

% Omegas

figure
subplot(2, 3, 1)
hold on
plot(sim_out.tout * dt, omega_1, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_omega_1, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 1', 'Test SS - Omega 1')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 2)
hold on
plot(sim_out.tout * dt, omega_2, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_omega_2, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 2', 'Test SS - Omega 2')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 3)
hold on
plot(sim_out.tout * dt, omega_3, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_omega_3, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 3', 'Test SS - Omega 3')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 4)
hold on
plot(sim_out.tout * dt, omega_4, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_omega_4, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 4', 'Test SS - Omega 4')
set(gca, 'FontSize', 15)
grid

subplot(2, 3, 5)
hold on
plot(sim_out.tout * dt, omega_5, 'LineWidth', 1.5)
plot(sim_out.tout * dt, ss_omega_5, 'LineStyle', '-.', 'LineWidth', 1)
legend('Multi-nodes scheme - Omega 5', 'Test SS - Omega 5')
set(gca, 'FontSize', 15)
grid

% figure
% plot(tt * dt, switch_signal.data)
% grid

%%

figure
hold on
grid
plot(sim_out.tout * dt, delta_1 * 180 / pi, 'LineWidth', 1.5)
plot(sim_out.tout * dt, delta_2 * 180 / pi, 'LineWidth', 1.5)
plot(sim_out.tout * dt, delta_3 * 180 / pi, 'LineWidth', 1.5)
plot(sim_out.tout * dt, delta_4 * 180 / pi, 'LineWidth', 1.5)
plot(sim_out.tout * dt, delta_5 * 180 / pi, 'LineWidth', 1.5)
legend('$\delta_{1}$', '$\delta_{2}$', '$\delta_{3}$', '$\delta_{4}$', '$\delta_{5}$','Interpreter','latex','FontSize',35)
% legend('Delta 1', 'Delta 2', 'Delta 3', 'Delta 4', 'Delta 5')

ylabel("Electrical angle (degrees)")
xlabel("Time (seconds)")
set(gca, 'FontSize', 35)