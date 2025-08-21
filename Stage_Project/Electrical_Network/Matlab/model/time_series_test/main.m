clc, clear, close all;

L = [0, 1;
    1, 0];

hi = [0.3244
    1.5886];

di = [0.5650;
    0.7844];

l_mat = [0 0.6511;
    0.6511 0];

dt = 0.2;

node_1_0 = [0; 0];
node_2_0 = [0; 0];

%% Compose state-space

% l_mat(1, 2) = 0;

A11 = [1, dt;
    -hi(1) * dt * (l_mat(1, 1) + l_mat(1, 2)), 1 - hi(1) * di(1) * dt];

A12 = [0, 0;
    hi(1) * l_mat(1, 2) * dt, 0];

A21 = [0, 0;
    hi(2) * l_mat(2, 1) * dt, 0];

A22 = [1, dt;
    -hi(2) * dt * (l_mat(2, 1) + l_mat(2, 2)), 1 - hi(2) * di(2) * dt];

A = [A11, A12;
    A21, A22];

B = [0;
    1];
B = blkdiag(B, B);


C = [1, 0;
    0, 1];
C = blkdiag(C, C);

%% 

n_max = 500;
tt = 1 : 500;
L_long = repmat(L, 1, 1, n_max);
L_long(:, :, 250 : end) = [0, -1; 0, 0] + L_long(:, :, 250 : end);
L = timeseries(L_long, tt);

u1 = double(tt >= 5);
u1 = timeseries(u1, tt);

u2 = 0 * u1;

%%

sim_out = sim("net_model_2_nodes.slx");

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