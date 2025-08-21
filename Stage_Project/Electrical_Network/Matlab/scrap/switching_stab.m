clc, clear, close all;

x0 = [7; -50];

i_max = 8;

A_vec = {diag([-2, -1]), [-1, 1; -1, -1]...
    [-0.1, 100; -100, -0.1], [-1, -2; 2, -1]...
    [-0.1, 100; -100, -0.1], [-0.2, -2; 2, -0.2]...
    [-0.1, 100; -100, -0.1], [-1, -2; 2, -1]};

x_ev = [x0];
t = [];

for i = 1 : i_max
    eig(A_vec{i})
    T = randn(2);
    ss_temp = ss(T * A_vec{i} * inv(T), zeros(2, 1), eye(2), 0);
    [x_ev_temp, t_temp] = initial(ss_temp, x_ev(:, end), 3);
    x_ev = [x_ev, x_ev_temp'];
    t = [t; t_temp];
end

figure

subplot(2, 1, 1)
% plot(t, x_ev(1, :))
plot(x_ev(1, :))
legend(['x1'])

subplot(2, 1, 2)
% plot(t, x_ev(2, :))
plot(x_ev(2, :))
legend(['x2'])