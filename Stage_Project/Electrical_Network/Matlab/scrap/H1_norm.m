clc, clear, close all;

% A = [-0.1, 100; -100, -0.1];
A = diag([-10, 190]);
T = randn(2);

A_p = T * A * inv(T);

tau = linspace(0.1, 1e5, 1e5)';

norm_val = zeros(1, 1e5);
for i = 1 : 1e5
    norm_val(i) = norm(eye(2) + tau(i) * A_p, 1);
end

plot(tau, norm_val)