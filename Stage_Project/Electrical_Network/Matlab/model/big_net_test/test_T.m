clc, clear, close all;

load('controllers.mat');
load('ss_matrices_A_conns.mat');
load('K_final.mat');

B = [0;
    1];
B = blkdiag(B, B, B, B, B);

C = [1, 0;
    0, 1];
C = blkdiag(C, C, C, C, C);

n = 10;
m = 5;
p = 10;

%% Compute K desired

% Uncomment those lines to use K - NRF
poles = 0.6 * [0.96, 0.96, 0.97, 0.99, 0.9, 0.96, 0.96, 0.97, 0.99, 0.9];
F_prime = -place(A_matrices(:, :, 9), B,  poles);
L_prime = -place(A_matrices(:, :, 9)', C, poles)';
K_final = ss(A_matrices(:, :, 9) + L_prime * C + B * F_prime, -L_prime, F_prime, 0, 1);

K_final.Ts = 1;

% Compute DCF for K desired

if ~all(abs(eig(K_final.a)) < 1)

    error('Controller is not stable')

end

U = K_final;
V = eye(p);

%

%% Compute T based on DCF

current_conf = 1;
J = C_ctrl(:, :, current_conf);
L = -B_ctrl(:, :, current_conf);
A_current = A_matrices(:, :, current_conf);

rez_fb = feedback(ss(A_current, B, C, 0, 1), K_final, 1);
if ~all(abs(eig(rez_fb.a)) < 1)

    error('Desired controller and system are not stable in feedback config')

end

% Compute DCF for current network model
r_cop = ss(A_current + B * J, [B -L], [J; C], eye(m + p), 1);
l_cop = ss(A_current + L * C, [-B L], [J; C], eye(m + p), 1);

M   = r_cop(1 : m, 1 : m);
Uz  = r_cop(1 : m, m + 1 : end);
N   = r_cop(m + 1 : end, 1 : m);
Vz  = r_cop(m + 1 : end, m + 1 : end);

Vtz =  l_cop(1 : m, 1 : m);
Utz = -l_cop(1 : m, m + 1 : end);
Nt  = -l_cop(m + 1 : end, 1 : m);
Mt  =  l_cop(m + 1 : end, m + 1 : end);

Z = Mt * V - Nt * U;
Q = M \ (U / Z - Uz);
Q = -ss(Q, 'min');

J = ss(A_current + B * J + L * C, [-L, B], [J; C], [zeros(m, p), eye(m); -eye(p), zeros(p, m)], 1);

J11 = J(1 : m, 1 : p);
J12 = J(1 : m, p + 1 : end);
J21 = J(m + 1 : end, 1 : p);
J22 = J(m + 1 : end, p + 1 : end);

T = lft(K_final, inv(J));
size(T)

T = ss(T, 'min')
size(T)

T_ref = ss(T_A(:, :, 9), T_B(:, :, 9), T_C(:, :, 9), T_D(:, :, 9), 1);