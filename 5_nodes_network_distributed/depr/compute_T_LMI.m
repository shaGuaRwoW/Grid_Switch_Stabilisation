clc, clear, close all;

% load('controllers.mat');
load('ss_matrices_A_conns.mat');
% load('K_final.mat');

load('distributed_sw_stable_ctrl_no_sparse_fast.mat');
L = -B_ctrl;
J = C_ctrl;

B = [0;
    1];
B = blkdiag(B, B, B, B, B);

C = [1, 0;
    0, 1];
C = blkdiag(C, C, C, C, C);

n = 10;
m = 5;
p = 10;

n_of_systems = 9;

%%

ni = [5; 5; 5; 5; 5];

At = sdpvar(repmat(sum(ni), 1, n_of_systems), repmat(sum(ni), 1, n_of_systems));
Bt = sdpvar(repmat(sum(ni), 1, n_of_systems), repmat(sum(n), 1, n_of_systems));
Ct = sdpvar(repmat(sum(m), 1, n_of_systems), repmat(sum(ni), 1, n_of_systems));

% P = [];
% for i = 1 : 5
%     temp_mat = 1 * randn(ni(i));
%     P = blkdiag(P, temp_mat' * temp_mat);
% end

P = eye(25);

%%

constr = [];
for i = 1 : n_of_systems

    constr = [constr;

        [P, (At{i} * P)';
        At{i} * P, P] >= 1e-3;

    ];

    constr = [constr;
    
        At{i}(1 : 5, 6 : end) == zeros(5, 20);

        At{i}(6 : 10, 1 : 5) == zeros(5, 5);
        At{i}(6 : 10, 11 : end) == zeros(5, 15);

        At{i}(11 : 15, 1 : 10) == zeros(5, 10);
        At{i}(11 : 15, 16 : end) == zeros(5, 10);

        At{i}(16 : 20, 1 : 15) == zeros(5, 15);
        At{i}(16 : 20, 21 : end) == zeros(5, 5);

        At{i}(21 : end, 1 : 20) == zeros(5, 20);

    ];

    obj = 0;
end

optimize(constr, obj)

%%

n_states = sum(ni);
T_A = zeros(n_states, n_states, n_of_systems);
T_B = zeros(n_states, 10, n_of_systems);
T_C = zeros(5, n_states, n_of_systems);
T_D = zeros(5, 10, n_of_systems);

for i = 1 : n_of_systems

    r_cop = ss(A_matrices(:, :, i) + B * J(:, :, i), [B -L(:, :, i)], [J(:, :, i); C], eye(m + p), 1);
    l_cop = ss(A_matrices(:, :, i) + L(:, :, i) * C, [-B L(:, :, i)], [J(:, :, i); C], eye(m + p), 1);
    
    M   = r_cop(1 : m, 1 : m);
    Uz  = r_cop(1 : m, m + 1 : end);
    N   = r_cop(m + 1 : end, 1 : m);
    Vz  = r_cop(m + 1 : end, m + 1 : end);
    
    Vtz =  l_cop(1 : m, 1 : m);
    Utz = -l_cop(1 : m, m + 1 : end);
    Nt  = -l_cop(m + 1 : end, 1 : m);
    Mt  =  l_cop(m + 1 : end, m + 1 : end);

    A_temp = evalfr(Vz * Nt, 1);
    B_temp = evalfr(N, 1);
    C_temp = evalfr(Nt, 1);

    % Q_val = sdpvar(m, n);
    B_sdpvar = sdpvar(n_states, n);

    mata = randn(m, n_states);

    constr_temp = [A_temp + B_temp * mata * inv(eye(n_states) - value(At{i})) * B_sdpvar * C_temp == zeros(10, 5);
        B_sdpvar(1 : ni(1), 7 : 8)                  == zeros(5, 2);
        B_sdpvar(ni(1) + 1 : 2 * ni(1), 5 : 6)      == zeros(5, 2);
        B_sdpvar(2 * ni(1) + 1 : 3 * ni(1), 3 : 4)  == zeros(5, 2);
        B_sdpvar(3 * ni(1) + 1 : 4 * ni(1), 1 : 2)  == zeros(5, 2);
    ];

    optimize(constr_temp, 0);

    value(norm(A_temp + B_temp * mata * inv(eye(n_states) - value(At{i})) * B_sdpvar * C_temp))

    T_A(:, :, i) = value(At{i});
    T_B(:, :, i) = value(B_sdpvar);
    T_C(:, :, i) = mata;
    % T_D = ;

    % T_B(1 : ni(1), 7 : 8, i)                  = zeros(5, 2);
    % T_B(ni(1) + 1 : 2 * ni(1), 5 : 6, i)      = zeros(5, 2);
    % T_B(2 * ni(1) + 1 : 3 * ni(1), 3 : 4, i)  = zeros(5, 2);
    % T_B(3 * ni(1) + 1 : 4 * ni(1), 1 : 2, i)  = zeros(5, 2);

end

save('T_LMI.mat', 'T_A', 'T_B', 'T_C');