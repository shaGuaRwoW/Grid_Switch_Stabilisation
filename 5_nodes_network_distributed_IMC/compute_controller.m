 clc, clear, close all;

%%

load('ss_matrices_A_conns.mat');
load('K_final.mat');

%%

n_of_systems = size(A_matrices, 3);

n = 10;
m = 5;
p = 10;

P = sdpvar(n, n);
U = sdpvar(repmat(m, 1, n_of_systems), repmat(n, 1, n_of_systems));
Y = sdpvar(repmat(n, 1, n_of_systems), repmat(p, 1, n_of_systems));
Q = sdpvar(n, n);

constr = [];

for i = 1 : n_of_systems

    constr = [constr; 
        [P, (A_matrices(:, :, i) * P + B * U{i})';
        A_matrices(:, :, i) * P + B * U{i}, P] >= 0; 
        
        [Q, (Q * A_matrices(:, :, i) + Y{i} * C)';
        Q * A_matrices(:, :, i) + Y{i} * C, Q] >= 0];

        % Q * A_matrices(:, :, i) + Y{i} * C == Q * A_db;
        % A_matrices(:, :, i) * P + B * U{i} == A_des * P];

end

optimize(constr)

P = value(P);
Q = value(Q);

%% Set up desired K - K_final

K_final.Ts = 1;

% Check if the obtained controller is stable 
if ~all(abs(eig(K_final.a)) < 1)

    error('Controller is not stable')

end

% Compute DCF for K desired
U_dcf_K = K_final;
V_dcf_K = eye(p);


%%

J = zeros(m, n, n_of_systems);
L = zeros(n, p, n_of_systems);

A_ctrl = zeros(n, n, n_of_systems);
D_ctrl = zeros(m, p, n_of_systems);

T = cell(n_of_systems, 1);

for i = 1 : n_of_systems

    % Check if desired controllers stabilizes the system

    rez_fb = feedback(ss(A_matrices(:, :, i), B, C, 0, 1), K_final, 1);
    if ~all(abs(eig(rez_fb.a)) < 1)
    
        error(['Desired controller and system ', num2str(i), ...
            ' are not stable in feedback configuration'])
    
    end

    J(:, :, i) = value(U{i}) / P;
    L(:, :, i) = Q \ value(Y{i});

    A_ctrl(:, :, i) = A_matrices(:, :, i) + L(:, :, i) * C + B * J(:, :, i);

    %% Compute T

    % LFT inverse way
    % Jg = ss(A_ctrl(:, :, i), [-L(:, :, i), B], [J(:, :, i); C], [zeros(m, p), eye(m); -eye(p), zeros(p, m)], 1);
    % 
    % T_temp = lft(K_final, inv(Jg));
    % T_temp = ss(T_temp, 'min');

    % DCF way

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
    
    Z = Mt * V_dcf_K - Nt * U_dcf_K;
    T_temp = M \ (U_dcf_K / Z - Uz);
    T_temp = -ss(T_temp, 'min');

    %%
    
    T{i} = T_temp;
    
    Kg = ss(A_ctrl(:, :, i), [-L(:, :, i), B], [J(:, :, i); C], [zeros(m, n), eye(m); -eye(n), zeros(n, m)], 1);
    Kd = lft(Kg, T{i});
    
    gapmetric(Kd, K_final)
end

B_ctrl = -L;
C_ctrl = J;

% save('controllers.mat', 'A_ctrl', 'B_ctrl', 'C_ctrl', 'D_ctrl', 'P', 'Q', 'J', 'L');

%% 

n_states = size(T{1}.A, 1);

prod_fin = 1;
for i = 1 : n_of_systems

    if size(T{i}.A, 1) ~= n_states

        error('The realizations do not have the same size')

    end

    prod_fin = prod_fin * all(abs(eig(T{i})) < 1);
    size(T{i})

    if ~all(abs(eig(T{i})) < 1)
    
        error(['Realization ', num2str(i), ' is not stable'])
    
    end

end


%%

T_A = zeros(n_states, n_states, n_of_systems);
T_B = zeros(n_states, 10, n_of_systems);
T_C = zeros(5, n_states, n_of_systems);
T_D = zeros(5, 10, n_of_systems);

for i = 1 : n_of_systems

    Pi = dlyap(T{i}.A', eye(n_states));
    
    sys_new = ss2ss(T{i}, sqrtm(Pi));
    % sys_new = ss2ss(T{i}, eye(n_states));
    
    T_A(:, :, i) = sys_new.A;
    T_B(:, :, i) = sys_new.B;
    T_C(:, :, i) = sys_new.C;
    T_D(:, :, i) = sys_new.D;
    
    norm(T_A(:, :, i), 'fro')
    
end

save('controllers.mat', 'A_ctrl', 'B_ctrl', 'C_ctrl', 'D_ctrl', 'T_A', 'T_B', 'T_C', 'T_D');