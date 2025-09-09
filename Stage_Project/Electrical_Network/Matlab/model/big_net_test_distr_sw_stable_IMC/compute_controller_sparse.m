clc, clear, close all;

%% Code based on the following article:

%    Sequential convex relaxation for convex optimization with bilinear
%    matrix equalities

%%

% load('ss_matrices_A_conns.mat');
load('ss_matrices_A_conns_acc.mat')


%%

% Depending on this constant, a certain structure is imposed for the base
% switching stabilizing controller. 
% super_sparse = 1 - Local controller i only uses locally estimated state
%                       to compute new estimate of state and control input
%                       (along with adjacent sensor readings)
%
% super_sparse = 0 - Local controller i uses state estimates from adjacent
%                       nodes in computing new local state estimate and
%                       control inputs
super_sparse = 0;

n_of_systems = size(A_matrices, 3);

n = size(A_matrices, 1);
m = size(B, 2);
p = size(C, 1);

%% Declare variables - 

% Xj - X matrix used in relaxing the bilinear constraint JP = U
% Yj - Y matrix used in relaxing the bilinear constraint JP = U

% Xl - X matrix used in relaxing the bilinear constraint QL = Y
% Yl - Y matrix used in relaxing the bilinear constraint QL = Y

% P - common Lyapunov function for the state feedback synthesis
% Q - common Lyapunov function for the state estimator synthesis

% J - state feedbacks
% Y - estimator feedbacks

% U - U = JP
% Y - Y = QL

Xj = repmat(eye(m, n), 1, 1, n_of_systems);
Yj = repmat(eye(n, n), 1, 1, n_of_systems);

Xl = repmat(eye(n, n), 1, 1, n_of_systems);
Yl = repmat(eye(n, p), 1, 1, n_of_systems);

P = sdpvar(n, n);
Q = sdpvar(n, n);

J = sdpvar(repmat(m, 1, n_of_systems), repmat(n, 1, n_of_systems));
L = sdpvar(repmat(n, 1, n_of_systems), repmat(p, 1, n_of_systems));

U = sdpvar(repmat(m, 1, n_of_systems), repmat(n, 1, n_of_systems));
Y = sdpvar(repmat(n, 1, n_of_systems), repmat(p, 1, n_of_systems));

gamma1 = 1;
gamma2 = 1;

% gamma1 = 1.4;
% gamma2 = 2;

%%

Nmax = 10;

for k = 1 : Nmax

    % By imposing this lower bound for P and Q we avoid cases where P and Q
    % are very close to 0
    constr = [P >= 5e-1 * eye(n); Q >= 1e-1 * eye(n)];
    % constr = [P >= 0; Q >= 0];
    
    obj = 0;
    
    % For each failure case
    for i = 1 : n_of_systems
    
        % Common Lyapunov function for state feedback and state estimator
        constr = [constr; 

            [P, gamma1 * (A_matrices(:, :, i) * P + B * U{i})';
            gamma1 * A_matrices(:, :, i) * P + gamma1 * B * U{i}, P] >= 0; 
            
            [Q, gamma2 * (Q * A_matrices(:, :, i) + Y{i} * C)';
            gamma2 * Q * A_matrices(:, :, i) + gamma2 * Y{i} * C, Q] >= 0

        ];
    
        % Compose the M matrices used in the convex relaxation 
        % Mj - used to relax JP = U
        % Ml - used to relax QL = Y
        Mj = [U{i} + Xj(:, :, i) * Yj(:, :, i) + J{i} * Yj(:, :, i) + Xj(:, :, i) * P, J{i} + Xj(:, :, i);
            P + Yj(:, :, i), eye(n)];
    
        Ml = [Y{i} + Xl(:, :, i) * Yl(:, :, i) + Q * Yl(:, :, i) + Xl(:, :, i) * L{i}, Q + Xl(:, :, i);
            L{i} + Yl(:, :, i), eye(n)];

        A_ctrl_i = A_matrices(:, :, i) + L{i} * C + B * J{i};

        % Impose structural constraints for each failure case based on the
        % topology of the network (the fixed one)
        constr = [constr;

            J{i}(1, 10 : 12) == zeros(1, 3);
            J{i}(2, 7 : 9) == zeros(1, 3);
            J{i}(3, 4 : 6) == zeros(1, 3);
            J{i}(4, 1 : 3) == zeros(1, 3);

            L{i}(1 : 3, 7 : 8) == zeros(3, 2);
            L{i}(4 : 6, 5 : 6) == zeros(3, 2);
            L{i}(7 : 9, 3 : 4) == zeros(3, 2);
            L{i}(10 : 12, 1 : 2) == zeros(3, 2);
            
        ];

        if super_sparse == 1
            % "Slow" behavior - Isolates the estimation dynamics for each node
            % Additionally, if super sparse structure is desired, control action
            % should not depend on adjacent estimated states
            constr = [constr;

                % Isolate estimator dynamics
                A_ctrl_i(1 : 3, 4 : end) == zeros(3, 12);

                A_ctrl_i(4 : 6, 1 : 3) == zeros(3, 3);
                A_ctrl_i(4 : 6, 7 : end) == zeros(3, 9);

                A_ctrl_i(7 : 9, 1 : 6) == zeros(3, 6);
                A_ctrl_i(7 : 9, 10 : end) == zeros(3, 6);

                A_ctrl_i(10 : 12, 1 : 9) == zeros(3, 9);
                A_ctrl_i(10 : 12, 13 : end) == zeros(3, 3);

                A_ctrl_i(13 : end, 1 : 12) == zeros(3, 12);

            ];
        else
            % "Fast" behavior - Only cuts down the connections which are not part of the topology
            constr = [constr;

                A_ctrl_i(1 : 3, 10 : 12) == zeros(3, 3);
                A_ctrl_i(4 : 6, 7 : 9) == zeros(3, 3);
                A_ctrl_i(7 : 9, 4 : 6) == zeros(3, 3);
                A_ctrl_i(10 : 12, 1 : 3) == zeros(3, 3);

            ];
        end
        
        obj = obj + norm(Mj, 'nuclear') + norm(Ml, 'nuclear');
    
    end
    
    % Impose structural constraints for each failure case based on the
    % connections which are being cut down
    constr = [constr;

        % Case 1 - line_12 and line_21 down
        J{1}(1, 4 : 6) == zeros(1, 3);
        L{1}(1 : 3, 3 : 4) == zeros(3, 2);

        J{1}(2, 1 : 3) == zeros(1, 3);
        L{1}(4 : 6, 1 : 2) == zeros(3, 2);

        % Case 2 - line_24 and line_42 down
        J{2}(2, 10 : 12) == zeros(1, 3);
        L{2}(4 : 6, 7 : 8) == zeros(3, 2);

        J{2}(4, 4 : 6) == zeros(1, 3);
        L{2}(10 : 12, 3 : 4) == zeros(3, 2);

        % Case 3 - line_34 and line_43 down
        J{3}(3, 10 : 12) == zeros(1, 3);
        L{3}(7 : 9, 7 : 8) == zeros(3, 2);

        J{3}(4, 7 : 9) == zeros(1, 3);
        L{3}(10 : 12, 5 : 6) == zeros(3, 2);

        % Case 4 - line_13 and line_31 down
        J{4}(1, 7 : 9) == zeros(1, 3);
        L{4}(1 : 3, 5 : 6) == zeros(3, 2);

        J{4}(3, 1 : 3) == zeros(1, 3);
        L{4}(7 : 9, 1 : 2) == zeros(3, 2);

        % Case 5 - line_15 and line_51 down
        J{5}(1, 13 : 15) == zeros(1, 3);
        L{5}(1 : 3, 9 : 10) == zeros(3, 2);

        J{5}(5, 1 : 3) == zeros(1, 3);
        L{5}(13 : 15, 1 : 2) == zeros(3, 2);

        % Case 6 - line_25 and line_52 down
        J{6}(2, 13 : 15) == zeros(1, 3);
        L{6}(4 : 6, 9 : 10) == zeros(3, 2);

        J{6}(5, 4 : 6) == zeros(1, 3);
        L{6}(13 : 15, 3 : 4) == zeros(3, 2);

        % Case 7 - line_45 and line_54 down
        J{7}(4, 13 : 15) == zeros(1, 3);
        L{7}(10 : 12, 9 : 10) == zeros(3, 2);

        J{7}(5, 10 : 12) == zeros(1, 3);
        L{7}(13 : 15, 7 : 8) == zeros(3, 2);

        % Case 8 - line_35 and line_53 down
        J{8}(3, 13 : 15) == zeros(1, 3);
        L{8}(7 : 9, 9 : 10) == zeros(3, 2);

        J{8}(5, 7 : 9) == zeros(1, 3);
        L{8}(13 : 15, 5 : 6) == zeros(3, 2);

    ];

    % If super sparse structure is not requested ("Fast" behavior case), some 
    % additional constraints have to be imposed to cut out
    % the remaining connections based on the failure case
    if super_sparse == 0
        A_ctrl_1 = A_matrices(:, :, 1) + L{1} * C + B * J{1};
        A_ctrl_2 = A_matrices(:, :, 2) + L{2} * C + B * J{2};
        A_ctrl_3 = A_matrices(:, :, 3) + L{3} * C + B * J{3};
        A_ctrl_4 = A_matrices(:, :, 4) + L{4} * C + B * J{4};
        A_ctrl_5 = A_matrices(:, :, 5) + L{5} * C + B * J{5};
        A_ctrl_6 = A_matrices(:, :, 6) + L{6} * C + B * J{6};
        A_ctrl_7 = A_matrices(:, :, 7) + L{7} * C + B * J{7};
        A_ctrl_8 = A_matrices(:, :, 8) + L{8} * C + B * J{8};
    
        constr = [constr;
    
            % Case 1 - line_12 and line_21 down
            A_ctrl_1(1 : 3, 4 : 6) == zeros(3, 3);
            A_ctrl_1(4 : 6, 1 : 3) == zeros(3, 3);
    
            % Case 2 - line_24 and line_42 down
            A_ctrl_2(4 : 6, 10 : 12) == zeros(3, 3);
            A_ctrl_2(10 : 12, 4 : 6) == zeros(3, 3);
    
            % Case 3 - line_34 and line_43 down
            A_ctrl_3(7 : 9, 10 : 12) == zeros(3, 3);
            A_ctrl_3(10 : 12, 7 : 9) == zeros(3, 3);
    
            % Case 4 - line_13 and line_31 down
            A_ctrl_4(1 : 3, 7 : 9) == zeros(3, 3);
            A_ctrl_4(7 : 9, 1 : 3) == zeros(3, 3);
    
            % Case 5 - line_15 and line_51 down
            A_ctrl_5(1 : 3, 13 : 15) == zeros(3, 3);
            A_ctrl_5(13 : 15, 1 : 3) == zeros(3, 3);
    
            % Case 6 - line_25 and line_52 down
            A_ctrl_6(4 : 6, 13 : 15) == zeros(3, 3);
            A_ctrl_6(13 : 15, 4 : 6) == zeros(3, 3);
    
            % Case 7 - line_45 and line_54 down
            A_ctrl_7(10 : 12, 13 : 15) == zeros(3, 3);
            A_ctrl_7(13 : 15, 10 : 12) == zeros(3, 3);
    
            % Case 8 - line_35 and line_53 down
            A_ctrl_8(7 : 9, 13 : 15) == zeros(3, 3);
            A_ctrl_8(13 : 15, 7 : 9) == zeros(3, 3);
    
        ];
    end

    optimize(constr, obj)
    
    for i = 1 : n_of_systems

        Xj(:, :, i) = -value(J{i});
        Yj(:, :, i) = -value(P);

        Xl(:, :, i) = -value(Q);
        Yl(:, :, i) = -value(L{i});

    end

end

%%

A_ctrl = zeros(n, n, n_of_systems);
J_mat = zeros(m, n, n_of_systems);
L_mat = zeros(n, p, n_of_systems);

for i = 1 : n_of_systems

    J_mat(:, :, i) = value(J{i});
    L_mat(:, :, i) = value(L{i});

    A_ctrl(:, :, i) = A_matrices(:, :, i) + L_mat(:, :, i) * C + B * J_mat(:, :, i);

end

B_ctrl = -L_mat;
C_ctrl = J_mat;
D_ctrl = zeros(m, p, n_of_systems);

n_states = 5;

T_A = zeros(n_states, n_states, n_of_systems);
T_B = zeros(n_states, 10, n_of_systems);
T_C = zeros(5, n_states, n_of_systems);
T_D = zeros(5, 10, n_of_systems);

if super_sparse == 1
    save('distributed_sw_stable_ctrl_sparse_slow.mat', 'A_ctrl', 'B_ctrl', 'C_ctrl', 'D_ctrl', 'T_A', 'T_B', 'T_C', 'T_D');
else
    save('distributed_sw_stable_ctrl_no_sparse_fast.mat', 'A_ctrl', 'B_ctrl', 'C_ctrl', 'D_ctrl', 'T_A', 'T_B', 'T_C', 'T_D');
end

