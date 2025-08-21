clc, clear, close all;

% load("ss_matrices_A_conns_acc.mat")
load('ss_matrices_A_conns.mat')
load('T_params.mat');
load('distributed_sw_stable_ctrl_no_sparse_fast.mat')

n_of_systems = size(A_matrices, 3);

n = size(A_matrices, 1);
m = size(B, 2);
p = size(C, 1);

z = tf('z', 1);
H = 1 / (z - 1) * eye(5);

n_states = 30;

T_A = zeros(5 * n_states, 5 * n_states, n_of_systems);
T_B = zeros(5 * n_states, 10, n_of_systems);
T_C = zeros(5, 5 * n_states, n_of_systems);
T_D = zeros(5, 10, n_of_systems);

%%

biggestmatever = cell(n_of_systems, 1);

for i = 1 : n_of_systems

    disp(['Failure case ', num2str(i), '###############################'])

    Pi = ss(A_matrices(:, :, i), B, C, 0, 1);

    Aki = A_ctrl(:, :, i);
    Li = -B_ctrl(:, :, i);
    Ji = C_ctrl(:, :, i);

    K11 = ss(Aki, -Li, Ji, 0, 1);
    K12 = ss(Aki, Bt, Ji, eye(m), 1);
    K21 = ss(Aki, -Li, Ct, eye(p), 1);
    K22 = ss(Aki, Bt, Ct, 0, 1);

    Tyw = feedback(Pi, series(K11, H), 1);
    Tyv = series(K12, feedback(series(H, Pi), K11, 1));
    Tew = series(Tyw, K21);
    Tev = series(Tyv, K21) + K22;

    Tyw = balreal(ss(Tyw, 'min'));
    Tyv = balreal(ss(Tyv, 'min'));
    Tew = balreal(ss(Tew, 'min'));
    Tev = balreal(ss(Tev, 'min'));

    Pgj = cell(5, 1);
    Tj = cell(5, 1);
    gammaj = zeros(5, 1);

    bigss = ss();
    bigmat = [];


    %% BLABLABLA


        Pg_temp = [Tyw, Tyv;
            Tew, Tev];
        Pg_temp = ss(Pg_temp, 'min');
        [K_temp, CLL, gamma_temp, inff] = hinfsyn(Pg_temp, 10, 5, [4 10]);
        K_temp_min = ss(K_temp, 'min');

    %%

    for j = 1 : 5

        Pgj{j} = [Tyw(2 * j - 1 : 2 * j, :), Tyv(2 * j - 1 : 2 * j, j);
            Tew(2 * j - 1 : 2 * j, :), Tev(2 * j - 1 : 2 * j, j)];

        Pgj{j} = ss(Pgj{j}, 'min');
        Pgj{j} = balreal(Pgj{j});

        [K_temp, Cll, gamma_temp, inff] = hinfsyn(Pgj{j}, 2, 1);

        if ~(all(abs(eig(K_temp)) < 1))
            error(['Controller of node ', num2str(j), ' is not stable'])
        end
        
        Tj{j} = K_temp;
        % Tj{j} = ss(K_temp, 'min');
        gammaj(j) = gamma_temp;

        bigmat = blkdiag(bigmat, Tj{j}.A);
        bigss = blkdiag(bigss, Tj{j});

        disp(['Node ', num2str(j), ', gamma = ', num2str(gamma_temp)])
        size(Tj{j})
    end

    biggestmatever{i} = bigmat;

    Pi = dlyap(bigss.A', eye(150));
    sys_new = ss2ss(bigss, sqrtm(Pi));

    for j = 1 : 5

        Tj_temp = ss(sys_new.A(n_states * (j - 1) + 1 : n_states * j, n_states * (j - 1) + 1: n_states * j), ...
            sys_new.B(n_states * (j - 1) + 1 : n_states * j, 2 * j - 1 : 2 * j), ...
            sys_new.C(j, n_states * (j - 1) + 1: n_states * j), ...
            sys_new.D(j, 2 * j - 1 : 2 * j), 1);

        % gapmetric(Tj{j}, Tj_temp)
        % norm(Tj{j} - Tj_temp)
        % gammaj{j} - norm(lft())
    end

    T_A(:, :, i) = sys_new.A;
    T_B(:, :, i) = sys_new.B;
    T_C(:, :, i) = sys_new.C;
    T_D(:, :, i) = sys_new.D;

    % norm(lft(Pg_temp, sys_new))
    % norm(lft(Pg_temp, 0 * sys_new))

    b = 1;

end

save('T_hinf.mat', 'T_A', 'T_B', 'T_C', 'T_D');

%%

% load('biggestmatever.mat')
% 
% n_states = 30;
% 
% P = sdpvar(n_states, n_states);
% Q = sdpvar(n_states, n_states);
% 
% constr = [P >= 1e-3 * eye(n_states);
%     Q >= 0];
% 
% for i = 1 : 2
% 
%     for j = 1 : 1
% 
%         constr = [constr;
%             biggestmatever{i}(n_states * (j - 1) + 1 : n_states * j, n_states * (j - 1) + 1: n_states * j)' * P * biggestmatever{i}(n_states * (j - 1) + 1 : n_states * j, n_states * (j - 1) + 1: n_states * j) - P + Q <= 0];
% 
%     end
% 
% end
% 
% optimize(constr, 0)

%% WORKS 