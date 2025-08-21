clc, clear, close all;

% load('controllers.mat');
% load('ss_matrices_A_conns.mat');
% load('K_final.mat');

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

X = repmat(eye(sum(ni)), 1, 1, n_of_systems);
Y = repmat(eye(sum(ni)), 1, 1, n_of_systems);

C = sdpvar(repmat(sum(ni), 1, n_of_systems), repmat(sum(ni), 1, n_of_systems)); % C = AP

P = sdpvar(sum(ni), sum(ni));

%%

Nmax = 10;

for k = 1 : Nmax

    constr = [P >= 0.1 * eye(sum(ni))];

    obj = 0;

    for i = 1 : n_of_systems
    
        constr = [constr;

            [P, C{i}';
            C{i}, P] >= 0;

        ];
    
        M = [C{i} + X(:, :, i) * Y(:, :, i) + At{i} * Y(:, :, i) + X(:, :, i) * P, At{i} + X(:, :, i);
            P + Y(:, :, i), eye(sum(ni))];

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

        obj = obj + norm(M, 'nuclear');
    end

    optimize(constr, obj)
    
    for i = 1 : n_of_systems

        X(:, :, i) = -value(At{i});
        Y(:, :, i) = -value(P);
    end

end