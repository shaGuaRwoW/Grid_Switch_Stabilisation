%% Preserve the individual node dynamics

load('rand_net_params.mat');

A_db = [];
for i=1:5
    A_db = blkdiag(A_db, [1 d_t;-1/d_t -1]);
end
L = A_db - A;
F = linsolve(B, A_des - A);

%% Form standard DCF over RH-inf of the network's model

r_cop = ss(A + B * F, [B -L], [F; C], eye(15), d_t);
l_cop = ss(A + L * C, [-B L], [F; C], eye(15), d_t);

%% Select component TFMs and apply transformation to facilitate sparsity

M = r_cop(1:5,1:5);
X = r_cop(1:5,6:15);
N = r_cop(6:15,1:5);
Y = r_cop(6:15,6:15);

Yt =  l_cop(1:5,1:5);
Xt = -l_cop(1:5,6:15);
Nt = -l_cop(6:15,1:5);
Mt =  l_cop(6:15,6:15);

%% Configure the desired model-matching problem

% A zero entry in the (i,j)'th position indicates that the j'th node cannot
% send information to the i'th one.
des_com_graph = [1 1 1 0 1;...
                 1 1 0 1 1;...
                 1 0 1 1 1;...
                 0 1 1 1 1;...
                 1 1 1 1 1;];
if (norm(des_com_graph - des_com_graph')>1e-6)
    error('The desired communication graph is not symmetrical!');
    % Recall that, in contrast to the platoon, we're assuming that the
    % communication pathways between the grid's nodes are bidirectional.
end

% We're modularizing the model-matching problem, so that we're splitting
% it according to the rows of the Q-parameter. The math involved is pretty
% gruesome, and the numerical part is *pure murder*, so we're gonna want to
% solve small problems, and then compose the solution out of partial answers.
M_vec=cell(1,5);
X_vec=cell(1,5);
for i=1:5
    for j=1:5
        if des_com_graph(i,j)<1
            M_vec{i} = ss([M_vec{i} Nt(:,j) Mt(:,[2*j-1 2*j])],'min');
            X_vec{i} = ss([X_vec{i} -Yt(i,j) -Xt(i,[2*j-1 2*j])],'min');
        end
    end
end

%% Solve model-matching problem (requires platform compatible with DSTools)

Q_sol_vec = cell(1,5);
N_vec = cell(1,5);
Q_vec = cell(1,5);
Q = tf(0,1) * ones(5,10);
info_sol = cell(1,5);
info_null = cell(1,5);

poles = 0; % We wanna have the poles as close as possible to the origin.
for i=1:5
    if isempty(X_vec{i}) || isempty(M_vec{i})
        Q_sol_vec{i} = tf(0,1) * ones(1,10);
        N_vec{i} = tf(0,1) * ones(1,10);
    else
        % Obtain a minimal degree solution
        [Q_sol_vec{i},info_sol{i}]=glsol(M_vec{i},X_vec{i},struct('mindeg',true));
        Q_sol_vec{i}=ss(tf(Q_sol_vec{i}),'min');
        if max(abs(eig(Q_sol_vec{i}))) >= 1 
            warning('Minimal degree solution is unstable! Resorting to regular solution.');
            [Q_sol_vec{i},info_sol{i}]=glsol(M_vec{i},X_vec{i},struct('poles',poles));
            Q_sol_vec{i}=ss(tf(Q_sol_vec{i}),'min');
             if max(abs(eig(Q_sol_vec{i}))) >= 1
                error('Regular solution is unstable!');
            end
        end
        if norm(ss(Q_sol_vec{i}*M_vec{i}-X_vec{i},'min'),inf) >= 1e-6
            error('Solution badly conditioned or inexistent!');
            disp(info_sol{i});
        end
        % Compute a basis for null space of the coefficient matrix, so we can
        % obtain the entire class of solutions for our model-matching (sub)problem.
        % Ideally, we'd like the basis to be simple (its degree is the sum of
        % the degrees of its rows), but that's much mode delicate (computationally).
        [N_vec{i},info_null{i}]=glnull(M_vec{i},struct('simple',false,'poles',poles));
        N_vec{i} = ss(N_vec{i}, 'min');
        if norm(ss(N_vec{i}*M_vec{i},'min'),inf) >= 1e-6
            error('Basis badly conditioned or inexistent!');
            disp(info_null{i});
        end
        if max(abs(eig(N_vec{i}))) >= 1
            error('Basis is unstable!');
        end
        % Use the bases to obtain solutions with zero feedthrough.
        Q_vec{i} = Q_sol_vec{i} + linsolve(N_vec{i}.d',-Q_sol_vec{i}.d')'*N_vec{i};
        Q(i,:) = Q_vec{i};
        Q = ss(tf(Q),'min');
    end
end

if norm(Q.d) >= 1e-6
    error('The Q-parameter is not strictly proper!');
else
    Q.d = 0 * Q.d; % The feedthrough is negligible, so we can just cancel
                   % it out, numerically.
end

%% Form the distributed controller

Yt_Q = ss(Yt + Q * Nt, 'min');
Xt_Q = ss(Xt + Q * Mt, 'min');
Yt_Q_diag = ss(Yt_Q(1,1),'min');
for i=2:5
    Yt_Q_diag = blkdiag(Yt_Q_diag,ss(Yt_Q(i,i),'min'));
end

if norm(Yt_Q_diag-eye(5),inf) < 1e-6
    Yt_Q_diag = eye(5);
end
Phi = eye(5) - Yt_Q_diag\Yt_Q;
Gamma = Yt_Q_diag\Xt_Q;

sparse_check_sum = 0;
for i=1:5
    sparse_check_sum = sparse_check_sum + norm(ss(Phi(i,i),'min'),inf);
    for j=1:5
        if des_com_graph(i,j)<1
            sparse_check_sum = sparse_check_sum + ...
            norm(ss(Phi(i,j),'min'),inf) + norm(ss(Gamma(i,[2*j-1 2*j]),'min'),inf);
        end
    end
end

if sparse_check_sum >= 1e6
    error('NRF does not possess asigned structure (check for bad conditioning)!');
end

K_D = ss([Phi Gamma], 'min');

%% Bring distributed subcontrollers to observable canonical form

K_D_sliced = cell(1,5);
decimal_places = 4; % We're gonna want to "brush away" all of the "numerical
                    % leftovers", in order to have a clean implementation.
for i=1:5
    K_D_sliced{i} = ss(K_D(i,:),'min');
    K_D_sliced{i}.a = round(K_D_sliced{i}.a,decimal_places);
    K_D_sliced{i}.b = round(K_D_sliced{i}.b,decimal_places);
    K_D_sliced{i} = ss2ocf(K_D_sliced{i});
    K_D_sliced{i}.a = round(K_D_sliced{i}.a,decimal_places);
    K_D_sliced{i}.b = round(K_D_sliced{i}.b,decimal_places);
end