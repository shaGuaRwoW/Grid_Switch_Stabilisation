clc, clear, close all;

% h = [0.3244; 1.5886];
% d = [0.5650; 0.7844];
% li1 = [0; 0.6511];
% li2 = [0.6511; 0];

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

% sim('net_model_2_nodes.slx')

%% Testing if configuration is OK

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