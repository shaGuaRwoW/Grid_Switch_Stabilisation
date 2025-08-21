clc, clear, close all;

hi = [0.3244
    1.5886];

di = [0.5650;
    0.7844];

l_mat = [0 0.6511;
    0.6511 0];

dt = 0.2;

node_1_0 = [0; 0];
node_2_0 = [0; 0];

%% Compose state-spaces

% Nominal state space

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

ss_nom = ss(A, B, C, 0, 1);

% Connection from 2 to 1 down

l_mat(1, 2) = 0;

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

ss_fault = ss(A, B, C, 0, 1);

%% Compute controllers

n = size(A, 1);
m = size(B, 2);
p = size(C, 1);

P = sdpvar(n, n);
Q = sdpvar(n, n);

U1 = sdpvar(m, n);
U2 = sdpvar(m, n);

Y1 = sdpvar(n, p);
Y2 = sdpvar(n, p);

lmi1 = [P, (ss_nom.A * P + ss_nom.B * U1)';
    ss_nom.A * P + ss_nom.B * U1, P] >= 0; % LMI Ctrb Sys nom
lmi2 = [P, (ss_fault.A * P + ss_fault.B * U2)';
    ss_fault.A * P + ss_fault.B * U2, P] >= 0; % LMI Ctrb Sys fault

lmi3 = [Q, (Q * ss_nom.A + Y1 * ss_nom.C)';
    Q * ss_nom.A + Y1 * ss_nom.C, Q] >= 0; % LMI Obsv Sys nom
lmi4 = [Q, (Q * ss_fault.A + Y2 * ss_fault.C)';
    Q * ss_fault.A + Y2 * ss_fault.C, Q] >= 0; % LMI Obsv Sys fault

final_constr = [lmi1, lmi2, lmi3, lmi4, P >= 0, Q >= 0];
objec = 0;

optimize(final_constr, objec);

P = value(P);
Q = value(Q);

U1 = value(U1);
U2 = value(U2);

Y1 = value(Y1);
Y2 = value(Y2);

L1 = Q \ Y1;
L2 = Q \ Y2;

J1 = U1 / P;
J2 = U2 / P;

ss_ctr_nom = ss(ss_nom.A + L1 * ss_nom.C + ss_nom.B * J1, -L1, J1, 0, 1);
ss_ctr_fault = ss(ss_fault.A + L2 * ss_fault.C + ss_fault.B * J2, -L2, J2, 0, 1);

% save("controller.mat", 'P', 'Q', 'U1', 'U2', 'Y1', 'Y2');
save("state_spaces.mat", 'ss_nom', 'ss_fault', 'ss_ctr_nom', 'ss_ctr_fault')