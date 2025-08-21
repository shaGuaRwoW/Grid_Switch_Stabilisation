%% Cleanup

clear; clc; close all;
load('first_layer_done.mat');

%% Choose variable bounds

U_lims = [0.10 0.10 0.10 0.10 0.10];
X_lims = [0.25 0.25 0.25 0.25 0.25;...
          0.10 0.10 0.10 0.10 0.10];

%% Select amplitude for noise generators

amp_n_x = 1e-3;
amp_n_w = 1e-3;
amp_n_u_f = 1e-3;
amp_n_u_s1 = 1e-3;
amp_n_u_s2 = 1e-3;

%% Configure network initial state

amp_x0 = diag(reshape(X_lims,length(A),1));
x0 = 2 * rand(length(A), 1) - 1;
x0 = amp_x0 * x0;

node_1_0 = [20, -30]' / 30;
node_2_0 = [-25, 31]' / 30;
node_3_0 = [10, 30]' / 30;
node_4_0 = [-10, -30]' / 30;
node_5_0 = [5, -7]' / 30;

x0 = [node_1_0; node_2_0; node_3_0; node_4_0; node_5_0];

%% Choose reference vales for network state vector

r_vec_1 =  0.2 * [1 1 1 1 1]';
r_vec_2 = -0.2 * [1 1 1 1 1]';

%% Impose prediction horizon and cost function weights

Q_cost = cell(5,1);
R_cost = cell(5,1);
for i = 1:5
    Q_cost{i} = 1e6;
    R_cost{i} = 1e0;
end