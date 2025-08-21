clc, clear, close all;

load('ss_matrices_A_conns.mat');

z = tf('z', 1);

H = 1 / (z - 1) * eye(5);

n_of_systems = size(A_matrices, 3);

A_new = zeros(15, 15, n_of_systems);
B_new = zeros(15, 5, n_of_systems);
C_new = zeros(10, 15, n_of_systems);

big_B = [eye(2); zeros(1, 2)];
small_B = [0; 0; 1];

T1 = blkdiag(big_B, big_B, big_B, big_B, big_B);
T2 = blkdiag(small_B, small_B, small_B, small_B, small_B);

T = [T1, T2];

for i = 1 : n_of_systems
    
    ss_temp = ss(A_matrices(:, :, i), B, C, 0, 1);
    ss_temp = series(H, ss_temp);

    ss_new = ss2ss(ss_temp, T);

    A_new(:, :, i) = ss_new.A;
    B_new(:, :, i) = ss_new.B;
    C_new(:, :, i) = ss_new.C;
    
end

%%

A_matrices = A_new;
B = B_new(:, :, 9);
C = C_new(:, :, 9);

save('ss_matrices_A_conns_acc.mat', 'A_matrices', 'B', 'C');