clc, clear, close all;

% Considered failure modes:
% 1 : line_12 and line_21 down
% 2 : line_24 and line_42 down
% 3 : line_34 and line_43 down
% 4 : line_13 and line_31 down
% 5 : line_15 and line_51 down
% 6 : line_25 and line_52 down
% 7 : line_45 and line_54 down
% 8 : line_35 and line_53 down

% Connection order: MSB to LSB
% line_12, line_13, line_15, line_34, line_51, line_52

% Based on the binary number formed by the above-mentioned sequence, we add
% one and we obtain the index of the matrix in the big vector defined below

B = [0;
    1];
B = blkdiag(B, B, B, B, B);

C = [1, 0;
    0, 1];
C = blkdiag(C, C, C, C, C);

A_matrices = zeros(10, 10, 9);

%% Compose state-spaces

% Nominal net
A_matrices(:, :, 1) = compute_A_conns([1, 2; 2, 1]);        % Mode 1 : Conns down between nodes 1 - 2
A_matrices(:, :, 2) = compute_A_conns([2, 4; 4, 2]);        % Mode 2 : Conns down between nodes 2 - 4
A_matrices(:, :, 3) = compute_A_conns([3, 4; 4, 3]);        % Mode 3 : Conns down between nodes 3 - 4
A_matrices(:, :, 4) = compute_A_conns([1, 3; 3, 1]);        % Mode 4 : Conns down between nodes 1 - 3
A_matrices(:, :, 5) = compute_A_conns([1, 5; 5, 1]);        % Mode 5 : Conns down between nodes 1 - 5
A_matrices(:, :, 6) = compute_A_conns([2, 5; 5, 2]);        % Mode 6 : Conns down between nodes 2 - 5
A_matrices(:, :, 7) = compute_A_conns([4, 5; 5, 4]);        % Mode 7 : Conns down between nodes 4 - 5
A_matrices(:, :, 8) = compute_A_conns([3, 5; 5, 3]);        % Mode 8 : Conns down between nodes 3 - 5
A_matrices(:, :, 9) = compute_A_conns([]);                  % Mode 9 : All good :)

save('ss_matrices_A_conns', 'A_matrices', 'B', 'C');