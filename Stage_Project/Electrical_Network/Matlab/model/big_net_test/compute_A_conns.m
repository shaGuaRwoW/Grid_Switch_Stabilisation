function [A] = compute_A_conns(faulty_lines)

    % Function receives a matrix with all the lines which are down as
    % follows: [i, j; k, z; ...] and return the state-space matrix A

    hi = [0.3244
        1.5886;
        0.6224;
        1.0571;
        0.3313];

    di = [0.5650;
        0.7844;
        0.7347;
        0.5060;
        0.6686];
    
    l_mat = [0 0.6511 0.7915 0 0.7513;
             0.6511 0 0 0.8433 0.6318;
             0.7915 0 0 0.8259 0.8620;
             0 0.8433 0.8259 0 0.8367;
             0.7513 0.6318 0.8620 0.8367 0];
    
    dt = 0.2;

    l_mat_temp = l_mat;

    if ~isempty(faulty_lines)
        n = size(faulty_lines, 1);
        for i = 1 : n
            l_mat_temp(faulty_lines(i, 1), faulty_lines(i, 2)) = 0;
        end
    end

    Aii = @(i) [1, dt;
        -hi(i) * dt * sum(l_mat_temp(i, :)), 1 - hi(i) * di(i) * dt];

    Aij = @(i, j) [0, 0;
                hi(i) * l_mat_temp(i, j) * dt, 0];

    A = [Aii(1), Aij(1, 2), Aij(1, 3), Aij(1, 4), Aij(1, 5);
        Aij(2, 1), Aii(2), Aij(2, 3), Aij(2, 4), Aij(2, 5);
        Aij(3, 1), Aij(3, 2), Aii(3), Aij(3, 4), Aij(3, 5);
        Aij(4, 1), Aij(4, 2), Aij(4, 3), Aii(4), Aij(4, 5);
        Aij(5, 1), Aij(5, 2), Aij(5, 3), Aij(5, 4), Aii(5)];

end

