
J = ss(A + B * F + L * C + L * D * F, [-L, B + L * D], [F; -(C + D * F)], ...
    [zeros(5, 10), eye(5); eye(10), -D], d_t);

K_final = lft(J, Q);
K_final = ss(K_final, 'min');

