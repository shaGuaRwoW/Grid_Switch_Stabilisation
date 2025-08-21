%% Select physical parameters and network topology

d_t = 0.2;

k_ij = rand(5)*0.5 + ones(5)*0.5;
k_ij = (k_ij + k_ij')/2;
d_i = rand(1,5)*0.5 + ones(1,5)*0.5;
inv_m_i = rand(1,5)*2;

%%

comm_graph = [0 1 1 0 1;...
              1 0 0 1 1;...
              1 0 0 1 1;...
              0 1 1 0 1;...
              1 1 1 1 0];
          
k_i = sum((comm_graph .* k_ij)')';
   
%% Form the state-space representation of the network

A = [];
A_des = [];
B = [];
A_dp = [];

for i = 1 : 5
    A_r = [];
    for j = 1 : 5
        if i==j
            A_ii = [1 d_t;...
                    -k_i(i)*inv_m_i(i)*d_t 1-d_i(i)*inv_m_i(i)*d_t];
            A_dp = blkdiag(A_dp,A_ii);
            A_r = [A_r A_ii];
            A_des = blkdiag(A_des, A_ii);
        else
            A_ij = [0 0;...
                    k_ij(i,j)*inv_m_i(i)*d_t 0]*comm_graph(i,j);
            A_r = [A_r A_ij];
        end
    end
    A = [A; A_r];
    B=blkdiag(B,[0;1]);
end

C=eye(10);
D=zeros(10,5);

abs(eig(A))
max(abs(eig(A_dp)))