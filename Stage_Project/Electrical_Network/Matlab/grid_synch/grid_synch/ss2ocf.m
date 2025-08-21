function [sys_o] = ss2ocf(sys)
    sys=ss(sys);
    [p,m]=size(sys.d);
    if p>1
        error('This function has been implemented only for MISO systems.')
    end
    tf_sys = tf(sys-sys.d);
    num_tf = poly(eig(sys.a));
    n_o=length(num_tf)-1;
    if n_o<1
        sys_o=ss(sys.d);
        return
    end
    A_o = zeros(n_o);
    A_o(:,1) = -(num_tf(2:end))';
    A_o(1:end-1,2:end)=eye(n_o-1);
    B_o=zeros(n_o,m);
    for i=1:m
        if norm(sys(i),inf)<1e-9
            temp = zeros(1,n_o+1);
        else
            temp = tf_sys(i).num{1};
        end
        B_o(:,i)=temp(2:end)';
        clear temp;
    end
    C_o=[1 zeros(1,n_o-1)];
    sys_o=ss(A_o,B_o,C_o,sys.d,sys.Ts);
end

