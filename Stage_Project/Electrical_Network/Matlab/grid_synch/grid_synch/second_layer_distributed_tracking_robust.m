function [uout] = second_layer_distributed_tracking_robust(currentx,currentr,currentu,ar_id,t)

if exist('Controller', 'var') < 1
    persistent Controller
end

if t == 0

    A = evalin('base','A');
    B = evalin('base','B');
    d_t = evalin('base','d_t');
    des_com_graph = evalin('base','des_com_graph');

    Q_cost = evalin('base','Q_cost');
    R_cost = evalin('base','R_cost');

    U_lims = evalin('base','U_lims');
    
    % Avoid explosion of internally defined variables in YALMIP
    yalmip('clear')
    
    % Setup the optimization problem

    eps_trk = sdpvar(1,1);
    eps_cmd = sdpvar(1,1);
    u_MPC = sdpvar(1,1);
    omega_next = sdpvar(1,1);
    r = sdpvar(1,1);
    ci = sdpvar(length(A),1);
    u_NRF = sdpvar(1,1);
    
    scal_mat = eye(10);
    for i=1:5
        if des_com_graph(ar_id,i)<1
            scal_mat([2*i-1 2*i],[2*i-1 2*i]) = zeros(2);
        end
    end
    
    objective = Q_cost{ar_id}*eps_trk + R_cost{ar_id}*eps_cmd ;
    constraints = [ omega_next == A(2*ar_id,:)*scal_mat*ci + B(2*ar_id,ar_id)*(u_MPC+u_NRF);...
                    eps_trk >= 0;...
                    eps_cmd >= 0;...
                    eps_cmd <= U_lims(ar_id);...
                    (u_MPC+u_NRF) <= eps_cmd;...
                    -eps_cmd <= (u_MPC+u_NRF);...
                    (r - ci(2*ar_id-1) - d_t*(ci(2*ar_id)+omega_next)) <= eps_trk;...
                    -eps_trk <= (r - ci(2*ar_id-1) - d_t*(ci(2*ar_id)+omega_next))];
    
    
    % Define an optimizer object which solves the problem for a particular
    % initial state and reference
    % Force optimizer to turn on dispay until you know things work
    % To see log you have to use debug breaks in code and run manually
    ops = sdpsettings('verbose',2);
    Controller{ar_id} = optimizer(constraints,objective,ops,{ci,r,u_NRF},u_MPC);
    
    % And use it here 
    [uout,problem] = Controller{ar_id}(currentx,currentr,currentu);
    uout = 0 * uout;
    if problem
       % Fix!
       error('Error ecountered just after initialization.');
    end
    
else    
    % Almost no overhead
    [uout,problem] = Controller{ar_id}(currentx,currentr,currentu);
    uout = 0 * uout;
    if problem
      % Debug, analyze, fix!
      error('Error ecountered during standard iteraition.');
    end 
end

end

