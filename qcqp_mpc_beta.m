function [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
            = qcqp_mpc_beta( mpc_prob, x0, x_ext0, beta_lb )
                                                                
    beta            = 0;
%     reset           = 0;
%     quadratic_con   = 1;
    
    [Q_o_tmp, c_o_tmp, alpha_o, A_tmp, b, lsense, lb_tmp, ub_tmp, Qc_tmp, q_tmp, beta_c_tmp]...
            = get_qcqp_mpc_matrices( mpc_prob, x0, x_ext0, beta );
        
    % Defining quadratic constraints
    [Qc, q, beta_c] = qc_matrices_mod(Qc_tmp, q_tmp, beta_c_tmp);
    
%     % Defining Q_o, c_o and alpha_o
%     Q_o = blkdiag(0, Q_o_tmp);                  % Mod for the beta problem
%     
%     c_o = [0, c_o_tmp];                         % Mod for the beta problem
    
    % Defining Q_o, c_o and alpha_o
    Q_o     = blkdiag(0, zeros(mpc_prob.N));  %  blkdiag(0, Q_o_tmp); %        % Mod for the beta problem
    
    c_o     = [1, zeros(1, mpc_prob.N)]; % [100, c_o_tmp]; % zeros(1,N+1); %                          % Mod for the beta problem
    
%     alpha_o = 0;
    
    % Linear constraints:  xN = 0
    A      = [zeros(size(A_tmp,1),1), A_tmp];   % Mod for the beta problem
    
    % Box constraints
    lb = [beta_lb; lb_tmp];
    ub = [inf; ub_tmp];

end

function [Qc, q, beta_c] = qc_matrices_mod(Qc_tmp, q_tmp, beta_c_tmp)  
    
    N       = size(Qc_tmp,1);

    Qc      = cell(N,1);
    q       = cell(N,1);
    beta_c  = cell(N,1);
    
    % Quadratic sopness constraints
    for i = 1:N
          
        Qz_sym       = blkdiag(0, Qc_tmp{i});       % Mod for the beta problem
        qz           = [-1; q_tmp{i}];              % Mod for the beta problem
        beta_        = beta_c_tmp{i};               % Mod for the beta problem
        
        Qc(i)        = {Qz_sym};
        q(i)         = {qz};  
        beta_c(i)    = {beta_};
    end
    
end
