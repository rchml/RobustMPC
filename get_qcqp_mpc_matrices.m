function [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c] ...
                    = get_qcqp_mpc_matrices( mpc_prob, x0, x_ext0, beta )
                
    Q_o     = mpc_prob.Q_o;
    c_o     = mpc_prob.c_o(x0);
    alpha_o = mpc_prob.alpha_o(x0); 
    A       = mpc_prob.A;
    b       = mpc_prob.b(x0);  
    lsense  = mpc_prob.lsense; 
    lb      = mpc_prob.lb; 
    ub      = mpc_prob.ub;
    
    Qc      = mpc_prob.Qc;
    q       = cell(mpc_prob.N+1,1);
    beta_c  = cell(mpc_prob.N+1,1);
    
    for i = 1:mpc_prob.N+1
        q(i)        = {mpc_prob.q{i}(x0)};
        beta_c(i)   = {mpc_prob.beta_c{i}(beta, x_ext0, x0)};
    end 
        
end