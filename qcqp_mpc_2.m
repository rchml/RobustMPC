function [Q_o, c_o, alpha_o, A_, b_, lsense, lb, ub, Qc, q, beta_c]...
            = qcqp_mpc_2 ( sys, x0, x_ext0, Q, P, R, N, gamma, beta, quadratic_con, reset)
           
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    
    n = size(A,2);
    p = size(B,2);
    
    tmp = cell(N,1);
    
    % State equality constraints
    tmp(1:end)  = {B};
    B_diag      = blkdiag(tmp{:});
    tmp(1:end)  = {A};
    A_blk       = [blkdiag(tmp{:}), zeros(n*N, n)];
    I_blk       = [zeros(n*N, n), eye(n*N)];
    A_x         = [B_diag, A_blk - I_blk, zeros(n*N, p*N)];
    b_x         = zeros(n*N, 1);
    
    % Initial state constraint
    A_x0        = [zeros(n, p*N), eye(n), zeros(n, N*n), zeros(n, p*N)];
    b_x0        = x0;
    
    % Final state constraint
    A_xf        = [zeros(n, p*N), zeros(n, N*n), eye(n), zeros(n, p*N)];
    b_xf        = zeros(n, 1);
    
    % Output equality constraints
    tmp(1:end)  = {D};
    D_diag      = blkdiag(tmp{:});
    tmp(1:end)  = {C};
    C_blk       = [blkdiag(tmp{:}), zeros(p*N, n)];
    A_y         = [D_diag, C_blk, -eye(p*N)]; 
    b_y         = zeros(p*N, 1);
    
    % Linear equality constraints
    A_ = [A_x; A_x0; A_xf; A_y];
    b_ = [b_x', b_x0', b_xf', b_y'];
    lsense = '=';
    
    % Box constraints
    lb = -inf * ones(N*(2*p+n)+n, 1);
    ub =  inf * ones(N*(2*p+n)+n,1);
    
    % Quadratic inequality constraints
    if quadratic_con
        [Qc, q, beta_c] = qc_matrices(sys, N, gamma, beta, x0, x_ext0); 
    else
    
        Qc     = {};
        q      = {};
        beta_c = {};
    end
    
    
    % Objective function    
    
    % Creating Q_
    tmp(1:end) = {Q};
    Q_ = blkdiag(tmp{:}, P);
    
    % Creating R_
    tmp(1:end) = {R};
    R_ = blkdiag(tmp{:});
    
    Q_o     = blkdiag(R_, Q_, zeros(p*N));
    c_o     = zeros(1, N*p + (N+1)*n + N*p);
    alpha_o = 0;
        
end

function [Qc, q, beta_c] = qc_matrices(sys, N, gamma, beta, x0, x_ext0)  
    
    Qc      = cell(N+1,1);
    q       = cell(N+1,1);
    beta_c  = cell(N+1,1);
    
    % Quadratic sopness constraints                     %MOD z0 <= beta
    for i = 0:N
        [Qz_sym, qz, beta_] = sopness_constraint(sys, N, i, beta, gamma, x0, x_ext0);
        Qc(i+1)        = {Qz_sym};
        q(i+1)         = {qz};  
        beta_c(i+1)    = {beta_};
    end
    
end

function [Qz_sym, qz, beta_] = sopness_constraint( sys, N, K, beta, gamma, x0, x_ext0 )

    n = size(sys.A,2);
    p = size(sys.B,2);
    
    % z0 <= beta
    if K == 0
        Qz_sym  = zeros(N*p + (N+1)*n + N*p);
        qz      = zeros(N*p + (N+1)*n + N*p, 1);
        beta_   = beta - x_ext0;
        return
    end
    
    % Creating Qz_sym
    I_blk   = blkdiag(eye(p*K), zeros(p*(N-K)));
    Qz      = [gamma*I_blk, zeros(p*N, (N+1)*n + N*p);
               zeros((N+1)*n, N*p + (N+1)*n + N*p);
               I_blk, zeros(p*N, (N+1)*n + N*p)];
    Qz_sym  = 0.5*(Qz + Qz'); 
    
    % Creating qz
    qz = zeros(N*p + (N+1)*n + N*p, 1);
    
    % Creating beta_
    beta_ = beta - x_ext0;
    
end



