function mpc_prob = qcqp_mpc_data( sys, Q, P, R, N, gamma, quadratic_con, z_fun)
                                   
    fprintf('qcqp_mpc_data: Computing problem data and structure...\n');

    n = size(sys.A,2);
    p = size(sys.B,2);
    
    tmp = cell(N,1);
    
    % Creating Q_
    tmp(1:end-1) = {Q};
    tmp(end) = {P};
    Q_ = blkdiag(tmp{:});
    
    % Creating R_
    tmp(1:end) = {R};
    R_ = blkdiag(tmp{:});
    
    % Creating T_
    tmp(1:end) = {sys.A};
    for i = 2:N
        tmp(i) = {tmp{i}*tmp{i-1}};
    end
    % Transpose each matrix in order to setting up T_ as a 'tall' matrix
    for i = 1:N
        tmp(i) = {tmp{i}'};
    end
    T_ = [tmp{:}]';
    
    % Creating S_
    tmp(1:end) = {sys.B};
    dB = blkdiag(tmp{:});
    v = [eye(n); T_(1:end-n,:)];
    c = [zeros(n,N*n); eye(n*(N-1),N*n)];
    tmp(1) = {v};
    for i = 2:N
        tmp(i) = {c*tmp{i-1}};
    end
    S_ = [tmp{:}]*dB;

    % Defining quadratic constraints    
        
    if quadratic_con

        [Qc_, q_fun, beta_c_fun, new_gamma] = qc_matrices(sys, N, gamma, z_fun);

        % Calculate vectors depending on initial states  
        Qc      = Qc_; 
        q       = q_fun; % cell(N,1);
        beta_c  = cell(N+1,1);
        for i = 1:N+1
            beta_c(i) = {@(beta, x_ext0, x0) beta + beta_c_fun{i}(x_ext0, x0)};
        end

%         [Qc, q, beta_c] = qc_matrices_old(sys, N, gamma, beta, x0, x_ext0); 
    else
    
        Qc     = {};
        q      = {};
        beta_c = {};
    end
    
    % Defining Q_o, c_o and alpha_o
    Q_o     = R_ + (S_')*Q_*S_;
    c_o     = @(in) 2*(in')*(T_')*Q_*S_;
    alpha_o = @(in) (in')*(Q + (T_')*Q_*T_)*in;
    
    % Linear constraints:  xN = 0
    I = eye(N*n);
    A = I((N-1)*n + 1:end,:)*S_;
    b = @(in) -I((N-1)*n + 1:end,:)*T_*in;
    lsense = '=';
    
    % Box constraints
    lb = -inf * ones(p*N,1);
    ub =  inf * ones(p*N,1);

    mpc_prob.sys        = sys;
    mpc_prob.Q          = Q;
    mpc_prob.P          = P;
    mpc_prob.R          = R;
    mpc_prob.N          = N;
    mpc_prob.gamma      = new_gamma;
    mpc_prob.z_fun      = z_fun;
    
    mpc_prob.Q_o        = Q_o;
    mpc_prob.c_o        = c_o;
    mpc_prob.alpha_o    = alpha_o;
    mpc_prob.A          = A;
    mpc_prob.b          = b;
    mpc_prob.lsense     = lsense;
    mpc_prob.lb         = lb;
    mpc_prob.ub         = ub;
    mpc_prob.Qc         = Qc;
    mpc_prob.q          = q;
    mpc_prob.beta_c     = beta_c;

    % Save problem structure
    save('mpc_prob.mat', '-struct', 'mpc_prob');
end

function [Qc, q, beta_c] = qc_matrices_old(sys, N, gamma, beta, x0, x_ext0)  
    
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
    
%     % Quadratic sopness constraints
%     for i = 1:N
%         [Qz_sym, qz, beta_] = sopness_constraint(sys, N, i, beta, gamma, x0, x_ext0);
%         Qc(i)        = {Qz_sym};
%         q(i)         = {qz};  
%         beta_c(i)    = {beta_};
%     end
%     
end

function [Qz_sym, qz, beta_] = sopness_constraint( sys, N, K, beta, gamma, x0, x_ext0 )

    n = size(sys.A,2);
    p = size(sys.B,2);
    tmp = cell(K,1);
    
    % z0 <= beta
    if K == 0
        Qz_sym  = zeros(N*p);
        qz      = zeros(N*p, 1);
        beta_   = beta - x_ext0;
        return
    end
    
    % Creating T_
    tmp(1:end) = {sys.A};
    for i = 2:K
        tmp(i) = {tmp{i}*tmp{i-1}};
    end
    % Transpose each matrix in order to setting up T_ as a 'tall' matrix
    for i = 1:K
        tmp(i) = {tmp{i}'};
    end
    T_ = [tmp{:}]';
    
    % Creating S_
    tmp(1:end) = {sys.B};
    dB = blkdiag(tmp{:});
    v = [eye(n); T_(1:end-n,:)];
    c = [zeros(n,K*n); eye(n*(K-1),K*n)];
    tmp(1) = {v};
    for i = 2:K
        tmp(i) = {c*tmp{i-1}};
    end
    S_ = [tmp{:}]*dB;
    
    % Creating Qz_sym
    tmp(1:end) = {sys.C};
    dC = blkdiag(tmp{:});
    M = dC*S_;                              %%%%%
    Mz = [zeros(p, K*p); M(1:(K-1)*p, :)];
    tmp(1:end) = {gamma*eye(p) + sys.D};    %%%%%
    dGD = blkdiag(tmp{:});
    Qz = Mz + dGD;
    Qz_sym_tmp = 0.5*(Qz + Qz');  
    Qz_sym = [Qz_sym_tmp, zeros(size(Qz_sym_tmp, 1), p*(N-K)); zeros(p*(N-K), p*N)];
    
    % Creating qz
    qz_tmp = dC*v*x0;                        %%%%%
    qz = [qz_tmp; zeros(N-K,1)];
    
    % Creating beta_
    beta_ = beta - x_ext0;
    
end

