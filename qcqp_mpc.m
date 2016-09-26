function [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
            = qcqp_mpc( sys, x0, x_ext0, Q, P, R, N, gamma, beta, quadratic_con, z_fun, reset )
                                   
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
    persistent Qc_ q_fun beta_c_fun
    
    if quadratic_con
        % Avoid to compute gradient and hessian at each iteration using reset
        if isempty(Qc_) || reset
            [Qc_, q_fun, beta_c_fun] = qc_matrices(sys, N, gamma, z_fun);
        end

        % Calculate vectors depending on initial states  
        Qc = Qc_; 
        q       = cell(N,1);
        beta_c  = cell(N,1);
        for i = 1:N
            q(i) = {q_fun{i}(x0)};
            beta_c(i) = {beta + beta_c_fun{i}(x_ext0, x0)};
        end

%         [Qc, q, beta_c] = qc_matrices_old(sys, N, gamma, beta, x0, x_ext0); 
    else
    
        Qc     = {};
        q      = {};
        beta_c = {};
    end
    
    % Defining Q_o, c_o and alpha_o
    Q_o = R_ + (S_')*Q_*S_;
    c_o = 2*(x0')*(T_')*Q_*S_;
    alpha_o = (x0')*(Q + (T_')*Q_*T_)*x0;
    
    % Linear constraints:  xN = 0
    I = eye(N*n);
    A = I((N-1)*n + 1:end,:)*S_;
    b = -I((N-1)*n + 1:end,:)*T_*x0;
    lsense = '=';
    
    % Box constraints
    lb = -inf * ones(p*N,1);
    ub =  inf * ones(p*N,1);

end

function [Qc, q, beta_c] = qc_matrices(sys, N, gamma, z_fun)  
    
    Qc      = cell(N,1);
    q       = cell(N,1);
    beta_c  = cell(N,1);
    
    % Quadratic sopness constraints
    for i = 1:N
        [Qz_sym, qz_fun, beta_fun] = sym_sopness_constraint(sys, N, i, gamma, z_fun);
        Qc(i)        = {Qz_sym};
        q(i)         = {qz_fun};  
        beta_c(i)    = {beta_fun};
    end
    
end


function [Qz, qz_fun, beta_fun] = sym_sopness_constraint( sys, N, K, gamma_, z_fun )

n = size(sys.B,1);
p = size(sys.B,2);

x0 = sym('x0', [n,1]);
u0 = sym('u0', [p,1]);
u  = sym('u', [p,N-1]);
assume(x0,'real')
assume(u0,'real');
assume(u,'real');

A  = sym('A', [n,n]);
B  = sym('B', [n,p]);
C  = sym('C', [p,n]);
D  = sym('D', [p,p]);

gamma = sym('gamma');
z0    = sym('z0');

% Choice of the extended state
if strcmp(z_fun, 'SOP')
    l = z_sop(K);
elseif strcmp(z_fun, 'L2')
    l = z_l2(K);
end

l = subs(l,A(:), sys.A(:));
l = subs(l,B(:), sys.B(:));
l = subs(l,C(:), sys.C(:));
l = subs(l,D(:), sys.D(:));
l = subs(l,gamma, gamma_);
%l = subs(l, x0, x0_);
u_t = u';
u_var = [u0(:); u_t(:)];

% Calculate Qz
h = 1/2*hessian(l, u_var);
Qz = double(h);

% Calculate qz
g = gradient(l, u_var);
g = subs(g, u_var, zeros(N,1));
qz_fun = matlabFunction(g, 'Vars', {x0});

% Calculate beta_
const = subs(l, u_var, zeros(N,1));
%const = subs(const, z0, x_ext0);
beta_fun = matlabFunction(-const, 'Vars', {z0, x0});
    

    function out = z_sop(k)

        if k == 1
            out = z0 + u0'*y(0) + gamma*(u0')*u0;                                                   %%%%%%%%%%%%
            return
        end
        out = z_sop(k-1) + transpose(u(:,k-1))*y(k-1) + gamma*transpose(u(:,k-1))*u(:,k-1);             %%%%%%%%%%%%
        
    end

    function out = z_l2(k)

        if k == 1
            out = z0 + u0'*u0 - gamma*(y(0)')*y(0);                                                   %%%%%%%%%%%%
            return
        end
        out = z_l2(k-1) + transpose(u(:,k-1))*u(:,k-1) - gamma*transpose(y(k-1))*y(k-1);             %%%%%%%%%%%%
        
    end

    function out = y(k)

        if k == 0
            out = C*x0 + D*u0;
            return
        end
        out = C*x(k) + D*u(:,k);
    end

    function out = x(k)

        if k == 1
            out = A*x0 + B*u0;
            return
        end
        out = A*x(k-1) + B*u(:,k-1);
        
    end

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

