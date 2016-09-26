function [Qc, q, beta_c, gamma] = qc_matrices(sys, N, gamma_max, z_fun)  
    
    Qc      = cell(N+1,1);
    q       = cell(N+1,1);
    beta_c  = cell(N+1,1);
    
    % Quadratic constraints
    for i = 0:N
        [Qz_sym, qz_fun, beta_fun] = sym_constraint(sys, N, i, z_fun);
        Qc(i+1)        = {Qz_sym};
        q(i+1)         = {qz_fun};  
        beta_c(i+1)    = {beta_fun};
    end
    
    % Selecting gamma in order to have an hessian matrix SDP and small gain condition
    if strcmp(z_fun, 'L2')
        gamma = get_gamma_l2_Qz_sdp(Qc, gamma_max);
        fprintf('Selected gamma: %f\n', gamma);
    elseif strcmp(z_fun, 'SOP')
        gamma = gamma_max;
    else
        error('qc_matrices: Invalid IQC selection (z_fun)');
    end
    
    % Use the choosed gamma
    for i = 0:N
        Qc(i+1)        = {Qc{i+1}(gamma)};
        q(i+1)         = {@(x0) q{i+1}(x0, gamma)};  
        beta_c(i+1)    = {@(z0, x0) beta_c{i+1}(z0, x0, gamma)};
    end
    
end

function gamma = get_gamma_l2_Qz_sdp(Qc, gamma0)

    gamma = gamma0;
    N = size(Qc, 1);
    i = 1;
    max_iter = 80;
    iter = max_iter;
    while i <= N
        if iter == 0
            error('get_gamma_l2_Qz_sdp: couldn''t find a suitable gamma, try again setting a greater max_iter value (current: %d)\n', max_iter);
        end
        if min(real(eig(Qc{i}(gamma)))) < 0
            gamma = gamma/2;
            iter = iter-1;
            continue
        end
        i = i+1;
    end
end


function [Qz_fun, qz_fun, beta_fun] = sym_constraint( sys, N, K, z_fun )

n = size(sys.B,1);
p = size(sys.B,2);

% z0 <= beta
if K == 0
    Qz_fun      = @(gamma) zeros(N*p);
    qz_fun      = @(x0, gamma) zeros(N*p, 1);
    beta_fun    = @(z0, x0, gamma) -z0;
    return
end

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
% l = subs(l,gamma, gamma_);

u_t = u';
u_var = [u0(:); u_t(:)];

% Calculate Qz
h = 1/2*hessian(l, u_var);
Qz_fun = matlabFunction(h, 'Vars', {gamma}); %double(h);

% Calculate qz
g = gradient(l, u_var);
g = subs(g, u_var, zeros(N,1));
qz_fun = matlabFunction(g, 'Vars', {x0, gamma});

% Calculate beta_
const = subs(l, u_var, zeros(N,1));
%const = subs(const, z0, x_ext0);
beta_fun = matlabFunction(-const, 'Vars', {z0, x0, gamma});                        % -const perche' in rhs della disuguaglianza
    

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
